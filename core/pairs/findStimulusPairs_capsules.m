function results = findStimulusPairs_capsules(folderPath, varargin)
% FINDSTIMULUSPAIRS_CAPSULES
% Find image pairs where one object's position is (nearly) identical
% (e.g., yellow) while the other object's position differs (e.g., purple).
%
% Assumptions:
% - Background is uniform grey
% - Two colored capsule-shaped objects: one yellow, one purple
% - Bitmaps named like '001.bmp' ... '328.bmp' (but code scans *.bmp anyway)
%
% Output:
% results struct with fields:
%   .pairs : Nx2 indices into file list
%   .files : file names
%   .feat  : per-image features
%   .params: settings used
%
% Usage:
%   results = findStimulusPairs_capsules('/path/to/folder');
%   results = findStimulusPairs_capsules('/path', 'Anchor','yellow');
%
% Optional name-value params:
%   'Anchor'              : 'yellow' or 'purple' (default 'yellow')
%   'SameCentroidTolPx'   : tolerance in pixels for "same position" (default 3)
%   'DifferentTolPx'      : minimum centroid difference in pixels for "different" (default 8)
%   'MinAreaPx'           : min area for object region (default 80)
%   'Verbose'             : true/false (default true)
%   'DebugPlot'           : true/false (default false)

% -------------------- Parameters --------------------
p = inputParser;
p.addRequired('folderPath', @(s)ischar(s)||isstring(s));
p.addParameter('Anchor','yellow', @(s)any(strcmpi(s,{'yellow','purple'})));
p.addParameter('SameCentroidTolPx', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('DifferentTolPx', 8, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinAreaPx', 80, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('DebugPlot', false, @(x)islogical(x)&&isscalar(x));
p.parse(folderPath, varargin{:});
S = p.Results;

folderPath = char(folderPath);

% -------------------- Gather files --------------------
files = dir(fullfile(folderPath, '*.bmp'));
if isempty(files)
    error('No .bmp files found in folder: %s', folderPath);
end

% Natural-sort if your filenames are numeric (001.bmp...)
names = {files.name};
[~, order] = sort_nat(names);
files = files(order);
names = {files.name};

n = numel(files);

% Preallocate feature arrays
feat(n) = struct( ...
    'yellowCentroid', [NaN NaN], ...
    'purpleCentroid', [NaN NaN], ...
    'yellowArea', NaN, ...
    'purpleArea', NaN, ...
    'yellowOrient', NaN, ...
    'purpleOrient', NaN ...
);

if S.Verbose
    fprintf('Scanning %d images in: %s\n', n, folderPath);
end

% -------------------- Extract features --------------------
for i = 1:n
    I = imread(fullfile(folderPath, files(i).name));
    if size(I,3) ~= 3
        error('Expected RGB image for %s', files(i).name);
    end


    
    
    
    % --- Robust background-based foreground mask (works for low-saturation objects) ---
    Iu = uint8(I);

    % Estimate background RGB from border median (very stable for uniform background)
    b = 10; % border width in pixels
    border = cat(1, ...
    reshape(Iu(1:b,:,:), [], 3), ...
    reshape(Iu(end-b+1:end,:,:), [], 3), ...
    reshape(Iu(:,1:b,:), [], 3), ...
    reshape(Iu(:,end-b+1:end,:), [], 3));
    bg = uint8(round(median(double(border), 1)));  % 1x3

    % Foreground = pixels different from background (use small tolerance for safety)
    tol = 2; % increase to ~5 if you have antialiasing/compression
    d = sqrt(sum((double(Iu) - double(reshape(bg,1,1,3))).^2, 3));
    fg = d > tol;

    % Now classify fg pixels into yellow vs purple using HSV hue (works even at low sat)
    Ihsv = rgb2hsv(Iu);
    H = Ihsv(:,:,1);
    V = Ihsv(:,:,3);

    % Broad hue ranges (based on your example: yellow ~0.19, purple/blue ~0.67)
    yellow = fg & (H > 0.06 & H < 0.28) & (V > 0.35);
    purple = fg & (H > 0.55 & H < 0.85) & (V > 0.35);
    
    % Cleanup
    yellow = bwareaopen(imfill(yellow,'holes'), S.MinAreaPx);
    purple = bwareaopen(imfill(purple,'holes'), S.MinAreaPx);

    % Keep the largest component (in case of tiny artifacts)
    yellow = keepLargestCC(yellow);
    purple = keepLargestCC(purple);

    % Regionprops
    ry = regionprops(yellow, 'Centroid','Area','Orientation');
    rp = regionprops(purple, 'Centroid','Area','Orientation');

    if isempty(ry) || isempty(rp)
        warning('Could not detect both objects in %s (yellow=%d, purple=%d).', ...
            files(i).name, ~isempty(ry), ~isempty(rp));
    else
        feat(i).yellowCentroid = ry(1).Centroid;
        feat(i).yellowArea     = ry(1).Area;
        feat(i).yellowOrient   = ry(1).Orientation;

        feat(i).purpleCentroid = rp(1).Centroid;
        feat(i).purpleArea     = rp(1).Area;
        feat(i).purpleOrient   = rp(1).Orientation;
    end

    if S.DebugPlot && mod(i,20)==1
        figure(1); clf;
        imshow(I); hold on;
        plot(feat(i).yellowCentroid(1), feat(i).yellowCentroid(2), 'y+', 'MarkerSize', 12, 'LineWidth',2);
        plot(feat(i).purpleCentroid(1), feat(i).purpleCentroid(2), 'm+', 'MarkerSize', 12, 'LineWidth',2);
        title(sprintf('%s', files(i).name), 'Interpreter','none');
        drawnow;
    end
end

% -------------------- Pairing logic --------------------
% thresholds (tune these)
anchorIoU_min = 0.95;     % "almost all pixels overlap"
otherIoU_max  = 0.60;     % "clearly different"
maxShiftPx    = 3;        % allow +/- a few px

pairs = [];
pairInfo = []; % store overlap scores and best shifts

for i = 1:n-1
    for j = i+1:n

        if anchorIsYellow
            A1 = masks(i).yellow;  A2 = masks(j).yellow;
            B1 = masks(i).purple;  B2 = masks(j).purple;
        else
            A1 = masks(i).purple;  A2 = masks(j).purple;
            B1 = masks(i).yellow;  B2 = masks(j).yellow;
        end

        if isempty(A1) || isempty(A2) || isempty(B1) || isempty(B2)
            continue;
        end

        % Best possible overlap under small translations
        [iouA, bestShiftA] = bestIoU_withShift(A1, A2, maxShiftPx);
        [iouB, bestShiftB] = bestIoU_withShift(B1, B2, maxShiftPx);

        if iouA >= anchorIoU_min && iouB <= otherIoU_max
            pairs(end+1,:) = [i j]; %#ok<AGROW>
            pairInfo(end+1) = struct( ...
                'anchorIoU', iouA, 'anchorShift', bestShiftA, ...
                'otherIoU',  iouB, 'otherShift',  bestShiftB ); %#ok<AGROW>
        end
    end
end

results.pairs = pairs;
results.pairInfo = pairInfo;

if S.Verbose
    fprintf('Found %d matching pairs (Anchor=%s, same<=%.1fpx, other>=%.1fpx)\n', ...
        size(pairs,1), S.Anchor, S.SameCentroidTolPx, S.DifferentTolPx);
end

% -------------------- Package results --------------------
results = struct();
results.pairs  = pairs;
results.files  = names;
results.feat   = feat;
results.params = S;

end

% ==================== Helper functions ====================

function BW = keepLargestCC(BW)
cc = bwconncomp(BW);
if cc.NumObjects <= 1
    return;
end
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, idx] = max(numPixels);
BW = false(size(BW));
BW(cc.PixelIdxList{idx}) = true;
end

function [sorted, idx] = sort_nat(c)
% Natural sort of filenames (001.bmp, 2.bmp, 10.bmp)
% Returns sorted cellstr and indices.
[sorted, idx] = sort(c);
try
    expr = '\d+';
    tokens = regexp(c, expr, 'match');
    nums = cellfun(@(t) iff(isempty(t), NaN, str2double(t{end})), tokens);
    [~, idx] = sortrows([isnan(nums(:)) nums(:)], [1 2]);
    sorted = c(idx);
catch
    % fall back to simple sort
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
