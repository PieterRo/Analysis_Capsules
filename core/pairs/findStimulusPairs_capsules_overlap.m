function results = findStimulusPairs_capsules_overlap(folderPath, varargin)
% FINDSTIMULUSPAIRS_CAPSULES_OVERLAP
% Find stimulus pairs where ONE object is (nearly) identical in pixel space
% (position + orientation), and the OTHER object differs.
%
% Designed for bitmaps with:
%   - homogeneous grey background
%   - exactly two flat-colored objects (yellow & purple) + background
%   - minimal/no antialiasing (but tolerates small artifacts via settings)
%
% Output:
%   results.pairs      : Nx2 indices into results.files
%   results.files      : file names (sorted)
%   results.pairInfo   : per-pair overlap diagnostics (IoU + best shift)
%   results.params     : parameters used
%   results.masks      : per-image masks (yellow/purple) (can be large)
%
% Usage:
%   R = findStimulusPairs_capsules_overlap("/path/to/folder");
%   R = findStimulusPairs_capsules_overlap("/path", 'Anchor','yellow', ...
%        'AnchorIoUmin',0.97,'OtherIoUmax',0.60,'MaxShiftPx',3);
%
%   % show filenames for first pairs
%   disp([string(R.files(R.pairs(:,1))) string(R.files(R.pairs(:,2)))])
%
% Notes:
% - "IoU" is intersection-over-union (Jaccard). 1.0 means perfect overlap.
% - Allowing small shifts makes it tolerant to a few-pixel translation.

% -------------------- Parameters --------------------
p = inputParser;
p.addRequired('folderPath', @(s)ischar(s)||isstring(s));
p.addParameter('Anchor','yellow', @(s)any(strcmpi(s,{'yellow','purple'})));
p.addParameter('AnchorIoUmin', 0.97, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('OtherIoUmax',  0.60, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MaxShiftPx',   3,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinAreaPx',    200,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose',      true, @(x)islogical(x)&&isscalar(x));
p.addParameter('KeepMasks',    true, @(x)islogical(x)&&isscalar(x)); % can be large memory
p.addParameter('DebugPair',    false,@(x)islogical(x)&&isscalar(x)); % plot found pairs
p.addParameter('DebugEvery',   0,    @(x)isnumeric(x)&&isscalar(x)&&x>=0); % plot every N images, 0=off
p.parse(folderPath, varargin{:});
S = p.Results;

folderPath = char(folderPath);
anchorIsYellow = strcmpi(S.Anchor,'yellow');

% -------------------- Gather files --------------------
files = dir(fullfile(folderPath, '*.bmp'));
if isempty(files)
    error('No .bmp files found in folder: %s', folderPath);
end

names = {files.name};
[~, order] = sort_nat(names);
files = files(order);
names = {files.name};
n = numel(files);

if S.Verbose
    fprintf('Scanning %d images in: %s\n', n, folderPath);
    fprintf('Pair criterion: Anchor=%s, AnchorIoU>=%.3f, OtherIoU<=%.3f, MaxShift=%d px\n', ...
        S.Anchor, S.AnchorIoUmin, S.OtherIoUmax, S.MaxShiftPx);
end

% -------------------- Extract masks per image --------------------
if S.KeepMasks
    masks(n) = struct('yellow', [], 'purple', []);
else
    masks = []; %#ok<NASGU>
end

valid = true(n,1);

for i = 1:n
    I = imread(fullfile(folderPath, files(i).name));
    if size(I,3) == 1
        I = repmat(I,1,1,3);
    end

    try
        [yellow, purple] = segmentByExactColors(I);

        % cleanup (should be minimal for your stimuli, but helps if artifacts exist)
        yellow = bwareaopen(yellow, S.MinAreaPx);
        purple = bwareaopen(purple, S.MinAreaPx);
        yellow = keepLargestCC(yellow);
        purple = keepLargestCC(purple);

        if ~any(yellow(:)) || ~any(purple(:))
            valid(i) = false;
            if S.Verbose
                warning('Empty mask after cleanup in %s (yellow=%d, purple=%d).', ...
                    files(i).name, any(yellow(:)), any(purple(:)));
            end
        end

        if S.KeepMasks
            masks(i).yellow = yellow;
            masks(i).purple = purple;
        end

        if S.DebugEvery > 0 && mod(i, S.DebugEvery) == 0
            figure(101); clf;
            subplot(1,3,1); imshow(I); title(files(i).name, 'Interpreter','none');
            subplot(1,3,2); imshow(yellow); title('Yellow mask');
            subplot(1,3,3); imshow(purple); title('Purple mask');
            drawnow;
        end

    catch ME
        valid(i) = false;
        if S.Verbose
            warning('Could not segment %s: %s', files(i).name, ME.message);
        end
    end
end

if S.Verbose
    fprintf('Valid images with both masks: %d / %d\n', nnz(valid), n);
end

if ~S.KeepMasks
    error(['KeepMasks=false is not supported in this version because overlap pairing ' ...
           'needs masks. Set KeepMasks=true.']);
end

% -------------------- Pair search using overlap --------------------
pairs = [];
pairInfo = struct('anchorIoU', {}, 'anchorShift', {}, 'otherIoU', {}, 'otherShift', {});

for i = 1:n-1
    if ~valid(i), continue; end
    for j = i+1:n
        if ~valid(j), continue; end

        if anchorIsYellow
            A1 = masks(i).yellow;  A2 = masks(j).yellow;
            B1 = masks(i).purple;  B2 = masks(j).purple;
        else
            A1 = masks(i).purple;  A2 = masks(j).purple;
            B1 = masks(i).yellow;  B2 = masks(j).yellow;
        end

        [iouA, shiftA] = bestIoU_withShift(A1, A2, S.MaxShiftPx);
        if iouA < S.AnchorIoUmin
            continue;
        end

        [iouB, shiftB] = bestIoU_withShift(B1, B2, S.MaxShiftPx);
        if iouB > S.OtherIoUmax
            continue;
        end

        pairs(end+1,:) = [i j]; %#ok<AGROW>
        pairInfo(end+1).anchorIoU = iouA; %#ok<AGROW>
        pairInfo(end).anchorShift = shiftA;
        pairInfo(end).otherIoU = iouB;
        pairInfo(end).otherShift = shiftB;

        if S.DebugPair
            I1 = imread(fullfile(folderPath, files(i).name));
            I2 = imread(fullfile(folderPath, files(j).name));
            figure(202); clf;
            subplot(2,2,1); imshow(I1); title("A: " + string(files(i).name), 'Interpreter','none');
            subplot(2,2,2); imshow(I2); title("B: " + string(files(j).name), 'Interpreter','none');

            subplot(2,2,3); imshow(overlayMasks(I1, A1, B1));
            title(sprintf('A masks (anchor=%s)', S.Anchor));

            subplot(2,2,4); imshow(overlayMasks(I2, A2, B2));
            title(sprintf('B masks | AnchorIoU=%.3f (shift [%d %d]), OtherIoU=%.3f (shift [%d %d])', ...
                iouA, shiftA(1), shiftA(2), iouB, shiftB(1), shiftB(2)));
            drawnow;
        end
    end
end

if S.Verbose
    fprintf('Found %d matching pairs.\n', size(pairs,1));
end

% -------------------- Package results --------------------
results = struct();
results.pairs = pairs;
results.files = names;
results.pairInfo = pairInfo;
results.params = S;
results.masks = masks;

end

% ==================== Helper functions ====================

function [yellowMask, purpleMask, bgRGB, objRGB] = segmentByExactColors(I)
% Segment yellow & purple objects by exact RGB colors.
% Assumes 3 flat colors total: background + two objects.
% Robustly identifies background via border median.

Iu = uint8(I);
if size(Iu,3) == 1
    Iu = repmat(Iu,1,1,3);
end

% Estimate background from border median
b = 10;
border = cat(1, ...
    reshape(Iu(1:b,:,:), [], 3), ...
    reshape(Iu(end-b+1:end,:,:), [], 3), ...
    reshape(Iu(:,1:b,:), [], 3), ...
    reshape(Iu(:,end-b+1:end,:), [], 3));
bgRGB = uint8(round(median(double(border), 1))); % 1x3

% Unique colors in image
cols = unique(reshape(Iu, [], 3), "rows");

% Remove background row (exact match). If not exact (rare), remove closest.
isBg = all(cols == bgRGB, 2);
if ~any(isBg)
    % pick closest color to bgRGB
    d = sum((double(cols) - double(bgRGB)).^2, 2);
    [~, k] = min(d);
    isBg = false(size(cols,1),1);
    isBg(k) = true;
end

objRGB = cols(~isBg, :);

if size(objRGB,1) ~= 2
    error("Expected exactly 2 object colors, found %d. (Total unique colors=%d)", ...
        size(objRGB,1), size(cols,1));
end

% Decide purple vs yellow:
% purple has higher Blue channel than yellow in your stimuli.
[~, idxPurple] = max(objRGB(:,3)); % max blue
idxYellow = 3 - idxPurple;

purpleRGB = objRGB(idxPurple,:);
yellowRGB = objRGB(idxYellow,:);

purpleMask = all(Iu == reshape(purpleRGB,1,1,3), 3);
yellowMask = all(Iu == reshape(yellowRGB,1,1,3), 3);
end

function BW = keepLargestCC(BW)
BW = logical(BW);
cc = bwconncomp(BW);
if cc.NumObjects <= 1
    return;
end
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, idx] = max(numPixels);
BW = false(size(BW));
BW(cc.PixelIdxList{idx}) = true;
end

function [bestIoU, bestShift] = bestIoU_withShift(A, B, maxShift)
% Compute maximum IoU between masks A and B allowing integer shifts in [-maxShift..maxShift].
% bestShift = [dx dy] applied to B (shift B to match A)

A = logical(A); B = logical(B);

bestIoU = -inf;
bestShift = [0 0];

for dy = -maxShift:maxShift
    for dx = -maxShift:maxShift
        Bs = shiftMask(B, dx, dy);

        inter = nnz(A & Bs);
        uni   = nnz(A | Bs);

        if uni == 0
            iou = 0;
        else
            iou = inter / uni;
        end

        if iou > bestIoU
            bestIoU = iou;
            bestShift = [dx dy];
        end
    end
end
end

function M = shiftMask(M0, dx, dy)
% Shift binary mask by (dx,dy) pixels. Positive dx shifts right, positive dy shifts down.
% Pads with false; no wrap-around.

M0 = logical(M0);
[h,w] = size(M0);
M = false(h,w);

x1 = max(1, 1+dx);
x2 = min(w, w+dx);
y1 = max(1, 1+dy);
y2 = min(h, h+dy);

x1s = max(1, 1-dx);
x2s = min(w, w-dx);
y1s = max(1, 1-dy);
y2s = min(h, h-dy);

if x1 <= x2 && y1 <= y2
    M(y1:y2, x1:x2) = M0(y1s:y2s, x1s:x2s);
end
end

function out = overlayMasks(I, A, B)
% Simple overlay: anchor in green, other in red on top of original image.
Iu = im2uint8(I);
if size(Iu,3)==1, Iu=repmat(Iu,1,1,3); end
out = Iu;

A = logical(A);
B = logical(B);

% add highlights (no fancy alpha; just boost channels)
out(:,:,2) = uint8(min(255, double(out(:,:,2)) + 120*double(A))); % green
out(:,:,1) = uint8(min(255, double(out(:,:,1)) + 120*double(B))); % red
end

function [sorted, idx] = sort_nat(c)
% Natural-ish sort for filenames containing digits (e.g., 001.bmp, 10.bmp)
% Falls back to normal sort if parsing fails.
[sorted, idx] = sort(c);
try
    expr = '\d+';
    tokens = regexp(c, expr, 'match');
    nums = cellfun(@(t) iff(isempty(t), NaN, str2double(t{end})), tokens);
    [~, idx] = sortrows([isnan(nums(:)) nums(:)], [1 2]);
    sorted = c(idx);
catch
    % fallback already done
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
