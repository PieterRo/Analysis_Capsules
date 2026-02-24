function results = findStimulusPairs_capsules_overlap_fast(folderPath, varargin)
% FINDSTIMULUSPAIRS_CAPSULES_OVERLAP_FAST
% Optimized pairing:
%  1) Segment masks by exact RGB colors (background + 2 objects)
%  2) Precompute regionprops (centroid + bounding box) per object per image
%  3) Pair search uses:
%       - cheap centroid pre-filter (reject most pairs fast)
%       - IoU computed only on small cropped ROIs (bounding-box IoU)
%       - small shift search in [-MaxShiftPx..MaxShiftPx]
%
% Output:
%   results.pairs      : Nx2 indices into results.files
%   results.files      : file names (sorted)
%   results.pairInfo   : per-pair diagnostics (IoU + best shifts)
%   results.params     : parameters used
%
% Usage:
%   folder = "/path/to/bitmaps";
%   R = findStimulusPairs_capsules_overlap_fast(folder, ...
%        'Anchor','yellow', 'AnchorIoUmin',0.99, 'OtherIoUmax',0.60, ...
%        'MaxShiftPx',3, 'MinAreaPx',200);
%
% Notes:
% - If your edges are pixel-sharp, set AnchorIoUmin closer to 0.995–1.0.
% - If there is any antialiasing/compression, use 0.95–0.98.

% -------------------- Parameters --------------------
p = inputParser;
p.addRequired('folderPath', @(s)ischar(s)||isstring(s));
p.addParameter('Anchor','yellow', @(s)any(strcmpi(s,{'yellow','purple'})));
p.addParameter('AnchorIoUmin', 0.98, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('OtherIoUmax',  0.60, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MaxShiftPx',   3,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinAreaPx',    200,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Verbose',      true, @(x)islogical(x)&&isscalar(x));
p.addParameter('KeepMasks',    false,@(x)islogical(x)&&isscalar(x)); % optional (memory)
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
    fprintf('Anchor=%s | AnchorIoU>=%.3f | OtherIoU<=%.3f | MaxShift=%d px\n', ...
        S.Anchor, S.AnchorIoUmin, S.OtherIoUmax, S.MaxShiftPx);
end

% -------------------- Precompute masks + props --------------------
% Store minimal info needed for fast pairing
% props(i).yellow / props(i).purple each has: Centroid [x y], Box [x1 y1 x2 y2]
props(n) = struct('yellow', [], 'purple', [], 'valid', false);

% Optional store masks (big). Not needed for final output.
if S.KeepMasks
    masks(n) = struct('yellow', [], 'purple', []);
else
    masks = []; %#ok<NASGU>
end

for i = 1:n
    I = imread(fullfile(folderPath, files(i).name));
    if size(I,3) == 1
        I = repmat(I,1,1,3);
    end

    try
        [yMask, pMask] = segmentByExactColors(I);

        % Cleanup
        yMask = bwareaopen(yMask, S.MinAreaPx);
        pMask = bwareaopen(pMask, S.MinAreaPx);
        yMask = keepLargestCC(yMask);
        pMask = keepLargestCC(pMask);

        if ~any(yMask(:)) || ~any(pMask(:))
            props(i).valid = false;
            continue;
        end

        ry = regionprops(yMask, 'Centroid','BoundingBox');
        rp = regionprops(pMask, 'Centroid','BoundingBox');
        if isempty(ry) || isempty(rp)
            props(i).valid = false;
            continue;
        end

        props(i).yellow = packProps(ry(1), size(yMask));
        props(i).purple = packProps(rp(1), size(pMask));
        props(i).valid = true;

        if S.KeepMasks
            masks(i).yellow = yMask;
            masks(i).purple = pMask;
        end

    catch
        props(i).valid = false;
    end
end

valid = [props.valid]';
if S.Verbose
    fprintf('Valid images with both objects: %d / %d\n', nnz(valid), n);
end

% -------------------- Pair search (fast) --------------------
pairs = zeros(0,2);
pairInfo = struct('anchorIoU', {}, 'anchorShift', {}, 'otherIoU', {}, 'otherShift', {});

% Open images on demand and segment again ONLY for candidate pairs.
% (This avoids storing 328 full masks in memory while still being fast, because
% only a small subset of pairs pass the centroid pre-filter.)
for i = 1:n-1
    if ~props(i).valid, continue; end
    for j = i+1:n
        if ~props(j).valid, continue; end

        if anchorIsYellow
            ci = props(i).yellow.Centroid; cj = props(j).yellow.Centroid;
        else
            ci = props(i).purple.Centroid; cj = props(j).purple.Centroid;
        end

        % ---- Cheap centroid pre-filter ----
        % If centroids are farther than MaxShiftPx+1, there is no way to get high IoU.
        if norm(ci - cj) > (S.MaxShiftPx + 1)
            continue;
        end

        % Load and segment only now (candidate)
        if S.KeepMasks
            y1 = masks(i).yellow; p1 = masks(i).purple;
            y2 = masks(j).yellow; p2 = masks(j).purple;
        else
            I1 = imread(fullfile(folderPath, files(i).name));
            I2 = imread(fullfile(folderPath, files(j).name));
            if size(I1,3)==1, I1 = repmat(I1,1,1,3); end
            if size(I2,3)==1, I2 = repmat(I2,1,1,3); end

            [y1,p1] = segmentByExactColors(I1);
            [y2,p2] = segmentByExactColors(I2);

            y1 = keepLargestCC(bwareaopen(y1, S.MinAreaPx));
            p1 = keepLargestCC(bwareaopen(p1, S.MinAreaPx));
            y2 = keepLargestCC(bwareaopen(y2, S.MinAreaPx));
            p2 = keepLargestCC(bwareaopen(p2, S.MinAreaPx));
        end

        if anchorIsYellow
            A1 = y1; A2 = y2;   B1 = p1; B2 = p2;
            boxA1 = props(i).yellow.Box; boxA2 = props(j).yellow.Box;
            boxB1 = props(i).purple.Box; boxB2 = props(j).purple.Box;
        else
            A1 = p1; A2 = p2;   B1 = y1; B2 = y2;
            boxA1 = props(i).purple.Box; boxA2 = props(j).purple.Box;
            boxB1 = props(i).yellow.Box; boxB2 = props(j).yellow.Box;
        end

        % ---- Fast IoU on cropped ROIs (anchor first) ----
        [iouA, shiftA] = bestIoU_withShift_bbox(A1, boxA1, A2, boxA2, S.MaxShiftPx);
        if iouA < S.AnchorIoUmin
            continue;
        end

        [iouB, shiftB] = bestIoU_withShift_bbox(B1, boxB1, B2, boxB2, S.MaxShiftPx);
        if iouB > S.OtherIoUmax
            continue;
        end

        pairs(end+1,:) = [i j]; %#ok<AGROW>
        pairInfo(end+1).anchorIoU = iouA; %#ok<AGROW>
        pairInfo(end).anchorShift = shiftA;
        pairInfo(end).otherIoU = iouB;
        pairInfo(end).otherShift = shiftB;
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

if S.KeepMasks
    results.masks = masks;
end
results.props = props; % handy for debugging / grouping

end

% ==================== Helper functions ====================

function P = packProps(r, imSize)
% r: regionprops struct with Centroid and BoundingBox [x y w h]
bb = r.BoundingBox;
x1 = max(1, floor(bb(1)));
y1 = max(1, floor(bb(2)));
x2 = min(imSize(2), ceil(bb(1) + bb(3)));
y2 = min(imSize(1), ceil(bb(2) + bb(4)));
P = struct('Centroid', r.Centroid, 'Box', [x1 y1 x2 y2]);
end

function [yellowMask, purpleMask, bgRGB, objRGB] = segmentByExactColors(I)
% Segment yellow & purple objects by exact RGB colors.
% Assumes: background + two object colors (usually exactly 3 unique colors).
Iu = uint8(I);
if size(Iu,3) == 1
    Iu = repmat(Iu,1,1,3);
end

% Background from border median
b = 10;
border = cat(1, ...
    reshape(Iu(1:b,:,:), [], 3), ...
    reshape(Iu(end-b+1:end,:,:), [], 3), ...
    reshape(Iu(:,1:b,:), [], 3), ...
    reshape(Iu(:,end-b+1:end,:), [], 3));
bgRGB = uint8(round(median(double(border), 1))); % 1x3

cols = unique(reshape(Iu, [], 3), "rows");

% Remove background: exact match if present, otherwise closest color
isBg = all(cols == bgRGB, 2);
if ~any(isBg)
    d = sum((double(cols) - double(bgRGB)).^2, 2);
    [~, k] = min(d);
    isBg = false(size(cols,1),1);
    isBg(k) = true;
end
objRGB = cols(~isBg, :);

if size(objRGB,1) ~= 2
    error("Expected exactly 2 object colors, found %d (unique colors total=%d).", ...
        size(objRGB,1), size(cols,1));
end

% Decide purple vs yellow by blue channel (works for your stimuli)
[~, idxPurple] = max(objRGB(:,3)); % max blue = purple
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

function [bestIoU, bestShift] = bestIoU_withShift_bbox(A, boxA, B, boxB, maxShift)
% Best IoU between A and B, allowing small shifts, computed on cropped ROI
% boxA/boxB: [x1 y1 x2 y2] in their respective images (full image coords).

A = logical(A); B = logical(B);
szA = size(A);

% Create a union ROI in A-coordinates big enough to hold both (plus padding).
pad = maxShift + 2;

% Expand boxes in their own images
roiA = expandBoxXY(boxA, szA, pad);
roiB = expandBoxXY(boxB, size(B), pad);

% Extract crops
Ac = A(roiA(2):roiA(4), roiA(1):roiA(3));
Bc0 = B(roiB(2):roiB(4), roiB(1):roiB(3));

% Resize/crop both to common canvas by embedding into same-size canvas
% (this avoids mismatched crop sizes when boxes differ)
h = max(size(Ac,1), size(Bc0,1)) + 2*pad;
w = max(size(Ac,2), size(Bc0,2)) + 2*pad;
Acan = false(h,w);
Bcan0 = false(h,w);

Acan(1:size(Ac,1), 1:size(Ac,2)) = Ac;
Bcan0(1:size(Bc0,1), 1:size(Bc0,2)) = Bc0;

bestIoU = -inf;
bestShift = [0 0];

for dy = -maxShift:maxShift
    for dx = -maxShift:maxShift
        Bcan = shiftMask(Bcan0, dx, dy);

        inter = nnz(Acan & Bcan);
        uni   = nnz(Acan | Bcan);
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

function box = expandBoxXY(b, sz, pad)
% b = [x1 y1 x2 y2] (already x1/y1/x2/y2)
x1 = max(1, b(1)-pad);
y1 = max(1, b(2)-pad);
x2 = min(sz(2), b(3)+pad);
y2 = min(sz(1), b(4)+pad);
box = [x1 y1 x2 y2];
end

function M = shiftMask(M0, dx, dy)
% Shift binary mask by (dx,dy). No wrap-around; pads with false.
M0 = logical(M0);
[h,w] = size(M0);
M = false(h,w);

x1 = max(1, 1+dx);  x2 = min(w, w+dx);
y1 = max(1, 1+dy);  y2 = min(h, h+dy);

x1s = max(1, 1-dx); x2s = min(w, w-dx);
y1s = max(1, 1-dy); y2s = min(h, h-dy);

if x1 <= x2 && y1 <= y2
    M(y1:y2, x1:x2) = M0(y1s:y2s, x1s:x2s);
end
end

function [sorted, idx] = sort_nat(c)
% Natural-ish sort for filenames containing digits (e.g., 001.bmp, 10.bmp)
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
