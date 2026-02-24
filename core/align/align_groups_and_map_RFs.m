function OUT = align_groups_and_map_RFs(folder, xRF, yRF, varargin)
% ALIGN_GROUPS_AND_MAP_RFS
%
% Bitmaps come in groups of 4: 001-004, 005-008, ..., 325-328.
% One group of 4 (default 061-064) is the template.
%
% For each other group:
%  - find best alignment to the corresponding template bitmap via:
%       optional mirroring (none / fliplr / flipud)
%       + rigid transform (rotation+translation)
%  - estimate transform from the 1st bitmap in the group, then CHECK the other 3.
%    If the check is poor (IoU < MismatchIoU) you can refit per-bitmap.
%
% RF centers:
%  - You provide ONE global RF set (xRF,yRF) defined in the *untransformed* bitmap space.
%  - The same transform found for a bitmap is applied to ALL RFs (vectorized).
%  - RFs are output in TEMPLATE pixel space, plus shifted coords (x-512,y-384).
%
% Classification for each (bitmap, RF):
%  - label: yellow / purple / background (tested against template masks)
%  - target/distractor: target = occluder = object with more visible pixels in template
%
% Polar coords in template:
%  - origin = stimulus-center (midpoint between closest boundaries of the two objects)
%  - angle=0 along target major axis
%  - theta_deg = RF angle relative to target axis
%  - rho = distance to stimulus-center
%
% OUTPUT:
%  OUT.templates(t)         : per-template info (masks, stimCenter, targetAxisAngleDeg, etc.)
%  OUT.perBitmap(b)         : transform + mirror used to map bitmap b into its template
%  OUT.RFmap                : table, one row per (bitmap, RFid)
%  OUT.template_distractorAnglesDeg : distractor direction estimate per template (check ~60°)
%
% REQUIREMENTS: Image Processing Toolbox.

% ---------------- parse options ----------------
p = inputParser;
p.addRequired('folder', @(s)ischar(s)||isstring(s));
p.addRequired('xRF', @(x)isnumeric(x)&&isvector(x));
p.addRequired('yRF', @(x)isnumeric(x)&&isvector(x)&&numel(x)==numel(xRF));

p.addParameter('TemplateFirst', 61, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('GroupSize', 4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('NumBitmaps', 328, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ImageSize', [768 1024], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('TryMirrorModes', {'none','fliplr','flipud'}, @(c)iscell(c)&&~isempty(c));
p.addParameter('TryRot90', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MinAreaPx', 200, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxRefineIter', 200, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter('CheckOther3', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('RefitIfMismatch', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MismatchIoU', 0.90, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

p.parse(folder, xRF, yRF, varargin{:});
S = p.Results;

folder = char(folder);
xRF = xRF(:); yRF = yRF(:);
N = numel(xRF);

H = S.ImageSize(1); W = S.ImageSize(2);

% ---------------- file mapping (###.bmp -> number) ----------------
bmpFiles = dir(fullfile(folder, '*.bmp'));
if isempty(bmpFiles), error('No .bmp files found in %s', folder); end

names = {bmpFiles.name};
[~, ord] = sort_nat(names);
bmpFiles = bmpFiles(ord);
names = {bmpFiles.name};

fileOf = strings(S.NumBitmaps,1);
for i = 1:numel(names)
    tok = regexp(names{i}, '(\d+)\.bmp$', 'tokens', 'once');
    if isempty(tok), continue; end
    n = str2double(tok{1});
    if n>=1 && n<=S.NumBitmaps
        fileOf(n) = string(names{i});
    end
end

missing = find(fileOf=="");
if ~isempty(missing) && S.Verbose
    fprintf('Note: %d bitmap numbers missing files (e.g. %s ...)\n', ...
        numel(missing), strjoin(string(missing(1:min(5,end))), ', '));
end

% ---------------- load templates ----------------
tplFirst = S.TemplateFirst;
G = S.GroupSize;
tplNums  = tplFirst:(tplFirst+G-1);

templates = struct([]);
for t = 1:G
    b = tplNums(t);
    assert(fileOf(b)~="", 'Template bitmap %d not found in folder.', b);

    I = imread(fullfile(folder, fileOf(b)));
    [yMask, pMask, bgRGB, yRGB, pRGB] = segmentByExactColors(I, S.MinAreaPx);

    templates(t).bmpNum = b;
    templates(t).file   = fileOf(b);
    templates(t).I      = I;
    templates(t).yMask  = yMask;
    templates(t).pMask  = pMask;
    templates(t).union  = yMask | pMask;
    templates(t).bgRGB  = bgRGB;
    templates(t).yRGB   = yRGB;
    templates(t).pRGB   = pRGB;

    templates(t).targetIsYellow = nnz(yMask) > nnz(pMask); % occluder (more visible pixels)
    [templates(t).stimCenter_xy, templates(t).targetAxisAngleDeg, templates(t).distrAngleDeg] = ...
        computeTemplateGeometry(templates(t));
end

if S.Verbose
    fprintf('Templates: %s\n', strjoin(string(tplNums), ', '));
    fprintf('Template distractor angles (deg, relative to target axis): %s\n', ...
        strjoin(string(round([templates.distrAngleDeg],1)), ', '));
end

% ---------------- registration config ----------------
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = S.MaxRefineIter;
optimizer.MinimumStepLength = 1e-6;

% ---------------- per-bitmap transforms ----------------
perBitmap = repmat(struct( ...
    'bmpNum',[], 'file',"", 'templateIdx',[], ...
    'mirrorMode',"", 'tform',[], 'iouToTemplate',NaN, 'notes',""), S.NumBitmaps, 1);

% Fill template group as identity
for k = 1:G
    b = tplNums(k);
    perBitmap(b).bmpNum = b;
    perBitmap(b).file = fileOf(b);
    perBitmap(b).templateIdx = k;
    perBitmap(b).mirrorMode = "none";
    perBitmap(b).tform = affine2d(eye(3));
    perBitmap(b).iouToTemplate = 1.0;
    perBitmap(b).notes = "template";
end

% Fit all other groups
nGroups = ceil(S.NumBitmaps / G);

for g = 1:nGroups
    grpStart = (g-1)*G + 1;
    grpNums = grpStart:(grpStart+G-1);
    grpNums = grpNums(grpNums<=S.NumBitmaps);

    % skip if this is the template group
    if ~isempty(intersect(grpNums, tplNums))
        continue;
    end

    b1 = grpNums(1);
    if fileOf(b1)==""
        continue;
    end

    % Fit from first bitmap in group to templateIdx=1
    I1 = imread(fullfile(folder, fileOf(b1)));
    [y1,p1] = segmentByExactColors(I1, S.MinAreaPx);
    u1 = y1 | p1;
    fixed = templates(1).union;

    best = findBestRigidTransform(u1, fixed, S.TryMirrorModes, S.TryRot90, optimizer, metric);

    if S.Verbose
        fprintf('Group %d (%d-%d): first=%d mirror=%s IoU=%.3f\n', ...
            g, grpNums(1), grpNums(end), b1, best.mirrorMode, best.iou);
    end

    for kk = 1:numel(grpNums)
        b = grpNums(kk);
        if fileOf(b)=="" || b < 1 || b > S.NumBitmaps
            continue;
        end
        tplIdx = kk;  % 1..4 alignment to 061..064 by position in group
        if tplIdx > G
            continue;
        end

        perBitmap(b).bmpNum = b;
        perBitmap(b).file = fileOf(b);
        perBitmap(b).templateIdx = tplIdx;

        Itpl = templates(tplIdx).union;

        I = imread(fullfile(folder, fileOf(b)));
        [yy,pp] = segmentByExactColors(I, S.MinAreaPx);
        u = yy | pp;

        % apply best-from-first
        uMoved = applyMirror(u, best.mirrorMode);
        uWarp  = imwarp(uMoved, best.tform, 'OutputView', imref2d(size(Itpl)));
        iou = IoU(uWarp, Itpl);

        tformUse = best.tform;
        mirrorUse = best.mirrorMode;
        note = "from-first";

        if S.CheckOther3 && S.RefitIfMismatch && (iou < S.MismatchIoU)
            best2 = findBestRigidTransform(u, Itpl, S.TryMirrorModes, S.TryRot90, optimizer, metric);
            uMoved2 = applyMirror(u, best2.mirrorMode);
            uWarp2  = imwarp(uMoved2, best2.tform, 'OutputView', imref2d(size(Itpl)));
            iou2 = IoU(uWarp2, Itpl);

            if iou2 > iou
                iou = iou2;
                tformUse = best2.tform;
                mirrorUse = best2.mirrorMode;
                note = sprintf("refit (was %.3f)", iou);
            end
        end

        perBitmap(b).tform = tformUse;
        perBitmap(b).mirrorMode = mirrorUse;
        perBitmap(b).iouToTemplate = iou;
        perBitmap(b).notes = note;
    end
end

% ---------------- map ALL RFs across ALL bitmaps (vectorized per bitmap) ----------------
RFrows = cell(S.NumBitmaps,1);

for b = 1:S.NumBitmaps
    if isempty(perBitmap(b).tform) || isempty(perBitmap(b).templateIdx) || perBitmap(b).file==""
        continue;
    end
    tplIdx = perBitmap(b).templateIdx;

    % mirror RFs in bitmap space
    [xm, ym] = mirrorPointVec(xRF, yRF, perBitmap(b).mirrorMode, [H W]);

    % affine transform to template space
    [xt, yt] = transformPointsForward(perBitmap(b).tform, xm, ym);

    % clamp to template bounds
    xi = round(min(max(xt,1), W));
    yi = round(min(max(yt,1), H));
    ind = sub2ind([H W], yi, xi);

    % classify using template masks (fast)
    yMask = templates(tplIdx).yMask;
    pMask = templates(tplIdx).pMask;

    label = repmat("background", N, 1);
    label(yMask(ind)) = "yellow";
    label(pMask(ind)) = "purple";

    % polar relative to stimulus center and target axis
    sc = templates(tplIdx).stimCenter_xy; % [x y]
    vx = xt - sc(1);
    vy = yt - sc(2);
    rho = hypot(vx, vy);
    theta_deg = wrapTo180(atan2d(vy, vx) - templates(tplIdx).targetAxisAngleDeg);

    % shift coords so (512,384) is screen center
    x_shift = xt - 512;
    y_shift = yt - 384;

    RFrows{b} = table( ...
        repmat(b,N,1), (1:N)', repmat(tplIdx,N,1), ...
        xt, yt, x_shift, y_shift, ...
        label, rho, theta_deg, ...
        'VariableNames', {'bmp','rfId','templateIdx','x_tpl','y_tpl','x_tpl_shift','y_tpl_shift','label','rho','theta_deg'});
end

RFmap = vertcat(RFrows{:});

% ---------------- pack output ----------------
OUT = struct();
OUT.templates = templates;
OUT.perBitmap = perBitmap;
OUT.RFmap = RFmap;
OUT.template_distractorAnglesDeg = [templates.distrAngleDeg];

end

% =====================================================================
% Helpers
% =====================================================================

function [yMask, pMask, bgRGB, yRGB, pRGB] = segmentByExactColors(I, minArea)
Iu = uint8(I);
if size(Iu,3)==1, Iu = repmat(Iu,1,1,3); end

% background from border median
b = 10;
border = cat(1, ...
    reshape(Iu(1:b,:,:), [], 3), ...
    reshape(Iu(end-b+1:end,:,:), [], 3), ...
    reshape(Iu(:,1:b,:), [], 3), ...
    reshape(Iu(:,end-b+1:end,:), [], 3));
bgRGB = uint8(round(median(double(border), 1)));

cols = unique(reshape(Iu, [], 3), "rows");

% remove background (exact or closest)
isBg = all(cols == bgRGB, 2);
if ~any(isBg)
    d = sum((double(cols) - double(bgRGB)).^2, 2);
    [~,k] = min(d);
    isBg = false(size(cols,1),1);
    isBg(k) = true;
end
objRGB = cols(~isBg,:);

if size(objRGB,1) ~= 2
    error('Expected 2 object colors; found %d (unique colors total=%d).', size(objRGB,1), size(cols,1));
end

% decide purple vs yellow by blue channel (purple has larger B)
[~, idxPurple] = max(objRGB(:,3));
idxYellow = 3 - idxPurple;
pRGB = objRGB(idxPurple,:);
yRGB = objRGB(idxYellow,:);

pMask = all(Iu == reshape(pRGB,1,1,3), 3);
yMask = all(Iu == reshape(yRGB,1,1,3), 3);

% cleanup
pMask = keepLargestCC(bwareaopen(pMask, minArea));
yMask = keepLargestCC(bwareaopen(yMask, minArea));
end

function BW = keepLargestCC(BW)
BW = logical(BW);
cc = bwconncomp(BW);
if cc.NumObjects <= 1, return; end
np = cellfun(@numel, cc.PixelIdxList);
[~,k] = max(np);
BW = false(size(BW));
BW(cc.PixelIdxList{k}) = true;
end

function best = findBestRigidTransform(movingMask, fixedMask, mirrorModes, tryRot90, optimizer, metric)
% imregtform expects numeric types (not logical) in many MATLAB versions
fixed  = single(logical(fixedMask));
moving0 = single(logical(movingMask));


rotInits = 0;
if tryRot90
    rotInits = [0 90 180 270];
end

best.iou = -inf;
best.tform = affine2d(eye(3));
best.mirrorMode = "none";

Rfixed = imref2d(size(fixed));

for mi = 1:numel(mirrorModes)
    mm = string(mirrorModes{mi});
    movingM = applyMirror(moving0, mm);

    for ang = rotInits
        movingR = imrotate(movingM, ang, 'nearest', 'crop');

        try
            init = rotAboutCenter(ang, size(movingR));
            tform = imregtform(movingR, Rfixed, fixed, Rfixed, 'rigid', optimizer, metric, ...
                'InitialTransformation', init, 'PyramidLevels', 3);
        catch
            tform = imregtform(movingR, fixed, 'rigid', optimizer, metric, 'PyramidLevels', 3);
        end

        moved = imwarp(movingR, tform, 'OutputView', Rfixed, 'Interp','nearest');

        iou = IoU(moved, fixed);

        if iou > best.iou
            best.iou = iou;
            best.tform = tform;
            best.mirrorMode = mm;
        end
    end
end
end

function moved = applyMirror(BW, mode)
switch string(mode)
    case "none"
        moved = BW;
    case "fliplr"
        moved = fliplr(BW);
    case "flipud"
        moved = flipud(BW);
    otherwise
        error('Unknown mirror mode: %s', mode);
end
end

function [xm, ym] = mirrorPointVec(x, y, mode, imSize)
H = imSize(1); W = imSize(2);
mode = string(mode);
switch mode
    case "none"
        xm = x; ym = y;
    case "fliplr"
        xm = (W + 1) - x; ym = y;
    case "flipud"
        xm = x; ym = (H + 1) - y;
    otherwise
        error("Unknown mirror mode: %s", mode);
end
end

function tform = rotAboutCenter(angleDeg, imSize)
H = imSize(1); W = imSize(2);
cx = (W+1)/2; cy = (H+1)/2;

a = deg2rad(angleDeg);
R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];

T1 = [1 0 0; 0 1 0; -cx -cy 1];
T2 = [1 0 0; 0 1 0;  cx  cy 1];

M = T1 * R * T2;
tform = affine2d(M);
end

function v = IoU(A, B)
A = logical(A); B = logical(B);
u = nnz(A | B);
if u == 0, v = 0; else, v = nnz(A & B)/u; end
end

function [stimCenter_xy, targetAxisAngleDeg, distractAngleDeg] = computeTemplateGeometry(T)
% stimulus center: midpoint between closest boundary points of the two masks
yB = bwperim(T.yMask);
pB = bwperim(T.pMask);

[Dp, idx] = bwdist(pB);
[yIdx, xIdx] = find(yB);

if isempty(yIdx) || isempty(find(pB,1))
    stimCenter_xy = [NaN NaN];
    targetAxisAngleDeg = NaN;
    distractAngleDeg = NaN;
    return;
end

dvals = Dp(sub2ind(size(Dp), yIdx, xIdx));
[~, k] = min(dvals);

y1 = yIdx(k); x1 = xIdx(k);
lin = idx(sub2ind(size(idx), y1, x1));
[y2, x2] = ind2sub(size(idx), lin);

stimCenter_xy = [(x1+x2)/2, (y1+y2)/2];

% target/distractor masks
if T.targetIsYellow
    targetMask = T.yMask; distractMask = T.pMask;
else
    targetMask = T.pMask; distractMask = T.yMask;
end

rt = regionprops(targetMask, 'Orientation');
rd = regionprops(distractMask, 'Centroid');

targetAxisAngleDeg = rt(1).Orientation;

vd = [rd(1).Centroid(1)-stimCenter_xy(1), rd(1).Centroid(2)-stimCenter_xy(2)];
distractAngleDeg = wrapTo180(atan2d(vd(2), vd(1)) - targetAxisAngleDeg);
end

function [sorted, idx] = sort_nat(c)
[sorted, idx] = sort(c);
try
    expr = '\d+';
    tokens = regexp(c, expr, 'match');
    nums = cellfun(@(t) iff(isempty(t), NaN, str2double(t{end})), tokens);
    [~, idx] = sortrows([isnan(nums(:)) nums(:)], [1 2]);
    sorted = c(idx);
catch
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
