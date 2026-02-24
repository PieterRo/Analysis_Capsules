function [img, meta] = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimNum, varargin)
%RENDER_STIM_FROM_ALLCOORDS Reconstruct stimulus bitmap from ALLCOORDS + RTAB384.
%
% Color rule (as requested):
%   colIdx == 1: colFig = purple,  colDist = yellowArm
%   colIdx == 2: colFig = yellowArm, colDist = purple
%
% Dots (3 circles) use a separate, brighter yellowDot.

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('Background', [128 128 128]/255, @(v) isnumeric(v) && (isscalar(v) || (isvector(v) && numel(v)==3)));
p.addParameter('DrawDots', true, @(v) islogical(v) && isscalar(v));
p.addParameter('DotDiameterPx', 16, @(v) isnumeric(v) && isscalar(v));

% Exact palette to reproduce your screenshot
defaultColors.purple    = [153 154 178]/255;  % lavender arm
defaultColors.yellowArm = [153 157 122]/255;  % pastel/olive arm
defaultColors.yellowDot = [255 253  56]/255;  % bright dot yellow

p.addParameter('Colors', defaultColors, @(s) isstruct(s));
p.parse(varargin{:});
opt = p.Results;

W = opt.ImageSize(1);
H = opt.ImageSize(2);

% Normalize background input to scalar gray or RGB
bg = opt.Background;
if isscalar(bg)
    bgRGB = [bg bg bg];
else
    bgRGB = bg(:)';
end

% Validate stimNum / RTAB
assert(stimNum >= 1 && stimNum <= size(RTAB384,1), 'stimNum out of range.');
assert(size(RTAB384,2) >= 8, 'RTAB384 must have at least 8 columns.');

% -----------------------
% Look up coords + params
% -----------------------
fieldName = sprintf('stim_%d', stimNum);
assert(isfield(ALLCOORDS, fieldName), 'ALLCOORDS missing field %s', fieldName);

s     = double(ALLCOORDS.(fieldName).s(:)');
tFig  = double(ALLCOORDS.(fieldName).t_fig(:)');
tBack = double(ALLCOORDS.(fieldName).t_back(:)');

widthPx = double(RTAB384(stimNum, 7));
colIdx  = double(RTAB384(stimNum, 8));

% -----------------------
% Color mapping (requested)
% -----------------------
switch colIdx
    case 1
        colFig  = opt.Colors.purple;
        colDist = opt.Colors.yellowArm;
    case 2
        colFig  = opt.Colors.yellowArm;
        colDist = opt.Colors.purple;
    otherwise
        error('Unexpected colIdx=%g (expected 1 or 2).', colIdx);
end

% -----------------------
% Initialize image
% -----------------------
img = zeros(H, W, 3);
img(:,:,1) = bgRGB(1);
img(:,:,2) = bgRGB(2);
img(:,:,3) = bgRGB(3);

% Centered coords -> image coords
toImg = @(p) [p(1) + W/2, H/2 - p(2)];

% Draw arms
img = draw_capsule(img, toImg(s), toImg(tBack), widthPx, colDist);
img = draw_capsule(img, toImg(s), toImg(tFig),  widthPx, colFig);

% Dots (bright yellow)
if opt.DrawDots
    dotR = opt.DotDiameterPx/2;
    yDot = opt.Colors.yellowDot;
    img = draw_filled_circle(img, toImg(s),     dotR, yDot);
    img = draw_filled_circle(img, toImg(tFig),  dotR, yDot);
    img = draw_filled_circle(img, toImg(tBack), dotR, yDot);
end

% Meta
meta = struct();
meta.stimNum   = stimNum;
meta.widthPx   = widthPx;
meta.colorIdx  = colIdx;
meta.s         = s;
meta.t_fig     = tFig;
meta.t_back    = tBack;
meta.col_fig   = colFig;
meta.col_dist  = colDist;
meta.colors    = opt.Colors;
meta.imageSize = [W H];
meta.background= bgRGB;

end

% =====================================================================
function img = draw_capsule(img, p1, p2, diameterPx, rgb)
[H,W,~] = size(img);
r = diameterPx/2;

v = p2 - p1;
L = hypot(v(1), v(2));
if L < 1e-6
    img = draw_filled_circle(img, p1, r, rgb);
    return
end
u = v / L;
n = [-u(2) u(1)];

A = p1 + n*r;
B = p1 - n*r;
C = p2 - n*r;
D = p2 + n*r;

polyX = [A(1) B(1) C(1) D(1)];
polyY = [A(2) B(2) C(2) D(2)];
mask = poly2mask(polyX, polyY, H, W);

[xx,yy] = meshgrid(1:W,1:H);
mask = mask | ((xx-p1(1)).^2 + (yy-p1(2)).^2 <= r^2);
mask = mask | ((xx-p2(1)).^2 + (yy-p2(2)).^2 <= r^2);

for c = 1:3
    chan = img(:,:,c);
    chan(mask) = rgb(c);
    img(:,:,c) = chan;
end
end

% =====================================================================
function img = draw_filled_circle(img, pc, r, rgb)
[H,W,~] = size(img);
[xx,yy] = meshgrid(1:W,1:H);
mask = ((xx-pc(1)).^2 + (yy-pc(2)).^2 <= r^2);

for c = 1:3
    chan = img(:,:,c);
    chan(mask) = rgb(c);
    img(:,:,c) = chan;
end
end
