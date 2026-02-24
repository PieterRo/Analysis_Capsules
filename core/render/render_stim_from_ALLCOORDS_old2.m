function [img, meta] = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimNum, varargin)
%RENDER_STIM_FROM_ALLCOORDS Reconstruct stimulus bitmap from ALLCOORDS + RTAB384.
%
%  [img, meta] = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimNum)
%  takes:
%     ALLCOORDS : struct with fields stim_1 ... stim_384
%                Each entry must contain:
%                   .s      : [x y] junction point (centered pixel coords)
%                   .t_fig  : [x y] endpoint of "figure/same" arm (centered px)
%                   .t_back : [x y] endpoint of "background/diff" arm (centered px)
%
%     RTAB384   : 384x8 numeric array, one row per stimulus, columns:
%                   1  stimulus id (should match row index)
%                   2  array/config index
%                   3  number of growth cones
%                   4  fg / class label (not needed for drawing)
%                   5  subtype label (not needed for drawing)
%                   6  angle (deg; usually already baked into coords)
%                   7  curve/line width (px)   <-- USED
%                   8  color assignment index  <-- USED (expects 1 or 2)
%
%     stimNum   : integer in [1..384]
%
%  Returns:
%     img  : HxWx3 double image in [0..1]
%     meta : struct with stimulus parameters used for rendering
%
%  This renderer draws two round-capped "capsule" segments:
%     s -> t_fig   and   s -> t_back
%  with diameter = RTAB384(stimNum,7).
%
%  Options (name-value pairs):
%     'ImageSize'     : [W H] pixels (default [1024 768])
%     'Background'    : gray level in [0..1] (default 0.5)
%     'DrawDots'      : true/false (default true) adds dots at endpoints
%     'DotDiameterPx' : diameter of endpoint dots (default 16)
%     'Colors'        : struct with fields:
%                       .green, .purple, .yellowDot, .yellowArm
%
%  Notes:
%   - Coordinates are assumed CENTERED pixels: (0,0) is screen center.
%   - If your ALLCOORDS uses another convention, adjust the toImg() mapping.
%
% Uses poly2mask (Image Processing Toolbox). If you don't have poly2mask,
% tell me and I'll provide a pure-MATLAB fallback.

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('Background', 0.5, @(v) isnumeric(v) && isscalar(v));
p.addParameter('DrawDots', true, @(v) islogical(v) && isscalar(v));
p.addParameter('DotDiameterPx', 16, @(v) isnumeric(v) && isscalar(v));

% Default palette
defaultColors.green  = [0.30 0.75 0.45];
defaultColors.purple = [0.70 0.60 0.85];

% Dots: saturated yellow
defaultColors.yellowDot = [1.00 0.90 0.00];

% Arms: lighter yellow (tint towards white), guaranteed lighter than yellowDot
armTint = 0.35; % 0 = identical to yellowDot, 1 = white
defaultColors.yellowArm = defaultColors.yellowDot + armTint*(1 - defaultColors.yellowDot);

p.addParameter('Colors', defaultColors, @(s) isstruct(s));
p.parse(varargin{:});
opt = p.Results;

W = opt.ImageSize(1);
H = opt.ImageSize(2);

% Validate stimNum
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

widthPx = double(RTAB384(stimNum, 7));   % curve width used by cgpenwid()
colIdx  = double(RTAB384(stimNum, 8));   % color assignment index (expects 1 or 2)

% -----------------------
% Define color mapping (your rule)
% -----------------------
% colIdx == 1: fig = purple, dist = light yellow
% colIdx == 2: fig = light yellow, dist = purple
switch colIdx
    case 1
        colFig  = opt.Colors.purple;
        colDist = opt.Colors.yellowArm;   % lighter than dot yellow
    case 2
        colFig  = opt.Colors.yellowArm;   % lighter than dot yellow
        colDist = opt.Colors.purple;
    otherwise
        error('Unexpected colIdx=%g (expected 1 or 2).', colIdx);
end

% -----------------------
% Initialize image
% -----------------------
img = repmat(opt.Background, [H W 3]);

% Helper: centered coords -> image pixel coords (float)
% Centered: x right, y up. Image: x right, y down.
toImg = @(p) [p(1) + W/2, H/2 - p(2)];

% Draw the two capsule strokes (round-capped thick segments)
% Distal/back arm first, then fig arm (fig can overwrite at the junction)
img = draw_capsule(img, toImg(s), toImg(tBack), widthPx, colDist);
img = draw_capsule(img, toImg(s), toImg(tFig),  widthPx, colFig);

% Optional endpoint dots (use saturated yellowDot)
if opt.DrawDots
    dotR = opt.DotDiameterPx/2;
    img = draw_filled_circle(img, toImg(s),     dotR, opt.Colors.yellowDot);
    img = draw_filled_circle(img, toImg(tFig),  dotR, opt.Colors.yellowDot);
    img = draw_filled_circle(img, toImg(tBack), dotR, opt.Colors.yellowDot);
end

% Pack metadata for your records
meta = struct();
meta.stimNum     = stimNum;
meta.widthPx     = widthPx;
meta.colorIdx    = colIdx;
meta.s           = s;
meta.t_fig       = tFig;
meta.t_back      = tBack;
meta.col_fig     = colFig;
meta.col_dist    = colDist;
meta.yellowDot   = opt.Colors.yellowDot;
meta.yellowArm   = opt.Colors.yellowArm;
meta.imageSize   = [W H];
meta.background  = opt.Background;

end

% =====================================================================
% Local helper: draw capsule (thick segment with round caps)
% =====================================================================
function img = draw_capsule(img, p1, p2, diameterPx, rgb)
% p1,p2 are IMAGE coords [x y] in pixels (float).
[H,W,~] = size(img);
r = diameterPx/2;  % radius

v = p2 - p1;
L = hypot(v(1), v(2));
if L < 1e-6
    img = draw_filled_circle(img, p1, r, rgb);
    return
end
u = v / L;
n = [-u(2) u(1)];

% Rectangle corners
A = p1 + n*r;
B = p1 - n*r;
C = p2 - n*r;
D = p2 + n*r;

% Polygon mask for the thick rectangle
polyX = [A(1) B(1) C(1) D(1)];
polyY = [A(2) B(2) C(2) D(2)];
mask = poly2mask(polyX, polyY, H, W);

% Add circular endcaps
[xx,yy] = meshgrid(1:W,1:H);
mask = mask | ((xx-p1(1)).^2 + (yy-p1(2)).^2 <= r^2);
mask = mask | ((xx-p2(1)).^2 + (yy-p2(2)).^2 <= r^2);

% Paint
for c = 1:3
    chan = img(:,:,c);
    chan(mask) = rgb(c);
    img(:,:,c) = chan;
end
end

% =====================================================================
% Local helper: filled circle
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

