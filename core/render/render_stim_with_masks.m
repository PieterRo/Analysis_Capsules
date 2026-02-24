function [img, masks, meta] = render_stim_with_masks(ALLCOORDS, RTAB384, stimNum, varargin)
% RENDER_STIM_WITH_MASKS
% Like render_stim_from_ALLCOORDS, but also returns masks for green/purple/background.
%
% Outputs:
%   img   : HxWx3 double in [0..1]
%   masks : struct with logical fields:
%           .green, .purple, .background, .figArm, .backArm
%   meta  : includes color assignment etc.

% -----------------------
% Parse inputs (same defaults as your renderer, but default DrawDots=false)
% -----------------------
p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('Background', 0.5, @(v) isnumeric(v) && isscalar(v));
p.addParameter('DrawDots', false, @(v) islogical(v) && isscalar(v)); % default FALSE for analysis
p.addParameter('DotDiameterPx', 16, @(v) isnumeric(v) && isscalar(v));
defaultColors.green  = [0.20 0.75 0.35];
defaultColors.purple = [0.70 0.60 0.85];
defaultColors.yellow = [1.00 0.90 0.00];
p.addParameter('Colors', defaultColors, @(s) isstruct(s));
p.parse(varargin{:});
opt = p.Results;

W = opt.ImageSize(1);
H = opt.ImageSize(2);

assert(stimNum >= 1 && stimNum <= size(RTAB384,1), 'stimNum out of range.');
assert(size(RTAB384,2) >= 8, 'RTAB384 must have at least 8 columns.');

fieldName = sprintf('stim_%d', stimNum);
assert(isfield(ALLCOORDS, fieldName), 'ALLCOORDS missing field %s', fieldName);

s     = double(ALLCOORDS.(fieldName).s(:)');
tFig  = double(ALLCOORDS.(fieldName).t_fig(:)');
tBack = double(ALLCOORDS.(fieldName).t_back(:)');

widthPx = double(RTAB384(stimNum, 7));
colIdx  = double(RTAB384(stimNum, 8));

% Color swap rule (same as your file)
if mod(colIdx,2)==1
    colFig  = opt.Colors.green;
    colBack = opt.Colors.purple;
    figColorName  = "green";
    backColorName = "purple";
else
    colFig  = opt.Colors.purple;
    colBack = opt.Colors.green;
    figColorName  = "purple";
    backColorName = "green";
end

% Initialize image + masks
img = repmat(opt.Background, [H W 3]);
maskFig  = false(H,W);
maskBack = false(H,W);

toImg = @(p) [p(1) + W/2, H/2 - p(2)];

% Build masks from geometry (no reading back pixels)
maskBack = capsule_mask(toImg(s), toImg(tBack), widthPx, H, W);
maskFig  = capsule_mask(toImg(s), toImg(tFig),  widthPx, H, W);

% Paint image (optional but useful for QC/plotting)
img = paint_mask(img, maskBack, colBack);
img = paint_mask(img, maskFig,  colFig);

% Optional endpoint dots (usually OFF for analysis)
if opt.DrawDots
    dotR = opt.DotDiameterPx/2;
    img = draw_filled_circle(img, toImg(s),     dotR, opt.Colors.yellow);
    img = draw_filled_circle(img, toImg(tFig),  dotR, opt.Colors.yellow);
    img = draw_filled_circle(img, toImg(tBack), dotR, opt.Colors.yellow);
end

label = zeros(H,W,'uint8');  % 0=background, 1=green, 2=purple

% Paint BACK arm first
if backColorName=="green"
    label(maskBack) = 1;
else
    label(maskBack) = 2;
end

% Paint FIG arm last (overwrites in overlap)
if figColorName=="green"
    label(maskFig) = 1;
else
    label(maskFig) = 2;
end

maskGreen  = (label == 1);
maskPurple = (label == 2);
maskBg     = (label == 0);

masks = struct();
masks.figArm     = maskFig;
masks.backArm    = maskBack;
masks.green      = maskGreen;
masks.purple     = maskPurple;
masks.background = maskBg;

meta = struct();
meta.stimNum   = stimNum;
meta.widthPx   = widthPx;
meta.colorIdx  = colIdx;
meta.s         = s;
meta.t_fig     = tFig;
meta.t_back    = tBack;
meta.col_fig   = colFig;
meta.col_back  = colBack;
meta.figColor  = figColorName;   % "green" or "purple"
meta.backColor = backColorName;  % "green" or "purple"
meta.imageSize = [W H];
meta.background= opt.Background;

end

% ---- helpers ----

function mask = capsule_mask(p1, p2, diameterPx, H, W)
r = diameterPx/2;

v = p2 - p1;
L = hypot(v(1), v(2));
[xx,yy] = meshgrid(1:W,1:H);

if L < 1e-6
    mask = ((xx-p1(1)).^2 + (yy-p1(2)).^2 <= r^2);
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

mask = mask | ((xx-p1(1)).^2 + (yy-p1(2)).^2 <= r^2);
mask = mask | ((xx-p2(1)).^2 + (yy-p2(2)).^2 <= r^2);
end

function img = paint_mask(img, mask, rgb)
for c = 1:3
    chan = img(:,:,c);
    chan(mask) = rgb(c);
    img(:,:,c) = chan;
end
end

function img = draw_filled_circle(img, pc, r, rgb)
[H,W,~] = size(img);
[xx,yy] = meshgrid(1:W,1:H);
mask = ((xx-pc(1)).^2 + (yy-pc(2)).^2 <= r^2);
img = paint_mask(img, mask, rgb);
end
