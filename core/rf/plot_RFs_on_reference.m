function outRF = plot_RFs_on_reference(StimDir, movingIdx, refIdx, x_rf, y_rf, varargin)
% Transform RF centers with the same bitmap alignment and plot on ref bitmap.

p = inputParser;
p.addParameter('flipY', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('rotTolDeg', 8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('centTolPx', 15, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('markerSize', 18, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('alpha', 0.35, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.parse(varargin{:});
flipY     = p.Results.flipY;
rotTolDeg = p.Results.rotTolDeg;
centTolPx = p.Results.centTolPx;
ms        = p.Results.markerSize;
a         = p.Results.alpha;

% --- Load bitmaps ---
movFile = fullfile(StimDir, sprintf('%03d.bmp', movingIdx));
refFile = fullfile(StimDir, sprintf('%03d.bmp', refIdx));
Imov = imread(movFile);
Iref = imread(refFile);

% --- Feature extraction (must exist on your path) ---
Fmov = extract_capsule_axes(Imov);
Fref = extract_capsule_axes(Iref);

% IMPORTANT: do NOT override stimCenter_xy here.
% Use the same center definition as in your working AlignBitmapsRFs.m.

% --- Compute alignment moving -> reference (robust retries) ---
[outAlign, Fmov_used] = robust_align(Fmov, Fref, rotTolDeg, centTolPx);

T = outAlign.tform;

% --- Convert RF offsets -> absolute pixels in MOVING bitmap coordinates ---
% Use the moving stimCenter that was used for alignment
Cmov = Fmov_used.stimCenter_xy(:)';   % [cx cy]

x_pix = Cmov(1) + x_rf(:);
if flipY
    y_pix = Cmov(2) - y_rf(:);       % visual +y up -> image y decreases
else
    y_pix = Cmov(2) + y_rf(:);
end

% --- Transform RF pixels into REFERENCE bitmap coordinates ---
[x_ref, y_ref] = transformPointsForward(T, x_pix, y_pix);

% --- Plot ONLY on reference bitmap ---
figure('Color','w');
ax = axes();
imshow(Iref, 'Parent', ax); hold(ax,'on'); axis(ax,'image');
title(ax, sprintf('RF centers onto %03d.bmp (from %03d.bmp)', refIdx, movingIdx));

scatter(ax, x_ref, y_ref, ms, 'filled', ...
    'MarkerFaceAlpha', a, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

% Mark reference stim center (detected)
Cref = Fref.stimCenter_xy(:)';
plot(ax, Cref(1), Cref(2), 'wx', 'MarkerSize', 12, 'LineWidth', 5);
plot(ax, Cref(1), Cref(2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);

% --- Pack outputs ---
outRF = struct();
outRF.align = outAlign;
outRF.tform = T;
outRF.x_pix_moving = x_pix;
outRF.y_pix_moving = y_pix;
outRF.x_pix_ref    = x_ref;
outRF.y_pix_ref    = y_ref;
outRF.Cmov = Cmov;
outRF.Cref = Cref;
end

% ===================== helper =====================

function [outAlign, Fmov_used] = robust_align(Fmov, Fref, rotTolDeg, centTolPx)
% Try normal; if it fails due to rot disagreement, try swapping Y/P labels.

Fmov_used = Fmov;

try
    outAlign = align_centroids_about_center(Fmov_used, Fref, ...
        'rotTolDeg', rotTolDeg, 'centTolPx', centTolPx);
    return
catch ME
    % fall through
end

% Retry 1: swap yellow/purple in moving (common failure mode if color ID flips)
Fswap = Fmov_used;
tmp = Fswap.yellow;
Fswap.yellow = Fswap.purple;
Fswap.purple = tmp;

try
    outAlign = align_centroids_about_center(Fswap, Fref, ...
        'rotTolDeg', rotTolDeg, 'centTolPx', centTolPx);
    Fmov_used = Fswap;
    return
catch ME2
    % If it still fails, rethrow original info but with extra context
    error("Alignment failed even after swapping Y/P in moving. Original error: %s | Swap retry error: %s", ...
        ME.message, ME2.message);
end
end
