% First RF_on_stim_V1


stimIdx = 3;                    % choose 1..328
cx = 512; cy = 384;             % bitmap pixel corresponding to RF (0,0)

% If RF y is positive upward (visual coordinates), set this true:
flipY = true;   % try true first; if it's mirrored vertically, set false

% --- Load bitmap ---
bmpFile = fullfile(StimDir, sprintf('%03d.bmp', stimIdx));
I = imread(bmpFile);
[H,W,~] = size(I);

% --- RF vectors (assumes you already have x,y for the current monkey) ---
x_rf = x(goodGlobal);   % RF centers in "pixel offsets" relative to (0,0)
y_rf = y(goodGlobal);

% Convert to bitmap pixels
x_pix = cx + x_rf;
if flipY
    y_pix = cy - y_rf;
else
    y_pix = cy + y_rf;
end

% Valid points: finite + inside image bounds
valid = isfinite(x_pix) & isfinite(y_pix) & x_pix>=1 & x_pix<=W & y_pix>=1 & y_pix<=H;

% ---- Create figure + axes ----
figure('Color','w');
ax = axes();                         % explicit axes handle
imshow(I, 'Parent', ax);
hold(ax, 'on');                      % IMPORTANT: after imshow, and on that axes
axis(ax, 'image');
title(ax, sprintf('Stim %03d with RF centers', stimIdx));


% Plot centers by area (use your existing colors cV1/cV4/cIT)
scatter(x_pix, y_pix, 18, cV1, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

% Mark the RF-origin / image center
plot(cx, cy, 'wx', 'MarkerSize', 12, 'LineWidth', 5); % white X
plot(cx, cy, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % black X on top

legend({'V1','(0,0) at (512,384)'}, 'Location','best');



% test rotation result first 
I1 = I;
I6 = imread(StimDir+"/061.bmp");

F1 = extract_capsule_axes(I1);
F6 = extract_capsule_axes(I6);

testje

figure; imshow(I1); hold on;
title('001.bmp with detected axes');

% yellow capsule
plot(F1.yellow.Centroid(1), F1.yellow.Centroid(2), 'y+', 'LineWidth',2, 'MarkerSize',12);
drawAxis(F1.yellow.Centroid, F1.yellow.Orientation, 120, 'y');

% purple capsule
plot(F1.purple.Centroid(1), F1.purple.Centroid(2), 'm+', 'LineWidth',2, 'MarkerSize',12);
drawAxis(F1.purple.Centroid, F1.purple.Orientation, 120, 'm');

plot(F1.stimCenter_xy(1), F1.stimCenter_xy(2), 'wo', 'MarkerFaceColor','w', 'MarkerSize',8);
hold off;

% next plot
figure; imshow(I6); hold on;
title('061.bmp with detected axes');

% yellow capsule
plot(F6.yellow.Centroid(1), F6.yellow.Centroid(2), 'y+', 'LineWidth',2, 'MarkerSize',12);
drawAxis(F6.yellow.Centroid, F6.yellow.Orientation, 120, 'y');

% purple capsule
plot(F6.purple.Centroid(1), F6.purple.Centroid(2), 'm+', 'LineWidth',2, 'MarkerSize',12);
drawAxis(F6.purple.Centroid, F6.purple.Orientation, 120, 'm');

plot(F6.stimCenter_xy(1), F6.stimCenter_xy(2), 'wo', 'MarkerFaceColor','w', 'MarkerSize',8);
hold off;


%plot RFs rotated

% do the alignment
out = align_centroids_about_center(F1, F6, 'rotTolDeg', 8, 'centTolPx', 15);
out.ok, out.mirror, out.rotDeg, out.errY, out.errP

% Example: plot V1 RFs from 001 onto reference 061

movingIdx = 1;
refIdx    = 61;

% choose channels (example: V1 first 512)

movingIdx = 1; refIdx = 61;

% example RF offsets (your x_pix/y_pix in RF_on_stim_V1 are offsets)
x_rf = x_pix(:);   % or whatever your variable names are in RF_on_stim_V1
y_rf = y_pix(:);

% outRF = plot_RFs_on_reference(StimDir, movingIdx, refIdx, x_rf, y_rf, 'centTolPx', 15);

% --- you already have this from your existing code ---
% out = align_centroids_about_center(Fmov, Fref, ...);
tform = out.tform;   % <-- THIS is the one and only transform used for the bitmap

% --- apply EXACT same transform as bitmap ---
[x_ref, y_ref] = transformPointsForward(tform, x_pix, y_pix);

% --- plot on the reference bitmap ---
figure; hold on;
Iref = imread(fullfile(StimDir, sprintf('%03d.bmp', refIdx)));
ax = axes(); imshow(Iref,'Parent',ax); hold(ax,'on'); axis(ax,'image');
scatter(ax, x_ref, y_ref, 14, 'filled');
hold off
