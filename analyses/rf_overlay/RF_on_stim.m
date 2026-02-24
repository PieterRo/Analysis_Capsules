% first run RFs

% RFs

stimIdx = 67;                    % choose 1..328
cx = 512; cy = 384;             % bitmap pixel corresponding to RF (0,0)

% If RF y is positive upward (visual coordinates), set this true:
flipY = true;   % try true first; if it's mirrored vertically, set false

%% --- Load bitmap ---
bmpFile = fullfile(StimDir, sprintf('%03d.bmp', stimIdx));
I = imread(bmpFile);
[H,W,~] = size(I);

%% --- RF vectors (assumes you already have x,y for the current monkey) ---
x_rf = x;   % RF centers in "pixel offsets" relative to (0,0)
y_rf = y;

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
idx = V1(valid(V1));
scatter(x_pix(idx), y_pix(idx), 18, cV1, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

idx = V4(valid(V4));
scatter(x_pix(idx), y_pix(idx), 22, cV4, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

idx = IT(valid(IT));
scatter(x_pix(idx), y_pix(idx), 22, cIT, 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.15);

% Mark the RF-origin / image center
plot(cx, cy, 'wx', 'MarkerSize', 12, 'LineWidth', 5); % white X
plot(cx, cy, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % black X on top

legend({'V1','V4','IT','(0,0) at (512,384)'}, 'Location','best');



