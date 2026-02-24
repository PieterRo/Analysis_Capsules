% First RF_on_stim
close all

NStim=184;
flipY = true;   % try true first; if it's mirrored vertically, set false
cx = 512; cy = 384;             % bitmap pixel corresponding to RF (0,0)

% --- RF vectors (assumes you already have x,y for the current monkey) ---
x_rf = x(goodGlobal);   % RF centers in "pixel offsets" relative to (0,0)
y_rf = y(goodGlobal);
NRF = numel(goodGlobal);

I6 = imread(StimDir+"/061.bmp");
F6 = extract_capsule_axes(I6);
refIdx    = 61;

% Convert to bitmap pixels
x_pix = cx + x_rf;
if flipY
    y_pix = cy - y_rf;
else
    y_pix = cy + y_rf;
end

N_RF_meas = NStim/2;
All_RFs_x = nan(NRF,N_RF_meas);
All_RFs_y = nan(NRF,N_RF_meas);

for StimLoop = 1:NStim/2
    stimIdx = StimLoop*2-1;
    bmpFile = fullfile(StimDir, sprintf('%03d.bmp', stimIdx));
    I = imread(bmpFile);
    [H,W,~] = size(I);
    
    % ---- Create figure + axes ----

    
    fig = figure('Color','w','WindowStyle','normal');
    set(fig,'Units','pixels');
    set(fig,'Position',[100 100 1800 1000]);   % explicit size
    tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    ax1 = subplot(1,2,1);
    imshow(I, 'Parent', ax1);
    hold(ax1, 'on');                      % IMPORTANT: after imshow, and on that axes
    axis(ax1, 'image');
    title(ax1, sprintf('Stim %03d with RF centers', stimIdx));

    % Plot centers by area (use your existing colors cV1/cV4/cIT)
    scatter(x_pix, y_pix, 18, cV1, 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

    % Mark the RF-origin / image center
    plot(cx, cy, 'wx', 'MarkerSize', 12, 'LineWidth', 5); % white X
    plot(cx, cy, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % black X on top
    legend({'V1','(0,0) at (512,384)'}, 'Location','best');

% test rotation result first 
    I1 = I;
    F1 = extract_capsule_axes(I1);

% do the alignment
    out = align_centroids_about_center(F1, F6, 'rotTolDeg', 14, 'centTolPx', 15);
    out.ok, out.mirror, out.rotDeg, out.errY, out.errP

% Example: plot V1 RFs from 001 onto reference 061
    movingIdx = stimIdx;
    tform = out.tform;   % <-- THIS is the one and only transform used for the bitmap

% --- apply EXACT same transform as bitmap ---
    [x_ref, y_ref] = transformPointsForward(tform, x_pix, y_pix);
    All_RFs_x(:,StimLoop) = x_ref;
    All_RFs_y(:,StimLoop) = y_ref;


% --- plot on the reference bitmap ---
    Iref = imread(fullfile(StimDir, sprintf('%03d.bmp', refIdx)));
    ax2 = subplot(1,2,2);
    imshow(Iref,'Parent',ax2); hold(ax2,'on'); axis(ax2,'image');
    scatter(ax2, x_ref, y_ref, 14, 'filled');
    hold off
end



fig = figure('Color','w','WindowStyle','normal');
hold on;
set(fig,'Units','pixels');
set(fig,'Position',[100 100 1800 1200]);   % explicit size
ax = axes();
imshow(Iref, 'Parent', ax);
hold(ax, 'on');                      % IMPORTANT: after imshow, and on that axes
axis(ax, 'image');
for SLoop = 1:NStim/2
    x_ref = All_RFs_x(:,SLoop);
    y_ref = All_RFs_y(:,SLoop);
    scatter(ax, x_ref, y_ref, 14, 'filled');
end
hold off


