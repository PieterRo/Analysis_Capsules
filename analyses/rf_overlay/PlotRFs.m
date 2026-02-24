% PlotRFs
% Plots RF center coordinates (x,y) on a 1024x768 pixel screen.

% --- USER INPUTS ---
W = 1024;                 % screen width in pixels
H = 768;                  % screen height in pixels
RFrange = 1:512;            % V1 RFs 
bgGray = 0.5;             % background gray level (0..1)
dotColor = [1 1 1];       % red
fpColor = [0 0 0];
dotSize = 8;              % marker size

% --- CHECKS ---
assert(exist('x','var')==1 && exist('y','var')==1, 'x and y must exist in the workspace.');
x2 = x(RFrange)+512; y2 = 381 - y(RFrange);
assert(numel(x2)==numel(y2), 'x and y must have the same length.');

% Optional: warn if points are outside screen bounds (image coordinates)
out = (x2 < 1 | x2 > W | y2 < 1 | y2 > H);
if any(out)
    fprintf('Warning: %d/%d points are outside [1..%d]x[1..%d].\n', sum(out), numel(x), W, H);
end

% --- PLOT ---
figure('Color','w');
ax = axes('Position',[0 0 1 1]);  % fill the figure

% Draw a uniform gray background as an image (avoids colormap pitfalls)
bg = bgGray * ones(H, W, 'double');
imshow(bg, 'Parent', ax, 'InitialMagnification', 'fit');

hold(ax, 'on');

% Make axes match pixel coordinates: (1,1) top-left, y increases downward
axis(ax, 'ij');                  % y-axis down (image coordinates)
axis(ax, [1 W 1 H]);             % exact pixel extents
axis(ax, 'image');               % equal aspect
axis(ax, 'off');

% Plot points
plot(ax, x2(1:512), y2(1:512), '.', 'Color', dotColor, 'MarkerSize', dotSize);

% plot center
plot(ax, 512, 384, '.', 'Color', fpColor, 'MarkerSize', 3 * dotSize);

title(ax, sprintf('RF centers (N=%d) on %dx%d screen', numel(x2), W, H));

% Keep title from shrinking the axes (optional; comment out if you prefer)
set(ax,'Position',[0 0 1 1]);
