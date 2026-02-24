function h = plot_stim_with_RFs(ALLCOORDS, RTAB384, stimNum, x_rf, y_rf)
% plot_stim_with_RFs
%
% Renders a stimulus using render_stim_from_ALLCOORDS and overlays RF centers.
%
% INPUTS
%   ALLCOORDS : struct with stimulus coordinates
%   RTAB384   : stimulus table
%   stimNum   : stimulus index (1..384)
%   x_rf      : vector of RF x-coordinates (centered, +right)
%   y_rf      : vector of RF y-coordinates (centered, +up)
%
% RF → image coordinates:
%   x_img = x_rf + 512
%   y_img = 384 - y_rf
%
% OUTPUT
%   h : handle struct (figure, axes, scatter)

% ------------------------------------------------------------
% Render stimulus (your existing routine)
% ------------------------------------------------------------
img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimNum);

[H, W, ~] = size(img);   % should be 768 x 1024

% ------------------------------------------------------------
% Convert RF coordinates to image pixel coordinates
% ------------------------------------------------------------
x_img = x_rf + W/2;      % = x + 512
y_img = H/2 - y_rf;      % = 384 - y

% ------------------------------------------------------------
% Plot
% ------------------------------------------------------------
h.fig = figure;
h.ax  = axes('Position',[0 0 1 1]);

imshow(img, 'InitialMagnification', 'fit');
axis image off
hold on

h.rf = scatter(x_img, y_img, ...
    20, ...              % marker size
    'w', ...             % face color
    'filled', ...
    'MarkerEdgeColor','k');

hold off

end
