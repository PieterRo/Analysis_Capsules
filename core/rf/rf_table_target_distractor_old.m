function T = rf_table_target_distractor(ALLCOORDS, RTAB384, stimNum, x_rf, y_rf, varargin)
% RF_TABLE_TARGET_DISTRACTOR
% For each RF center: index, x/y in pixels, color at center (yellowArm/purple/gray),
% and whether it is on target arm (t_fig), distractor arm (t_back), or background.
%
% Uses analytic masks from render_stim_with_masks2 (no pixel color read-back).
%
% Options:
%   'ClipToImage' : true (default). If false, out-of-bounds stays NaN/out_of_bounds.

p = inputParser;
p.addParameter('ClipToImage', true, @(b) islogical(b) && isscalar(b));
p.parse(varargin{:});
opt = p.Results;

% Render + masks (dots OFF by default inside renderer)
[~, masks, meta] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimNum);

W = meta.imageSize(1);
H = meta.imageSize(2);

% RF -> pixel coords
x_px = x_rf(:) + W/2;   % e.g. +512
y_px = H/2 - y_rf(:);   % e.g. 384 - y
N = numel(x_px);

% Bounds and sampling indices
inBounds = x_px>=1 & x_px<=W & y_px>=1 & y_px<=H;
x_s = round(x_px);
y_s = round(y_px);

if opt.ClipToImage
    x_s(~inBounds) = max(1, min(W, x_s(~inBounds)));
    y_s(~inBounds) = max(1, min(H, y_s(~inBounds)));
end

idx = sub2ind([H W], y_s, x_s);

% (3) Color at RF center
isYellow = masks.yellowArm(idx);
isPurple = masks.purple(idx);

center_color = strings(N,1);
center_color(:) = "gray";
center_color(isYellow) = "yellowArm";
center_color(isPurple) = "purple";

if ~opt.ClipToImage
    center_color(~inBounds) = "out_of_bounds";
end

% (4) Target/distractor/background from arm masks
onTarget     = masks.figArm(idx);   % t_fig arm
onDistractor = masks.backArm(idx);  % t_back arm

assignment = strings(N,1);
assignment(:) = "background";
assignment(onDistractor) = "distractor";
assignment(onTarget)     = "target";

if ~opt.ClipToImage
    assignment(~inBounds) = "out_of_bounds";
end

% Helpful flag: overlaps (junction region might be in both if they overlap)
overlap = onTarget & onDistractor;

RF = (1:N).';
T = table(RF, x_px, y_px, center_color, assignment, overlap, inBounds, ...
    'VariableNames', {'RF','x_px','y_px','center_color','assignment','overlap','inBounds'});

end
