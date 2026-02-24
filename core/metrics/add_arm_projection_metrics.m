function T = add_arm_projection_metrics(T, ALLCOORDS, stimNum, varargin)
% ADD_ARM_PROJECTION_METRICS
% Adds arm-projection metrics for RFs on target/distractor objects.
%
% Requires T contains:
%   T.x_px, T.y_px (pixel coords)
%   T.assignment   (string/categorical: "target", "distractor", ...)
%
% Adds columns:
%   T.along        : signed distance from s to projection point along arm axis (px)
%                   (NOT clamped; endcaps may yield <0 or > arm length)
%   T.proj_x/y     : projection point on infinite arm axis (px)
%   T.perp_signed  : signed perpendicular distance (px), positive toward other arm
%
% Options (Name-Value):
%   'ImageSize' : [W H] (default [1024 768])

p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.parse(varargin{:});
W = p.Results.ImageSize(1);
H = p.Results.ImageSize(2);

% ---- get geometry for this stimulus ----
fieldName = sprintf('stim_%d', stimNum);
assert(isfield(ALLCOORDS, fieldName), 'ALLCOORDS missing field %s', fieldName);

s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

toPx = @(q) [q(1) + W/2, H/2 - q(2)];

s_px     = toPx(s);
tFig_px  = toPx(tFig);
tBack_px = toPx(tBack);

% Unit directions for each arm axis
uT = unit_vec(tFig_px - s_px);   % target axis direction
uD = unit_vec(tBack_px - s_px);  % distractor axis direction

% Normals pointing TOWARD the other arm
% Target: + toward distractor
nT = toward_other_normal(s_px, tFig_px, tBack_px);
% Distractor: + toward target
nD = toward_other_normal(s_px, tBack_px, tFig_px);

% ---- allocate outputs ----
N = height(T);
along       = nan(N,1);
proj_x      = nan(N,1);
proj_y      = nan(N,1);
perp_signed = nan(N,1);

% Make sure assignment is string for comparisons
assign = T.assignment;
if iscategorical(assign)
    assign = string(assign);
elseif ischar(assign)
    assign = string(assign);
end

isTarget = (assign == "target");
isDistr  = (assign == "distractor");

pRF = [T.x_px, T.y_px];

% ---- target RFs ----
if any(isTarget)
    pt = pRF(isTarget,:);
    a  = pt - s_px;                  % Nx2
    t  = a * uT(:);                  % Nx1 along-axis coordinate (unclamped)
    proj = s_px + t .* uT;           % Nx2 projection point on infinite axis
    d  = sum((pt - proj) .* nT, 2);  % signed perp distance (+ toward distractor)

    along(isTarget)       = t;
    proj_x(isTarget)      = proj(:,1);
    proj_y(isTarget)      = proj(:,2);
    perp_signed(isTarget) = d;
end

% ---- distractor RFs ----
if any(isDistr)
    pd = pRF(isDistr,:);
    a  = pd - s_px;
    t  = a * uD(:);
    proj = s_px + t .* uD;
    d  = sum((pd - proj) .* nD, 2);  % signed perp distance (+ toward target)

    along(isDistr)       = t;
    proj_x(isDistr)      = proj(:,1);
    proj_y(isDistr)      = proj(:,2);
    perp_signed(isDistr) = d;
end

% Attach
T.along       = along;
T.proj_x      = proj_x;
T.proj_y      = proj_y;
T.perp_signed = perp_signed;

end

% ---------- helpers ----------
function u = unit_vec(v)
L = hypot(v(1), v(2));
if L < 1e-12
    u = [0 0];
else
    u = v / L;
end
end

function n = toward_other_normal(s_px, t_arm_px, t_other_px)
% Returns a unit normal to the arm axis pointing toward the other arm.
v = t_arm_px - s_px;
Lv = hypot(v(1), v(2));
if Lv < 1e-12
    n = [0 0];
    return
end
u = v / Lv;

w = t_other_px - s_px;
w_perp = w - dot(w, u) * u;   % component perpendicular to arm axis
Lw = hypot(w_perp(1), w_perp(2));

if Lw < 1e-12
    % arms nearly collinear -> no well-defined "toward other" side
    n = [0 0];
else
    n = w_perp / Lw;
end
end
