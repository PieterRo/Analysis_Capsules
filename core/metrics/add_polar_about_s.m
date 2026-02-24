function T = add_polar_about_s(T, ALLCOORDS, stimNum, varargin)
% ADD_POLAR_ABOUT_S
% Adds polar coordinates (r, theta) of RFs relative to stimulus junction s.
%
% Convention:
%   - Origin at s (in pixel coordinates).
%   - 0 degrees points along the INTERNAL bisector between the two arm axes:
%         b = unit(u_target + u_distractor)
%   - Angles increase toward the TARGET arm (u_target).
%   - theta in [0, 360).
%   - r is Euclidean distance from s in pixels.
%
% Requires T has:
%   T.x_px, T.y_px  (pixel coords)
%
% Options (Name-Value):
%   'ImageSize'         : [W H] (default [1024 768])
%   'OnlyBackground'    : true/false (default false)
%   'BackgroundLabel'   : label used in T.assignment for background (default "background")
%
% Adds columns:
%   T.r_s
%   T.theta_s_deg

p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('OnlyBackground', false, @(b) islogical(b) && isscalar(b));
p.addParameter('BackgroundLabel', "background", @(s) ischar(s) || isstring(s));
p.parse(varargin{:});

W = p.Results.ImageSize(1);
H = p.Results.ImageSize(2);
onlyBg = p.Results.OnlyBackground;
bgLabel = string(p.Results.BackgroundLabel);

% --- stimulus geometry in pixel coords ---
fieldName = sprintf('stim_%d', stimNum);
assert(isfield(ALLCOORDS, fieldName), 'ALLCOORDS missing field %s', fieldName);

s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

toPx = @(q) [q(1) + W/2, H/2 - q(2)];

s_px     = toPx(s);
tFig_px  = toPx(tFig);
tBack_px = toPx(tBack);

uT = unit_vec(tFig_px  - s_px);  % target direction
uD = unit_vec(tBack_px - s_px);  % distractor direction

% Internal bisector (between arms)
b_raw = uT + uD;
if hypot(b_raw(1), b_raw(2)) < 1e-10
    % Arms nearly opposite -> internal bisector undefined
    % Fallback: choose b perpendicular to uT (arbitrary but stable), and warn via NaNs
    b = [NaN NaN];
else
    b = unit_vec(b_raw);
end

% Determine sign so that TARGET is at positive angle from bisector
if any(isnan(b))
    thetaT = NaN;
else
    thetaT = atan2(cross2(b, uT), dot(b, uT)); % radians in [-pi, pi]
end
flipSign = ~isnan(thetaT) && (thetaT < 0);

% --- compute r/theta for RFs ---
N = height(T);
r_s = nan(N,1);
theta_deg = nan(N,1);

% Optional restriction to background
idx = true(N,1);
if onlyBg
    if ismember("assignment", string(T.Properties.VariableNames))
        assign = T.assignment;
        if iscategorical(assign), assign = string(assign); end
        idx = (assign == bgLabel);
    else
        error('OnlyBackground=true requires T.assignment column.');
    end
end

pRF = [T.x_px, T.y_px];
v = pRF(idx,:) - s_px;               % Nx2 vectors from s to RF
r = sqrt(sum(v.^2, 2));
vhat = v ./ max(r, eps);             % avoid div-by-0

if any(isnan(b))
    % If bisector undefined, keep angles NaN
    theta = nan(sum(idx),1);
else
    theta = atan2(cross2_row(b, vhat), dot_row(b, vhat)); % signed radians [-pi, pi]
    if flipSign
        theta = -theta;
    end
end

r_s(idx) = r;
theta_deg(idx) = mod(rad2deg(theta), 360);

T.r_s = r_s;
T.theta_s_deg = theta_deg;

end

% -------- helpers --------
function u = unit_vec(v)
L = hypot(v(1), v(2));
if L < 1e-12
    u = [0 0];
else
    u = v / L;
end
end

function z = cross2(a, b)
% scalar 2D cross product a x b
z = a(1)*b(2) - a(2)*b(1);
end

function z = cross2_row(a, B)
% cross2(a, B(i,:)) for each row i
z = a(1)*B(:,2) - a(2)*B(:,1);
end

function d = dot_row(a, B)
% dot(a, B(i,:)) for each row i
d = a(1)*B(:,1) + a(2)*B(:,2);
end
