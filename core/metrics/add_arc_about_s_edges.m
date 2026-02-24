function T = add_arc_about_s_edges(T, ALLCOORDS, RTAB384, stimNum, varargin)
% ADD_ARC_ABOUT_S_EDGES
% Edge-based angular normalization for background space around s.
% Normalizes along the circle BETWEEN OBJECT EDGES (not centerlines).
%
% Adds:
%   T.r_s
%   T.arc_isInner_edge
%   T.arc_frac_edge

p = inputParser;
p.addParameter('ImageSize', [1024 768], @(v) isnumeric(v) && numel(v)==2);
p.parse(varargin{:});
W = p.Results.ImageSize(1);
H = p.Results.ImageSize(2);

% width in px (GC)
w = double(RTAB384(stimNum, 7));
if ~(isfinite(w) && w > 0)
    error('Invalid widthPx for stim %d', stimNum);
end
rad = w/2;

% geometry
fieldName = sprintf('stim_%d', stimNum);
s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

toPx = @(q) [q(1) + W/2, H/2 - q(2)];
s_px  = toPx(s);
tT_px = toPx(tFig);   % target axis
tD_px = toPx(tBack);  % distractor axis

uT = unit_vec(tT_px - s_px);
uD = unit_vec(tD_px - s_px);

% signed inner angle from distractor axis to target axis
alpha = signed_angle(uD, uT);                 % (-pi, pi]
if abs(alpha) < 1e-10
    error('Arms nearly collinear; edge-arc undefined for stim %d', stimNum);
end
sgn = sign(alpha); if sgn==0, sgn = 1; end

% long way angle (from uD to uT)
alpha_long = alpha - sgn*2*pi;                % opposite sign, magnitude 2pi-|alpha|

% RF vectors
pRF = [T.x_px, T.y_px];
v = pRF - s_px;
r = sqrt(sum(v.^2,2));

% store r_s
T.r_s = r;

% unit directions
vhat = v ./ max(r, eps);

% signed shortest angle from uD axis to RF direction
beta = zeros(height(T),1);
for i = 1:height(T)
    beta(i) = signed_angle(uD, vhat(i,:));    % (-pi,pi]
end

% edge half-angles as function of radius
% clamp r to at least rad to avoid asin>1
r_eff = max(r, rad + 1e-6);
Delta = asin(min(1, rad ./ r_eff));           % in [0, pi/2)

DeltaD = Delta;   % same width for both arms
DeltaT = Delta;

% ---- INNER free-space (between edges) ----
alpha_free = alpha - sgn*(DeltaD + DeltaT);   % shrinks inner arc by both edge wedges

% beta measured from distractor EDGE toward target
beta_free = beta - sgn*DeltaD;

% inner membership: point lies between edges along inner wedge
% i.e., beta has same sign as alpha and lies within [DeltaD, alpha-DeltaT]
isInner = (sgn*beta >= DeltaD) & (sgn*beta <= (sgn*alpha - DeltaT)) & (abs(alpha_free) > 1e-6);

% fraction along inner free-space
arc_frac = nan(height(T),1);
arc_frac(isInner) = beta_free(isInner) ./ alpha_free(isInner);

% ---- OUTER free-space (long way around) ----
% Convert beta to long-path equivalent, then remove distractor edge wedge from start.
idxOuter = ~isInner;

% map beta onto the long traversal direction (same sign as alpha_long)
beta_long = beta;
if alpha_long < 0
    % want beta_long in (-2pi, 0]
    tmp = beta(idxOuter);
    tmp(tmp > 0) = tmp(tmp > 0) - 2*pi;
    beta_long(idxOuter) = tmp;
else
    % want beta_long in [0, 2pi)
    tmp = beta(idxOuter);
    tmp(tmp < 0) = tmp(tmp < 0) + 2*pi;
    beta_long(idxOuter) = tmp;
end

% outer free-space gets larger by removing wedges from inner region
alpha_long_free = alpha_long + sgn*(DeltaD + DeltaT);

% OUTER arc starts at OUTER edge of distractor: -sgn*DeltaD
beta_long_free = beta_long + sgn*DeltaD;  % subtract(-sgn*DeltaD)
arc_frac(idxOuter) = beta_long_free(idxOuter) ./ alpha_long_free(idxOuter);

% clamp
arc_frac = max(min(arc_frac, 1), 0);

T.arc_isInner_edge = isInner;
T.arc_frac_edge    = arc_frac;

end

% ---- helpers ----
function u = unit_vec(v)
L = hypot(v(1), v(2));
u = v / max(L, eps);
end

function ang = signed_angle(a, b)
ang = atan2(a(1)*b(2) - a(2)*b(1), a(1)*b(1) + a(2)*b(2));
end
