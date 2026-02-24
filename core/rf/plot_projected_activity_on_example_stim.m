function h = plot_projected_activity_on_example_stim(Tall, ALLCOORDS, RTAB384, exampleStimNum, R, SNR, varargin)
% PLOT_PROJECTED_ACTIVITY_ON_EXAMPLE_STIM
%
% Plots normalized V1 activity (false color) at RF locations projected into the
% coordinate frame of an example stimulus, using the same GC-based projection
% as plot_projected_RFs_on_example_stim.
%
% Complementary stimuli: within each block of 8, pairs are (1,5), (2,6), (3,7), (4,8)
% and this repeats for 9-16, 17-24, ... up to 384.
%
% Inputs:
%   Tall: struct array where Tall(stimNum).T is a table with RF location fields:
%         assignment, along_GC, perp_signed_GC, r_s_GC, arc_frac_edge, arc_isInner_edge
%   ALLCOORDS, RTAB384: stimulus geometry tables/structs used for rendering
%   exampleStimNum: stimulus number that defines canonical frame
%   R.meanAct: [1024 x 384 x 70] double
%   R.timeWindows: [70 x 2] (ms)
%   SNR: struct with fields muSpont, muYellowEarly, muYellowLate, muPurpleEarly, muPurpleLate
%
% Options (name/value):
%   'TimeBin'        (default 1)  : which of the 70 bins to use
%   'UseOnlyV1'      (default true): use sites 1:512 only
%   'SiteIdx'        (default []) : optional subset (indices within chosen site range)
%   'StimIdx'        (default []) : optional subset of stimuli to include (1..384)
%   'OnlyOnObjects'  (default true): only plot target/distractor assigned sites
%   'MarkerSize'     (default 12)
%   'AlphaMax'       (default 0.85)
%   'AlphaThresh'    (default 0.10) : normalized magnitude below this becomes fully transparent
%   'ClipRange'      (default [-1.0 2.0]) : clip normalized values to this range before coloring

p = inputParser;
p.addParameter('TimeBin', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('UseOnlyV1', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SiteIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('StimIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('OnlyOnObjects', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MarkerSize', 12, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AlphaMax', 0.95, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AlphaThresh', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ClipRange', [-1.0 2.0], @(x) isnumeric(x) && numel(x)==2);
p.parse(varargin{:});
opt = p.Results;

W = 1024; H = 768;
toPx = @(q) [q(1) + W/2, H/2 - q(2)];

% ---- Example stimulus geometry (same as your RF routine) ----
fieldName = sprintf('stim_%d', exampleStimNum);
s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

s_px = toPx(s);
tT   = toPx(tFig);
tD   = toPx(tBack);

uT = (tT - s_px) / norm(tT - s_px);
uD = (tD - s_px) / norm(tD - s_px);

nT = perpTowardOther(s_px, tT, tD);
nD = perpTowardOther(s_px, tD, tT);

widthEx = double(RTAB384(exampleStimNum,7));
radEx   = widthEx/2;

alphaEx = signedAngle(uD, uT);
sgn = sign(alphaEx); if sgn==0, sgn=1; end
alphaLongEx = alphaEx - sgn*2*pi;

% ---- Determine site set ----
if opt.UseOnlyV1
    baseSites = 1:512;
else
    baseSites = 1:size(R.meanAct,1);
end
if ~isempty(opt.SiteIdx)
    sites = baseSites(opt.SiteIdx);
else
    sites = baseSites;
end

% ---- Stimulus set ----
if isempty(opt.StimIdx)
    stimList = 1:384;
else
    stimList = opt.StimIdx(:)';
end

% ---- Normalization vectors (per site) ----
muSpont = double(SNR.muSpont(:));
muTop   = max( [double(SNR.muYellowEarly(:)), double(SNR.muYellowLate(:)), ...
                double(SNR.muPurpleEarly(:)), double(SNR.muPurpleLate(:))], [], 2 );

% Scale: (top - spont) is typically what you want if muTop includes baseline
scale = muTop - muSpont;
scale(scale <= 1e-9) = 1e-9;

% ---- Accumulators for scatter ----
X = []; Y = []; C = []; A = [];

% To avoid double counting complementary pairs, only take the "lower" member of each pair.
seenPair = false(384,1);

for stimNum = stimList
    comp = complementaryStim(stimNum);
    pairKey = min(stimNum, comp);

    if seenPair(pairKey), continue; end
    seenPair(pairKey) = true;

    % Table with RF locations projected into example frame (assumed by your pipeline)
    if stimNum > numel(Tall) || ~isfield(Tall(stimNum),'T')
        continue;
    end
    T = Tall(stimNum).T;

    % Restrict rows to requested sites (assumes row order matches site index)
    % If your table has a site-id column instead, swap this selection accordingly.
    T = T(sites, :);

    assign = string(T.assignment);

    % Compute positions for target/distractor/background entries (same formulas)
    x = nan(height(T),1);
    y = nan(height(T),1);

    % Target
    idxT = (assign=="target") & ~isnan(T.along_GC);
    if any(idxT)
        along = T.along_GC(idxT) * widthEx;
        perp  = T.perp_signed_GC(idxT) * widthEx;
        pT = s_px + along.*uT + perp.*nT;
        x(idxT) = pT(:,1); y(idxT) = pT(:,2);
    end

    % Distractor
    idxD = (assign=="distractor") & ~isnan(T.along_GC);
    if any(idxD)
        along = T.along_GC(idxD) * widthEx;
        perp  = T.perp_signed_GC(idxD) * widthEx;
        pD = s_px + along.*uD + perp.*nD;
        x(idxD) = pD(:,1); y(idxD) = pD(:,2);
    end

    % Background (edge-based arc)
    idxB = (assign=="background") & ~isnan(T.r_s_GC);
    if any(idxB)
        r_px = T.r_s_GC(idxB) * widthEx;
        frac = T.arc_frac_edge(idxB);
        isIn = T.arc_isInner_edge(idxB);

        r_eff = max(r_px, radEx + 1e-6);
        Delta = asin(min(1, radEx ./ r_eff));

        alphaFree     = alphaEx     - sgn*(Delta+Delta);
        alphaLongFree = alphaLongEx + sgn*(Delta+Delta);

        beta = zeros(size(frac));
        beta(isIn)  =  sgn*Delta(isIn) + frac(isIn).*alphaFree(isIn);
        beta(~isIn) = -sgn*Delta(~isIn) + frac(~isIn).*alphaLongFree(~isIn);

        cb = cos(beta); sb = sin(beta);
        vx = cb*uD(1) - sb*uD(2);
        vy = sb*uD(1) + cb*uD(2);

        pB = s_px + [r_px.*vx, r_px.*vy];
        x(idxB) = pB(:,1); y(idxB) = pB(:,2);
    end

    % Optional: keep only object-assigned points
    if opt.OnlyOnObjects
        keep = (assign=="target") | (assign=="distractor");
    else
        keep = ~isnan(x) & ~isnan(y);
    end

    if ~any(keep), continue; end

    % ---- Activity: average across complementary stimuli ----
    tb = opt.TimeBin;
    act1 = squeeze(double(R.meanAct(sites, stimNum, tb)));
    act2 = squeeze(double(R.meanAct(sites, comp,    tb)));
    act  = 0.5*(act1 + act2);

    % baseline subtract + normalize (per site)
    z = (act - muSpont(sites)) ./ scale(sites);

    % Clip range for stable coloring
    z = min(max(z, opt.ClipRange(1)), opt.ClipRange(2));

    % Map to RGB + alpha
    [rgb, alpha] = valueToColorAlpha(z, opt.AlphaMax, opt.AlphaThresh, opt.ClipRange);

    % Keep only plotted points (same indexing as sites)
    X = [X; x(keep)];
    Y = [Y; y(keep)];
    C = [C; rgb(keep,:)];
    A = [A; alpha(keep)];
end

% ---- Plotting ----
figure('Color',[0.5 0.5 0.5]);
ax = axes('Position',[0 0 1 1]); hold(ax,'on');

img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, exampleStimNum);
imshow(img,'Parent',ax,'InitialMagnification','fit');
set(ax,'Position',[0 0 1 1]);
set(ax,'Color',[0.5 0.5 0.5]);
axis(ax,'ij');

hSc = scatter(ax, X, Y, opt.MarkerSize, C, 'filled');

% ---- Transparency depends ONLY on activity magnitude (alpha from valueToColorAlpha) ----
% A is per-point alpha computed from z (near 0 -> transparent; large |z| -> opaque).
alphaValues = A(:);
alphaValues(~isfinite(alphaValues)) = 0;
alphaValues = max(0, min(1, alphaValues));

% Apply (version-safe): prefer per-point alpha, fallback to uniform alpha.
appliedPerPointAlpha = false;
if isprop(hSc,'AlphaData') && isprop(hSc,'MarkerFaceAlpha')
    try
        hSc.MarkerFaceAlpha = 'flat';
        hSc.AlphaData       = alphaValues;
        if isprop(hSc,'AlphaDataMapping')
            hSc.AlphaDataMapping = 'none';
        end
        if isprop(hSc,'MarkerEdgeAlpha')
            hSc.MarkerEdgeAlpha = 'flat';
        end
        appliedPerPointAlpha = true;
    catch
        appliedPerPointAlpha = false;
    end
end

if ~appliedPerPointAlpha
    % Older MATLAB fallback: uniform alpha only
    aMean = mean(alphaValues);
    if isprop(hSc,'MarkerFaceAlpha')
        hSc.MarkerFaceAlpha = aMean;
    end
    if isprop(hSc,'MarkerEdgeAlpha')
        hSc.MarkerEdgeAlpha = aMean;
    end
end

% Expand axes a bit
xAll = [X; 1; W];
yAll = [Y; 1; H];
margin = 20;
xlim(ax, [min(xAll)-margin, max(xAll)+margin]);
ylim(ax, [min(yAll)-margin, max(yAll)+margin]);
axis(ax,'equal');
set(ax,'YDir','reverse');

hFrame = rectangle(ax,'Position',[0.5 0.5 W H], 'EdgeColor',[0.85 0.85 0.85], 'LineWidth',1);
uistack(hFrame,'top');

h = struct();
h.fig = gcf;
h.ax  = ax;
h.nPoints = numel(X);
h.timeBin = opt.TimeBin;
h.timeWindow = R.timeWindows(opt.TimeBin,:);

fprintf('Plotted %d activity points (time bin %d: %g-%g ms)\n', ...
    h.nPoints, opt.TimeBin, h.timeWindow(1), h.timeWindow(2));

end

% -------------------- Helper: complementary pairing --------------------
function comp = complementaryStim(i)
% i in 1..384
block = floor((i-1)/8);          % 0..47
pos   = mod(i-1,8) + 1;          % 1..8
if pos <= 4
    comp = block*8 + (pos+4);
else
    comp = block*8 + (pos-4);
end
end

% -------------------- Helper: color/alpha mapping --------------------
function [rgb, alpha] = valueToColorAlpha(z, alphaMax, alphaThresh, clipRange)
%VALUE TO COLOR + ALPHA (perceptually more balanced)
% z should already be clipped to [clipRange(1) clipRange(2)].

zmin = clipRange(1);
zmax = clipRange(2);

% --- Alpha from magnitude with dead-zone around 0 ---
mag = abs(z);
magMax = max(abs([zmin zmax]));
% Gamma shaping makes small magnitudes fade more smoothly (near-zero more transparent)
alphaGamma = 1.0;
alphaMinNonZero = 0.08;  % keep weak-but-nonzero activity slightly visible

t = (mag - alphaThresh) ./ max(1e-9, (magMax - alphaThresh));  % 0..1
t = min(max(t, 0), 1);
t = t .^ alphaGamma;

alpha = alphaMinNonZero + (alphaMax - alphaMinNonZero) .* t;
alpha(mag <= 0) = 0;

rgb = zeros(numel(z),3);

neg = z < 0;
pos = z > 0;
zer = ~neg & ~pos;

% --- Colors ---
% Baseline should be light/neutral (lets stimulus show through when alpha small)
c0_neg = [0.88 0.93 1.00];   % very light blue (near 0-)
c1_neg = [0.00 0.70 1.00];   % bright sky-blue (strong negative)  <-- key change

cR = [1.00 0.10 0.10];       % red
cO = [1.00 0.60 0.05];       % orange
cY = [1.00 1.00 0.15];       % yellow (high)

% Negative: interpolate 0 -> zmin using high-luminance blue
if any(neg)
    % t = 0 at 0, t = 1 at zmin (zmin is negative)
    t = min(max(z(neg) / zmin, 0), 1);  % dividing by negative flips sign as desired
    rgb(neg,:) = (1-t).*c0_neg + t.*c1_neg;
end

% Positive: red -> orange -> yellow
if any(pos)
    t = min(max(z(pos) / zmax, 0), 1); % 0..1
    rgbPos = zeros(sum(pos),3);

    % 0..0.6: red -> orange
    t1 = min(t/0.6, 1);
    rgbPos = (1-t1).*cR + t1.*cO;

    % 0.6..1: orange -> yellow
    hi = t > 0.6;
    if any(hi)
        t2 = (t(hi)-0.6)/0.4;
        rgbPos(hi,:) = (1-t2).*cO + t2.*cY;
    end

    rgb(pos,:) = rgbPos;
end

% Exactly zero: neutral (won’t matter much because alpha ~ 0)
rgb(zer,:) = repmat([1 1 1], sum(zer), 1);
end

% -------------------- Geometry helpers (same as your file) --------------------
function n = perpTowardOther(s_px, t_arm, t_other)
v = t_arm - s_px;
u = v / norm(v);
w = t_other - s_px;
w_perp = w - dot(w,u)*u;
n = w_perp / norm(w_perp);
end

function ang = signedAngle(a,b)
ang = atan2(a(1)*b(2)-a(2)*b(1), dot(a,b));
end
