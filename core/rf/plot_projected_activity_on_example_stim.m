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
%   'CoverageSiteIdx' (default []): optional RF-coverage subset; defaults to SiteIdx
%   'StimIdx'        (default []) : optional subset of stimuli to include (1..384)
%   'OnlyOnObjects'  (default true): only plot target/distractor assigned sites
%   'MarkerSize'     (default 12)
%   'AlphaMax'       (default 0.85)
%   'AlphaThresh'    (default 0.10) : normalized magnitude below this becomes fully transparent
%   'ClipRange'      (default [-1.0 2.0]) : clip normalized values to this range before coloring
%   'PlotMode'       (default 'scatter'): 'scatter' or density-normalized 'smoothed'
%   'SmoothSigmaPx'  (default 15): Gaussian sigma in pixels for smoothed mode
%   'SmoothSupportFraction' (default 0.01): minimum local weight relative to its peak
%   'SmoothActivityAlphaGain' (default 1): display-only alpha gain for activity
%   'SmoothFramePaddingPx' (default 0): white margin around the stimulus frame
%   'CoverageColor'  (default [0.70 0.70 0.70]): neutral RF-coverage color
%   'CoverageAlphaMax' (default 0.55): maximum opacity of RF coverage
%   'ContourLineWidth' (default 2): target/distractor outline width in smoothed mode

p = inputParser;
p.addParameter('TimeBin', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('UseOnlyV1', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SiteIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('CoverageSiteIdx', [], ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('StimIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('OnlyOnObjects', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MarkerSize', 12, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AlphaMax', 0.95, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AlphaThresh', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ClipRange', [-1.0 2.0], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('PlotMode', 'scatter', @(x) ischar(x) || isstring(x));
p.addParameter('SmoothSigmaPx', 15, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('SmoothSupportFraction', 0.01, ...
    @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('SmoothActivityAlphaGain', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('SmoothFramePaddingPx', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
p.addParameter('CoverageColor', [0.70 0.70 0.70], ...
    @(x) isnumeric(x) && numel(x)==3 && all(x>=0) && all(x<=1));
p.addParameter('CoverageAlphaMax', 0.55, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('ContourLineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
opt = p.Results;
plotMode = validatestring(char(opt.PlotMode), {'scatter', 'smoothed'});

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
if ~isempty(opt.CoverageSiteIdx)
    coverageSites = baseSites(opt.CoverageSiteIdx);
else
    coverageSites = sites;
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

% ---- Accumulators for projected samples ----
X = []; Y = []; V = []; C = []; A = [];
XCoverage = []; YCoverage = [];

geom = struct('s_px', s_px, 'uT', uT, 'uD', uD, 'nT', nT, 'nD', nD, ...
    'widthEx', widthEx, 'radEx', radEx, 'alphaEx', alphaEx, ...
    'alphaLongEx', alphaLongEx, 'sgn', sgn);

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
    TAll = Tall(stimNum).T;

    % Restrict response values to selected sites, while coverage may use all RFs.
    [x, y, keep] = project_table_rows( ...
        TAll(sites,:), geom, opt.OnlyOnObjects);
    if strcmp(plotMode, 'smoothed')
        [xCoverage, yCoverage, keepCoverage] = project_table_rows( ...
            TAll(coverageSites,:), geom, opt.OnlyOnObjects);
        XCoverage = [XCoverage; xCoverage(keepCoverage)];
        YCoverage = [YCoverage; yCoverage(keepCoverage)];
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
    V = [V; z(keep)];
    C = [C; rgb(keep,:)];
    A = [A; alpha(keep)];
end

% ---- Plotting ----
if strcmp(plotMode, 'smoothed')
    figColor = [1 1 1];
else
    figColor = [0.5 0.5 0.5];
end
figure('Color', figColor);
ax = axes('Position',[0 0 1 1]); hold(ax,'on');
set(ax,'Position',[0 0 1 1]);
set(ax,'Color', figColor);
axis(ax,'ij');

smoothField = [];
smoothWeight = [];
coverageWeight = [];
nRasterPoints = 0;
nCoverageRasterPoints = 0;

if strcmp(plotMode, 'scatter')
    img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, exampleStimNum);
    imshow(img,'Parent',ax,'InitialMagnification','fit');
    set(ax,'Position',[0 0 1 1]);
    hSc = scatter(ax, X, Y, opt.MarkerSize, C, 'filled');

    % Transparency depends only on activity magnitude.
    alphaValues = A(:);
    alphaValues(~isfinite(alphaValues)) = 0;
    alphaValues = max(0, min(1, alphaValues));

    % Apply per-point alpha when supported; otherwise use uniform alpha.
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
        aMean = mean(alphaValues);
        if isprop(hSc,'MarkerFaceAlpha')
            hSc.MarkerFaceAlpha = aMean;
        end
        if isprop(hSc,'MarkerEdgeAlpha')
            hSc.MarkerEdgeAlpha = aMean;
        end
    end

    % Expand axes a bit to retain all projected scatter points.
    xAll = [X; 1; W];
    yAll = [Y; 1; H];
    margin = 20;
    xlim(ax, [min(xAll)-margin, max(xAll)+margin]);
    ylim(ax, [min(yAll)-margin, max(yAll)+margin]);
    axis(ax,'equal');
    set(ax,'YDir','reverse');

    hFrame = rectangle(ax,'Position',[0.5 0.5 W H], ...
        'EdgeColor',[0.85 0.85 0.85], 'LineWidth',1);
    uistack(hFrame,'top');
else
    [smoothField, smoothWeight, coverageWeight, rgbField, alphaField, ...
        coverageAlpha, nRasterPoints, nCoverageRasterPoints, ...
        displayXLim, displayYLim] = ...
        smooth_projected_field(X, Y, V, XCoverage, YCoverage, W, H, ...
            opt.SmoothSigmaPx, ...
            opt.SmoothSupportFraction, opt.AlphaMax, opt.AlphaThresh, ...
            opt.ClipRange, opt.SmoothActivityAlphaGain, ...
            opt.CoverageAlphaMax, opt.SmoothFramePaddingPx);

    coverageRGB = repmat(reshape(opt.CoverageColor, 1, 1, 3), ...
        [size(coverageAlpha,1) size(coverageAlpha,2) 1]);
    coverageRGB = uint8(round(255 * coverageRGB));
    hCoverage = image(ax, displayXLim, displayYLim, coverageRGB);
    set(hCoverage, 'AlphaData', single(coverageAlpha));

    rgbDisplay = uint8(round(255 * rgbField));
    hIm = image(ax, displayXLim, displayYLim, rgbDisplay);
    set(hIm, 'AlphaData', single(alphaField));

    [~, masks] = render_stim_with_masks2( ...
        ALLCOORDS, RTAB384, exampleStimNum, 'Background', [1 1 1]);
    contour(ax, double(masks.figArm), [0.5 0.5], '-', ...
        'Color', [0.10 0.10 0.10], 'LineWidth', opt.ContourLineWidth);
    contour(ax, double(masks.backArm), [0.5 0.5], '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', opt.ContourLineWidth);

    xlim(ax, displayXLim);
    ylim(ax, displayYLim);
    axis(ax, 'image');
    axis(ax, 'off');
    set(ax, 'YDir', 'reverse');
end

h = struct();
h.fig = gcf;
h.ax  = ax;
h.nPoints = numel(X);
h.timeBin = opt.TimeBin;
h.timeWindow = R.timeWindows(opt.TimeBin,:);
h.plotMode = plotMode;
h.smoothSigmaPx = opt.SmoothSigmaPx;
h.smoothFramePaddingPx = round(opt.SmoothFramePaddingPx);
h.nRasterPoints = nRasterPoints;
h.nCoveragePoints = numel(XCoverage);
h.nCoverageRasterPoints = nCoverageRasterPoints;
h.smoothField = smoothField;
h.smoothWeight = smoothWeight;
h.coverageWeight = coverageWeight;
h.colorScaleValues = [];
h.colorScaleRGB = [];
if strcmp(plotMode, 'smoothed')
    h.xLimits = displayXLim;
    h.yLimits = displayYLim;
    h.colorScaleValues = linspace(opt.ClipRange(1), opt.ClipRange(2), 301);
    [scaleRGB, scaleAlpha] = valueToColorAlpha( ...
        h.colorScaleValues(:), opt.AlphaMax, opt.AlphaThresh, opt.ClipRange);
    scaleRGB = emphasize_smoothed_positive_red( ...
        scaleRGB, h.colorScaleValues(:), opt.ClipRange);
    scaleAlpha = min(opt.SmoothActivityAlphaGain .* scaleAlpha, 1);
    scaleAlpha(abs(h.colorScaleValues(:)) <= opt.AlphaThresh) = 0;
    scaleRGB = scaleAlpha .* scaleRGB + (1-scaleAlpha) .* ones(size(scaleRGB));
    h.colorScaleRGB = reshape(uint8(round(255 .* scaleRGB)), ...
        [1 numel(h.colorScaleValues) 3]);
else
    h.xLimits = xlim(ax);
    h.yLimits = ylim(ax);
end

fprintf('Plotted %d activity points (time bin %d: %g-%g ms)\n', ...
    h.nPoints, opt.TimeBin, h.timeWindow(1), h.timeWindow(2));
if strcmp(plotMode, 'smoothed')
    fprintf(['Smoothed %d in-frame points with Gaussian sigma %.1f px, ' ...
             'support fraction %.3g, and %d px frame padding.\n'], ...
        nRasterPoints, opt.SmoothSigmaPx, opt.SmoothSupportFraction, ...
        round(opt.SmoothFramePaddingPx));
    fprintf('RF coverage: %d projected points, %d inside the display frame.\n', ...
        h.nCoveragePoints, h.nCoverageRasterPoints);
end

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

% -------------------- Helper: project RF table rows --------------------
function [x, y, keep] = project_table_rows(T, geom, onlyOnObjects)
assign = string(T.assignment);
x = nan(height(T),1);
y = nan(height(T),1);

idxT = (assign=="target") & ~isnan(T.along_GC);
if any(idxT)
    along = T.along_GC(idxT) * geom.widthEx;
    perp  = T.perp_signed_GC(idxT) * geom.widthEx;
    pT = geom.s_px + along.*geom.uT + perp.*geom.nT;
    x(idxT) = pT(:,1);
    y(idxT) = pT(:,2);
end

idxD = (assign=="distractor") & ~isnan(T.along_GC);
if any(idxD)
    along = T.along_GC(idxD) * geom.widthEx;
    perp  = T.perp_signed_GC(idxD) * geom.widthEx;
    pD = geom.s_px + along.*geom.uD + perp.*geom.nD;
    x(idxD) = pD(:,1);
    y(idxD) = pD(:,2);
end

idxB = (assign=="background") & ~isnan(T.r_s_GC);
if any(idxB)
    rPx = T.r_s_GC(idxB) * geom.widthEx;
    frac = T.arc_frac_edge(idxB);
    isInner = T.arc_isInner_edge(idxB);
    rEff = max(rPx, geom.radEx + 1e-6);
    delta = asin(min(1, geom.radEx ./ rEff));
    alphaFree = geom.alphaEx - geom.sgn*(delta+delta);
    alphaLongFree = geom.alphaLongEx + geom.sgn*(delta+delta);

    beta = zeros(size(frac));
    beta(isInner) = geom.sgn*delta(isInner) + ...
        frac(isInner).*alphaFree(isInner);
    beta(~isInner) = -geom.sgn*delta(~isInner) + ...
        frac(~isInner).*alphaLongFree(~isInner);

    cb = cos(beta);
    sb = sin(beta);
    vx = cb*geom.uD(1) - sb*geom.uD(2);
    vy = sb*geom.uD(1) + cb*geom.uD(2);
    pB = geom.s_px + [rPx.*vx, rPx.*vy];
    x(idxB) = pB(:,1);
    y(idxB) = pB(:,2);
end

if onlyOnObjects
    keep = (assign=="target") | (assign=="distractor");
else
    keep = isfinite(x) & isfinite(y);
end
end

% -------------------- Helper: smoothed image-space average --------------------
function [zField, weightField, coverageWeight, rgbField, alphaField, ...
        coverageAlpha, nValid, nCoverageValid, xLimits, yLimits] = ...
        smooth_projected_field(X, Y, V, XCoverage, YCoverage, W, H, ...
            sigmaPx, supportFraction, alphaMax, alphaThresh, clipRange, ...
            activityAlphaGain, coverageAlphaMax, framePaddingPx)
framePaddingPx = round(framePaddingPx);
xLimits = [1-framePaddingPx, W+framePaddingPx];
yLimits = [1-framePaddingPx, H+framePaddingPx];
canvasW = W + 2*framePaddingPx;
canvasH = H + 2*framePaddingPx;

xi = round(double(X(:))) - xLimits(1) + 1;
yi = round(double(Y(:))) - yLimits(1) + 1;
v = double(V(:));
valid = isfinite(xi) & isfinite(yi) & isfinite(v) & ...
    xi >= 1 & xi <= canvasW & yi >= 1 & yi <= canvasH;
nValid = nnz(valid);
assert(nValid > 0, 'No finite projected samples fall inside the image frame.');

xiCoverage = round(double(XCoverage(:))) - xLimits(1) + 1;
yiCoverage = round(double(YCoverage(:))) - yLimits(1) + 1;
validCoverage = isfinite(xiCoverage) & isfinite(yiCoverage) & ...
    xiCoverage >= 1 & xiCoverage <= canvasW & ...
    yiCoverage >= 1 & yiCoverage <= canvasH;
nCoverageValid = nnz(validCoverage);
assert(nCoverageValid > 0, ...
    'No finite RF-coverage samples fall inside the image frame.');

lin = sub2ind([canvasH canvasW], yi(valid), xi(valid));
sumMap = reshape(accumarray( ...
    lin, v(valid), [canvasH*canvasW 1], @sum, 0), [canvasH canvasW]);
countMap = reshape(accumarray( ...
    lin, 1, [canvasH*canvasW 1], @sum, 0), [canvasH canvasW]);
linCoverage = sub2ind([canvasH canvasW], ...
    yiCoverage(validCoverage), xiCoverage(validCoverage));
coverageMap = reshape(accumarray( ...
    linCoverage, 1, [canvasH*canvasW 1], @sum, 0), ...
    [canvasH canvasW]);

radius = max(1, ceil(3*sigmaPx));
kx = -radius:radius;
g = exp(-0.5 * (kx./sigmaPx).^2);
g = g / sum(g);

sumBlur = conv2(conv2(sumMap, g, 'same'), g', 'same');
weightField = conv2(conv2(countMap, g, 'same'), g', 'same');
coverageWeight = conv2(conv2(coverageMap, g, 'same'), g', 'same');
zField = sumBlur ./ max(weightField, eps);
zField = min(max(zField, clipRange(1)), clipRange(2));

peakWeight = max(weightField(:));
minWeight = supportFraction * peakWeight;
support = weightField >= minWeight;
zField(~support) = NaN;

zForColor = zField;
zForColor(~isfinite(zForColor)) = 0;
[rgb, valueAlpha] = valueToColorAlpha( ...
    zForColor(:), alphaMax, alphaThresh, clipRange);
rgb = emphasize_smoothed_positive_red(rgb, zForColor(:), clipRange);
rgbField = reshape(rgb, [canvasH canvasW 3]);
valueAlpha = reshape(valueAlpha, [canvasH canvasW]);

supportAlpha = (weightField - minWeight) ./ max(4*minWeight, eps);
supportAlpha = min(max(supportAlpha, 0), 1);
alphaField = activityAlphaGain .* valueAlpha .* supportAlpha;
alphaField(abs(zForColor) <= alphaThresh) = 0;
alphaField = min(alphaField, 1);
alphaField(~support | ~isfinite(alphaField)) = 0;

peakCoverageWeight = max(coverageWeight(:));
minCoverageWeight = supportFraction * peakCoverageWeight;
coverageSupport = coverageWeight >= minCoverageWeight;
coverageAlpha = (coverageWeight - minCoverageWeight) ./ ...
    max(4*minCoverageWeight, eps);
coverageAlpha = coverageAlphaMax * min(max(coverageAlpha, 0), 1);
coverageAlpha(~coverageSupport | ~isfinite(coverageAlpha)) = 0;
end

function rgb = emphasize_smoothed_positive_red(rgb, z, clipRange)
% Keep the smoothed positive-response map red except near its upper limit.
positive = z > 0;
if ~any(positive)
    return;
end

t = min(max(z(positive) ./ clipRange(2), 0), 1);
deepRed = [0.72 0.00 0.00];
brightRed = [1.00 0.10 0.04];
yellow = [1.00 0.92 0.12];

rgbPositive = zeros(nnz(positive), 3);
redRange = t <= 0.8;
redMix = t(redRange) ./ 0.8;
rgbPositive(redRange,:) = ...
    (1-redMix).*deepRed + redMix.*brightRed;

highRange = ~redRange;
if any(highRange)
    highMix = (t(highRange)-0.8) ./ 0.2;
    rgbPositive(highRange,:) = ...
        (1-highMix).*brightRed + highMix.*yellow;
end
rgb(positive,:) = rgbPositive;
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
