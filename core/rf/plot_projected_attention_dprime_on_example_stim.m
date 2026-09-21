function h = plot_projected_attention_dprime_on_example_stim( ...
    Tall, ALLCOORDS, RTAB384, exampleStimNum, R, SNR, varargin)
% PLOT_PROJECTED_ATTENTION_DPRIME_ON_EXAMPLE_STIM
% Project a spatial T-D d-prime map, including RFs on the gray background.
%
% Stimulus pairs [1 6], [2 5], [3 8], and [4 7] (repeated per block of 8)
% exchange target and distractor while preserving the matched visual layout.
% Both directions are projected into the canonical target/distractor frame.
% Local Gaussian-weighted first and second moments then define d-prime.

p = inputParser;
p.addParameter('TimeBin', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('SiteIdx', [], ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('CoverageSiteIdx', [], ...
    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('ExcludeOverlap', true, ...
    @(x) islogical(x) && isscalar(x));
p.addParameter('SmoothSigmaPx', 15, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('SmoothSupportFraction', 0.01, ...
    @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('FramePaddingPx', 100, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
p.addParameter('ClipRange', [-0.5 0.5], ...
    @(x) isnumeric(x) && numel(x)==2 && x(1)<0 && x(2)>0);
p.addParameter('AlphaMax', 0.95, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('AlphaThresh', 0.02, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('AlphaGain', 2.5, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('CoverageColor', [0.70 0.70 0.70], ...
    @(x) isnumeric(x) && numel(x)==3 && all(x>=0) && all(x<=1));
p.addParameter('CoverageAlphaMax', 0.55, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('ContourLineWidth', 4, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
opt = p.Results;

nSites = size(R.meanAct, 1);
assert(size(R.meanAct,2) == 384 && size(R.meanSqAct,2) == 384, ...
    'Expected response moments for 384 stimuli.');
assert(size(R.meanAct,3) >= opt.TimeBin && ...
       size(R.meanSqAct,3) >= opt.TimeBin, ...
    'TimeBin exceeds the available response windows.');
assert(numel(Tall) == 384, 'Tall must contain 384 stimulus tables.');

if isempty(opt.SiteIdx)
    sites = (1:nSites).';
else
    sites = unique(double(opt.SiteIdx(:)), 'stable');
end
if isempty(opt.CoverageSiteIdx)
    coverageSites = (1:nSites).';
else
    coverageSites = unique(double(opt.CoverageSiteIdx(:)), 'stable');
end
assert(all(sites>=1 & sites<=nSites & sites==floor(sites)), ...
    'SiteIdx contains invalid local site indices.');
assert(all(coverageSites>=1 & coverageSites<=nSites & ...
    coverageSites==floor(coverageSites)), ...
    'CoverageSiteIdx contains invalid local site indices.');

W = 1024;
H = 768;
geom = canonical_geometry(ALLCOORDS, RTAB384, exampleStimNum, W, H);

b = double(SNR.muSpont(:));
topMat = [double(SNR.muYellowEarly(:)), double(SNR.muYellowLate(:)), ...
    double(SNR.muPurpleEarly(:)), double(SNR.muPurpleLate(:))];
assert(numel(b)>=nSites && size(topMat,1)>=nSites, ...
    'SNR normalization vectors are smaller than the response site range.');
scale = max(topMat, [], 2) - b;
scale(~isfinite(scale) | scale<=0) = NaN;

[nTrials, perSiteTrials] = normalize_trial_counts(R.nTrials, nSites);
pairs = attention_pairs();

X = [];
Y = [];
MuT = [];
M2T = [];
MuD = [];
M2D = [];
Weight = [];
Group = [];
XCoverage = [];
YCoverage = [];

requiredVars = ["assignment","overlap","center_color", ...
    "along_GC","perp_signed_GC","r_s_GC", ...
    "arc_frac_edge","arc_isInner_edge"];

for pairIdx = 1:size(pairs,1)
    for direction = 1:2
        if direction == 1
            stimT = pairs(pairIdx,1);
            stimD = pairs(pairIdx,2);
        else
            stimT = pairs(pairIdx,2);
            stimD = pairs(pairIdx,1);
        end

        TT = Tall(stimT).T;
        TD = Tall(stimD).T;
        if ~istable(TT) || ~istable(TD) || ...
                height(TT)<nSites || height(TD)<nSites || ...
                any(~ismember(requiredVars, string(TT.Properties.VariableNames))) || ...
                any(~ismember(requiredVars, string(TD.Properties.VariableNames)))
            continue;
        end

        [xCoverage, yCoverage, keepCoverage] = ...
            project_rows(TT(coverageSites,:), geom, false);
        XCoverage = [XCoverage; xCoverage(keepCoverage)]; %#ok<AGROW>
        YCoverage = [YCoverage; yCoverage(keepCoverage)]; %#ok<AGROW>

        TS = TT(sites,:);
        DS = TD(sites,:);
        [x, y, keepGeometry, group] = project_rows( ...
            TS, geom, opt.ExcludeOverlap);

        assignT = string(TS.assignment);
        assignD = string(DS.assignment);
        compatible = ...
            (assignT=="target" & assignD=="distractor") | ...
            (assignT=="distractor" & assignD=="target") | ...
            (assignT=="background" & assignD=="background");
        compatible = compatible & ...
            (string(TS.center_color) == string(DS.center_color));
        if opt.ExcludeOverlap
            compatible = compatible & ~logical(TS.overlap) & ~logical(DS.overlap);
        end

        [muT, m2T] = normalized_moments( ...
            R, sites, stimT, opt.TimeBin, b, scale);
        [muD, m2D] = normalized_moments( ...
            R, sites, stimD, opt.TimeBin, b, scale);
        w = matched_trial_weight(nTrials, perSiteTrials, sites, stimT, stimD);

        keep = keepGeometry & compatible & isfinite(muT) & isfinite(m2T) & ...
            isfinite(muD) & isfinite(m2D) & isfinite(w) & w>0;
        X = [X; x(keep)]; %#ok<AGROW>
        Y = [Y; y(keep)]; %#ok<AGROW>
        MuT = [MuT; muT(keep)]; %#ok<AGROW>
        M2T = [M2T; m2T(keep)]; %#ok<AGROW>
        MuD = [MuD; muD(keep)]; %#ok<AGROW>
        M2D = [M2D; m2D(keep)]; %#ok<AGROW>
        Weight = [Weight; w(keep)]; %#ok<AGROW>
        Group = [Group; group(keep)]; %#ok<AGROW>
    end
end

assert(~isempty(X), 'No matched attention samples were available for plotting.');

[dprimeField, localWeight, coverageWeight, rgbField, alphaField, ...
    coverageAlpha, nRasterPoints, nCoverageRasterPoints, xLimits, yLimits] = ...
    rasterize_dprime(X, Y, MuT, M2T, MuD, M2D, Weight, ...
        XCoverage, YCoverage, W, H, opt);

fig = figure('Color', 'w');
ax = axes('Position', [0 0 1 1]);
hold(ax, 'on');
set(ax, 'Color', 'w', 'YDir', 'reverse');

coverageRGB = repmat(reshape(opt.CoverageColor, 1, 1, 3), ...
    [size(coverageAlpha,1) size(coverageAlpha,2) 1]);
hCoverage = image(ax, xLimits, yLimits, ...
    uint8(round(255 .* coverageRGB)));
set(hCoverage, 'AlphaData', single(coverageAlpha));

hActivity = image(ax, xLimits, yLimits, ...
    uint8(round(255 .* rgbField)));
set(hActivity, 'AlphaData', single(alphaField));

[~, masks] = render_stim_with_masks2( ...
    ALLCOORDS, RTAB384, exampleStimNum, 'Background', [1 1 1]);
contour(ax, double(masks.figArm), [0.5 0.5], '-', ...
    'Color', [0.05 0.05 0.05], 'LineWidth', opt.ContourLineWidth);
contour(ax, double(masks.backArm), [0.5 0.5], '--', ...
    'Color', [0.18 0.18 0.18], 'LineWidth', opt.ContourLineWidth);

xlim(ax, xLimits);
ylim(ax, yLimits);
axis(ax, 'image');
axis(ax, 'off');

[scaleValues, scaleRGB] = make_color_scale(opt);
finiteField = dprimeField(isfinite(dprimeField));

h = struct();
h.fig = fig;
h.ax = ax;
h.timeBin = opt.TimeBin;
h.timeWindow = double(R.timeWindows(opt.TimeBin,:));
h.nSites = numel(sites);
h.nSamples = numel(X);
h.nTargetSamples = nnz(Group==1);
h.nDistractorSamples = nnz(Group==2);
h.nBackgroundSamples = nnz(Group==3);
h.nRasterPoints = nRasterPoints;
h.nCoveragePoints = numel(XCoverage);
h.nCoverageRasterPoints = nCoverageRasterPoints;
h.dprimeField = dprimeField;
h.localWeight = localWeight;
h.coverageWeight = coverageWeight;
h.xLimits = xLimits;
h.yLimits = yLimits;
h.clipRange = opt.ClipRange;
h.colorScaleValues = scaleValues;
h.colorScaleRGB = scaleRGB;
h.fieldRange = [min(finiteField) max(finiteField)];
h.fractionClipped = mean(finiteField<opt.ClipRange(1) | ...
    finiteField>opt.ClipRange(2));

fprintf(['Attention d-prime samples: target=%d, distractor=%d, ' ...
    'background=%d (total=%d).\n'], ...
    h.nTargetSamples, h.nDistractorSamples, h.nBackgroundSamples, h.nSamples);
fprintf(['D-prime field range=[%.4g %.4g], %.2f%% of supported pixels ' ...
    'outside display range [%.3g %.3g].\n'], ...
    h.fieldRange, 100*h.fractionClipped, opt.ClipRange);
end

function geom = canonical_geometry(ALLCOORDS, RTAB384, stimNum, W, H)
toPx = @(q) [q(1)+W/2, H/2-q(2)];
f = sprintf('stim_%d', stimNum);
sPx = toPx(double(ALLCOORDS.(f).s(:))');
tT = toPx(double(ALLCOORDS.(f).t_fig(:))');
tD = toPx(double(ALLCOORDS.(f).t_back(:))');
uT = (tT-sPx) ./ norm(tT-sPx);
uD = (tD-sPx) ./ norm(tD-sPx);
nT = perpendicular_toward(sPx, tT, tD);
nD = perpendicular_toward(sPx, tD, tT);
if istable(RTAB384)
    widthEx = double(RTAB384{stimNum,7});
else
    widthEx = double(RTAB384(stimNum,7));
end
radEx = widthEx/2;
alphaEx = signed_angle(uD, uT);
sgn = sign(alphaEx);
if sgn==0, sgn=1; end
geom = struct('sPx',sPx,'uT',uT,'uD',uD,'nT',nT,'nD',nD, ...
    'widthEx',widthEx,'radEx',radEx,'alphaEx',alphaEx, ...
    'alphaLongEx',alphaEx-sgn*2*pi,'sgn',sgn);
end

function pairs = attention_pairs()
basePairs = [1 6; 2 5; 3 8; 4 7];
pairs = zeros(192,2);
row = 0;
for block = 0:47
    pairs(row+(1:4),:) = basePairs + 8*block;
    row = row+4;
end
end

function [nTrials, perSite] = normalize_trial_counts(nTrials, nSites)
if isvector(nTrials)
    nTrials = double(nTrials(:)');
    assert(numel(nTrials)==384, 'Trial-count vector must contain 384 values.');
    perSite = false;
else
    assert(size(nTrials,1)>=nSites && size(nTrials,2)==384, ...
        'Per-site trial counts must be [nSites x 384].');
    nTrials = double(nTrials);
    perSite = true;
end
end

function w = matched_trial_weight(nTrials, perSite, sites, stimT, stimD)
if perSite
    w = min(nTrials(sites,stimT), nTrials(sites,stimD));
else
    w = repmat(min(nTrials(stimT),nTrials(stimD)), numel(sites), 1);
end
end

function [mu, m2] = normalized_moments(R, sites, stim, tb, b, scale)
ex = squeeze(double(R.meanAct(sites,stim,tb)));
ex2 = squeeze(double(R.meanSqAct(sites,stim,tb)));
bs = b(sites);
ss = scale(sites);
mu = (ex-bs) ./ ss;
m2 = (ex2-2.*bs.*ex+bs.^2) ./ ss.^2;
end

function [x, y, keep, group] = project_rows(T, geom, excludeOverlap)
assign = string(T.assignment);
x = nan(height(T),1);
y = nan(height(T),1);
group = zeros(height(T),1);

idxT = assign=="target" & isfinite(T.along_GC);
idxD = assign=="distractor" & isfinite(T.along_GC);
idxB = assign=="background" & isfinite(T.r_s_GC);
group(idxT) = 1;
group(idxD) = 2;
group(idxB) = 3;

if any(idxT)
    along = double(T.along_GC(idxT))*geom.widthEx;
    perp = double(T.perp_signed_GC(idxT))*geom.widthEx;
    p = geom.sPx + along.*geom.uT + perp.*geom.nT;
    x(idxT) = p(:,1);
    y(idxT) = p(:,2);
end
if any(idxD)
    along = double(T.along_GC(idxD))*geom.widthEx;
    perp = double(T.perp_signed_GC(idxD))*geom.widthEx;
    p = geom.sPx + along.*geom.uD + perp.*geom.nD;
    x(idxD) = p(:,1);
    y(idxD) = p(:,2);
end
if any(idxB)
    rPx = double(T.r_s_GC(idxB))*geom.widthEx;
    frac = double(T.arc_frac_edge(idxB));
    isInner = logical(T.arc_isInner_edge(idxB));
    rEff = max(rPx, geom.radEx+1e-6);
    delta = asin(min(1,geom.radEx./rEff));
    alphaFree = geom.alphaEx-geom.sgn*(delta+delta);
    alphaLongFree = geom.alphaLongEx+geom.sgn*(delta+delta);
    beta = zeros(size(frac));
    beta(isInner) = geom.sgn*delta(isInner) + ...
        frac(isInner).*alphaFree(isInner);
    beta(~isInner) = -geom.sgn*delta(~isInner) + ...
        frac(~isInner).*alphaLongFree(~isInner);
    cb = cos(beta);
    sb = sin(beta);
    vx = cb*geom.uD(1)-sb*geom.uD(2);
    vy = sb*geom.uD(1)+cb*geom.uD(2);
    p = geom.sPx + [rPx.*vx,rPx.*vy];
    x(idxB) = p(:,1);
    y(idxB) = p(:,2);
end

keep = isfinite(x) & isfinite(y) & group>0;
if excludeOverlap && ismember('overlap',T.Properties.VariableNames)
    keep = keep & ~logical(T.overlap);
end
end

function [dprimeField, weightField, coverageWeight, rgbField, alphaField, ...
        coverageAlpha, nValid, nCoverageValid, xLimits, yLimits] = ...
        rasterize_dprime(X,Y,MuT,M2T,MuD,M2D,Weight, ...
            XCoverage,YCoverage,W,H,opt)
pad = round(opt.FramePaddingPx);
xLimits = [1-pad,W+pad];
yLimits = [1-pad,H+pad];
canvasW = W+2*pad;
canvasH = H+2*pad;

xi = round(double(X(:)))-xLimits(1)+1;
yi = round(double(Y(:)))-yLimits(1)+1;
w = double(Weight(:));
valid = isfinite(xi) & isfinite(yi) & isfinite(w) & w>0 & ...
    isfinite(MuT) & isfinite(M2T) & isfinite(MuD) & isfinite(M2D) & ...
    xi>=1 & xi<=canvasW & yi>=1 & yi<=canvasH;
nValid = nnz(valid);
assert(nValid>0, 'No attention samples fall inside the display frame.');
lin = sub2ind([canvasH canvasW],yi(valid),xi(valid));

weightMap = map_sum(lin,w(valid),canvasH,canvasW);
sumTMap = map_sum(lin,w(valid).*MuT(valid),canvasH,canvasW);
sumT2Map = map_sum(lin,w(valid).*M2T(valid),canvasH,canvasW);
sumDMap = map_sum(lin,w(valid).*MuD(valid),canvasH,canvasW);
sumD2Map = map_sum(lin,w(valid).*M2D(valid),canvasH,canvasW);

xc = round(double(XCoverage(:)))-xLimits(1)+1;
yc = round(double(YCoverage(:)))-yLimits(1)+1;
validCoverage = isfinite(xc) & isfinite(yc) & ...
    xc>=1 & xc<=canvasW & yc>=1 & yc<=canvasH;
nCoverageValid = nnz(validCoverage);
linCoverage = sub2ind([canvasH canvasW], ...
    yc(validCoverage),xc(validCoverage));
coverageMap = map_sum(linCoverage,ones(nCoverageValid,1),canvasH,canvasW);

radius = max(1,ceil(3*opt.SmoothSigmaPx));
kx = -radius:radius;
g = exp(-0.5*(kx./opt.SmoothSigmaPx).^2);
g = g./sum(g);
blur = @(m) conv2(conv2(m,g,'same'),g','same');

weightField = blur(weightMap);
sumT = blur(sumTMap);
sumT2 = blur(sumT2Map);
sumD = blur(sumDMap);
sumD2 = blur(sumD2Map);
coverageWeight = blur(coverageMap);

muT = sumT./max(weightField,eps);
muD = sumD./max(weightField,eps);
varT = max(0,sumT2./max(weightField,eps)-muT.^2);
varD = max(0,sumD2./max(weightField,eps)-muD.^2);
denom = sqrt(0.5*(varT+varD));
dprimeField = (muT-muD)./denom;

minWeight = opt.SmoothSupportFraction*max(weightField(:));
support = weightField>=minWeight & isfinite(dprimeField) & denom>1e-6;
dprimeField(~support) = NaN;

displayValue = min(max(dprimeField,opt.ClipRange(1)),opt.ClipRange(2));
displayValue(~isfinite(displayValue)) = 0;
[rgbField,alphaField] = signed_display(displayValue,opt);
supportAlpha = (weightField-minWeight)./max(4*minWeight,eps);
supportAlpha = min(max(supportAlpha,0),1);
alphaField = alphaField.*supportAlpha;
alphaField(~support) = 0;

minCoverage = opt.SmoothSupportFraction*max(coverageWeight(:));
coverageAlpha = (coverageWeight-minCoverage)./max(4*minCoverage,eps);
coverageAlpha = opt.CoverageAlphaMax*min(max(coverageAlpha,0),1);
coverageAlpha(coverageWeight<minCoverage | ~isfinite(coverageAlpha)) = 0;
end

function m = map_sum(lin,v,H,W)
m = reshape(accumarray(lin,v,[H*W 1],@sum,0),[H W]);
end

function [rgbField,alphaField] = signed_display(z,opt)
z = double(z);
rgbField = ones([size(z),3]);
neg = z<0;
pos = z>0;
if any(neg(:))
    t = min(abs(z(neg))./abs(opt.ClipRange(1)),1);
    t = t(:);
    c0 = [0.35 0.68 0.95];
    c1 = [0.03 0.22 0.72];
    c = (1-t).*c0+t.*c1;
    for j=1:3
        layer = rgbField(:,:,j);
        layer(neg) = c(:,j);
        rgbField(:,:,j) = layer;
    end
end
if any(pos(:))
    t = min(z(pos)./opt.ClipRange(2),1);
    t = t(:);
    c0 = [1.00 0.38 0.38];
    c1 = [0.82 0.00 0.02];
    c = (1-t).*c0+t.*c1;
    for j=1:3
        layer = rgbField(:,:,j);
        layer(pos) = c(:,j);
        rgbField(:,:,j) = layer;
    end
end
magMax = max(abs(opt.ClipRange));
tAlpha = (abs(z)-opt.AlphaThresh)./max(magMax-opt.AlphaThresh,eps);
tAlpha = min(max(tAlpha,0),1);
alphaField = min(1,opt.AlphaGain*opt.AlphaMax*tAlpha);
alphaField(abs(z)<=opt.AlphaThresh) = 0;
end

function [values,rgb] = make_color_scale(opt)
values = linspace(opt.ClipRange(1),opt.ClipRange(2),301);
[color,alpha] = signed_display(values,opt);
alpha = reshape(alpha,[1 numel(values) 1]);
rgb = alpha.*color+(1-alpha).*ones(size(color));
rgb = uint8(round(255.*rgb));
end

function n = perpendicular_toward(sPx,tArm,tOther)
u = (tArm-sPx)./norm(tArm-sPx);
w = tOther-sPx;
wPerp = w-dot(w,u)*u;
n = wPerp./norm(wPerp);
end

function a = signed_angle(u,v)
a = atan2(u(1)*v(2)-u(2)*v(1),dot(u,v));
end
