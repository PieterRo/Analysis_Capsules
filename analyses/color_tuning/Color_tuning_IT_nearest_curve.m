function OUT = Color_tuning_IT_nearest_curve(varargin)
%COLOR_TUNING_IT_NEAREST_CURVE IT color tuning from nearest curve distance.
%
% For each RF center and stimulus, distance is measured to the nearest
% visible yellow and purple stimulus pixel. Complementary color-swap pairs
% are then compared with a weighted regression:
%
%   normalized response(a) - normalized response(b)
%       = intercept + beta * signed nearest-color distance + error
%
% Positive beta means a stronger response when yellow is nearer to the RF
% center; negative beta means a stronger response when purple is nearer.

p = inputParser;
p.addParameter('Monkey', 'N', @(x) any(strcmpi(string(x),["N" "F"])));
p.addParameter('Config', struct(), @isstruct);
p.addParameter('ResponseSNRThreshold', 0.7, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ResponseTopKQuartets', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('MinPairs', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 3 && x == floor(x));
p.addParameter('PixelsPerDegree', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
p.addParameter('Alpha', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(exist('bwdist', 'file') == 2, ...
    'Color_tuning_IT_nearest_curve requires Image Processing Toolbox.');

cfg = config();
if ~isempty(fieldnames(opt.Config)); cfg = opt.Config; end
suffix = upper(char(opt.Monkey));
monkeyName = 'Mr Nilson'; if suffix=='F'; monkeyName='Figaro'; end
geometryFile = fullfile(cfg.matDir, ['Tall_IT_lines_' suffix '.mat']);
responseFile = fullfile(cfg.matDir, ['SNR_capsules_' suffix '_d12.mat']);
assert(isfile(geometryFile), 'Missing %s. Run Line_Stimuli_IT first.', geometryFile);
assert(isfile(responseFile), 'Missing %s.', responseFile);

G = load(geometryFile, 'Tall_IT', 'ALLCOORDS', 'RTAB384', 'RFrange');
S = load(responseFile, 'R');
validateInputs(G, S.R);

if isempty(opt.PixelsPerDegree)
    assert(suffix=='N', 'Figaro requires explicit PixelsPerDegree; use 1 for native-pixel selection tests.');
    calibrationFile = fullfile(cfg.extrasRoot, 'monkeyN', 'RFs', ...
        'BarMap_Nilson.mat');
    C = load(calibrationFile, 'pixperdeg');
    pixelsPerDegree = double(C.pixperdeg);
else
    calibrationFile = '';
    pixelsPerDegree = double(opt.PixelsPerDegree);
end

globalSites = double(G.RFrange(:));
nSites = numel(globalSites);
nStimuli = numel(G.Tall_IT);
nQuartets = nStimuli / 4;

snr = compute_snr_per_color_region(S.R, G.Tall_IT, globalSites);
baseline = double(snr.muSpont(:));
noiseSD = double(snr.sdSpont(:));
validNoise = isfinite(baseline) & isfinite(noiseSD) & noiseSD > 0;
[visualResponseScore, visuallyResponsive] = computeVisualResponseSelection( ...
    S.R, globalSites, baseline, noiseSD, opt.ResponseTopKQuartets, ...
    opt.ResponseSNRThreshold, nQuartets);

T = G.Tall_IT(1).T;
xPx = double(T.x_px(:));
yPx = double(T.y_px(:));
hasRFcenter = isfinite(xPx) & isfinite(yPx) & ...
    xPx >= 1 & xPx <= 1024 & yPx >= 1 & yPx <= 768;

fprintf('Building nearest yellow/purple curve distances for %d stimuli.\n', nStimuli);
[distanceYellowPx, distancePurplePx, geometryQC] = buildNearestDistances( ...
    G.ALLCOORDS, G.RTAB384, xPx, yPx, hasRFcenter, nStimuli);

[pairsA, pairsB] = complementaryPairs(nStimuli);
signedDistancePx = 0.5 * ((distancePurplePx(:, pairsA) - distanceYellowPx(:, pairsA)) - ...
    (distancePurplePx(:, pairsB) - distanceYellowPx(:, pairsB)));
meanNearestDistancePx = 0.5 * ( ...
    min(distanceYellowPx(:, pairsA), distancePurplePx(:, pairsA)) + ...
    min(distanceYellowPx(:, pairsB), distancePurplePx(:, pairsB)));

eligibleSite = validNoise & visuallyResponsive & hasRFcenter;
assert(any(eligibleSite), ...
    'No IT sites passed the RF-center, normalization, and response criteria.');
fprintf('IT nearest-curve color tuning (%s)\n',monkeyName);
fprintf('  RF centers available: %d / %d sites.\n', nnz(hasRFcenter), nSites);
fprintf('  Valid spontaneous normalization: %d / %d sites.\n', nnz(validNoise), nSites);
fprintf('  Visual criterion: top-%d quartet |SNR| > %.2f.\n', ...
    opt.ResponseTopKQuartets, opt.ResponseSNRThreshold);
fprintf('  Included sites: %d / %d.\n', nnz(eligibleSite), nSites);

early = fitWindow(S.R, globalSites, baseline, noiseSD, eligibleSite, ...
    signedDistancePx, meanNearestDistancePx, pairsA, pairsB, 2, opt);
late = fitWindow(S.R, globalSites, baseline, noiseSD, eligibleSite, ...
    signedDistancePx, meanNearestDistancePx, pairsA, pairsB, 3, opt);

OUT = struct();
OUT.description = ['IT color tuning from the minimum RF-center distance to ' ...
    'visible yellow versus purple capsule pixels.'];
OUT.monkey = monkeyName;
OUT.distanceDefinition = ['For each stimulus, bwdist is applied separately to ' ...
    'the rendered yellow and purple masks; values are sampled at the RF center.'];
OUT.predictorDefinition = ['Within complementary pair a/b: 0.5 * ' ...
    '[(dPurple-dYellow)_a - (dPurple-dYellow)_b]. Positive means yellow is ' ...
    'nearer in stimulus a than in b. Predictor is standardized within site.'];
OUT.responseDefinition = ['Complementary-pair difference of stimulus mean ' ...
    'responses normalized by spontaneous SD.'];
OUT.weightDefinition = ['Inverse estimated variance of the complementary-pair ' ...
    'response difference, using the balanced trial count min(nA,nB).'];
OUT.siteInclusion = ['Finite in-bounds RF center, valid spontaneous SD, and ' ...
    'historical top-K absolute quartet-response SNR criterion.'];
OUT.responseSNRThreshold = opt.ResponseSNRThreshold;
OUT.responseTopKQuartets = opt.ResponseTopKQuartets;
OUT.minPairs = opt.MinPairs;
OUT.alpha = opt.Alpha;
OUT.pixelsPerDegree = pixelsPerDegree;
OUT.calibrationFile = calibrationFile;
OUT.geometryFile = geometryFile;
OUT.responseFile = responseFile;
OUT.siteIndexLocal = find(eligibleSite);
OUT.siteIndexGlobal = globalSites(eligibleSite);
OUT.hasRFcenter = hasRFcenter;
OUT.visualResponseScoreAllSites = visualResponseScore;
OUT.distanceYellowPx = distanceYellowPx(eligibleSite, :);
OUT.distancePurplePx = distancePurplePx(eligibleSite, :);
OUT.signedDistancePx = signedDistancePx(eligibleSite, :);
OUT.signedDistanceDeg = OUT.signedDistancePx / pixelsPerDegree;
OUT.meanNearestDistancePx = meanNearestDistancePx(eligibleSite, :);
OUT.meanNearestDistanceDeg = OUT.meanNearestDistancePx / pixelsPerDegree;
OUT.pairsA = pairsA;
OUT.pairsB = pairsB;
OUT.geometryQC = geometryQC;
OUT.early = early;
OUT.late = late;

stem = sprintf('Color_tuning_IT_nearest_curve_respSNR%.2g_%s', ...
    opt.ResponseSNRThreshold, suffix);
stem = strrep(stem, '.', 'p');
OUT.resultFile = fullfile(cfg.resultsDir, [stem '.mat']);
OUT.figureFile = fullfile(cfg.resultsDir, [stem '.png']);

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end
if opt.SaveOutputs
    save(OUT.resultFile, 'OUT', '-v7.3');
    if ~isempty(fig)
        print(fig, OUT.figureFile, '-dpng', '-r180');
    end
    fprintf('Saved %s\n', OUT.resultFile);
    if ~isempty(fig)
        fprintf('Saved %s\n', OUT.figureFile);
    end
end
if ~isempty(fig)
    OUT.figure = fig;
end

printWindowSummary('Early (40-240 ms)', early, opt.Alpha);
printWindowSummary('Late (300-500 ms)', late, opt.Alpha);
end

function validateInputs(G, R)
requiredGeometry = {'Tall_IT','ALLCOORDS','RTAB384','RFrange'};
for i = 1:numel(requiredGeometry)
    assert(isfield(G, requiredGeometry{i}), ...
        'Geometry file lacks %s.', requiredGeometry{i});
end
assert(mod(numel(G.Tall_IT),8) == 0, 'Expected complete eight-stimulus blocks.');
assert(size(G.RTAB384, 1) == numel(G.Tall_IT), 'Expected matching stimulus-table rows.');
assert(mod(numel(G.Tall_IT), 8) == 0, ...
    'Stimulus count must be divisible by eight.');
assert(all(size(R.meanAct) == size(R.meanSqAct)), ...
    'R.meanAct and R.meanSqAct must have equal dimensions.');
assert(size(R.meanAct, 2) == numel(G.Tall_IT) && size(R.meanAct, 3) >= 3, ...
    'Expected matching stimuli and three windows.');
assert(max(G.RFrange) <= size(R.meanAct, 1), ...
    'R does not contain every IT site in RFrange.');
end

function [distanceY, distanceP, QC] = buildNearestDistances( ...
        ALLCOORDS, RTAB384, xPx, yPx, hasCenter, nStimuli)
nSites = numel(xPx);
distanceY = nan(nSites, nStimuli);
distanceP = nan(nSites, nStimuli);
sampleIndex = sub2ind([768 1024], round(yPx(hasCenter)), round(xPx(hasCenter)));

for stimulus = 1:nStimuli
    [~, masks] = render_stim_with_masks2(ALLCOORDS, RTAB384, stimulus, ...
        'ImageSize', [1024 768], 'DrawDots', false);
    mapY = bwdist(masks.yellowArm);
    mapP = bwdist(masks.purple);
    distanceY(hasCenter, stimulus) = mapY(sampleIndex);
    distanceP(hasCenter, stimulus) = mapP(sampleIndex);
end

[pairsA, pairsB] = complementaryPairs(nStimuli);
swapErrorY = abs(distanceY(:, pairsA) - distanceP(:, pairsB));
swapErrorP = abs(distanceP(:, pairsA) - distanceY(:, pairsB));
QC.maxComplementSwapErrorPx = max([swapErrorY(:); swapErrorP(:)], [], 'omitnan');
QC.meanComplementSwapErrorPx = mean([swapErrorY(:); swapErrorP(:)], 'omitnan');
QC.nSitesWithRFcenter = nnz(hasCenter);
QC.nSitesWithoutRFcenter = nnz(~hasCenter);
assert(QC.maxComplementSwapErrorPx < 1e-9, ...
    'Complementary yellow/purple masks do not swap exactly (max error %.3g px).', ...
    QC.maxComplementSwapErrorPx);
end

function W = fitWindow(R, globalSites, baseline, noiseSD, eligibleSite, ...
        signedDistancePx, meanNearestDistancePx, pairsA, pairsB, windowIndex, opt)
nSites = numel(globalSites);
nPairs = numel(pairsA);
nTrials = double(R.nTrials);

fields = {'beta','intercept','tStat','pValue','rSquaredPct','nPairs', ...
    'predictorMeanPx','predictorScalePx','effectivePairN','medianNearestDistancePx'};
for i = 1:numel(fields)
    W.(fields{i}) = nan(nSites, 1);
end
W.responsePairDifference = nan(nSites, nPairs);
W.predictorZ = nan(nSites, nPairs);
W.pairWeight = nan(nSites, nPairs);

for site = find(eligibleSite(:))'
    globalSite = globalSites(site);
    mu = squeeze(double(R.meanAct(globalSite, :, windowIndex)))';
    meanSq = squeeze(double(R.meanSqAct(globalSite, :, windowIndex)))';
    if isvector(nTrials)
        nTrialSite = double(nTrials(:));
    else
        nTrialSite = double(nTrials(globalSite, :))';
    end

    response = (mu - baseline(site)) / noiseSD(site);
    pairResponse = response(pairsA) - response(pairsB);
    nEff = min(nTrialSite(pairsA), nTrialSite(pairsB));
    sampleVar = nan(size(mu));
    enoughTrials = nTrialSite > 1 & isfinite(mu) & isfinite(meanSq);
    sampleVar(enoughTrials) = max(0, meanSq(enoughTrials) - mu(enoughTrials).^2) .* ...
        nTrialSite(enoughTrials) ./ (nTrialSite(enoughTrials) - 1);
    pairVariance = (sampleVar(pairsA) + sampleVar(pairsB)) ./ ...
        nEff / (noiseSD(site)^2);

    x = signedDistancePx(site, :)';
    nearestDistance = meanNearestDistancePx(site, :)';
    valid = isfinite(x) & isfinite(pairResponse) & isfinite(pairVariance) & ...
        pairVariance >= 0 & isfinite(nEff) & nEff > 1;
    if nnz(valid) < opt.MinPairs
        continue;
    end

    positiveVariance = pairVariance(valid & pairVariance > 0);
    if isempty(positiveVariance)
        varianceFloor = 1e-6;
    else
        varianceFloor = max(median(positiveVariance) * 1e-3, 1e-6);
    end
    weight = 1 ./ max(pairVariance(valid), varianceFloor);
    xv = x(valid);
    yv = pairResponse(valid);
    xMean = sum(weight .* xv) / sum(weight);
    xScale = sqrt(sum(weight .* (xv - xMean).^2) / sum(weight));
    if ~(isfinite(xScale) && xScale > 0)
        continue;
    end
    xz = (xv - xMean) / xScale;

    X = [ones(numel(xz), 1), xz];
    sqrtWeight = sqrt(weight);
    Xw = X .* sqrtWeight;
    yw = yv .* sqrtWeight;
    coefficient = Xw \ yw;
    fitted = X * coefficient;
    residual = yv - fitted;
    dof = numel(yv) - 2;
    rss = sum(weight .* residual.^2);
    weightedMean = sum(weight .* yv) / sum(weight);
    tss = sum(weight .* (yv - weightedMean).^2);
    covariance = (rss / dof) * inv(Xw' * Xw); %#ok<MINV>
    slopeSE = sqrt(max(covariance(2, 2), 0));
    if slopeSE > 0
        tStat = coefficient(2) / slopeSE;
        pValue = 2 * (1 - tcdf(abs(tStat), dof));
    elseif coefficient(2) ~= 0
        tStat = sign(coefficient(2)) * Inf;
        pValue = 0;
    else
        tStat = 0;
        pValue = 1;
    end

    W.beta(site) = coefficient(2);
    W.intercept(site) = coefficient(1);
    W.tStat(site) = tStat;
    W.pValue(site) = pValue;
    W.rSquaredPct(site) = 100 * (1 - rss / tss);
    W.nPairs(site) = nnz(valid);
    W.predictorMeanPx(site) = xMean;
    W.predictorScalePx(site) = xScale;
    W.effectivePairN(site) = sum(weight)^2 / sum(weight.^2);
    W.medianNearestDistancePx(site) = median(nearestDistance(valid), 'omitnan');
    W.responsePairDifference(site, valid) = yv;
    W.predictorZ(site, valid) = xz;
    W.pairWeight(site, valid) = weight;
end

W.includedSite = eligibleSite & isfinite(W.beta);
W.significant = W.includedSite & isfinite(W.pValue) & W.pValue < opt.Alpha;
W.significantFDR = false(nSites, 1);
W.significantFDR(W.includedSite) = benjaminiHochberg(W.pValue(W.includedSite), opt.Alpha);
W.yellowPreferring = W.includedSite & W.beta > 0;
W.purplePreferring = W.includedSite & W.beta < 0;
W.windowIndex = windowIndex;
W.windowMs = double(R.timeWindows(windowIndex, :));
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Name', 'IT nearest-curve color tuning', ...
    'NumberTitle', 'off', 'Tag', 'IT_nearest_curve_color_tuning', ...
    'Position', [80 80 1500 920]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
drawBinnedTuning(OUT.early, 'Early: 40-240 ms', true);
nexttile;
drawBinnedTuning(OUT.late, 'Late: 300-500 ms', false);

nexttile;
hold on;
edges = linspace(min([OUT.early.beta; OUT.late.beta], [], 'omitnan'), ...
    max([OUT.early.beta; OUT.late.beta], [], 'omitnan'), 25);
histogram(OUT.early.beta(OUT.early.includedSite), edges, ...
    'FaceColor', [0.95 0.72 0.18], 'FaceAlpha', 0.55, 'EdgeColor', 'none');
histogram(OUT.late.beta(OUT.late.includedSite), edges, ...
    'FaceColor', [0.48 0.24 0.62], 'FaceAlpha', 0.50, 'EdgeColor', 'none');
xline(0, 'k:');
xlabel('Nearest-color slope (response SD per predictor SD)');
ylabel('Sites');
title('Site tuning coefficients');
legend({'Early','Late'}, 'Location', 'best');
grid on;

nexttile;
valid = OUT.early.includedSite & OUT.late.includedSite;
scatter(OUT.early.beta(valid), OUT.late.beta(valid), 26, ...
    [0.25 0.45 0.70], 'filled', 'MarkerFaceAlpha', 0.65);
hold on;
limits = [min([OUT.early.beta(valid); OUT.late.beta(valid)]), ...
    max([OUT.early.beta(valid); OUT.late.beta(valid)])];
if limits(1) == limits(2)
    limits = limits + [-1 1];
end
plot(limits, limits, 'k--', 'LineWidth', 1);
xline(0, 'k:');
yline(0, 'k:');
xlim(limits);
ylim(limits);
axis square;
grid on;
xlabel('Early slope');
ylabel('Late slope');
if nnz(valid) >= 3
    r = corr(OUT.early.beta(valid), OUT.late.beta(valid), 'Rows', 'complete');
else
    r = NaN;
end
title(sprintf('Early versus late | N=%d, r=%.2f', nnz(valid), r));

sgtitle(sprintf(['IT color tuning from nearest curve | N=%d sites | ' ...
    'positive = yellow preference'], numel(OUT.siteIndexLocal)));
end

function drawBinnedTuning(W, panelTitle, showYLabel)
edges = [-inf -1.5 -0.75 -0.25 0.25 0.75 1.5 inf];
centers = [-2 -1.125 -0.5 0 0.5 1.125 2];
nSites = numel(W.beta);
siteBinMean = nan(nSites, numel(centers));
for site = find(W.includedSite(:))'
    x = W.predictorZ(site, :);
    y = W.responsePairDifference(site, :);
    for bin = 1:numel(centers)
        use = isfinite(x) & isfinite(y) & x > edges(bin) & x <= edges(bin + 1);
        if any(use)
            siteBinMean(site, bin) = mean(y(use));
        end
    end
end
mu = mean(siteBinMean, 1, 'omitnan');
n = sum(isfinite(siteBinMean), 1);
sem = std(siteBinMean, 0, 1, 'omitnan') ./ sqrt(n);
errorbar(centers, mu, sem, 'o-', 'Color', [0.18 0.18 0.18], ...
    'MarkerFaceColor', [0.95 0.72 0.18], 'LineWidth', 1.4, 'CapSize', 4);
hold on;
xline(0, 'k:');
yline(0, 'k:');
grid on;
xlabel('Yellow-nearer distance contrast (site SD)');
if showYLabel
    ylabel('Complementary-pair response difference (spontaneous SD)');
end
title(sprintf('%s | p<0.05: %d/%d', panelTitle, ...
    nnz(W.significant), nnz(W.includedSite)));
end

function printWindowSummary(label, W, alpha)
included = W.includedSite;
significant = W.significant;
fprintf('%s\n', label);
fprintf('  Fitted sites: %d. Median pairs/site: %.0f.\n', ...
    nnz(included), median(W.nPairs(included), 'omitnan'));
fprintf('  Median beta: %.3f response SD per predictor SD.\n', ...
    median(W.beta(included), 'omitnan'));
fprintf('  Uncorrected p < %.2f: %d (%d yellow, %d purple).\n', ...
    alpha, nnz(significant), nnz(significant & W.beta > 0), ...
    nnz(significant & W.beta < 0));
fprintf('  Benjamini-Hochberg FDR q < %.2f: %d.\n', ...
    alpha, nnz(W.significantFDR));
end

function significant = benjaminiHochberg(pValue, q)
pValue = pValue(:);
significant = false(size(pValue));
valid = isfinite(pValue);
[pSorted, order] = sort(pValue(valid));
n = numel(pSorted);
if n == 0
    return;
end
pass = find(pSorted <= (1:n)' * q / n, 1, 'last');
if isempty(pass)
    return;
end
validIndex = find(valid);
significant(validIndex(order(1:pass))) = true;
end

function [pairsA, pairsB] = complementaryPairs(nStimuli)
pairsA = find(mod((1:nStimuli) - 1, 8) < 4)';
pairsB = pairsA + 4;
assert(numel(pairsA) == nStimuli / 2, ...
    'Expected one complementary partner for every stimulus.');
end

function [score, responsive] = computeVisualResponseSelection( ...
        R, globalSites, baseline, noiseSD, topK, threshold, nQuartets)
nSites = numel(globalSites);
quartetEarly = nan(nSites, nQuartets);
quartetLate = nan(nSites, nQuartets);
nTrials = double(R.nTrials);

for site = 1:nSites
    globalSite = globalSites(site);
    if isvector(nTrials)
        nTrialSite = double(nTrials(:));
    else
        nTrialSite = double(nTrials(globalSite, :))';
    end
    early = squeeze(double(R.meanAct(globalSite, :, 2)))';
    late = squeeze(double(R.meanAct(globalSite, :, 3)))';
    for quartet = 1:nQuartets
        members = quartetMembers(quartet);
        quartetEarly(site, quartet) = weightedStimulusMean( ...
            early(members), nTrialSite(members));
        quartetLate(site, quartet) = weightedStimulusMean( ...
            late(members), nTrialSite(members));
    end
end

quartetEarly = (quartetEarly - baseline) ./ noiseSD;
quartetLate = (quartetLate - baseline) ./ noiseSD;
score = max([topKAbsoluteMean(quartetEarly, topK), ...
    topKAbsoluteMean(quartetLate, topK)], [], 2, 'omitnan');
responsive = isfinite(score) & score > threshold;
end

function members = quartetMembers(quartet)
block = floor((quartet - 1) / 2);
if mod(quartet, 2) == 1
    offsets = [1 2 5 6];
else
    offsets = [3 4 7 8];
end
members = 8 * block + offsets;
end

function mu = weightedStimulusMean(response, nTrials)
valid = isfinite(response) & isfinite(nTrials) & nTrials > 0;
if any(valid)
    mu = sum(response(valid) .* nTrials(valid)) / sum(nTrials(valid));
else
    mu = NaN;
end
end

function value = topKAbsoluteMean(responseByQuartet, topK)
nSites = size(responseByQuartet, 1);
value = nan(nSites, 1);
for site = 1:nSites
    x = abs(responseByQuartet(site, :));
    x = sort(x(isfinite(x)), 'descend');
    if ~isempty(x)
        value(site) = mean(x(1:min(topK, numel(x))));
    end
end
end
