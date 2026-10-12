function OUT = Color_direction_tuning_IT(varargin)
%COLOR_DIRECTION_TUNING_IT Test IT tuning for purple-to-yellow direction.
%
% Each visible color configuration occurs once with yellow and once with
% purple as target. Those two conditions are averaged first, so target and
% distractor status do not define the effect. For every IT site, a weighted
% regression then asks whether cos(theta) and sin(theta), where theta points
% from the purple capsule center to the yellow capsule center, explain the
% target-balanced response.

p = inputParser;
p.addParameter('Monkey', 'N', @(x) any(strcmpi(string(x),["N" "F"])));
p.addParameter('Config', struct(), @isstruct);
p.addParameter('ResponseSNRThreshold', 0.7, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ResponseTopKQuartets', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('MinPairs', 80, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 8 && x == floor(x));
p.addParameter('Alpha', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();
if ~isempty(fieldnames(opt.Config)); cfg = opt.Config; end
suffix = upper(char(opt.Monkey));
monkeyName = 'Mr Nilson'; if suffix=='F'; monkeyName='Figaro'; end
geometryFile = fullfile(cfg.matDir, ['Tall_IT_lines_' suffix '.mat']);
responseFile = fullfile(cfg.matDir, ['SNR_capsules_' suffix '_d12.mat']);
assert(isfile(geometryFile), 'Missing %s.', geometryFile);
assert(isfile(responseFile), 'Missing %s.', responseFile);

G = load(geometryFile, 'Tall_IT', 'ALLCOORDS', 'RTAB384', 'RFrange');
S = load(responseFile, 'R');
validateInputs(G, S.R);

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
eligibleSite = validNoise & visuallyResponsive;
assert(any(eligibleSite), 'No IT sites passed the response criterion.');

[configurationStimuli, complementaryConfigurations] = ...
    buildConfigurationMap(nStimuli);
[stimulusThetaDeg, configurationThetaDeg, geometryQC] = ...
    buildAndValidateColorGeometry(G.ALLCOORDS, G.RTAB384, ...
    configurationStimuli, complementaryConfigurations);

[configurationResponse, configurationVariance] = ...
    buildBalancedConfigurationResponses(S.R, globalSites, baseline, ...
    noiseSD, configurationStimuli);

early = fitDirectionWindow(configurationResponse(:, :, 2), ...
    configurationVariance(:, :, 2), configurationThetaDeg, ...
    complementaryConfigurations, eligibleSite, opt);
late = fitDirectionWindow(configurationResponse(:, :, 3), ...
    configurationVariance(:, :, 3), configurationThetaDeg, ...
    complementaryConfigurations, eligibleSite, opt);

OUT = struct();
OUT.description = ['IT tuning for the direction from the purple capsule ' ...
    'center to the yellow capsule center, after target-status balancing ' ...
    'and without occupancy correction.'];
OUT.monkey = monkeyName;
OUT.directionConvention = ['0 deg: yellow right of purple; 90 deg: yellow ' ...
    'above purple; 180 deg: yellow left; 270 deg: yellow below.'];
OUT.targetBalancing = ['Responses are averaged equally across two stimuli ' ...
    'with identical colored curves outside their shared overlap and opposite ' ...
    'target and overlap color.'];
OUT.model = ['Within each complementary configuration pair, response(A) - ' ...
    'response(B) = intercept + 2*betaCos*cos(color direction A) + ' ...
    '2*betaSin*sin(color direction A).'];
OUT.directionTest = ['Two-degree-of-freedom F-test for adding cosine and ' ...
    'sine to the intercept-only model.'];
OUT.responseSNRThreshold = opt.ResponseSNRThreshold;
OUT.responseTopKQuartets = opt.ResponseTopKQuartets;
OUT.minPairs = opt.MinPairs;
OUT.alpha = opt.Alpha;
OUT.siteIndexLocal = find(eligibleSite);
OUT.siteIndexGlobal = globalSites(eligibleSite);
OUT.visualResponseScoreAllSites = visualResponseScore;
OUT.configurationStimuli = configurationStimuli;
OUT.complementaryConfigurations = complementaryConfigurations;
OUT.stimulusThetaDeg = stimulusThetaDeg;
OUT.configurationThetaDeg = configurationThetaDeg;
OUT.configurationResponse = configurationResponse(eligibleSite, :, :);
OUT.configurationVariance = configurationVariance(eligibleSite, :, :);
OUT.geometryQC = geometryQC;
OUT.early = early;
OUT.late = late;

stem = sprintf('Color_direction_tuning_IT_respSNR%.2g_%s', ...
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

fprintf('\nIT color-direction tuning: %d / %d responsive sites included.\n', ...
    nnz(eligibleSite), nSites);
printWindowSummary('Early (40-240 ms)', early, opt.Alpha);
printWindowSummary('Late (300-500 ms)', late, opt.Alpha);
end

function validateInputs(G, R)
assert(mod(numel(G.Tall_IT),8) == 0 && size(G.RTAB384, 1) == numel(G.Tall_IT), ...
    'Expected matching geometry for complete eight-stimulus blocks.');
assert(all(size(R.meanAct) == size(R.meanSqAct)), ...
    'R.meanAct and R.meanSqAct dimensions differ.');
assert(size(R.meanAct, 2) == numel(G.Tall_IT) && size(R.meanAct, 3) >= 3, ...
    'Expected matching stimuli and three response windows.');
end

function [configurationStimuli, complementary] = buildConfigurationMap(nStimuli)
assert(mod(nStimuli, 8) == 0, 'Stimulus count must be divisible by eight.');
nBlocks = nStimuli / 8;
configurationStimuli = zeros(4 * nBlocks, 2);
complementary = zeros(2 * nBlocks, 2);
for block = 1:nBlocks
    baseStimulus = 8 * (block - 1);
    rows = 4 * (block - 1) + (1:4);
    configurationStimuli(rows, :) = baseStimulus + ...
        [1 6; 2 5; 3 8; 4 7];
    pairRows = 2 * (block - 1) + (1:2);
    complementary(pairRows, :) = rows([1 2; 3 4]);
end
assert(numel(unique(configurationStimuli(:))) == nStimuli, ...
    'Every stimulus must occur in exactly one target-balanced pair.');
end

function [thetaStimulus, thetaConfiguration, QC] = ...
        buildAndValidateColorGeometry(ALLCOORDS, RTAB384, ...
        configurationStimuli, complementary)
nStimuli = size(RTAB384, 1);
thetaStimulus = nan(nStimuli, 1);
for stimulus = 1:nStimuli
    field = sprintf('stim_%d', stimulus);
    coords = ALLCOORDS.(field);
    junction = double(coords.s(:))';
    centerFig = 0.5 * (junction + double(coords.t_fig(:))');
    centerBack = 0.5 * (junction + double(coords.t_back(:))');
    colorIndex = double(RTAB384(stimulus, 8));
    if colorIndex == 1
        centerPurple = centerFig;
        centerYellow = centerBack;
    elseif colorIndex == 2
        centerYellow = centerFig;
        centerPurple = centerBack;
    else
        error('Unexpected color index %g for stimulus %d.', ...
            colorIndex, stimulus);
    end
    difference = centerYellow - centerPurple;
    thetaStimulus(stimulus) = mod(atan2d(difference(2), difference(1)), 360);
end

nConfigurations = size(configurationStimuli, 1);
thetaConfiguration = nan(nConfigurations, 1);
maxDuplicateAngleError = 0;
maxYellowMaskMismatch = 0;
maxPurpleMaskMismatch = 0;
maxYellowMismatchOutsideOverlap = 0;
maxPurpleMismatchOutsideOverlap = 0;
maxUnionMaskMismatch = 0;
for configuration = 1:nConfigurations
    stimuli = configurationStimuli(configuration, :);
    angleError = abs(angleDifference( ...
        thetaStimulus(stimuli(1)), thetaStimulus(stimuli(2))));
    maxDuplicateAngleError = max(maxDuplicateAngleError, angleError);
    thetaConfiguration(configuration) = thetaStimulus(stimuli(1));

    [~, masksA] = render_stim_with_masks2( ...
        ALLCOORDS, RTAB384, stimuli(1), 'DrawDots', false);
    [~, masksB] = render_stim_with_masks2( ...
        ALLCOORDS, RTAB384, stimuli(2), 'DrawDots', false);
    maxYellowMaskMismatch = max(maxYellowMaskMismatch, ...
        nnz(xor(masksA.yellowArm, masksB.yellowArm)));
    maxPurpleMaskMismatch = max(maxPurpleMaskMismatch, ...
        nnz(xor(masksA.purple, masksB.purple)));
    sharedOverlap = (masksA.figArm & masksA.backArm) | ...
        (masksB.figArm & masksB.backArm);
    outsideOverlap = ~sharedOverlap;
    maxYellowMismatchOutsideOverlap = max( ...
        maxYellowMismatchOutsideOverlap, nnz(outsideOverlap & ...
        xor(masksA.yellowArm, masksB.yellowArm)));
    maxPurpleMismatchOutsideOverlap = max( ...
        maxPurpleMismatchOutsideOverlap, nnz(outsideOverlap & ...
        xor(masksA.purple, masksB.purple)));
    unionA = masksA.figArm | masksA.backArm;
    unionB = masksB.figArm | masksB.backArm;
    maxUnionMaskMismatch = max(maxUnionMaskMismatch, nnz(xor(unionA, unionB)));
end

oppositeError = abs(abs(angleDifference( ...
    thetaConfiguration(complementary(:, 1)), ...
    thetaConfiguration(complementary(:, 2)))) - 180);

QC = struct();
QC.nConfigurations = nConfigurations;
QC.maxDuplicateAngleErrorDeg = maxDuplicateAngleError;
QC.maxComplementDirectionErrorDeg = max(oppositeError);
QC.maxDuplicateYellowMaskMismatchPixels = maxYellowMaskMismatch;
QC.maxDuplicatePurpleMaskMismatchPixels = maxPurpleMaskMismatch;
QC.maxYellowMismatchOutsideOverlapPixels = maxYellowMismatchOutsideOverlap;
QC.maxPurpleMismatchOutsideOverlapPixels = maxPurpleMismatchOutsideOverlap;
QC.maxUnionMaskMismatchPixels = maxUnionMaskMismatch;
assert(maxDuplicateAngleError < 1e-9, ...
    'Target-balanced duplicate directions differ by %.3g deg.', ...
    maxDuplicateAngleError);
assert(QC.maxComplementDirectionErrorDeg < 1e-9, ...
    'Complementary color directions are not exactly opposite.');
assert(maxYellowMismatchOutsideOverlap == 0 && ...
    maxPurpleMismatchOutsideOverlap == 0, ...
    'Target-balanced color curves differ outside their shared overlap.');
assert(maxUnionMaskMismatch == 0, ...
    'Target-balanced stimulus pairs have different total curve geometry.');
fprintf(['Color-geometry validation: %d target-balanced configurations; ' ...
    'colored curves identical outside overlap; complement angle error ' ...
    '%.3g deg.\n'], nConfigurations, QC.maxComplementDirectionErrorDeg);
end

function [response, variance] = buildBalancedConfigurationResponses( ...
        R, globalSites, baseline, noiseSD, configurationStimuli)
nSites = numel(globalSites);
nConfigurations = size(configurationStimuli, 1);
nWindows = size(R.meanAct, 3);
response = nan(nSites, nConfigurations, nWindows);
variance = nan(nSites, nConfigurations, nWindows);
nTrials = double(R.nTrials);

for site = 1:nSites
    globalSite = globalSites(site);
    if isvector(nTrials)
        nTrialSite = double(nTrials(:));
    else
        nTrialSite = double(nTrials(globalSite, :))';
    end
    for window = 1:nWindows
        mu = squeeze(double(R.meanAct(globalSite, :, window)))';
        meanSq = squeeze(double(R.meanSqAct(globalSite, :, window)))';
        varianceOfMean = nan(size(mu));
        valid = nTrialSite > 1 & isfinite(mu) & isfinite(meanSq);
        varianceOfMean(valid) = max(0, meanSq(valid) - mu(valid).^2) ./ ...
            (nTrialSite(valid) - 1);
        normalizedMu = (mu - baseline(site)) / noiseSD(site);
        normalizedVariance = varianceOfMean / (noiseSD(site)^2);
        for configuration = 1:nConfigurations
            stimuli = configurationStimuli(configuration, :);
            % Equal weighting ensures target color is balanced even if trial
            % counts differ between the two otherwise identical stimuli.
            response(site, configuration, window) = ...
                0.5 * sum(normalizedMu(stimuli));
            variance(site, configuration, window) = ...
                0.25 * sum(normalizedVariance(stimuli));
        end
    end
end
end

function W = fitDirectionWindow(response, responseVariance, thetaDeg, ...
        complementary, eligibleSite, opt)
nSites = size(response, 1);
fields = {'betaCos','betaSin','directionAmplitude','preferredDirectionDeg', ...
    'pValue','fStat','partialRSquaredPct','deltaRSquaredPct','fullRSquaredPct', ...
    'nPairs','effectivePairN'};
for field = fields
    W.(field{1}) = nan(nSites, 1);
end
pairA = complementary(:, 1);
pairB = complementary(:, 2);
pairTheta = thetaDeg(pairA);
cosTheta = cosd(pairTheta(:));
sinTheta = sind(pairTheta(:));

for site = find(eligibleSite(:))'
    responseA = double(response(site, pairA))';
    responseB = double(response(site, pairB))';
    y = responseA - responseB;
    variance = double(responseVariance(site, pairA) + ...
        responseVariance(site, pairB))';
    valid = isfinite(y) & isfinite(variance) & variance >= 0 & ...
        isfinite(cosTheta) & isfinite(sinTheta);
    if nnz(valid) < opt.MinPairs
        continue;
    end

    positiveVariance = variance(valid & variance > 0);
    if isempty(positiveVariance)
        varianceFloor = 1e-6;
    else
        varianceFloor = max(median(positiveVariance) * 1e-3, 1e-6);
    end
    weight = 1 ./ max(variance(valid), varianceFloor);
    yv = y(valid);
    direction = [cosTheta(valid), sinTheta(valid)];
    reducedDesign = ones(nnz(valid), 1);
    fullDesign = [reducedDesign, direction];
    fit = nestedWeightedFit(yv, weight, reducedDesign, fullDesign, 2);

    % The complementary response difference contains twice the amplitude
    % of the underlying single-configuration first harmonic.
    betaCos = 0.5 * fit.fullCoefficient(end - 1);
    betaSin = 0.5 * fit.fullCoefficient(end);
    W.betaCos(site) = betaCos;
    W.betaSin(site) = betaSin;
    W.directionAmplitude(site) = hypot(betaCos, betaSin);
    W.preferredDirectionDeg(site) = mod(atan2d(betaSin, betaCos), 360);
    W.pValue(site) = fit.pValue;
    W.fStat(site) = fit.fStat;
    W.partialRSquaredPct(site) = fit.partialRSquaredPct;
    W.deltaRSquaredPct(site) = fit.deltaRSquaredPct;
    W.fullRSquaredPct(site) = fit.fullRSquaredPct;
    W.nPairs(site) = nnz(valid);
    W.effectivePairN(site) = sum(weight)^2 / sum(weight.^2);
end

W.includedSite = eligibleSite & isfinite(W.pValue);
W.significant = W.includedSite & W.pValue < opt.Alpha;
W.significantFDR = false(nSites, 1);
W.significantFDR(W.includedSite) = benjaminiHochberg( ...
    W.pValue(W.includedSite), opt.Alpha);
W.directionConvention = 'purple-to-yellow';
end

function fit = nestedWeightedFit(y, weight, reducedDesign, fullDesign, addedDF)
sqrtWeight = sqrt(weight);
reducedCoefficient = (reducedDesign .* sqrtWeight) \ (y .* sqrtWeight);
fullCoefficient = (fullDesign .* sqrtWeight) \ (y .* sqrtWeight);
reducedResidual = y - reducedDesign * reducedCoefficient;
fullResidual = y - fullDesign * fullCoefficient;
rssReduced = sum(weight .* reducedResidual.^2);
rssFull = sum(weight .* fullResidual.^2);
weightedMean = sum(weight .* y) / sum(weight);
rssNull = sum(weight .* (y - weightedMean).^2);
dfError = numel(y) - size(fullDesign, 2);

improvement = max(0, rssReduced - rssFull);
if dfError > 0 && rssFull > 0
    fStat = (improvement / addedDF) / (rssFull / dfError);
    pValue = 1 - fcdf(fStat, addedDF, dfError);
else
    fStat = NaN;
    pValue = NaN;
end

fit = struct();
fit.reducedCoefficient = reducedCoefficient;
fit.fullCoefficient = fullCoefficient;
fit.rssReduced = rssReduced;
fit.rssFull = rssFull;
fit.fStat = fStat;
fit.pValue = pValue;
fit.partialRSquaredPct = 100 * improvement / rssReduced;
fit.deltaRSquaredPct = 100 * improvement / rssNull;
fit.fullRSquaredPct = 100 * (1 - rssFull / rssNull);
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Position', [80 80 1450 900], ...
    'Name', 'IT color-direction tuning', 'NumberTitle', 'off');

earlyPref = OUT.early.preferredDirectionDeg(OUT.early.significant);
latePref = OUT.late.preferredDirectionDeg(OUT.late.significant);
axEarly = polaraxes('Parent', fig, 'Position', [0.07 0.55 0.38 0.34]);
polarhistogram(axEarly, deg2rad(earlyPref), 12, ...
    'FaceColor', [0.94 0.66 0.16], 'FaceAlpha', 0.8);
axEarly.ThetaZeroLocation = 'right';
axEarly.ThetaDir = 'counterclockwise';
title(axEarly, sprintf('Early 40-240 ms | p<0.05: %d, FDR: %d', ...
    nnz(OUT.early.significant), nnz(OUT.early.significantFDR)));

axLate = polaraxes('Parent', fig, 'Position', [0.55 0.55 0.38 0.34]);
polarhistogram(axLate, deg2rad(latePref), 12, ...
    'FaceColor', [0.48 0.31 0.68], 'FaceAlpha', 0.8);
axLate.ThetaZeroLocation = 'right';
axLate.ThetaDir = 'counterclockwise';
title(axLate, sprintf('Late 300-500 ms | p<0.05: %d, FDR: %d', ...
    nnz(OUT.late.significant), nnz(OUT.late.significantFDR)));

axEffect = axes('Parent', fig, 'Position', [0.09 0.10 0.37 0.32]);
hold(axEffect, 'on');
edges = linspace(0, max([OUT.early.partialRSquaredPct; ...
    OUT.late.partialRSquaredPct], [], 'omitnan'), 25);
histogram(axEffect, OUT.early.partialRSquaredPct(OUT.early.includedSite), ...
    edges, 'FaceColor', [0.94 0.66 0.16], 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none');
histogram(axEffect, OUT.late.partialRSquaredPct(OUT.late.includedSite), ...
    edges, 'FaceColor', [0.48 0.31 0.68], 'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');
xlabel(axEffect, 'Variance explained by color direction (%)');
ylabel(axEffect, 'Sites');
title(axEffect, 'Color-direction effect size');
legend(axEffect, {'Early','Late'}, 'Location', 'best');
formatAxes(axEffect);

axCount = axes('Parent', fig, 'Position', [0.58 0.10 0.32 0.32]);
counts = [nnz(OUT.early.significant), nnz(OUT.early.significantFDR); ...
    nnz(OUT.late.significant), nnz(OUT.late.significantFDR)];
bar(axCount, counts);
set(axCount, 'XTickLabel', {'Early','Late'});
ylabel(axCount, 'Direction-selective sites');
legend(axCount, {'p < 0.05','FDR q < 0.05'}, 'Location', 'best');
title(axCount, sprintf('Color-direction tests | N=%d', ...
    numel(OUT.siteIndexLocal)));
formatAxes(axCount);

sgtitle(fig, ['IT color-direction tuning | angle points from purple to yellow ' ...
    '| 0 deg right, 90 deg up'], 'FontWeight', 'bold');
end

function printWindowSummary(label, W, alpha)
included = W.includedSite;
significant = W.significant;
fprintf('%s\n', label);
fprintf('  Fitted sites: %d; median complementary pairs/site: %.0f.\n', ...
    nnz(included), median(W.nPairs(included), 'omitnan'));
fprintf(['  Color direction: p < %.2f for %d sites; ' ...
    'FDR q < %.2f for %d sites.\n'], alpha, nnz(significant), ...
    alpha, nnz(W.significantFDR));
fprintf('  Median variance explained: %.2f%% (significant sites %.2f%%).\n', ...
    median(W.partialRSquaredPct(included), 'omitnan'), ...
    median(W.partialRSquaredPct(significant), 'omitnan'));
if any(significant)
    theta = W.preferredDirectionDeg(significant);
    resultant = mean(exp(1i * deg2rad(theta)));
    fprintf('  Significant-site circular mean: %.1f deg; resultant length %.2f.\n', ...
        mod(rad2deg(angle(resultant)), 360), abs(resultant));
end
end

function difference = angleDifference(a, b)
difference = mod(a - b + 180, 360) - 180;
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

function formatAxes(ax)
grid(ax, 'on');
box(ax, 'off');
set(ax, 'FontSize', 10, 'Layer', 'top');
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
