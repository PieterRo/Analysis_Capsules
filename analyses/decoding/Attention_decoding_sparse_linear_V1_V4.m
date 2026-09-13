function OUT = Attention_decoding_sparse_linear_V1_V4(varargin)
%ATTENTION_DECODING_SPARSE_LINEAR_V1_V4 Compare sparse attention decoders.
%
% Within each quartet, the first matched stimulus pair is class +1 and the
% complementary pair is class -1. Site polarity is fixed to that arbitrary
% orientation: target sites in the +1 pair are positive and distractor
% sites are negative. Responses remain centered per site and quartet.

p = inputParser;
p.addParameter('ElasticNetAlpha', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('NumLambda', 60, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 10 && x == floor(x));
p.addParameter('LambdaRatio', 1e-4, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('SvmLambdaRange', [1e-6 1], ...
    @(x) isnumeric(x) && isequal(size(x), [1 2]) && ...
    all(isfinite(x)) && all(x > 0) && x(1) < x(2));
p.addParameter('CoefficientTolerance', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ChunkTrials', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();

fprintf('Recomputing the matched full-Fisher reference decoders.\n');
fisherV1 = Attention_decoding_V1_proof( ...
    'ResponseCentering', 'quartet', 'CurveMarginDeg', 0, ...
    'MinSitesPerQuartet', 20, 'SaveOutputs', false, 'MakeFigure', false);
fisherV4 = Attention_decoding_V4_proof( ...
    'ResponseCentering', 'quartet', 'MinSitesPerQuartet', 20, ...
    'SaveOutputs', false, 'MakeFigure', false);

OUT = struct();
OUT.description = ['In-sample comparison of quartet-centered full-Fisher, ' ...
    'logistic lasso, elastic-net logistic, and L1-SVM attention decoders.'];
OUT.monkey = 'Mr Nilson';
OUT.classDefinition = ['Within each quartet, the first geometry-matched ' ...
    'stimulus pair is +1 and the complementary pair is -1.'];
OUT.targetAlignedScoreDefinition = ...
    'Class decision score multiplied by the trial class; values above zero are correct.';
OUT.responseCentering = 'quartet';
OUT.responseWindow = [300 500];
OUT.curveMarginDeg = 0;
OUT.minSitesPerQuartet = 20;
OUT.elasticNetAlpha = opt.ElasticNetAlpha;
OUT.noCrossValidation = true;
OUT.selectionRule = ['Maximum in-sample accuracy; ties use the fewest ' ...
    'nonzero coefficients and then the largest lambda.'];

OUT.V1 = analyzeRegion('V1', fisherV1, cfg, opt);
OUT.V4 = analyzeRegion('V4', fisherV4, cfg, opt);

methodNames = {'Full Fisher', 'Logistic lasso', ...
    sprintf('Elastic net (alpha = %.2g)', opt.ElasticNetAlpha), 'L1 SVM'};
OUT.methodNames = methodNames;
OUT.accuracy = [OUT.V1.accuracy; OUT.V4.accuracy];
OUT.nNonzero = [OUT.V1.nNonzero; OUT.V4.nNonzero];

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

resultFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.png');
OUT.resultFile = resultFile;
OUT.figureFile = figureFile;

if opt.SaveOutputs
    save(resultFile, 'OUT', '-v7.3');
    if ~isempty(fig)
        print(fig, figureFile, '-dpng', '-r180');
    end
    fprintf('Saved %s\n', resultFile);
    if ~isempty(fig)
        fprintf('Saved %s\n', figureFile);
    end
end

if ~isempty(fig)
    OUT.figure = fig;
end

printSummary(OUT);

end

function R = analyzeRegion(region, fisher, cfg, opt)
fprintf('\nBuilding %s quartet-pair feature matrix.\n', region);
D = buildFeatures(region, fisher, cfg, opt.ChunkTrials);

weights = classBalancingWeights(D.classLabel);
lasso = fitLogisticPath(D.featureMatrix, D.classLabel, weights, 1, opt);
elastic = fitLogisticPath(D.featureMatrix, D.classLabel, weights, ...
    opt.ElasticNetAlpha, opt);
svm = fitSvmPath(D.featureMatrix, D.classLabel, weights, opt);

fisherClassScore = D.classLabel .* double(fisher.SFullFisher(:));
fisherAccuracy = mean(signClass(fisherClassScore) == D.classLabel);

R = D;
R.fullFisher = struct();
R.fullFisher.classScore = fisherClassScore;
R.fullFisher.targetAlignedScore = D.classLabel .* fisherClassScore;
R.fullFisher.accuracy = fisherAccuracy;
R.fullFisher.nNonzero = size(D.featureMatrix, 2);
R.fullFisher.note = ['Fisher uses stimulus-specific covariance weights; ' ...
    'nNonzero is the candidate site pool, not one global coefficient vector.'];
R.logisticLasso = addTargetAlignedScores(lasso, D);
R.elasticNetLogistic = addTargetAlignedScores(elastic, D);
R.l1Svm = addTargetAlignedScores(svm, D);
R.accuracy = [R.fullFisher.accuracy, R.logisticLasso.accuracy, ...
    R.elasticNetLogistic.accuracy, R.l1Svm.accuracy];
R.nNonzero = [R.fullFisher.nNonzero, R.logisticLasso.nNonzero, ...
    R.elasticNetLogistic.nNonzero, R.l1Svm.nNonzero];

fprintf('%s: %d trials, %d quartets, %d candidate sites.\n', ...
    region, numel(D.classLabel), numel(unique(D.quartet)), ...
    size(D.featureMatrix, 2));
end

function D = buildFeatures(region, fisher, cfg, chunkTrials)
dataDir = fullfile(cfg.dataRoot, 'Mr Nilson');
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
tb = double(m2.tb);
tb = tb(:)';

switch region
    case 'V1'
        globalSites = (1:512)';
        geom = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'), 'Tall_V1');
        Tall = geom.Tall_V1;
        snrData = load(fullfile(cfg.matDir, ...
            'SNR_V1_byColor_byWindow.mat'), 'SNR');
        baseline = double(snrData.SNR.muSpont(:));
    case 'V4'
        globalSites = (513:768)';
        geom = load(fullfile(cfg.matDir, 'Tall_V4_lines_N.mat'), 'Tall_V4');
        responseData = load(fullfile(cfg.matDir, ...
            'SNR_capsules_N_d12.mat'), 'R');
        Tall = geom.Tall_V4;
        snr = compute_snr_per_color_region(responseData.R, Tall, globalSites);
        baseline = double(snr.muSpont(:));
    otherwise
        error('Unknown region %s.', region);
end

assert(strcmp(fisher.responseCentering, 'quartet'), ...
    '%s Fisher reference is not quartet-centered.', region);
assert(isequal(double(fisher.responseWindowRequested), [300 500]), ...
    '%s Fisher reference uses an unexpected response window.', region);
assert(fisher.minSitesPerQuartet == 20, ...
    '%s Fisher reference uses an unexpected quartet threshold.', region);

trialIndex = double(fisher.trialIndex(:));
stimulus = double(fisher.stimulus(:));
quartet = double(fisher.quartet(:));
classLabel = ones(numel(stimulus), 1);
classLabel(mod(stimulus, 2) == 0) = -1;

eligibleMask = logical(fisher.eligibleSiteMask(:));
eligibleSite = find(eligibleMask);
responseScale = double(fisher.responseScale(:));
quartetMean = double(fisher.quartetMeanResponse);
timeMask = tb >= fisher.responseWindowRequested(1) & ...
    tb <= fisher.responseWindowRequested(2);
timeSamples = find(timeMask);
assert(~isempty(timeSamples), '%s response window does not overlap tb.', region);

[canonicalSign, activeByStimulus] = buildCanonicalSigns(Tall, eligibleMask);
validateQuartetPairing(canonicalSign, activeByStimulus, ...
    fisher.includedQuartet, region);

nTrials = numel(trialIndex);
nEligible = numel(eligibleSite);
featureMatrix = zeros(nTrials, nEligible);
availability = false(nTrials, nEligible);
normalizedResponseByTrial = nan(nTrials, nEligible);
canonicalSignByTrial = zeros(nTrials, nEligible);
[~, nAllTrials, ~] = size(m1, 'normMUA');
rowByGlobalTrial = zeros(nAllTrials, 1);
rowByGlobalTrial(trialIndex) = 1:nTrials;

for firstGlobalTrial = 1:chunkTrials:nAllTrials
    lastGlobalTrial = min(nAllTrials, firstGlobalTrial + chunkTrials - 1);
    globalRange = firstGlobalTrial:lastGlobalTrial;
    selected = rowByGlobalTrial(globalRange) > 0;
    if ~any(selected)
        continue;
    end
    rows = rowByGlobalTrial(globalRange(selected));
    raw = double(m1.normMUA(globalSites, globalRange, timeSamples));
    trialResponse = reshape(mean(raw, 3, 'omitnan'), ...
        numel(globalSites), numel(globalRange));
    trialResponse = trialResponse(:, selected);
    normalized = (trialResponse - baseline) ./ responseScale;
    centered = normalized - quartetMean(:, quartet(rows));

    stim = stimulus(rows);
    signs = canonicalSign(:, stim);
    active = activeByStimulus(:, stim) & isfinite(centered);
    features = centered .* signs;
    features(~active) = 0;
    featureMatrix(rows, :) = features(eligibleSite, :)';
    availability(rows, :) = active(eligibleSite, :)';
    normalizedResponseByTrial(rows, :) = normalized(eligibleSite, :)';
    canonicalSignByTrial(rows, :) = signs(eligibleSite, :)';
end

nSitesUsed = sum(availability, 2);
assert(isequal(nSitesUsed, double(fisher.nSitesUsed(:))), ...
    '%s sparse features do not reproduce the Fisher site counts.', region);

variablePredictor = any(availability, 1) & ...
    std(featureMatrix, 0, 1) > 0;
featureMatrix = featureMatrix(:, variablePredictor);
availability = availability(:, variablePredictor);
normalizedResponseByTrial = normalizedResponseByTrial(:, variablePredictor);
canonicalSignByTrial = canonicalSignByTrial(:, variablePredictor);
eligibleSite = eligibleSite(variablePredictor);

assert(all(isfinite(featureMatrix(:))), ...
    '%s feature matrix contains non-finite values.', region);
assert(all(ismember(classLabel, [-1 1])), ...
    '%s class labels must be -1 or +1.', region);

D = struct();
D.region = region;
D.trialIndex = trialIndex;
D.stimulus = stimulus;
D.quartet = quartet;
D.classLabel = classLabel;
D.positiveClassStimulusParity = 'odd stimulus number within each 8-stimulus block';
D.siteIndexLocal = eligibleSite;
D.siteIndexGlobal = globalSites(eligibleSite);
D.featureMatrix = featureMatrix;
D.availability = availability;
D.normalizedResponse = normalizedResponseByTrial;
D.canonicalSign = canonicalSignByTrial;
D.nSitesUsed = sum(availability, 2);
D.quartetMeanResponse = quartetMean;
D.responseScale = responseScale;
D.baseline = baseline;
D.responseWindowSampled = [tb(timeSamples(1)), tb(timeSamples(end))];
D.pairingValidation = 'All non-overlap site polarities matched within stimulus pairs.';
D.featureDefinition = ['Quartet-centered normalized response multiplied by ' ...
    'the site polarity of the arbitrarily chosen positive stimulus pair; ' ...
    'background and overlap entries are zero.'];
D.ALLMATColumnsUsedViaFisherReference = [1 9 11];
D.ALLMATNote = ['The Fisher reference supplies ALLMAT-based trial/stimulus, ' ...
    'correctness, and day selection; site assignments come from Tall.'];
end

function [canonicalSign, active] = buildCanonicalSigns(Tall, eligibleMask)
nSites = numel(eligibleMask);
nStimuli = numel(Tall);
canonicalSign = zeros(nSites, nStimuli);
active = false(nSites, nStimuli);

for stim = 1:nStimuli
    T = Tall(stim).T;
    assignment = string(T.assignment(1:nSites));
    role = zeros(nSites, 1);
    role(assignment == "target") = 1;
    role(assignment == "distractor") = -1;
    if ismember('overlap', T.Properties.VariableNames)
        overlap = logical(T.overlap(1:nSites));
    else
        overlap = false(nSites, 1);
    end
    active(:, stim) = eligibleMask & role ~= 0 & ~overlap;
    classSign = 1;
    if mod(stim, 2) == 0
        classSign = -1;
    end
    canonicalSign(:, stim) = classSign * role;
    canonicalSign(~active(:, stim), stim) = 0;
end
end

function validateQuartetPairing(signByStimulus, activeByStimulus, ...
    includedQuartet, region)
nMismatch = 0;
nCompared = 0;
for quartet = includedQuartet(:)'
    members = quartetMembers(quartet);
    reference = signByStimulus(:, members(1));
    for memberIdx = 2:numel(members)
        comparison = signByStimulus(:, members(memberIdx));
        compare = activeByStimulus(:, members(1)) & ...
            activeByStimulus(:, members(memberIdx));
        nMismatch = nMismatch + nnz(reference(compare) ~= comparison(compare));
        nCompared = nCompared + nnz(compare);
    end
end
assert(nMismatch == 0, ...
    '%s has %d inconsistent site polarities across quartet pairs.', ...
    region, nMismatch);
fprintf('%s pairing validation: %d site-stimulus comparisons, 0 mismatches.\n', ...
    region, nCompared);
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

function weights = classBalancingWeights(classLabel)
n = numel(classLabel);
nPositive = nnz(classLabel == 1);
nNegative = nnz(classLabel == -1);
assert(nPositive > 0 && nNegative > 0, 'Both classes must contain trials.');
weights = zeros(n, 1);
weights(classLabel == 1) = n / (2 * nPositive);
weights(classLabel == -1) = n / (2 * nNegative);
end

function model = fitLogisticPath(X, classLabel, weights, alpha, opt)
[B, stats] = lassoglm(X, classLabel == 1, 'binomial', ...
    'Alpha', alpha, 'NumLambda', opt.NumLambda, ...
    'LambdaRatio', opt.LambdaRatio, 'Standardize', true, ...
    'Weights', weights, 'CV', 'resubstitution');
score = X * B + stats.Intercept;
accuracy = mean(signClass(score) == classLabel, 1);
nNonzero = countNonzero(B, opt.CoefficientTolerance);
bestIndex = chooseBestModel(accuracy, nNonzero, stats.Lambda);

model = packageModel(B, stats.Intercept, stats.Lambda, accuracy, ...
    nNonzero, bestIndex, X, 'logistic', alpha);
model.deviancePath = stats.Deviance;
end

function model = fitSvmPath(X, classLabel, weights, opt)
lambda = logspace(log10(opt.SvmLambdaRange(1)), ...
    log10(opt.SvmLambdaRange(2)), opt.NumLambda);
predictorMean = mean(X, 1);
predictorScale = std(X, 0, 1);
predictorScale(~isfinite(predictorScale) | predictorScale <= 0) = 1;
standardizedX = (X - predictorMean) ./ predictorScale;
fitted = fitclinear(standardizedX, classLabel, 'Learner', 'svm', ...
    'Regularization', 'lasso', 'Solver', 'sparsa', ...
    'Lambda', lambda, 'Weights', weights);
B = double(fitted.Beta) ./ predictorScale(:);
bias = double(fitted.Bias) - predictorMean * B;
score = X * B + bias;

for modelIdx = 1:size(score, 2)
    accuracyForward = mean(signClass(score(:, modelIdx)) == classLabel);
    accuracyReverse = mean(signClass(-score(:, modelIdx)) == classLabel);
    if accuracyReverse > accuracyForward
        B(:, modelIdx) = -B(:, modelIdx);
        bias(modelIdx) = -bias(modelIdx);
        score(:, modelIdx) = -score(:, modelIdx);
    end
end

accuracy = mean(signClass(score) == classLabel, 1);
nNonzero = countNonzero(B, opt.CoefficientTolerance);
bestIndex = chooseBestModel(accuracy, nNonzero, lambda);
model = packageModel(B, bias, lambda, accuracy, nNonzero, ...
    bestIndex, X, 'svm', 1);
end

function model = packageModel(B, bias, lambda, accuracy, nNonzero, ...
    bestIndex, X, learner, alpha)
coefficient = B(:, bestIndex);
intercept = bias(bestIndex);
model = struct();
model.learner = learner;
model.alpha = alpha;
model.lambda = lambda(bestIndex);
model.coefficient = coefficient;
model.intercept = intercept;
model.classScore = X * coefficient + intercept;
model.accuracy = accuracy(bestIndex);
model.nNonzero = nNonzero(bestIndex);
model.bestPathIndex = bestIndex;
model.lambdaPath = lambda;
model.accuracyPath = accuracy;
model.nNonzeroPath = nNonzero;
model.coefficientPath = B;
model.interceptPath = bias;
end

function model = addTargetAlignedScores(model, D)
model.targetAlignedScore = D.classLabel .* model.classScore;
denominator = D.availability * abs(model.coefficient);
weightedResponse = D.featureMatrix * model.coefficient;
model.classScoreNormalized = nan(size(model.classScore));
hasWeight = denominator > 0;
model.classScoreNormalized(hasWeight) = ...
    weightedResponse(hasWeight) ./ denominator(hasWeight);
model.targetAlignedScoreNormalized = ...
    D.classLabel .* model.classScoreNormalized;
model.sumAbsAvailableWeight = denominator;
end

function bestIndex = chooseBestModel(accuracy, nNonzero, lambda)
valid = isfinite(accuracy) & nNonzero > 0;
assert(any(valid), 'Regularization path has no nonzero finite model.');
bestAccuracy = max(accuracy(valid));
candidate = find(valid & abs(accuracy - bestAccuracy) < 1e-12);
fewestSites = min(nNonzero(candidate));
candidate = candidate(nNonzero(candidate) == fewestSites);
[~, largestLambda] = max(lambda(candidate));
bestIndex = candidate(largestLambda);
end

function nNonzero = countNonzero(B, tolerance)
scale = max(abs(B), [], 1);
threshold = tolerance .* max(scale, 1);
nNonzero = sum(abs(B) > threshold, 1);
end

function predicted = signClass(score)
predicted = ones(size(score));
predicted(score < 0) = -1;
end

function fig = makeFigure(OUT)
colors = [0.36 0.24 0.55; 0.18 0.49 0.42; ...
    0.79 0.45 0.15; 0.69 0.20 0.24];
fig = figure('Color', 'w', 'Position', [80 80 1380 800]);

subplot(2, 2, 1);
b = bar(100 * OUT.accuracy, 'grouped');
for idx = 1:numel(b)
    b(idx).FaceColor = colors(idx, :);
end
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('In-sample accuracy (%)');
ylim([50 100]);
legend(OUT.methodNames, 'Location', 'southoutside', ...
    'Orientation', 'horizontal');
title('Same trials within each region');
formatAxes(gca);

subplot(2, 2, 2);
b = bar(OUT.nNonzero, 'grouped');
for idx = 1:numel(b)
    b(idx).FaceColor = colors(idx, :);
end
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('Sites with nonzero weight');
title('Selected population size');
formatAxes(gca);

subplot(2, 2, 3);
plotPaths(OUT.V1, colors);
title(sprintf('V1 regularization paths (N = %d trials)', ...
    numel(OUT.V1.classLabel)));

subplot(2, 2, 4);
plotPaths(OUT.V4, colors);
title(sprintf('V4 regularization paths (N = %d trials)', ...
    numel(OUT.V4.classLabel)));

sgtitle(['Nilson attention decoding: quartet-centered sparse linear ' ...
    'models, exact curve, minimum 20 sites'], 'FontWeight', 'bold');
end

function plotPaths(region, colors)
models = {region.logisticLasso, region.elasticNetLogistic, region.l1Svm};
for idx = 1:numel(models)
    model = models{idx};
    plot(model.nNonzeroPath, 100 * model.accuracyPath, '-o', ...
        'Color', colors(idx + 1, :), 'LineWidth', 1.3, ...
        'MarkerSize', 3, 'DisplayName', modelLabel(model));
    hold on;
    plot(model.nNonzero, 100 * model.accuracy, 'o', ...
        'Color', colors(idx + 1, :), 'MarkerFaceColor', colors(idx + 1, :), ...
        'MarkerSize', 8, 'HandleVisibility', 'off');
end
plot(region.fullFisher.nNonzero, 100 * region.fullFisher.accuracy, 'd', ...
    'Color', colors(1, :), 'MarkerFaceColor', colors(1, :), ...
    'MarkerSize', 8, 'DisplayName', 'Full Fisher');
xlabel('Sites with nonzero weight');
ylabel('In-sample accuracy (%)');
ylim([50 100]);
legend('Location', 'southeast');
formatAxes(gca);
end

function label = modelLabel(model)
if strcmp(model.learner, 'svm')
    label = 'L1 SVM';
elseif model.alpha == 1
    label = 'Logistic lasso';
else
    label = sprintf('Elastic net (alpha = %.2g)', model.alpha);
end
end

function formatAxes(ax)
grid(ax, 'on');
box(ax, 'off');
set(ax, 'FontSize', 11, 'Layer', 'top');
end

function printSummary(OUT)
fprintf('\nSparse linear attention decoder summary (in-sample)\n');
fprintf('  Method                         V1 accuracy/sites   V4 accuracy/sites\n');
for methodIdx = 1:numel(OUT.methodNames)
    fprintf('  %-30s %6.2f%% / %3d     %6.2f%% / %3d\n', ...
        OUT.methodNames{methodIdx}, ...
        100 * OUT.accuracy(1, methodIdx), OUT.nNonzero(1, methodIdx), ...
        100 * OUT.accuracy(2, methodIdx), OUT.nNonzero(2, methodIdx));
end
end
