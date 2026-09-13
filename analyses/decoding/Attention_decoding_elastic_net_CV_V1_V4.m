function OUT = Attention_decoding_elastic_net_CV_V1_V4(varargin)
%ATTENTION_DECODING_ELASTIC_NET_CV_V1_V4 Tune alpha and lambda per area.
%
% Uses repeated, quartet- and class-balanced 10-fold cross-validation.
% Quartet means are estimated from training trials within each fold and
% applied to both training and held-out trials. V1 and V4 are tuned
% independently using the same grid and one-standard-error selection rule.

p = inputParser;
p.addParameter('AlphaGrid', [0.1 0.25 0.5 0.75 1], ...
    @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && ...
    all(x > 0) && all(x <= 1));
p.addParameter('NumLambda', 24, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 10 && x == floor(x));
p.addParameter('LambdaRatio', 1e-3, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('NumFolds', 10, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && x == floor(x));
p.addParameter('NumRepeats', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('RandomSeed', 130926, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('CoefficientTolerance', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
opt.AlphaGrid = unique(double(opt.AlphaGrid(:)'));

cfg = config();
baseFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');
base = loadFeatureSource(baseFile);

OUT = struct();
OUT.description = ['Repeated quartet-balanced cross-validation for ' ...
    'elastic-net logistic attention decoders.'];
OUT.monkey = 'Mr Nilson';
OUT.responseCentering = 'training-fold quartet means';
OUT.responseWindow = base.responseWindow;
OUT.curveMarginDeg = base.curveMarginDeg;
OUT.minSitesPerQuartet = base.minSitesPerQuartet;
OUT.alphaGrid = opt.AlphaGrid;
OUT.numLambda = opt.NumLambda;
OUT.lambdaRatio = opt.LambdaRatio;
OUT.numFolds = opt.NumFolds;
OUT.numRepeats = opt.NumRepeats;
OUT.randomSeed = opt.RandomSeed;
OUT.selectionRule = ['Sparsest alpha/lambda combination within one standard ' ...
    'error of the highest mean cross-validated accuracy.'];

OUT.V1 = crossValidateRegion(base.V1, opt, opt.RandomSeed + 1000);
OUT.V4 = crossValidateRegion(base.V4, opt, opt.RandomSeed + 2000);

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

resultFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.png');
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

function base = loadFeatureSource(baseFile)
needsRebuild = true;
if isfile(baseFile)
    data = load(baseFile, 'OUT');
    base = data.OUT;
    needsRebuild = ~isfield(base.V1, 'normalizedResponse') || ...
        ~isfield(base.V1, 'canonicalSign') || ...
        ~isfield(base.V4, 'normalizedResponse') || ...
        ~isfield(base.V4, 'canonicalSign');
end

if needsRebuild
    fprintf('Rebuilding the feature source with uncentered trial responses.\n');
    base = Attention_decoding_sparse_linear_V1_V4( ...
        'MakeFigure', false, 'SaveOutputs', true);
end
end

function R = crossValidateRegion(D, opt, seed)
region = D.region;
fprintf('\n%s elastic-net cross-validation: %d trials, %d sites.\n', ...
    region, numel(D.classLabel), size(D.normalizedResponse, 2));

fullTrainingMask = true(numel(D.classLabel), 1);
fullX = makeFoldFeatures(D, fullTrainingMask);
maxDifference = max(abs(fullX(:) - D.featureMatrix(:)));
assert(maxDifference < 1e-10, ...
    '%s full-data quartet centering differs from the source by %.3g.', ...
    region, maxDifference);

nAlpha = numel(opt.AlphaGrid);
nLambda = opt.NumLambda;
nTrials = numel(D.classLabel);
lambdaFraction = logspace(log10(opt.LambdaRatio), 0, nLambda);
lambdaGrid = zeros(nAlpha, nLambda);
fullWeights = classBalancingWeights(D.classLabel);

for alphaIdx = 1:nAlpha
    [~, pilot] = lassoglm(fullX, D.classLabel == 1, 'binomial', ...
        'Alpha', opt.AlphaGrid(alphaIdx), 'NumLambda', nLambda, ...
        'LambdaRatio', opt.LambdaRatio, 'Standardize', true, ...
        'Weights', fullWeights, 'CV', 'resubstitution');
    lambdaGrid(alphaIdx, :) = max(pilot.Lambda) * lambdaFraction;
end

scoreByTrial = nan(nTrials, nLambda, nAlpha, opt.NumRepeats);
nNonzeroByFold = nan(nAlpha, nLambda, opt.NumRepeats, opt.NumFolds);
foldIndex = zeros(nTrials, opt.NumRepeats);

for repeatIdx = 1:opt.NumRepeats
    foldIndex(:, repeatIdx) = balancedFoldIndex( ...
        D.quartet, D.classLabel, opt.NumFolds, seed + repeatIdx);
    fprintf('  repeat %d/%d\n', repeatIdx, opt.NumRepeats);
    for fold = 1:opt.NumFolds
        isTest = foldIndex(:, repeatIdx) == fold;
        isTrain = ~isTest;
        X = makeFoldFeatures(D, isTrain);
        trainWeights = classBalancingWeights(D.classLabel(isTrain));

        for alphaIdx = 1:nAlpha
            [B, stats] = lassoglm(X(isTrain, :), ...
                D.classLabel(isTrain) == 1, 'binomial', ...
                'Alpha', opt.AlphaGrid(alphaIdx), ...
                'Lambda', lambdaGrid(alphaIdx, :), ...
                'Standardize', true, 'Weights', trainWeights, ...
                'CV', 'resubstitution');
            scoreByTrial(isTest, :, alphaIdx, repeatIdx) = ...
                X(isTest, :) * B + stats.Intercept;
            nNonzeroByFold(alphaIdx, :, repeatIdx, fold) = ...
                countNonzero(B, opt.CoefficientTolerance);
        end
        fprintf('    fold %d/%d complete\n', fold, opt.NumFolds);
    end
end

accuracyByRepeat = nan(nAlpha, nLambda, opt.NumRepeats);
accuracyByFold = nan(nAlpha, nLambda, opt.NumRepeats, opt.NumFolds);
for repeatIdx = 1:opt.NumRepeats
    for alphaIdx = 1:nAlpha
        scores = scoreByTrial(:, :, alphaIdx, repeatIdx);
        accuracyByRepeat(alphaIdx, :, repeatIdx) = ...
            mean(signClass(scores) == D.classLabel, 1);
        for fold = 1:opt.NumFolds
            use = foldIndex(:, repeatIdx) == fold;
            accuracyByFold(alphaIdx, :, repeatIdx, fold) = ...
                mean(signClass(scores(use, :)) == D.classLabel(use), 1);
        end
    end
end

meanAccuracy = mean(accuracyByRepeat, 3);
foldSamples = reshape(accuracyByFold, nAlpha, nLambda, []);
accuracySE = std(foldSamples, 0, 3) / sqrt(size(foldSamples, 3));
meanNonzero = mean(mean(nNonzeroByFold, 4), 3);
[alphaIdx, lambdaIdx, bestAccuracy, oneSEThreshold] = ...
    selectHyperparameters(meanAccuracy, accuracySE, meanNonzero);

selectedScores = squeeze(scoreByTrial(:, lambdaIdx, alphaIdx, :));
if opt.NumRepeats == 1
    selectedScores = selectedScores(:);
end
crossValidatedScore = mean(selectedScores, 2);
crossValidatedAccuracy = mean( ...
    signClass(crossValidatedScore) == D.classLabel);

selectedAlpha = opt.AlphaGrid(alphaIdx);
selectedLambda = lambdaGrid(alphaIdx, lambdaIdx);
[finalB, finalStats] = lassoglm(fullX, D.classLabel == 1, 'binomial', ...
    'Alpha', selectedAlpha, 'Lambda', selectedLambda, ...
    'Standardize', true, 'Weights', fullWeights, 'CV', 'resubstitution');
finalClassScore = fullX * finalB + finalStats.Intercept;
finalNonzero = countNonzero(finalB, opt.CoefficientTolerance);

denominator = D.availability * abs(finalB);
normalizedScore = nan(size(finalClassScore));
hasWeight = denominator > 0;
normalizedScore(hasWeight) = ...
    (fullX(hasWeight, :) * finalB) ./ denominator(hasWeight);

R = struct();
R.region = region;
R.trialIndex = D.trialIndex;
R.stimulus = D.stimulus;
R.quartet = D.quartet;
R.classLabel = D.classLabel;
R.siteIndexLocal = D.siteIndexLocal;
R.siteIndexGlobal = D.siteIndexGlobal;
R.foldIndex = foldIndex;
R.alphaGrid = opt.AlphaGrid;
R.lambdaGrid = lambdaGrid;
R.lambdaFraction = lambdaFraction;
R.accuracyByRepeat = accuracyByRepeat;
R.accuracyByFold = accuracyByFold;
R.meanAccuracy = meanAccuracy;
R.accuracySE = accuracySE;
R.meanNonzero = meanNonzero;
R.nNonzeroByFold = nNonzeroByFold;
R.bestMeanAccuracy = bestAccuracy;
R.oneSEThreshold = oneSEThreshold;
R.selectedAlphaIndex = alphaIdx;
R.selectedLambdaIndex = lambdaIdx;
R.selectedAlpha = selectedAlpha;
R.selectedLambda = selectedLambda;
R.selectedLambdaFraction = lambdaFraction(lambdaIdx);
R.selectedMeanAccuracy = meanAccuracy(alphaIdx, lambdaIdx);
R.selectedAccuracySE = accuracySE(alphaIdx, lambdaIdx);
R.selectedMeanNonzero = meanNonzero(alphaIdx, lambdaIdx);
R.crossValidatedClassScoreByRepeat = selectedScores;
R.crossValidatedClassScore = crossValidatedScore;
R.crossValidatedTargetAlignedScore = ...
    D.classLabel .* crossValidatedScore;
R.crossValidatedAccuracyOfMeanScore = crossValidatedAccuracy;
R.finalCoefficient = finalB;
R.finalIntercept = finalStats.Intercept;
R.finalNonzero = finalNonzero;
R.finalClassScore = finalClassScore;
R.finalTargetAlignedScore = D.classLabel .* finalClassScore;
R.finalInSampleAccuracy = mean( ...
    signClass(finalClassScore) == D.classLabel);
R.finalClassScoreNormalized = normalizedScore;
R.finalTargetAlignedScoreNormalized = D.classLabel .* normalizedScore;
R.finalSumAbsAvailableWeight = denominator;
R.selectionNote = ['Alpha and lambda were selected independently for this ' ...
    'area with the common one-standard-error rule.'];

fprintf(['%s selected alpha %.3g, lambda %.4g (fraction %.4g): ' ...
    'CV %.2f%% +/- %.2f%%, mean %.1f sites; final %d sites.\n'], ...
    region, selectedAlpha, selectedLambda, ...
    R.selectedLambdaFraction, 100 * R.selectedMeanAccuracy, ...
    100 * R.selectedAccuracySE, R.selectedMeanNonzero, R.finalNonzero);
end

function X = makeFoldFeatures(D, trainingMask)
nQuartets = max(D.quartet);
nSites = size(D.normalizedResponse, 2);
trainingMean = nan(nQuartets, nSites);
for quartet = unique(D.quartet(:))'
    use = trainingMask & D.quartet == quartet;
    assert(any(use), 'Training fold lacks quartet %d.', quartet);
    trainingMean(quartet, :) = ...
        mean(D.normalizedResponse(use, :), 1, 'omitnan');
end

centered = D.normalizedResponse - trainingMean(D.quartet, :);
X = centered .* D.canonicalSign;
X(~D.availability | ~isfinite(X)) = 0;
end

function foldIndex = balancedFoldIndex(quartet, classLabel, nFolds, seed)
rng(seed, 'twister');
foldIndex = zeros(numel(classLabel), 1);
for quartetValue = unique(quartet(:))'
    for classValue = [-1 1]
        trials = find(quartet == quartetValue & classLabel == classValue);
        assert(numel(trials) >= 2, ...
            'Quartet %d class %d has fewer than two trials.', ...
            quartetValue, classValue);
        trials = trials(randperm(numel(trials)));
        foldOrder = randperm(nFolds);
        labels = repmat(foldOrder, 1, ceil(numel(trials) / nFolds));
        foldIndex(trials) = labels(1:numel(trials));
    end
end
assert(all(foldIndex > 0), 'Some trials were not assigned to a fold.');
for fold = 1:nFolds
    assert(all(ismember([-1 1], unique(classLabel(foldIndex == fold))')), ...
        'Fold %d does not contain both classes.', fold);
end
end

function weights = classBalancingWeights(classLabel)
n = numel(classLabel);
nPositive = nnz(classLabel == 1);
nNegative = nnz(classLabel == -1);
weights = zeros(n, 1);
weights(classLabel == 1) = n / (2 * nPositive);
weights(classLabel == -1) = n / (2 * nNegative);
end

function [alphaIdx, lambdaIdx, bestAccuracy, threshold] = ...
    selectHyperparameters(meanAccuracy, accuracySE, meanNonzero)
[bestAccuracy, bestLinearIndex] = max(meanAccuracy(:));
threshold = bestAccuracy - accuracySE(bestLinearIndex);
candidate = isfinite(meanAccuracy) & meanAccuracy >= threshold & ...
    meanNonzero > 0;
assert(any(candidate(:)), 'No model satisfies the one-SE selection rule.');

fewestSites = min(meanNonzero(candidate));
candidate = candidate & abs(meanNonzero - fewestSites) < 1e-12;
candidateIndex = find(candidate);
[~, bestCandidate] = max(meanAccuracy(candidateIndex));
[alphaIdx, lambdaIdx] = ind2sub(size(meanAccuracy), ...
    candidateIndex(bestCandidate));
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
fig = figure('Color', 'w', 'Position', [80 80 1380 800]);
colors = lines(numel(OUT.alphaGrid));

subplot(2, 2, 1);
plotRegionPaths(OUT.V1, colors);
title(sprintf('V1: selected alpha %.2g, lambda fraction %.3g', ...
    OUT.V1.selectedAlpha, OUT.V1.selectedLambdaFraction));

subplot(2, 2, 2);
plotRegionPaths(OUT.V4, colors);
title(sprintf('V4: selected alpha %.2g, lambda fraction %.3g', ...
    OUT.V4.selectedAlpha, OUT.V4.selectedLambdaFraction));

subplot(2, 2, 3);
accuracy = 100 * [OUT.V1.selectedMeanAccuracy, OUT.V4.selectedMeanAccuracy];
accuracySE = 100 * [OUT.V1.selectedAccuracySE, OUT.V4.selectedAccuracySE];
bar(accuracy, 'FaceColor', [0.24 0.48 0.58]);
hold on;
errorbar(1:2, accuracy, accuracySE, '.k', 'LineWidth', 1.4);
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('Cross-validated accuracy (%)');
ylim([50 100]);
title('One-SE selected models');
formatAxes(gca);

subplot(2, 2, 4);
siteCount = [OUT.V1.finalNonzero, OUT.V4.finalNonzero];
bar(siteCount, 'FaceColor', [0.72 0.39 0.19]);
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('Nonzero sites in final fit');
title('Final decoder population');
formatAxes(gca);

sgtitle(sprintf(['Nilson elastic-net attention decoder: %dx%d-fold ' ...
    'quartet-balanced cross-validation'], OUT.numRepeats, OUT.numFolds), ...
    'FontWeight', 'bold');
end

function plotRegionPaths(R, colors)
for alphaIdx = 1:numel(R.alphaGrid)
    plot(R.meanNonzero(alphaIdx, :), ...
        100 * R.meanAccuracy(alphaIdx, :), '-o', ...
        'Color', colors(alphaIdx, :), 'LineWidth', 1.3, ...
        'MarkerSize', 3, ...
        'DisplayName', sprintf('alpha = %.2g', R.alphaGrid(alphaIdx)));
    hold on;
end
plot(R.selectedMeanNonzero, 100 * R.selectedMeanAccuracy, 'kp', ...
    'MarkerFaceColor', [1 0.82 0.18], 'MarkerSize', 11, ...
    'DisplayName', 'one-SE selection');
chanceLine = yline(50, ':', 'Color', [0.4 0.4 0.4]);
chanceLine.HandleVisibility = 'off';
xlabel('Mean sites with nonzero weight');
ylabel('Cross-validated accuracy (%)');
ylim([50 100]);
legend('Location', 'southeast');
formatAxes(gca);
end

function formatAxes(ax)
grid(ax, 'on');
box(ax, 'off');
set(ax, 'FontSize', 11, 'Layer', 'top');
end

function printSummary(OUT)
fprintf('\nElastic-net cross-validation summary\n');
regions = {'V1', 'V4'};
for idx = 1:numel(regions)
    R = OUT.(regions{idx});
    fprintf(['  %s: alpha %.3g, lambda %.4g; CV %.2f%% +/- %.2f%%; ' ...
        'mean fold sites %.1f; final sites %d; final in-sample %.2f%%.\n'], ...
        regions{idx}, R.selectedAlpha, R.selectedLambda, ...
        100 * R.selectedMeanAccuracy, 100 * R.selectedAccuracySE, ...
        R.selectedMeanNonzero, R.finalNonzero, ...
        100 * R.finalInSampleAccuracy);
end
end
