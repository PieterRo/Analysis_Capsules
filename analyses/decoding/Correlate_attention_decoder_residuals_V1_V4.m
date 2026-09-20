function OUT = Correlate_attention_decoder_residuals_V1_V4(varargin)
%CORRELATE_ATTENTION_DECODER_RESIDUALS_V1_V4 Correlate final decoder scores.
%
% Matches V1 and V4 trials, uses target-aligned scores from the final
% cross-validation-selected decoders, subtracts each area's mean score per
% quartet, and computes the across-trial Pearson correlation.

p = inputParser;
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();
sourceFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat');
source = load(sourceFile, 'OUT');
decoder = source.OUT;

[trialIndex, v1Index, v4Index] = intersect( ...
    decoder.V1.trialIndex, decoder.V4.trialIndex);
quartet = decoder.V1.quartet(v1Index);
classLabel = decoder.V1.classLabel(v1Index);
assert(isequal(quartet, decoder.V4.quartet(v4Index)), ...
    'V1 and V4 quartet labels differ for matched trials.');
assert(isequal(classLabel, decoder.V4.classLabel(v4Index)), ...
    'V1 and V4 class labels differ for matched trials.');

scoreV1 = double(decoder.V1.finalTargetAlignedScore(v1Index));
scoreV4 = double(decoder.V4.finalTargetAlignedScore(v4Index));
assert(all(isfinite(scoreV1)) && all(isfinite(scoreV4)), ...
    'Matched final decoder scores must be finite.');

[residualV1, meanByQuartetV1] = subtractQuartetMean(scoreV1, quartet);
[residualV4, meanByQuartetV4] = subtractQuartetMean(scoreV4, quartet);
[correlation, pValue] = corr(residualV1, residualV4, 'Type', 'Pearson');

nTrials = numel(trialIndex);
fisherZ = atanh(max(-1 + eps, min(1 - eps, correlation)));
zHalfWidth = 1.96 / sqrt(nTrials - 3);
confidenceInterval = tanh(fisherZ + [-1 1] * zHalfWidth);
regression = [ones(nTrials, 1), residualV1] \ residualV4;

OUT = struct();
OUT.description = ['Pearson correlation between final target-aligned V1 ' ...
    'and V4 decoder scores after subtracting each area''s mean per quartet.'];
OUT.monkey = decoder.monkey;
OUT.sourceFile = sourceFile;
OUT.scoreSource = 'finalTargetAlignedScore';
OUT.residualization = 'Separate V1 and V4 means over shared trials per quartet.';
OUT.trialIndex = trialIndex;
OUT.quartet = quartet;
OUT.classLabel = classLabel;
OUT.scoreV1 = scoreV1;
OUT.scoreV4 = scoreV4;
OUT.residualV1 = residualV1;
OUT.residualV4 = residualV4;
OUT.meanByQuartetV1 = meanByQuartetV1;
OUT.meanByQuartetV4 = meanByQuartetV4;
OUT.nTrials = nTrials;
OUT.nQuartets = numel(unique(quartet));
OUT.pearsonR = correlation;
OUT.pValue = pValue;
OUT.confidenceInterval95 = confidenceInterval;
OUT.regressionIntercept = regression(1);
OUT.regressionSlope = regression(2);
OUT.decoderV1 = struct('alpha', decoder.V1.selectedAlpha, ...
    'lambda', decoder.V1.selectedLambda, ...
    'nSites', decoder.V1.finalNonzero);
OUT.decoderV4 = struct('alpha', decoder.V4.selectedAlpha, ...
    'lambda', decoder.V4.selectedLambda, ...
    'nSites', decoder.V4.finalNonzero);

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

resultFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_residual_correlation_V1_V4_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_residual_correlation_V1_V4_N.png');
OUT.resultFile = resultFile;
OUT.figureFile = figureFile;

if opt.SaveOutputs
    save(resultFile, 'OUT');
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

fprintf(['V1-V4 final decoder residual correlation: r = %.4f, ' ...
    '95%% CI [%.4f %.4f], p = %.4g, N = %d trials, %d quartets.\n'], ...
    OUT.pearsonR, OUT.confidenceInterval95(1), ...
    OUT.confidenceInterval95(2), OUT.pValue, OUT.nTrials, OUT.nQuartets);
end

function [residual, meanByQuartet] = subtractQuartetMean(score, quartet)
nQuartets = max(quartet);
meanByQuartet = nan(nQuartets, 1);
residual = nan(size(score));
for quartetValue = unique(quartet(:))'
    use = quartet == quartetValue;
    meanByQuartet(quartetValue) = mean(score(use));
    residual(use) = score(use) - meanByQuartet(quartetValue);
end
assert(all(isfinite(residual)), 'Quartet residualization produced non-finite values.');
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Position', [100 100 880 760]);
scatter(OUT.residualV1, OUT.residualV4, 14, [0.24 0.43 0.56], ...
    'filled', 'MarkerFaceAlpha', 0.28, 'MarkerEdgeAlpha', 0.12);
hold on;
xline(0, ':', 'Color', [0.55 0.55 0.55]);
yline(0, ':', 'Color', [0.55 0.55 0.55]);

xFit = linspace(min(OUT.residualV1), max(OUT.residualV1), 200);
yFit = OUT.regressionIntercept + OUT.regressionSlope * xFit;
plot(xFit, yFit, '-', 'Color', [0.78 0.18 0.16], 'LineWidth', 2.2);

xlabel('V1 decoder-score residual');
ylabel('V4 decoder-score residual');
title({'Trial-wise V1-V4 decoder correlation after quartet-mean subtraction', ...
    sprintf('Pearson r = %.3f, 95%% CI [%.3f, %.3f], p = %.2g', ...
    OUT.pearsonR, OUT.confidenceInterval95(1), ...
    OUT.confidenceInterval95(2), OUT.pValue)});
text(0.97, 0.05, sprintf('N = %d trials\nN = %d quartets', ...
    OUT.nTrials, OUT.nQuartets), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 5);
grid on;
box off;
set(gca, 'FontSize', 12, 'Layer', 'top');
end
