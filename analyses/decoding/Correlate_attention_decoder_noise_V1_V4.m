function OUT = Correlate_attention_decoder_noise_V1_V4(varargin)
%CORRELATE_ATTENTION_DECODER_NOISE_V1_V4 Correlate within-stimulus scores.
%
% Matches V1 and V4 trials, uses target-aligned scores from the final
% cross-validation-selected decoders, subtracts each area's mean score for
% every stimulus, and correlates the pooled trial residuals.

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
stimulus = decoder.V1.stimulus(v1Index);
quartet = decoder.V1.quartet(v1Index);
classLabel = decoder.V1.classLabel(v1Index);
assert(isequal(stimulus, decoder.V4.stimulus(v4Index)), ...
    'V1 and V4 stimulus labels differ for matched trials.');
assert(isequal(quartet, decoder.V4.quartet(v4Index)), ...
    'V1 and V4 quartet labels differ for matched trials.');
assert(isequal(classLabel, decoder.V4.classLabel(v4Index)), ...
    'V1 and V4 class labels differ for matched trials.');

scoreV1 = double(decoder.V1.finalTargetAlignedScore(v1Index));
scoreV4 = double(decoder.V4.finalTargetAlignedScore(v4Index));
assert(all(isfinite(scoreV1)) && all(isfinite(scoreV4)), ...
    'Matched final decoder scores must be finite.');

[residualV1, stimulusValues, meanByStimulusV1, trialsPerStimulus] = ...
    subtractStimulusMean(scoreV1, stimulus);
[residualV4, stimulusValuesV4, meanByStimulusV4, trialsPerStimulusV4] = ...
    subtractStimulusMean(scoreV4, stimulus);
assert(isequal(stimulusValues, stimulusValuesV4) && ...
    isequal(trialsPerStimulus, trialsPerStimulusV4), ...
    'V1 and V4 stimulus groupings differ.');

[correlation, pValue] = corr(residualV1, residualV4, 'Type', 'Pearson');
nTrials = numel(trialIndex);
fisherZ = atanh(max(-1 + eps, min(1 - eps, correlation)));
zHalfWidth = 1.96 / sqrt(nTrials - 3);
confidenceInterval = tanh(fisherZ + [-1 1] * zHalfWidth);
regression = [ones(nTrials, 1), residualV1] \ residualV4;

OUT = struct();
OUT.description = ['Pearson correlation between final target-aligned V1 ' ...
    'and V4 decoder scores after subtracting each area''s mean per stimulus.'];
OUT.monkey = decoder.monkey;
OUT.sourceFile = sourceFile;
OUT.scoreSource = 'finalTargetAlignedScore';
OUT.residualization = 'Separate V1 and V4 means over shared trials per stimulus.';
OUT.trialIndex = trialIndex;
OUT.stimulus = stimulus;
OUT.quartet = quartet;
OUT.classLabel = classLabel;
OUT.scoreV1 = scoreV1;
OUT.scoreV4 = scoreV4;
OUT.residualV1 = residualV1;
OUT.residualV4 = residualV4;
OUT.stimulusValues = stimulusValues;
OUT.meanByStimulusV1 = meanByStimulusV1;
OUT.meanByStimulusV4 = meanByStimulusV4;
OUT.trialsPerStimulus = trialsPerStimulus;
OUT.nTrials = nTrials;
OUT.nStimuli = numel(stimulusValues);
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
    'Attention_decoder_noise_correlation_V1_V4_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_noise_correlation_V1_V4_N.png');
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

fprintf(['V1-V4 decoder noise correlation: r = %.4f, ' ...
    '95%% CI [%.4f %.4f], p = %.4g, N = %d trials, %d stimuli.\n'], ...
    OUT.pearsonR, OUT.confidenceInterval95(1), ...
    OUT.confidenceInterval95(2), OUT.pValue, OUT.nTrials, OUT.nStimuli);
end

function [residual, stimulusValues, meanByStimulus, trialsPerStimulus] = ...
        subtractStimulusMean(score, stimulus)
[stimulusValues, ~, stimulusIndex] = unique(stimulus(:), 'sorted');
nStimuli = numel(stimulusValues);
meanByStimulus = nan(nStimuli, 1);
trialsPerStimulus = zeros(nStimuli, 1);
residual = nan(size(score));
for stimulusIdx = 1:nStimuli
    use = stimulusIndex == stimulusIdx;
    trialsPerStimulus(stimulusIdx) = nnz(use);
    meanByStimulus(stimulusIdx) = mean(score(use));
    residual(use) = score(use) - meanByStimulus(stimulusIdx);
end
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Position', [100 100 760 700]);
ax = axes(fig);
hold(ax, 'on');

scatter(ax, OUT.residualV1, OUT.residualV4, 20, ...
    [0.25 0.55 0.82], 'filled', ...
    'MarkerFaceAlpha', 0.30, 'MarkerEdgeAlpha', 0.15);

xLimits = paddedLimits(OUT.residualV1);
yLimits = paddedLimits(OUT.residualV4);
plot(ax, xLimits, OUT.regressionIntercept + ...
    OUT.regressionSlope * xLimits, '-', ...
    'Color', [0.82 0.12 0.12], 'LineWidth', 2.5);
plot(ax, xLimits, [0 0], '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);
plot(ax, [0 0], yLimits, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);

xlim(ax, xLimits);
ylim(ax, yLimits);
xlabel(ax, 'V1 decoder residual (within stimulus)');
ylabel(ax, 'V4 decoder residual (within stimulus)');
title(ax, {'Noise correlation of V1 and V4 decoder scores', ...
    sprintf('r = %.3f, 95%% CI [%.3f, %.3f], p = %.2g', ...
    OUT.pearsonR, OUT.confidenceInterval95(1), ...
    OUT.confidenceInterval95(2), OUT.pValue)});
text(ax, 0.97, 0.05, sprintf('N = %d trials\nN = %d stimuli', ...
    OUT.nTrials, OUT.nStimuli), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontSize', 11, 'Color', [0.20 0.20 0.20]);
axis(ax, 'square');
box(ax, 'off');
set(ax, 'FontSize', 12, 'LineWidth', 1, 'TickDir', 'out');
end

function limits = paddedLimits(values)
limits = [min(values), max(values)];
padding = 0.05 * diff(limits);
if padding == 0
    padding = max(1, abs(limits(1)) * 0.05);
end
limits = limits + [-padding padding];
end
