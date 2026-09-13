function OUT = Correlate_attention_color_decoder_noise_V1_V4(varargin)
%CORRELATE_ATTENTION_COLOR_DECODER_NOISE_V1_V4 Compare decoder residuals.
%
% Uses final target-aligned attention and color decoder scores on the same
% V1/V4 trials. Each score is residualized per stimulus before computing
% four cross-area correlations: attention-attention, color-color, and the
% two attention-color combinations.

p = inputParser;
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();
attentionFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat');
colorFile = fullfile(cfg.resultsDir, ...
    'Color_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat');
attentionData = load(attentionFile, 'OUT');
colorData = load(colorFile, 'OUT');
attention = attentionData.OUT;
color = colorData.OUT;

assert(isequal(attention.V1.trialIndex, color.V1.trialIndex), ...
    'Attention and color V1 trial selections differ.');
assert(isequal(attention.V4.trialIndex, color.V4.trialIndex), ...
    'Attention and color V4 trial selections differ.');

[trialIndex, v1Index, v4Index] = intersect( ...
    attention.V1.trialIndex, attention.V4.trialIndex);
stimulus = attention.V1.stimulus(v1Index);
quartet = attention.V1.quartet(v1Index);
assert(isequal(stimulus, attention.V4.stimulus(v4Index)), ...
    'V1 and V4 stimulus labels differ for matched trials.');
assert(isequal(quartet, attention.V4.quartet(v4Index)), ...
    'V1 and V4 quartet labels differ for matched trials.');
assert(isequal(stimulus, color.V1.stimulus(v1Index)) && ...
    isequal(stimulus, color.V4.stimulus(v4Index)), ...
    'Attention and color stimulus labels differ.');

score = struct();
score.attentionV1 = double(attention.V1.finalTargetAlignedScore(v1Index));
score.attentionV4 = double(attention.V4.finalTargetAlignedScore(v4Index));
score.colorV1 = double(color.V1.finalTargetAlignedScore(v1Index));
score.colorV4 = double(color.V4.finalTargetAlignedScore(v4Index));

scoreNames = fieldnames(score);
residual = struct();
meanByStimulus = struct();
stimulusValues = [];
trialsPerStimulus = [];
for scoreIdx = 1:numel(scoreNames)
    name = scoreNames{scoreIdx};
    assert(all(isfinite(score.(name))), ...
        '%s contains non-finite decoder scores.', name);
    [residual.(name), values, means, counts] = ...
        subtractStimulusMean(score.(name), stimulus);
    meanByStimulus.(name) = means;
    if isempty(stimulusValues)
        stimulusValues = values;
        trialsPerStimulus = counts;
    else
        assert(isequal(stimulusValues, values) && ...
            isequal(trialsPerStimulus, counts), ...
            'Stimulus grouping differs between decoder scores.');
    end
end

comparison = struct( ...
    'name', {'attentionV1_attentionV4', 'colorV1_colorV4', ...
        'attentionV1_colorV4', 'colorV1_attentionV4'}, ...
    'xField', {'attentionV1', 'colorV1', 'attentionV1', 'colorV1'}, ...
    'yField', {'attentionV4', 'colorV4', 'colorV4', 'attentionV4'}, ...
    'xLabel', {'V1 attention', 'V1 color', 'V1 attention', 'V1 color'}, ...
    'yLabel', {'V4 attention', 'V4 color', 'V4 color', 'V4 attention'});

nTrials = numel(trialIndex);
statistics = struct();
for comparisonIdx = 1:numel(comparison)
    item = comparison(comparisonIdx);
    x = residual.(item.xField);
    y = residual.(item.yField);
    statistics.(item.name) = correlationStatistics(x, y, nTrials);
end

OUT = struct();
OUT.description = ['Cross-area Pearson correlations between final ' ...
    'target-aligned attention and color decoder scores after separate ' ...
    'within-stimulus mean subtraction.'];
OUT.monkey = attention.monkey;
OUT.attentionSourceFile = attentionFile;
OUT.colorSourceFile = colorFile;
OUT.scoreSource = 'finalTargetAlignedScore';
OUT.residualization = ['Each of four decoder scores separately centered ' ...
    'over shared trials within each exact stimulus, then pooled.'];
OUT.trialIndex = trialIndex;
OUT.stimulus = stimulus;
OUT.quartet = quartet;
OUT.stimulusValues = stimulusValues;
OUT.trialsPerStimulus = trialsPerStimulus;
OUT.score = score;
OUT.residual = residual;
OUT.meanByStimulus = meanByStimulus;
OUT.comparison = comparison;
OUT.statistics = statistics;
OUT.nTrials = nTrials;
OUT.nStimuli = numel(stimulusValues);
OUT.nQuartets = numel(unique(quartet));

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

resultFile = fullfile(cfg.resultsDir, ...
    'Attention_color_decoder_noise_correlations_V1_V4_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_color_decoder_noise_correlations_V1_V4_N.png');
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

printSummary(OUT);
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
assert(all(isfinite(residual)), ...
    'Within-stimulus residualization produced non-finite values.');
end

function stats = correlationStatistics(x, y, nTrials)
[r, pValue] = corr(x, y, 'Type', 'Pearson');
fisherZ = atanh(max(-1 + eps, min(1 - eps, r)));
zHalfWidth = 1.96 / sqrt(nTrials - 3);
stats = struct();
stats.pearsonR = r;
stats.pValue = pValue;
stats.confidenceInterval95 = tanh(fisherZ + [-1 1] * zHalfWidth);
regression = [ones(nTrials, 1), x] \ y;
stats.regressionIntercept = regression(1);
stats.regressionSlope = regression(2);
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Position', [80 80 1320 1050]);
panelColors = [0.20 0.48 0.72; 0.48 0.32 0.62; ...
    0.18 0.55 0.45; 0.78 0.42 0.16];

for comparisonIdx = 1:numel(OUT.comparison)
    item = OUT.comparison(comparisonIdx);
    stats = OUT.statistics.(item.name);
    x = OUT.residual.(item.xField);
    y = OUT.residual.(item.yField);
    ax = subplot(2, 2, comparisonIdx);
    hold(ax, 'on');
    scatter(ax, x, y, 16, panelColors(comparisonIdx, :), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.10);
    xLimits = paddedLimits(x);
    yLimits = paddedLimits(y);
    plot(ax, xLimits, stats.regressionIntercept + ...
        stats.regressionSlope * xLimits, '-', ...
        'Color', [0.78 0.12 0.12], 'LineWidth', 2.2);
    plot(ax, xLimits, [0 0], '-', ...
        'Color', [0.70 0.70 0.70], 'LineWidth', 0.7);
    plot(ax, [0 0], yLimits, '-', ...
        'Color', [0.70 0.70 0.70], 'LineWidth', 0.7);
    xlim(ax, xLimits);
    ylim(ax, yLimits);
    xlabel(ax, sprintf('%s decoder residual', item.xLabel));
    ylabel(ax, sprintf('%s decoder residual', item.yLabel));
    title(ax, {sprintf('%s vs %s', item.xLabel, item.yLabel), ...
        sprintf('r = %.3f, 95%% CI [%.3f, %.3f], p = %.2g', ...
        stats.pearsonR, stats.confidenceInterval95(1), ...
        stats.confidenceInterval95(2), stats.pValue)});
    axis(ax, 'square');
    box(ax, 'off');
    set(ax, 'FontSize', 10, 'LineWidth', 1, 'TickDir', 'out');
end

sgtitle(sprintf(['V1-V4 decoder noise correlations after stimulus-mean ' ...
    'subtraction; N = %d trials, %d stimuli'], ...
    OUT.nTrials, OUT.nStimuli), 'FontWeight', 'bold');
end

function limits = paddedLimits(values)
limits = [min(values), max(values)];
padding = 0.05 * diff(limits);
if padding == 0
    padding = max(1, abs(limits(1)) * 0.05);
end
limits = limits + [-padding padding];
end

function printSummary(OUT)
fprintf('\nV1-V4 decoder noise correlations\n');
for comparisonIdx = 1:numel(OUT.comparison)
    item = OUT.comparison(comparisonIdx);
    stats = OUT.statistics.(item.name);
    fprintf('  %-28s r = % .4f, 95%% CI [% .4f % .4f], p = %.4g\n', ...
        item.name, stats.pearsonR, stats.confidenceInterval95(1), ...
        stats.confidenceInterval95(2), stats.pValue);
end
fprintf('  N = %d trials, %d stimuli, %d quartets.\n', ...
    OUT.nTrials, OUT.nStimuli, OUT.nQuartets);
end
