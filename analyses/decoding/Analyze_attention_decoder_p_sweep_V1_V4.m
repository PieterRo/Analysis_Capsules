function SUMMARY = Analyze_attention_decoder_p_sweep_V1_V4()
%ANALYZE_ATTENTION_DECODER_P_SWEEP_V1_V4 Summarize site-level p selection.

cfg = config();
labels = {'all', 'p < 0.20', 'p < 0.10', 'p < 0.05'};
suffix = {'', '_p0p2', '_p0p1', '_p0p05'};
nLevels = numel(labels);
V1 = cell(nLevels, 1);
V4 = cell(nLevels, 1);

for level = 1:nLevels
    v1File = fullfile(cfg.resultsDir, sprintf( ...
        ['Attention_decoder_V1_proof_margin0deg_minSites20_' ...
        'quartetCentered%s_N.mat'], suffix{level}));
    v4File = fullfile(cfg.resultsDir, sprintf( ...
        ['Attention_decoder_V4_fullFisher_minSites20_' ...
        'quartetCentered%s_N.mat'], suffix{level}));
    v1Data = load(v1File, 'OUT');
    v4Data = load(v4File, 'OUT');
    V1{level} = v1Data.OUT;
    V4{level} = v4Data.OUT;
end

SUMMARY = struct();
SUMMARY.description = ['Site-level attention p-value sensitivity analysis ' ...
    'for quartet-centered full-Fisher V1 and V4 decoders.'];
SUMMARY.labels = labels;
SUMMARY.pThreshold = [Inf 0.20 0.10 0.05];
SUMMARY.nSitesV1 = zeros(1, nLevels);
SUMMARY.nSitesV4 = zeros(1, nLevels);
SUMMARY.nQuartetsV1 = zeros(1, nLevels);
SUMMARY.nQuartetsV4 = zeros(1, nLevels);
SUMMARY.nTrialsV1 = zeros(1, nLevels);
SUMMARY.nTrialsV4 = zeros(1, nLevels);
SUMMARY.accuracyV1 = zeros(1, nLevels);
SUMMARY.accuracyV4 = zeros(1, nLevels);
SUMMARY.commonTrialIndex = cell(nLevels, 1);
SUMMARY.commonScoreV1 = cell(nLevels, 1);
SUMMARY.commonScoreV4 = cell(nLevels, 1);
SUMMARY.commonAccuracyV1 = zeros(1, nLevels);
SUMMARY.commonAccuracyV4 = zeros(1, nLevels);

for level = 1:nLevels
    SUMMARY.nSitesV1(level) = nnz(V1{level}.eligibleSiteMask);
    SUMMARY.nSitesV4(level) = nnz(V4{level}.eligibleSiteMask);
    SUMMARY.nQuartetsV1(level) = numel(V1{level}.includedQuartet);
    SUMMARY.nQuartetsV4(level) = numel(V4{level}.includedQuartet);
    SUMMARY.nTrialsV1(level) = V1{level}.nScored;
    SUMMARY.nTrialsV4(level) = V4{level}.nScored;
    SUMMARY.accuracyV1(level) = V1{level}.accuracyFullFisher;
    SUMMARY.accuracyV4(level) = V4{level}.accuracyFullFisher;

    [commonTrial, v1Idx, v4Idx] = intersect( ...
        V1{level}.trialIndex, V4{level}.trialIndex);
    SUMMARY.commonTrialIndex{level} = commonTrial;
    SUMMARY.commonScoreV1{level} = V1{level}.SFullFisher(v1Idx);
    SUMMARY.commonScoreV4{level} = V4{level}.SFullFisher(v4Idx);
    SUMMARY.commonAccuracyV1(level) = ...
        mean(SUMMARY.commonScoreV1{level} > 0);
    SUMMARY.commonAccuracyV4(level) = ...
        mean(SUMMARY.commonScoreV4{level} > 0);
end

fixedTrialV1 = V1{1}.trialIndex;
fixedTrialV4 = V4{1}.trialIndex;
for level = 2:nLevels
    fixedTrialV1 = intersect(fixedTrialV1, V1{level}.trialIndex);
    fixedTrialV4 = intersect(fixedTrialV4, V4{level}.trialIndex);
end
SUMMARY.fixedTrialIndexV1 = fixedTrialV1;
SUMMARY.fixedTrialIndexV4 = fixedTrialV4;
SUMMARY.fixedAccuracyV1 = zeros(1, nLevels);
SUMMARY.fixedAccuracyV4 = zeros(1, nLevels);

for level = 1:nLevels
    [presentV1, v1Idx] = ismember(fixedTrialV1, V1{level}.trialIndex);
    [presentV4, v4Idx] = ismember(fixedTrialV4, V4{level}.trialIndex);
    assert(all(presentV1) && all(presentV4), ...
        'Fixed trial set is missing from one of the decoder results.');
    SUMMARY.fixedAccuracyV1(level) = ...
        mean(V1{level}.SFullFisher(v1Idx) > 0);
    SUMMARY.fixedAccuracyV4(level) = ...
        mean(V4{level}.SFullFisher(v4Idx) > 0);
end

fig = makeFigure(SUMMARY);
resultFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_pThreshold_sweep_V1_V4_quartetCentered_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_pThreshold_sweep_V1_V4_quartetCentered_N.png');
SUMMARY.resultFile = resultFile;
SUMMARY.figureFile = figureFile;
save(resultFile, 'SUMMARY');
print(fig, figureFile, '-dpng', '-r180');
SUMMARY.figure = fig;

fprintf('Saved %s\n', resultFile);
fprintf('Saved %s\n', figureFile);

end

function fig = makeFigure(S)
x = 1:numel(S.labels);
cV1 = [0.48 0.32 0.62];
cV4 = [0.75 0.35 0.18];
fig = figure('Color', 'w', 'Position', [100 100 1300 760]);

subplot(2, 2, 1);
plotPair(x, 100 * S.accuracyV1, 100 * S.accuracyV4, cV1, cV4);
ylabel('Accuracy (%)');
title('All available trials');

subplot(2, 2, 2);
plotPair(x, 100 * S.fixedAccuracyV1, 100 * S.fixedAccuracyV4, cV1, cV4);
ylabel('Accuracy (%)');
title(sprintf('Fixed trial sets: V1 N = %d; V4 N = %d', ...
    numel(S.fixedTrialIndexV1), numel(S.fixedTrialIndexV4)));

subplot(2, 2, 3);
plotPair(x, S.nSitesV1, S.nSitesV4, cV1, cV4);
ylabel('Selected sites');
title('Site-level selection');

subplot(2, 2, 4);
plot(x, cellfun(@numel, S.commonTrialIndex), '-o', ...
    'Color', [0.20 0.45 0.38], 'LineWidth', 1.8, 'MarkerFaceColor', 'w');
ylabel('Common V1-V4 trials');
title('Trials available for later correlation');
formatAxes(gca, x, S.labels);

sgtitle(['Nilson quartet-centered full-Fisher decoder: ' ...
    'site-level attention p-value sensitivity'], 'FontWeight', 'bold');
end

function plotPair(x, yV1, yV4, cV1, cV4)
plot(x, yV1, '-o', 'Color', cV1, 'LineWidth', 1.8, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'V1');
hold on;
plot(x, yV4, '-o', 'Color', cV4, 'LineWidth', 1.8, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'V4');
legend('Location', 'best');
formatAxes(gca, x, {'all', 'p < 0.20', 'p < 0.10', 'p < 0.05'});
end

function formatAxes(ax, x, labels)
set(ax, 'XTick', x, 'XTickLabel', labels, 'FontSize', 11, 'Layer', 'top');
xlim(ax, [0.7, numel(x) + 0.3]);
grid(ax, 'on');
box(ax, 'off');
end
