function H = Plot_attention_decoding_quartet_centered_V1_V4()
%PLOT_ATTENTION_DECODING_QUARTET_CENTERED_V1_V4 Compare full-Fisher scores.

cfg = config();
v1File = fullfile(cfg.resultsDir, ...
    'Attention_decoder_V1_proof_margin0deg_minSites20_quartetCentered_N.mat');
v4File = fullfile(cfg.resultsDir, ...
    'Attention_decoder_V4_fullFisher_minSites20_quartetCentered_N.mat');
v1Data = load(v1File, 'OUT');
v4Data = load(v4File, 'OUT');
V1 = v1Data.OUT;
V4 = v4Data.OUT;

assert(strcmp(V1.responseCentering, 'quartet'), ...
    'The V1 result is not quartet-centered.');
assert(strcmp(V4.responseCentering, 'quartet'), ...
    'The V4 result is not quartet-centered.');
assert(V1.minSitesPerQuartet == 20 && V4.minSitesPerQuartet == 20, ...
    'Expected a minimum of 20 sites per quartet.');

allScores = [V1.SFullFisher(:); V4.SFullFisher(:)];
[~, edges] = histcounts(allScores, 'BinMethod', 'fd', ...
    'Normalization', 'probability');
probV1 = histcounts(V1.SFullFisher, edges, 'Normalization', 'probability');
probV4 = histcounts(V4.SFullFisher, edges, 'Normalization', 'probability');
yMax = 1.08 * max([probV1(:); probV4(:)]);
xLimits = [edges(1), edges(end)];

fig = figure('Color', 'w', 'Position', [100 100 1250 500]);
axV1 = subplot(1, 2, 1);
plotPanel(axV1, V1, edges, xLimits, yMax, [0.48 0.32 0.62]);
title(axV1, sprintf('V1: %.1f%% (%d/%d), N = %d quartets', ...
    100 * V1.accuracyFullFisher, V1.nCorrectFullFisher, V1.nScored, ...
    numel(V1.includedQuartet)));

axV4 = subplot(1, 2, 2);
plotPanel(axV4, V4, edges, xLimits, yMax, [0.75 0.35 0.18]);
title(axV4, sprintf('V4: %.1f%% (%d/%d), N = %d quartets', ...
    100 * V4.accuracyFullFisher, V4.nCorrectFullFisher, V4.nScored, ...
    numel(V4.includedQuartet)));

sgtitle(['Nilson quartet-centered full Fisher attention decoder, 300-500 ms; ' ...
    'RF centers on curves, at least 20 sites per quartet'], ...
    'FontWeight', 'bold');

figureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_V1_V4_fullFisher_quartetCentered_minSites20_N.png');
print(fig, figureFile, '-dpng', '-r180');

H = struct('figure', fig, 'axV1', axV1, 'axV4', axV4, ...
    'figureFile', figureFile, 'v1ResultFile', v1File, ...
    'v4ResultFile', v4File);
fprintf('Saved %s\n', figureFile);

end

function plotPanel(ax, result, edges, xLimits, yMax, color)
histogram(ax, result.SFullFisher, edges, 'Normalization', 'probability', ...
    'FaceColor', color, 'EdgeColor', 'w', 'FaceAlpha', 0.9);
hold(ax, 'on');
xline(ax, 0, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.8);
xlabel(ax, 'Attention score S');
ylabel(ax, 'Probability');
xlim(ax, xLimits);
ylim(ax, [0 yMax]);
box(ax, 'off');
grid(ax, 'on');
set(ax, 'FontSize', 12, 'Layer', 'top');
end
