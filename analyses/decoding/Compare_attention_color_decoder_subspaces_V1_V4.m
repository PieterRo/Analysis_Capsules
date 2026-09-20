function OUT = Compare_attention_color_decoder_subspaces_V1_V4(varargin)
%COMPARE_ATTENTION_COLOR_DECODER_SUBSPACES_V1_V4 Compare decoder axes.
%
% Reconstructs each decoder's effective weight vector in the original
% neural-response coordinates for every stimulus. The absolute cosine and
% principal angle quantify overlap between the one-dimensional attention
% and color decoder subspaces. Both decoders use 300-500 ms activity.

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
featureFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');
attentionData = load(attentionFile, 'OUT');
colorData = load(colorFile, 'OUT');
featureData = load(featureFile, 'OUT');
attention = attentionData.OUT;
color = colorData.OUT;
feature = featureData.OUT;

assert(isequal(double(attention.responseWindow), [300 500]) && ...
    isequal(double(color.responseWindow), [300 500]), ...
    'Attention and color decoders must both use 300-500 ms activity.');

OUT = struct();
OUT.description = ['Stimulus-specific principal angles between final ' ...
    'attention and color decoder axes in original neural-response space.'];
OUT.monkey = attention.monkey;
OUT.responseWindow = [300 500];
OUT.attentionSourceFile = attentionFile;
OUT.colorSourceFile = colorFile;
OUT.featureSourceFile = featureFile;
OUT.angleDefinition = ['acosd(abs(cosine)) after applying each decoder''s ' ...
    'stimulus-specific sign and RF availability mask; 90 degrees is orthogonal.'];
OUT.V1 = analyzeRegion('V1', attention.V1, color.V1, feature.V1);
OUT.V4 = analyzeRegion('V4', attention.V4, color.V4, feature.V4);

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

resultFile = fullfile(cfg.resultsDir, ...
    'Attention_color_decoder_subspace_angles_V1_V4_N.mat');
figureFile = fullfile(cfg.resultsDir, ...
    'Attention_color_decoder_subspace_angles_V1_V4_N.png');
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

function R = analyzeRegion(region, attention, color, feature)
assert(isequal(attention.trialIndex, color.trialIndex) && ...
    isequal(attention.trialIndex, feature.trialIndex), ...
    '%s trial ordering differs between sources.', region);
assert(isequal(attention.stimulus, color.stimulus) && ...
    isequal(attention.stimulus, feature.stimulus), ...
    '%s stimulus ordering differs between sources.', region);
assert(isequal(attention.siteIndexGlobal, color.siteIndexGlobal) && ...
    isequal(attention.siteIndexGlobal, feature.siteIndexGlobal), ...
    '%s site ordering differs between sources.', region);

attentionCoefficient = double(attention.finalCoefficient(:));
colorCoefficient = double(color.finalCoefficient(:));
coefficientCosine = cosineSimilarity( ...
    attentionCoefficient, colorCoefficient);

stimulusValues = unique(attention.stimulus(:), 'sorted');
nStimuli = numel(stimulusValues);
signedCosine = nan(nStimuli, 1);
absoluteCosine = nan(nStimuli, 1);
principalAngleDeg = nan(nStimuli, 1);
nAvailableSites = zeros(nStimuli, 1);
nJointNonzeroSites = zeros(nStimuli, 1);

for stimulusIdx = 1:nStimuli
    stimulus = stimulusValues(stimulusIdx);
    row = find(attention.stimulus == stimulus, 1);
    attentionSign = double(feature.canonicalSign(row, :))';
    colorSign = double(color.classLabel(row) .* color.colorRole(row, :))';
    attentionWeight = attentionCoefficient .* attentionSign;
    colorWeight = colorCoefficient .* colorSign;
    nAvailableSites(stimulusIdx) = nnz( ...
        attentionSign ~= 0 & colorSign ~= 0);
    nJointNonzeroSites(stimulusIdx) = nnz( ...
        attentionWeight ~= 0 & colorWeight ~= 0);
    signedCosine(stimulusIdx) = cosineSimilarity( ...
        attentionWeight, colorWeight);
    absoluteCosine(stimulusIdx) = abs(signedCosine(stimulusIdx));
    principalAngleDeg(stimulusIdx) = acosd(absoluteCosine(stimulusIdx));
end

assert(all(isfinite(principalAngleDeg)), ...
    '%s contains undefined stimulus-specific angles.', region);

R = struct();
R.region = region;
R.stimulus = stimulusValues;
R.attentionCoefficient = attentionCoefficient;
R.colorCoefficient = colorCoefficient;
R.coefficientCosine = coefficientCosine;
R.coefficientPrincipalAngleDeg = acosd(abs(coefficientCosine));
R.signedCosine = signedCosine;
R.absoluteCosine = absoluteCosine;
R.principalAngleDeg = principalAngleDeg;
R.nAvailableSites = nAvailableSites;
R.nJointNonzeroSites = nJointNonzeroSites;
R.nCandidateSites = numel(attentionCoefficient);
R.nAttentionSites = nnz(attentionCoefficient);
R.nColorSites = nnz(colorCoefficient);
R.nOverlappingNonzeroSites = nnz( ...
    attentionCoefficient ~= 0 & colorCoefficient ~= 0);
R.medianAbsoluteCosine = median(absoluteCosine);
R.meanAbsoluteCosine = mean(absoluteCosine);
R.medianSharedVariance = median(absoluteCosine .^ 2);
R.medianPrincipalAngleDeg = median(principalAngleDeg);
R.principalAngleIqrDeg = prctile(principalAngleDeg, [25 75]);
R.principalAngleRangeDeg = [min(principalAngleDeg), max(principalAngleDeg)];
end

function cosine = cosineSimilarity(x, y)
denominator = norm(x) * norm(y);
assert(denominator > 0, 'Cannot compare a zero-length decoder axis.');
cosine = dot(x, y) / denominator;
cosine = max(-1, min(1, cosine));
end

function fig = makeFigure(OUT)
fig = figure('Color', 'w', 'Position', [100 100 1160 520]);
regions = {'V1', 'V4'};
colors = [0.20 0.48 0.72; 0.78 0.42 0.16];
edges = 65:2:91;

for regionIdx = 1:numel(regions)
    R = OUT.(regions{regionIdx});
    ax = subplot(1, 2, regionIdx);
    histogram(ax, R.principalAngleDeg, edges, ...
        'FaceColor', colors(regionIdx, :), 'EdgeColor', 'w', ...
        'FaceAlpha', 0.88);
    hold(ax, 'on');
    medianLine = xline(ax, R.medianPrincipalAngleDeg, '-', ...
        'Color', [0.78 0.12 0.12], 'LineWidth', 2.2);
    orthogonalLine = xline(ax, 90, ':', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1.5);
    xlim(ax, [65 91]);
    xlabel(ax, 'Attention-color principal angle (degrees)');
    ylabel(ax, 'Number of stimuli');
    title(ax, sprintf(['%s: median %.1f deg, median cos^2 = %.1f%%\n' ...
        '%d overlapping nonzero sites'], regionIdxLabel(regions{regionIdx}), ...
        R.medianPrincipalAngleDeg, 100 * R.medianSharedVariance, ...
        R.nOverlappingNonzeroSites));
    legend(ax, [medianLine orthogonalLine], {'Median', 'Orthogonal'}, ...
        'Location', 'northwest');
    box(ax, 'off');
    set(ax, 'FontSize', 11, 'LineWidth', 1, 'TickDir', 'out');
end

sgtitle(['Attention versus color decoder subspaces, 300-500 ms; ' ...
    'stimulus-specific effective weights'], 'FontWeight', 'bold');
end

function label = regionIdxLabel(region)
label = region;
end

function printSummary(OUT)
fprintf('\nAttention-color decoder subspace comparison, 300-500 ms\n');
for region = {'V1', 'V4'}
    R = OUT.(region{1});
    fprintf(['  %s: median angle %.2f deg, IQR [%.2f %.2f], ' ...
        'median cos^2 %.2f%%; coefficient angle %.2f deg; ' ...
        '%d shared nonzero sites.\n'], region{1}, ...
        R.medianPrincipalAngleDeg, R.principalAngleIqrDeg(1), ...
        R.principalAngleIqrDeg(2), 100 * R.medianSharedVariance, ...
        R.coefficientPrincipalAngleDeg, R.nOverlappingNonzeroSites);
end
end
