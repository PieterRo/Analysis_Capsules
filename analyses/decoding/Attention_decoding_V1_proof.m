function OUT = Attention_decoding_V1_proof(varargin)
%ATTENTION_DECODING_V1_PROOF Trial-wise V1 attention read-out for Nilson.
%
% This is an in-sample proof-of-principle. Attention d-prime weights and
% trial scores are estimated from the same day-1/2 data. For each trial,
% sites whose RF centers lie on or near the target receive +d-prime and
% sites on or near the distractor receive -d-prime. Distant background and
% overlap sites are omitted.

p = inputParser;
p.addParameter('SNRThreshold', 0.7, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ResponseWindow', [300 500], ...
    @(x) isnumeric(x) && isequal(size(x), [1 2]) && x(1) < x(2));
p.addParameter('CovarianceWindow', [-200 0], ...
    @(x) isnumeric(x) && isequal(size(x), [1 2]) && x(1) < x(2));
p.addParameter('CovarianceShrinkage', 0.1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
p.addParameter('MinSitesPerQuartet', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == floor(x));
p.addParameter('CurveMarginDeg', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('PixelsPerDegree', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
p.addParameter('Days', [1 2], @(x) isnumeric(x) && isvector(x) && ~isempty(x));
p.addParameter('OnlyCorrect', true, @(x) islogical(x) && isscalar(x));
p.addParameter('ExcludeOverlap', true, @(x) islogical(x) && isscalar(x));
p.addParameter('ChunkTrials', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();

geomData = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'), 'Tall_V1');
respData = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'), 'R');
snrData = load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat'), 'SNR');
Tall_V1 = geomData.Tall_V1;
R3 = respData.R;
SNR = snrData.SNR;

if isempty(opt.PixelsPerDegree)
    calibrationFile = fullfile(cfg.extrasRoot, 'monkeyN', 'RFs', ...
        'BarMap_Nilson.mat');
    calibrationData = load(calibrationFile, 'pixperdeg');
    pixelsPerDegree = double(calibrationData.pixperdeg);
else
    calibrationFile = '';
    pixelsPerDegree = double(opt.PixelsPerDegree);
end
curveMarginPx = opt.CurveMarginDeg * pixelsPerDegree;

if opt.CurveMarginDeg > 0
    coordData = load(fullfile(cfg.logsDir, ...
        'ObjAtt_lines_monkeyN_20220201_B1.mat'), 'ALLCOORDS');
    ALLCOORDS = coordData.ALLCOORDS;
else
    ALLCOORDS = [];
end

assert(numel(Tall_V1) == 384, 'Expected geometry for 384 stimuli.');
assert(size(R3.meanAct, 2) == 384, 'Expected responses for 384 stimuli.');

timeIdx = find(all(abs(double(R3.timeWindows) - opt.ResponseWindow) < 1e-9, 2), 1);
if isempty(timeIdx)
    error('Requested response window %s is absent from SNR_capsules_N_d12.mat.', ...
        mat2str(opt.ResponseWindow));
end

attOpts = struct('v1Sites', 1:512, 'timeIdx', timeIdx, ...
    'excludeOverlap', opt.ExcludeOverlap, 'verbose', false);
attention = attention_modulation_V1_3bin(R3, Tall_V1, SNR, attOpts);

snrMatrix = [SNR.yellowEarly, SNR.yellowLate, ...
    SNR.purpleEarly, SNR.purpleLate];
bestSNR = max(snrMatrix, [], 2, 'omitnan');
visuallyDriven = isfinite(bestSNR) & bestSNR > opt.SNRThreshold;

dprime = double(attention.dprime(:));
responseMidpoint = 0.5 * (double(attention.muT(:)) + double(attention.muD(:)));
pooledAttentionSD = sqrt(0.5 * ...
    (double(attention.varT(:)) + double(attention.varD(:))));
baseline = double(SNR.muSpont(:));
topResponse = max([SNR.muYellowEarly(:), SNR.muYellowLate(:), ...
    SNR.muPurpleEarly(:), SNR.muPurpleLate(:)], [], 2);
responseScale = double(topResponse - baseline);
responseScale(~isfinite(responseScale) | responseScale <= 0) = NaN;

hasWeight = attention.validSite(:) & isfinite(dprime) & ...
    isfinite(responseMidpoint) & isfinite(responseScale) & ...
    isfinite(pooledAttentionSD) & pooledAttentionSD > 0;
eligibleSite = visuallyDriven & hasWeight;
fisherWeight = dprime ./ pooledAttentionSD;
eligibleSites = find(eligibleSite);
eligibleRank = zeros(512, 1);
eligibleRank(eligibleSites) = 1:numel(eligibleSites);

targetByStimulus = false(512, 384);
distractorByStimulus = false(512, 384);
overlapByStimulus = false(512, 384);
for stimNum = 1:384
    [targetByStimulus(:, stimNum), distractorByStimulus(:, stimNum), ...
        overlapByStimulus(:, stimNum)] = curveAssignments( ...
        Tall_V1(stimNum), ALLCOORDS, curveMarginPx);
end

dataDir = fullfile(cfg.dataRoot, 'Mr Nilson');
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:)';
[nChannels, nTrials, nTimes] = size(m1, 'normMUA');

assert(nChannels >= 512, 'normMUA has fewer than 512 V1 channels.');
assert(size(ALLMAT, 1) == nTrials, 'ALLMAT and normMUA trial counts differ.');
assert(numel(tb) == nTimes, 'tb and normMUA time dimensions differ.');
assert(size(ALLMAT, 2) >= 11, 'Expected the 11-column Nilson ALLMAT format.');

stimPerTrial = ALLMAT(:, 1);
includeTrial = ismember(ALLMAT(:, 11), opt.Days(:));
if opt.OnlyCorrect
    includeTrial = includeTrial & ALLMAT(:, 9) == 1;
end
includeTrial = includeTrial & isfinite(stimPerTrial) & ...
    stimPerTrial >= 1 & stimPerTrial <= 384 & stimPerTrial == floor(stimPerTrial);

timeMask = tb >= opt.ResponseWindow(1) & tb <= opt.ResponseWindow(2);
if ~any(timeMask)
    error('Response window %s does not overlap tb.', mat2str(opt.ResponseWindow));
end
timeSamples = find(timeMask);
covarianceMask = tb >= opt.CovarianceWindow(1) & tb <= opt.CovarianceWindow(2);
if ~any(covarianceMask)
    error('Covariance window %s does not overlap tb.', mat2str(opt.CovarianceWindow));
end
covarianceSamples = find(covarianceMask);

S = nan(nTrials, 1);
SFisher = nan(nTrials, 1);
SFullFisher = nan(nTrials, 1);
nSitesUsed = zeros(nTrials, 1);
nTargetSites = zeros(nTrials, 1);
nDistractorSites = zeros(nTrials, 1);
sumAbsWeight = nan(nTrials, 1);
sumAbsFisherWeight = nan(nTrials, 1);
centeredResponseByTrial = nan(512, nTrials);
prestimResponseByTrial = nan(numel(eligibleSites), nTrials);

fprintf(['V1 attention proof-of-principle: reading %d trials in chunks of %d ' ...
    '(%g-%g ms).\n'], nTrials, opt.ChunkTrials, tb(timeSamples(1)), tb(timeSamples(end)));

for firstTrial = 1:opt.ChunkTrials:nTrials
    lastTrial = min(nTrials, firstTrial + opt.ChunkTrials - 1);
    trialRange = firstTrial:lastTrial;
    localInclude = includeTrial(trialRange);
    if ~any(localInclude)
        continue;
    end

    raw = double(m1.normMUA(1:512, trialRange, timeSamples));
    trialResponse = mean(raw, 3, 'omitnan');
    rawPrestim = double(m1.normMUA(1:512, trialRange, covarianceSamples));
    prestimResponse = mean(rawPrestim, 3, 'omitnan');
    prestimResponse = prestimResponse(eligibleSites, :);
    prestimResponse = (prestimResponse - baseline(eligibleSites)) ./ ...
        responseScale(eligibleSites);
    prestimResponseByTrial(:, trialRange(localInclude)) = ...
        prestimResponse(:, localInclude);

    localTrials = find(localInclude);
    for jj = localTrials(:)'
        globalTrial = trialRange(jj);
        stimNum = stimPerTrial(globalTrial);
        isTarget = targetByStimulus(:, stimNum);
        isDistractor = distractorByStimulus(:, stimNum);
        onCurve = isTarget | isDistractor;

        if opt.ExcludeOverlap
            notOverlap = ~overlapByStimulus(:, stimNum);
        else
            notOverlap = true(512, 1);
        end

        normalizedResponse = (trialResponse(:, jj) - baseline) ./ responseScale;
        centeredResponse = normalizedResponse - responseMidpoint;
        centeredResponseByTrial(:, globalTrial) = centeredResponse;
        use = eligibleSite & onCurve & notOverlap & isfinite(centeredResponse);

        role = zeros(512, 1);
        role(isTarget) = 1;
        role(isDistractor) = -1;
        effectiveWeight = dprime .* role;
        denominator = sum(abs(dprime(use)));
        effectiveFisherWeight = fisherWeight .* role;
        fisherDenominator = sum(abs(fisherWeight(use)));

        if denominator > 0 && fisherDenominator > 0
            S(globalTrial) = sum(effectiveWeight(use) .* centeredResponse(use)) / denominator;
            SFisher(globalTrial) = ...
                sum(effectiveFisherWeight(use) .* centeredResponse(use)) / fisherDenominator;
            nSitesUsed(globalTrial) = nnz(use);
            nTargetSites(globalTrial) = nnz(use & isTarget);
            nDistractorSites(globalTrial) = nnz(use & isDistractor);
            sumAbsWeight(globalTrial) = denominator;
            sumAbsFisherWeight(globalTrial) = fisherDenominator;
        end
    end

    fprintf('  trials %d-%d of %d\n', firstTrial, lastTrial, nTrials);
end

covarianceTrial = includeTrial & all(isfinite(prestimResponseByTrial), 1)';
nCovarianceTrials = nnz(covarianceTrial);
if nCovarianceTrials <= numel(eligibleSites)
    error(['Only %d complete prestimulus trials for %d sites; cannot estimate ' ...
        'the full covariance robustly.'], nCovarianceTrials, numel(eligibleSites));
end

prestimCovariance = cov(prestimResponseByTrial(:, covarianceTrial)');
covarianceDiagonal = diag(diag(prestimCovariance));
regularizedCovariance = ...
    (1 - opt.CovarianceShrinkage) * prestimCovariance + ...
    opt.CovarianceShrinkage * covarianceDiagonal;
covarianceCondition = cond(prestimCovariance);
regularizedCovarianceCondition = cond(regularizedCovariance);

deltaAttention = double(attention.muT(:) - attention.muD(:));
fullWeightByStimulus = cell(384, 1);
fullSiteByStimulus = cell(384, 1);
fullDenominatorByStimulus = nan(384, 1);
siteCountByStimulus = zeros(384, 1);

for stimNum = 1:384
    isTarget = targetByStimulus(:, stimNum);
    isDistractor = distractorByStimulus(:, stimNum);
    onCurve = isTarget | isDistractor;

    if opt.ExcludeOverlap
        notOverlap = ~overlapByStimulus(:, stimNum);
    else
        notOverlap = true(512, 1);
    end

    use = eligibleSite & onCurve & notOverlap;
    siteIdx = find(use);
    if isempty(siteIdx)
        continue;
    end

    covarianceIdx = eligibleRank(siteIdx);
    role = ones(numel(siteIdx), 1);
    role(isDistractor(siteIdx)) = -1;
    alignedCovariance = regularizedCovariance(covarianceIdx, covarianceIdx) .* ...
        (role * role');
    alignedWeight = alignedCovariance \ deltaAttention(siteIdx);

    fullSiteByStimulus{stimNum} = siteIdx;
    fullWeightByStimulus{stimNum} = alignedWeight .* role;
    fullDenominatorByStimulus(stimNum) = sum(abs(alignedWeight));
    siteCountByStimulus(stimNum) = numel(siteIdx);
end

quartetSiteCount = nan(96, 1);
for quartetIdx = 1:96
    members = quartetMembers(quartetIdx);
    memberCounts = siteCountByStimulus(members);
    if any(memberCounts ~= memberCounts(1))
        error('Curve-site count is not constant within quartet %d.', quartetIdx);
    end
    quartetSiteCount(quartetIdx) = memberCounts(1);
end

for globalTrial = find(includeTrial)'
    stimNum = stimPerTrial(globalTrial);
    siteIdx = fullSiteByStimulus{stimNum};
    weight = fullWeightByStimulus{stimNum};
    denominator = fullDenominatorByStimulus(stimNum);
    if isempty(siteIdx) || ~isfinite(denominator) || denominator <= 0
        continue;
    end

    response = centeredResponseByTrial(siteIdx, globalTrial);
    finiteResponse = isfinite(response) & isfinite(weight);
    if ~any(finiteResponse)
        continue;
    end
    SFullFisher(globalTrial) = ...
        sum(weight(finiteResponse) .* response(finiteResponse)) / ...
        sum(abs(weight(finiteResponse)));
end

finiteScoreTrial = includeTrial & isfinite(S) & isfinite(SFisher) & ...
    isfinite(SFullFisher);
trialQuartet = stimulusToQuartet(stimPerTrial);
coverageTrial = false(nTrials, 1);
validQuartet = trialQuartet >= 1 & trialQuartet <= 96;
coverageTrial(validQuartet) = ...
    quartetSiteCount(trialQuartet(validQuartet)) >= opt.MinSitesPerQuartet;
scoredTrial = finiteScoreTrial & coverageTrial;
nScored = nnz(scoredTrial);
if nScored == 0
    error('No trials received a finite decoder score.');
end
unscoredTrial = includeTrial & ~scoredTrial;
noFiniteScoreTrial = includeTrial & ~finiteScoreTrial;
belowCoverageTrial = finiteScoreTrial & ~coverageTrial;

nCorrect = nnz(S(scoredTrial) > 0);
accuracy = nCorrect / nScored;
nCorrectFisher = nnz(SFisher(scoredTrial) > 0);
accuracyFisher = nCorrectFisher / nScored;
nCorrectFullFisher = nnz(SFullFisher(scoredTrial) > 0);
accuracyFullFisher = nCorrectFullFisher / nScored;

stimulus = stimPerTrial(scoredTrial);
quartet = stimulusToQuartet(stimulus);

OUT = struct();
OUT.description = ['In-sample proof-of-principle comparing d-prime, diagonal ' ...
    'Fisher, and prestimulus full-covariance Fisher weights.'];
OUT.monkey = 'Mr Nilson';
OUT.region = 'V1';
OUT.responseWindowRequested = opt.ResponseWindow;
OUT.responseWindowSampled = [tb(timeSamples(1)), tb(timeSamples(end))];
OUT.covarianceWindowRequested = opt.CovarianceWindow;
OUT.covarianceWindowSampled = ...
    [tb(covarianceSamples(1)), tb(covarianceSamples(end))];
OUT.covarianceShrinkage = opt.CovarianceShrinkage;
OUT.curveMarginDeg = opt.CurveMarginDeg;
OUT.curveMarginPx = curveMarginPx;
OUT.pixelsPerDegree = pixelsPerDegree;
OUT.pixelCalibrationFile = calibrationFile;
OUT.minSitesPerQuartet = opt.MinSitesPerQuartet;
OUT.quartetSiteCount = quartetSiteCount;
OUT.includedQuartet = find(quartetSiteCount >= opt.MinSitesPerQuartet);
OUT.nCovarianceTrials = nCovarianceTrials;
OUT.prestimCovariance = prestimCovariance;
OUT.regularizedCovariance = regularizedCovariance;
OUT.covarianceCondition = covarianceCondition;
OUT.regularizedCovarianceCondition = regularizedCovarianceCondition;
OUT.days = opt.Days(:)';
OUT.onlyCorrect = opt.OnlyCorrect;
OUT.snrThreshold = opt.SNRThreshold;
OUT.excludeOverlap = opt.ExcludeOverlap;
OUT.visualSiteMask = visuallyDriven;
OUT.eligibleSiteMask = eligibleSite;
OUT.bestSNR = bestSNR;
OUT.attentionDprime = dprime;
OUT.pooledAttentionSD = pooledAttentionSD;
OUT.fisherWeight = fisherWeight;
OUT.responseMidpoint = responseMidpoint;
OUT.responseScale = responseScale;
OUT.nIncluded = nnz(includeTrial);
OUT.trialIndex = find(scoredTrial);
OUT.stimulus = stimulus;
OUT.quartet = quartet;
OUT.S = S(scoredTrial);
OUT.SFisher = SFisher(scoredTrial);
OUT.SFullFisher = SFullFisher(scoredTrial);
OUT.nSitesUsed = nSitesUsed(scoredTrial);
OUT.nTargetSites = nTargetSites(scoredTrial);
OUT.nDistractorSites = nDistractorSites(scoredTrial);
OUT.sumAbsWeight = sumAbsWeight(scoredTrial);
OUT.sumAbsFisherWeight = sumAbsFisherWeight(scoredTrial);
OUT.nScored = nScored;
OUT.nUnscored = nnz(unscoredTrial);
OUT.nWithoutFiniteScore = nnz(noFiniteScoreTrial);
OUT.nBelowCoverageThreshold = nnz(belowCoverageTrial);
OUT.belowCoverageTrialIndex = find(belowCoverageTrial);
OUT.belowCoverageStimulus = stimPerTrial(belowCoverageTrial);
OUT.belowCoverageQuartet = stimulusToQuartet(OUT.belowCoverageStimulus);
OUT.unscoredTrialIndex = find(unscoredTrial);
OUT.unscoredStimulus = stimPerTrial(unscoredTrial);
OUT.unscoredQuartet = stimulusToQuartet(OUT.unscoredStimulus);
OUT.nCorrect = nCorrect;
OUT.accuracy = accuracy;
OUT.nCorrectFisher = nCorrectFisher;
OUT.accuracyFisher = accuracyFisher;
OUT.nCorrectFullFisher = nCorrectFullFisher;
OUT.accuracyFullFisher = accuracyFullFisher;

fprintf('\nV1 attention proof-of-principle\n');
fprintf('  Visually driven sites (bestSNR > %.2f): %d / 512\n', ...
    opt.SNRThreshold, nnz(visuallyDriven));
fprintf('  Sites with a finite attention weight: %d / 512\n', nnz(eligibleSite));
fprintf('  Curve margin: %.2f deg (%.2f px beyond the capsule edge)\n', ...
    opt.CurveMarginDeg, curveMarginPx);
fprintf('  Scored trials: %d / %d included (%d without eligible curve sites)\n', ...
    nScored, OUT.nIncluded, OUT.nWithoutFiniteScore);
fprintf('  Minimum sites per quartet: %d (%d trials below threshold)\n', ...
    opt.MinSitesPerQuartet, OUT.nBelowCoverageThreshold);
if OUT.nUnscored > 0
    fprintf('  Unscored quartets: %s\n', mat2str(unique(OUT.unscoredQuartet(:))'));
end
fprintf('  Median sites per trial: %.1f (target %.1f, distractor %.1f)\n', ...
    median(OUT.nSitesUsed), median(OUT.nTargetSites), median(OUT.nDistractorSites));
fprintf('  d-prime correct (S > 0): %d / %d = %.2f%%\n', ...
    nCorrect, nScored, 100 * accuracy);
fprintf('  Fisher correct (S > 0): %d / %d = %.2f%%\n', ...
    nCorrectFisher, nScored, 100 * accuracyFisher);
fprintf(['  Full-covariance Fisher correct (S > 0): %d / %d = %.2f%% ' ...
    '(shrinkage %.2f)\n'], nCorrectFullFisher, nScored, ...
    100 * accuracyFullFisher, opt.CovarianceShrinkage);
fprintf('  Prestimulus covariance trials: %d; condition %.3g -> %.3g\n', ...
    nCovarianceTrials, covarianceCondition, regularizedCovarianceCondition);

if opt.MakeFigure
    fig = figure('Color', 'w', 'Position', [100 100 1580 500]);
    subplot(1, 3, 1);
    histogram(OUT.S, 'BinMethod', 'fd', 'Normalization', 'probability', ...
        'FaceColor', [0.35 0.43 0.52], 'EdgeColor', 'w', 'FaceAlpha', 0.9);
    hold on;
    xline(0, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.8);
    xlabel('Attention score S');
    ylabel('Probability');
    title(sprintf('d-prime: P(S > 0) = %.1f%% (%d/%d)', ...
        100 * accuracy, nCorrect, nScored));
    box off;
    grid on;
    set(gca, 'FontSize', 12, 'Layer', 'top');

    subplot(1, 3, 2);
    histogram(OUT.SFisher, 'BinMethod', 'fd', 'Normalization', 'probability', ...
        'FaceColor', [0.24 0.52 0.38], 'EdgeColor', 'w', 'FaceAlpha', 0.9);
    hold on;
    xline(0, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.8);
    xlabel('Attention score S');
    ylabel('Probability');
    title(sprintf('Fisher: P(S > 0) = %.1f%% (%d/%d)', ...
        100 * accuracyFisher, nCorrectFisher, nScored));
    box off;
    grid on;
    set(gca, 'FontSize', 12, 'Layer', 'top');

    subplot(1, 3, 3);
    histogram(OUT.SFullFisher, 'BinMethod', 'fd', 'Normalization', 'probability', ...
        'FaceColor', [0.48 0.32 0.62], 'EdgeColor', 'w', 'FaceAlpha', 0.9);
    hold on;
    xline(0, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.8);
    xlabel('Attention score S');
    ylabel('Probability');
    title(sprintf('Full Fisher: P(S > 0) = %.1f%% (%d/%d)', ...
        100 * accuracyFullFisher, nCorrectFullFisher, nScored));
    box off;
    grid on;
    set(gca, 'FontSize', 12, 'Layer', 'top');

    sgtitle(sprintf(['Nilson V1 attention read-out, %g-%g ms; %.1f deg margin; ' ...
        'quartets with at least %d sites (N = %d)'], ...
        opt.ResponseWindow(1), opt.ResponseWindow(2), ...
        opt.CurveMarginDeg, opt.MinSitesPerQuartet, ...
        numel(OUT.includedQuartet)), ...
        'FontWeight', 'bold');
    OUT.figure = fig;
else
    fig = [];
end

if opt.SaveOutputs
    if exist(cfg.resultsDir, 'dir') ~= 7
        mkdir(cfg.resultsDir);
    end
    thresholdLabel = sprintf('minSites%d', opt.MinSitesPerQuartet);
    marginLabel = strrep(sprintf('%g', opt.CurveMarginDeg), '.', 'p');
    analysisLabel = sprintf('margin%sdeg_%s', marginLabel, thresholdLabel);
    resultFile = fullfile(cfg.resultsDir, sprintf( ...
        'Attention_decoder_V1_proof_%s_N.mat', analysisLabel));
    figureFile = fullfile(cfg.resultsDir, sprintf( ...
        'Attention_decoder_V1_fisher_comparison_%s_N.png', analysisLabel));
    OUT.resultFile = resultFile;
    OUT.figureFile = figureFile;

    if isfield(OUT, 'figure')
        figureHandle = OUT.figure;
        OUT = rmfield(OUT, 'figure');
    else
        figureHandle = [];
    end
    save(resultFile, 'OUT', '-v7.3');
    if ~isempty(figureHandle)
        OUT.figure = figureHandle;
    end
    if ~isempty(fig)
        print(fig, figureFile, '-dpng', '-r180');
    end
    fprintf('  Saved %s\n', resultFile);
    if ~isempty(fig)
        fprintf('  Saved %s\n', figureFile);
    end
end

end

function quartet = stimulusToQuartet(stimulus)
block = floor((stimulus - 1) / 8);
position = mod(stimulus - 1, 8) + 1;
quartet = 2 * block + 1 + ismember(position, [3 4 7 8]);
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

function [isTarget, isDistractor, overlap] = curveAssignments( ...
    stimulusGeometry, ALLCOORDS, marginPx)
T = stimulusGeometry.T;
assignment = string(T.assignment(1:512));
isTarget = assignment == "target";
isDistractor = assignment == "distractor";
overlap = logical(T.overlap(1:512));

if marginPx == 0
    return;
end

fieldName = sprintf('stim_%d', stimulusGeometry.stimNum);
coords = ALLCOORDS.(fieldName);
toPx = @(point) [point(1) + 512, 384 - point(2)];
s = toPx(double(coords.s(:))');
tTarget = toPx(double(coords.t_fig(:))');
tDistractor = toPx(double(coords.t_back(:))');
radiusPx = double(stimulusGeometry.widthPx) / 2 + marginPx;
points = [double(T.x_px(1:512)), double(T.y_px(1:512))];

isTarget = pointSegmentDistance(points, s, tTarget) <= radiusPx;
isDistractor = pointSegmentDistance(points, s, tDistractor) <= radiusPx;
overlap = isTarget & isDistractor;
end

function distance = pointSegmentDistance(points, segmentStart, segmentEnd)
segment = segmentEnd - segmentStart;
segmentLengthSquared = sum(segment .^ 2);
if segmentLengthSquared == 0
    projection = repmat(segmentStart, size(points, 1), 1);
else
    fraction = ((points - segmentStart) * segment') / segmentLengthSquared;
    fraction = max(0, min(1, fraction));
    projection = segmentStart + fraction .* segment;
end
distance = hypot(points(:, 1) - projection(:, 1), ...
    points(:, 2) - projection(:, 2));
end
