function OUT = Attention_decoding_V4_proof(varargin)
%ATTENTION_DECODING_V4_PROOF In-sample full-Fisher attention read-out.
%
% Uses visually driven Nilson V4 sites whose RF centers lie on the target
% or distractor curve. The full covariance is estimated from prestimulus
% activity. Attention weights and trial scores use the same day-1/2 data.

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
globalSites = (513:768).';
nSites = numel(globalSites);
nStimuli = 384;
nQuartets = 96;

geomData = load(fullfile(cfg.matDir, 'Tall_V4_lines_N.mat'), 'Tall_V4');
respData = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'), 'R');
Tall_V4 = geomData.Tall_V4;
R3 = respData.R;

assert(numel(Tall_V4) == nStimuli, 'Expected geometry for 384 stimuli.');
assert(height(Tall_V4(1).T) == nSites, ...
    'Tall_V4 rows must correspond to global channels 513:768.');
assert(size(R3.meanAct, 1) >= globalSites(end), ...
    'Response summary has fewer than 768 channels.');
assert(size(R3.meanAct, 2) == nStimuli, ...
    'Expected responses for 384 stimuli.');

timeIdx = find(all(abs(double(R3.timeWindows) - opt.ResponseWindow) < 1e-9, 2), 1);
if isempty(timeIdx)
    error('Requested response window %s is absent from SNR_capsules_N_d12.mat.', ...
        mat2str(opt.ResponseWindow));
end

SNR = compute_snr_per_color_region(R3, Tall_V4, globalSites);
R3V4 = subsetResponse(R3, globalSites);
attOpts = struct('v1Sites', 1:nSites, 'timeIdx', timeIdx, ...
    'excludeOverlap', opt.ExcludeOverlap, 'verbose', false);
attention = attention_modulation_V1_3bin(R3V4, Tall_V4, SNR, attOpts);

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
eligibleSites = find(eligibleSite);
eligibleRank = zeros(nSites, 1);
eligibleRank(eligibleSites) = 1:numel(eligibleSites);

targetByStimulus = false(nSites, nStimuli);
distractorByStimulus = false(nSites, nStimuli);
overlapByStimulus = false(nSites, nStimuli);
for stimNum = 1:nStimuli
    T = Tall_V4(stimNum).T;
    assignment = string(T.assignment(1:nSites));
    targetByStimulus(:, stimNum) = assignment == "target";
    distractorByStimulus(:, stimNum) = assignment == "distractor";
    if ismember('overlap', T.Properties.VariableNames)
        overlapByStimulus(:, stimNum) = logical(T.overlap(1:nSites));
    end
end

dataDir = fullfile(cfg.dataRoot, 'Mr Nilson');
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:)';
[nChannels, nTrials, nTimes] = size(m1, 'normMUA');

assert(nChannels >= globalSites(end), 'normMUA has fewer than 768 channels.');
assert(size(ALLMAT, 1) == nTrials, 'ALLMAT and normMUA trial counts differ.');
assert(numel(tb) == nTimes, 'tb and normMUA time dimensions differ.');
assert(size(ALLMAT, 2) >= 11, 'Expected the 11-column Nilson ALLMAT format.');

stimPerTrial = ALLMAT(:, 1);
includeTrial = ismember(ALLMAT(:, 11), opt.Days(:));
if opt.OnlyCorrect
    includeTrial = includeTrial & ALLMAT(:, 9) == 1;
end
includeTrial = includeTrial & isfinite(stimPerTrial) & ...
    stimPerTrial >= 1 & stimPerTrial <= nStimuli & ...
    stimPerTrial == floor(stimPerTrial);

timeMask = tb >= opt.ResponseWindow(1) & tb <= opt.ResponseWindow(2);
if ~any(timeMask)
    error('Response window %s does not overlap tb.', mat2str(opt.ResponseWindow));
end
timeSamples = find(timeMask);
covarianceMask = tb >= opt.CovarianceWindow(1) & tb <= opt.CovarianceWindow(2);
if ~any(covarianceMask)
    error('Covariance window %s does not overlap tb.', ...
        mat2str(opt.CovarianceWindow));
end
covarianceSamples = find(covarianceMask);

centeredResponseByTrial = nan(nSites, nTrials);
prestimResponseByTrial = nan(numel(eligibleSites), nTrials);

fprintf(['V4 full-Fisher attention decoder: reading %d trials in chunks of %d ' ...
    '(%g-%g ms).\n'], nTrials, opt.ChunkTrials, ...
    tb(timeSamples(1)), tb(timeSamples(end)));

for firstTrial = 1:opt.ChunkTrials:nTrials
    lastTrial = min(nTrials, firstTrial + opt.ChunkTrials - 1);
    trialRange = firstTrial:lastTrial;
    localInclude = includeTrial(trialRange);
    if ~any(localInclude)
        continue;
    end

    raw = double(m1.normMUA(globalSites, trialRange, timeSamples));
    trialResponse = mean(raw, 3, 'omitnan');
    normalizedResponse = (trialResponse - baseline) ./ responseScale;
    centeredResponseByTrial(:, trialRange(localInclude)) = ...
        normalizedResponse(:, localInclude) - responseMidpoint;

    rawPrestim = double(m1.normMUA(globalSites, trialRange, covarianceSamples));
    prestimResponse = mean(rawPrestim, 3, 'omitnan');
    prestimResponse = prestimResponse(eligibleSites, :);
    prestimResponse = (prestimResponse - baseline(eligibleSites)) ./ ...
        responseScale(eligibleSites);
    prestimResponseByTrial(:, trialRange(localInclude)) = ...
        prestimResponse(:, localInclude);

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
siteByStimulus = cell(nStimuli, 1);
weightByStimulus = cell(nStimuli, 1);
denominatorByStimulus = nan(nStimuli, 1);
siteCountByStimulus = zeros(nStimuli, 1);

for stimNum = 1:nStimuli
    isTarget = targetByStimulus(:, stimNum);
    isDistractor = distractorByStimulus(:, stimNum);
    onCurve = isTarget | isDistractor;
    if opt.ExcludeOverlap
        notOverlap = ~overlapByStimulus(:, stimNum);
    else
        notOverlap = true(nSites, 1);
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

    siteByStimulus{stimNum} = siteIdx;
    weightByStimulus{stimNum} = alignedWeight .* role;
    denominatorByStimulus(stimNum) = sum(abs(alignedWeight));
    siteCountByStimulus(stimNum) = numel(siteIdx);
end

quartetSiteCount = nan(nQuartets, 1);
for quartetIdx = 1:nQuartets
    members = quartetMembers(quartetIdx);
    memberCounts = siteCountByStimulus(members);
    if any(memberCounts ~= memberCounts(1))
        error('Curve-site count is not constant within quartet %d.', quartetIdx);
    end
    quartetSiteCount(quartetIdx) = memberCounts(1);
end

S = nan(nTrials, 1);
nSitesUsed = zeros(nTrials, 1);
nTargetSites = zeros(nTrials, 1);
nDistractorSites = zeros(nTrials, 1);
sumAbsWeight = nan(nTrials, 1);

for globalTrial = find(includeTrial)'
    stimNum = stimPerTrial(globalTrial);
    siteIdx = siteByStimulus{stimNum};
    weight = weightByStimulus{stimNum};
    denominator = denominatorByStimulus(stimNum);
    if isempty(siteIdx) || ~isfinite(denominator) || denominator <= 0
        continue;
    end

    response = centeredResponseByTrial(siteIdx, globalTrial);
    finiteResponse = isfinite(response) & isfinite(weight);
    if ~any(finiteResponse)
        continue;
    end

    usedSites = siteIdx(finiteResponse);
    S(globalTrial) = sum(weight(finiteResponse) .* response(finiteResponse)) / ...
        sum(abs(weight(finiteResponse)));
    nSitesUsed(globalTrial) = numel(usedSites);
    nTargetSites(globalTrial) = ...
        nnz(targetByStimulus(usedSites, stimNum));
    nDistractorSites(globalTrial) = ...
        nnz(distractorByStimulus(usedSites, stimNum));
    sumAbsWeight(globalTrial) = sum(abs(weight(finiteResponse)));
end

finiteScoreTrial = includeTrial & isfinite(S);
trialQuartet = stimulusToQuartet(stimPerTrial);
coverageTrial = false(nTrials, 1);
validQuartet = trialQuartet >= 1 & trialQuartet <= nQuartets;
coverageTrial(validQuartet) = ...
    quartetSiteCount(trialQuartet(validQuartet)) >= opt.MinSitesPerQuartet;
scoredTrial = finiteScoreTrial & coverageTrial;
nScored = nnz(scoredTrial);
if nScored == 0
    error('No trials received a finite decoder score.');
end

noFiniteScoreTrial = includeTrial & ~finiteScoreTrial;
belowCoverageTrial = finiteScoreTrial & ~coverageTrial;
unscoredTrial = includeTrial & ~scoredTrial;
nCorrect = nnz(S(scoredTrial) > 0);
accuracy = nCorrect / nScored;
stimulus = stimPerTrial(scoredTrial);
quartet = stimulusToQuartet(stimulus);

OUT = struct();
OUT.description = ['In-sample V4 attention decoder using prestimulus ' ...
    'full-covariance Fisher weights.'];
OUT.monkey = 'Mr Nilson';
OUT.region = 'V4';
OUT.globalSites = globalSites;
OUT.responseWindowRequested = opt.ResponseWindow;
OUT.responseWindowSampled = [tb(timeSamples(1)), tb(timeSamples(end))];
OUT.covarianceWindowRequested = opt.CovarianceWindow;
OUT.covarianceWindowSampled = ...
    [tb(covarianceSamples(1)), tb(covarianceSamples(end))];
OUT.covarianceShrinkage = opt.CovarianceShrinkage;
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
OUT.responseMidpoint = responseMidpoint;
OUT.responseScale = responseScale;
OUT.nIncluded = nnz(includeTrial);
OUT.trialIndex = find(scoredTrial);
OUT.stimulus = stimulus;
OUT.quartet = quartet;
OUT.S = S(scoredTrial);
OUT.SFullFisher = OUT.S;
OUT.nSitesUsed = nSitesUsed(scoredTrial);
OUT.nTargetSites = nTargetSites(scoredTrial);
OUT.nDistractorSites = nDistractorSites(scoredTrial);
OUT.sumAbsWeight = sumAbsWeight(scoredTrial);
OUT.nScored = nScored;
OUT.nUnscored = nnz(unscoredTrial);
OUT.nWithoutFiniteScore = nnz(noFiniteScoreTrial);
OUT.nBelowCoverageThreshold = nnz(belowCoverageTrial);
OUT.unscoredTrialIndex = find(unscoredTrial);
OUT.unscoredStimulus = stimPerTrial(unscoredTrial);
OUT.unscoredQuartet = stimulusToQuartet(OUT.unscoredStimulus);
OUT.nCorrect = nCorrect;
OUT.nCorrectFullFisher = nCorrect;
OUT.accuracy = accuracy;
OUT.accuracyFullFisher = accuracy;

fprintf('\nV4 full-Fisher attention decoder\n');
fprintf('  Visually driven sites (bestSNR > %.2f): %d / %d\n', ...
    opt.SNRThreshold, nnz(visuallyDriven), nSites);
fprintf('  Sites with a finite attention weight: %d / %d\n', ...
    nnz(eligibleSite), nSites);
fprintf('  Scored trials: %d / %d included (%d without eligible curve sites)\n', ...
    nScored, OUT.nIncluded, OUT.nWithoutFiniteScore);
fprintf('  Minimum sites per quartet: %d (%d trials below threshold)\n', ...
    opt.MinSitesPerQuartet, OUT.nBelowCoverageThreshold);
fprintf('  Included quartets: %d / %d\n', numel(OUT.includedQuartet), nQuartets);
fprintf('  Median sites per trial: %.1f (target %.1f, distractor %.1f)\n', ...
    median(OUT.nSitesUsed), median(OUT.nTargetSites), ...
    median(OUT.nDistractorSites));
fprintf('  Full-covariance Fisher correct (S > 0): %d / %d = %.2f%%\n', ...
    nCorrect, nScored, 100 * accuracy);
fprintf('  Prestimulus covariance trials: %d; condition %.3g -> %.3g\n', ...
    nCovarianceTrials, covarianceCondition, regularizedCovarianceCondition);

if opt.MakeFigure
    fig = figure('Color', 'w', 'Position', [100 100 650 500]);
    histogram(OUT.S, 'BinMethod', 'fd', 'Normalization', 'probability', ...
        'FaceColor', [0.75 0.35 0.18], 'EdgeColor', 'w', 'FaceAlpha', 0.9);
    hold on;
    xline(0, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.8);
    xlabel('Attention score S');
    ylabel('Probability');
    title(sprintf(['Nilson V4 full Fisher, %g-%g ms: P(S > 0) = ' ...
        '%.1f%% (%d/%d); N = %d quartets'], opt.ResponseWindow(1), ...
        opt.ResponseWindow(2), 100 * accuracy, nCorrect, nScored, ...
        numel(OUT.includedQuartet)));
    box off;
    grid on;
    set(gca, 'FontSize', 12, 'Layer', 'top');
    OUT.figure = fig;
else
    fig = [];
end

if opt.SaveOutputs
    if exist(cfg.resultsDir, 'dir') ~= 7
        mkdir(cfg.resultsDir);
    end
    thresholdLabel = sprintf('minSites%d', opt.MinSitesPerQuartet);
    resultFile = fullfile(cfg.resultsDir, sprintf( ...
        'Attention_decoder_V4_fullFisher_%s_N.mat', thresholdLabel));
    figureFile = fullfile(cfg.resultsDir, sprintf( ...
        'Attention_decoder_V4_fullFisher_%s_N.png', thresholdLabel));
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

function Rlocal = subsetResponse(R, globalSites)
Rlocal = R;
Rlocal.meanAct = R.meanAct(globalSites, :, :);
Rlocal.meanSqAct = R.meanSqAct(globalSites, :, :);
if ~isvector(R.nTrials)
    Rlocal.nTrials = R.nTrials(globalSites, :);
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
