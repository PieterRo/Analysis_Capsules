function OUT = Compare_within_area_capsule_noise_N(varargin)
%COMPARE_WITHIN_AREA_CAPSULE_NOISE_N Capsule-conditioned V1/V4 noise correlation.
%
% Uses unique unordered non-self pairs within V1 and within V4. For every
% site and trial, the exact-stimulus mean has already been subtracted in the
% source residual file. Pairwise correlations are accumulated separately
% when both RF centers lie on the same capsule, on different capsules, on
% target-target, or on distractor-distractor.

p = inputParser;
p.addParameter('AttentionAlpha', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
p.addParameter('MinTrialsPerCondition', 20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 5 && x == floor(x));
p.addParameter('NumBootstraps', 10000, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 100 && x == floor(x));
p.addParameter('NumDistanceBins', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 3 && x == floor(x));
p.addParameter('RandomSeed', 190926, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Visible', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();
sourceFile = fullfile(cfg.repoRoot, 'results', ...
    'Trial_noise_coupling_V1_V4_IT_N.mat');
attentionFiles = { ...
    fullfile(cfg.resultsDir, 'OUT_attention_modulation_3bin_timeIdx3.mat'), ...
    fullfile(cfg.resultsDir, ...
    'OUT_attention_modulation_V4_N_3bin_timeIdx3.mat')};
geometryFiles = {fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'), ...
    fullfile(cfg.matDir, 'Tall_V4_lines_N.mat')};
calibrationFile = fullfile(cfg.extrasRoot, 'monkeyN', 'RFs', ...
    'BarMap_Nilson.mat');
required = [{sourceFile}, attentionFiles, geometryFiles, {calibrationFile}];
for file = 1:numel(required)
    assert(isfile(required{file}), 'Missing source file %s.', required{file});
end

sourceData = load(sourceFile, 'OUT');
S = sourceData.OUT;
assert(isequal(double(S.windowSampled(2, :)), [300 499]), ...
    'Expected residual samples from 300 through 499 ms.');
stimulus = double(S.stimulus(:));
residual = double(S.stimulusResidual(:, :, 2));
assert(size(residual, 2) == numel(stimulus), ...
    'Residual trials and stimulus labels disagree.');
assert(all(stimulus == floor(stimulus)) && min(stimulus) >= 1, ...
    'Stimulus labels must be positive integer geometry indices.');

attentionV1 = load(attentionFiles{1}, 'OUT');
attentionV4 = load(attentionFiles{2}, 'OUT');
geometryV1 = load(geometryFiles{1}, 'Tall_V1');
geometryV4 = load(geometryFiles{2}, 'Tall_V4');
calibration = load(calibrationFile, 'pixperdeg');
pixelsPerDegree = double(calibration.pixperdeg);
assert(isfinite(pixelsPerDegree) && pixelsPerDegree > 0, ...
    'Invalid pixels-per-degree calibration.');

areaInput(1) = areaDefinition('V1', 'V1-V1', 1, 0, ...
    geometryV1.Tall_V1, attentionV1.OUT);
areaInput(2) = areaDefinition('V4', 'V4-V4', 2, 512, ...
    geometryV4.Tall_V4, attentionV4.OUT);

areaCell = cell(2, 1);
for areaIndex = 1:2
    definition = areaInput(areaIndex);
    allSites = double(S.areaSites{definition.sourceArea}(:));
    allPositions = double(S.areaPosition{definition.sourceArea}(:));
    localSites = allSites - definition.globalOffset;
    pValue = double(definition.attention.pValueTD(localSites));
    dprime = double(definition.attention.dprime(localSites));
    selected = isfinite(pValue) & pValue < opt.AttentionAlpha & ...
        isfinite(dprime) & dprime > 0;
    sites = allSites(selected);
    positions = allPositions(selected);
    pValue = pValue(selected);
    dprime = dprime(selected);
    localSites = localSites(selected);
    assert(numel(unique(sites)) == numel(sites) && ...
        all(diff(sites) > 0), 'Selected %s sites are not unique and sorted.', ...
        definition.name);
    assert(all(isfinite(residual(positions, :)), 'all'), ...
        'Selected %s residuals contain non-finite values.', definition.name);

    roles = capsuleRoles(definition.geometry, localSites);
    centersPx = rfCenters(definition.geometry, localSites);
    fprintf(['%s within-area capsule analysis: %d/%d sites selected ' ...
        '(p<%.2f, T>D); %d possible unique pairs.\n'], ...
        definition.name, numel(sites), numel(allSites), ...
        opt.AttentionAlpha, nchoosek(numel(sites), 2));
    areaCell{areaIndex} = analyzeArea(definition, sites, positions, pValue, ...
        dprime, centersPx, pixelsPerDegree, residual, stimulus, roles, ...
        opt, opt.RandomSeed + areaIndex);
end
area = vertcat(areaCell{:});

summary = makeSummary(area);
distanceQC = vertcat(area.distanceBins);

OUT = struct();
OUT.monkey = 'Mr Nilson';
OUT.responseWindow = [300 500];
OUT.sampledWindow = [300 499];
OUT.sourceFile = sourceFile;
OUT.attentionFiles = attentionFiles;
OUT.geometryFiles = geometryFiles;
OUT.calibrationFile = calibrationFile;
OUT.pixelsPerDegree = pixelsPerDegree;
OUT.area = area;
OUT.summary = summary;
OUT.distanceQC = distanceQC;
OUT.options = opt;
OUT.residualDefinition = ['For every site and trial, subtract the mean ' ...
    '300-500 ms response for that exact stimulus.'];
OUT.siteSelectionDefinition = sprintf(['Both sites have attention p<%.2f ' ...
    'and positive target-versus-distractor d-prime.'], opt.AttentionAlpha);
OUT.conditionDefinition = ['Target=+1 and distractor=-1 only for RF ' ...
    'centers assigned to one non-overlapping capsule. Overlap, background, ' ...
    'or unassigned RF centers contribute to no condition.'];
OUT.pairDefinition = ['Unique unordered physical site pairs only (i<j); ' ...
    'self-pairs and reversed duplicates are excluded.'];
OUT.inferenceDefinition = sprintf(['Bars are back-transformed means of ' ...
    'pairwise Fisher-z correlations. A comparison uses one common pair ' ...
    'set with at least %d trials in every displayed condition. Each ' ...
    'bootstrap resamples the single within-area site set once and weights ' ...
    'pair i<j by count_i*count_j. Error bars are bootstrap SEM; p-values ' ...
    'are two-sided bootstrap-Wald tests in Fisher-z space, uncorrected.'], ...
    opt.MinTrialsPerCondition);
OUT.distanceCaveat = ['RF separation is not experimentally matched by ' ...
    'capsule condition. Pair tables and distance-bin diagnostics are saved; ' ...
    'the primary figure is the requested direct matched-pair comparison.'];

stem = sprintf(['Within_area_capsule_noise_V1_V4_300to500ms_' ...
    'attentionSigP%g_minTrials%d_N'], opt.AttentionAlpha, ...
    opt.MinTrialsPerCondition);
OUT.resultFile = fullfile(cfg.resultsDir, [stem '.mat']);
OUT.summaryFile = fullfile(cfg.resultsDir, [stem '_summary.csv']);
OUT.distanceQCFile = fullfile(cfg.resultsDir, ...
    [stem '_distance_qc.csv']);
OUT.pairFiles = {fullfile(cfg.resultsDir, [stem '_V1_pairs.csv']), ...
    fullfile(cfg.resultsDir, [stem '_V4_pairs.csv'])};
OUT.figureFile = fullfile(cfg.resultsDir, [stem '.png']);
OUT.figurePdf = fullfile(cfg.resultsDir, [stem '.pdf']);

fig = makeFigure(OUT, opt);
if opt.SaveOutputs
    save(OUT.resultFile, 'OUT', '-v7.3');
    writetable(summary, OUT.summaryFile);
    writetable(distanceQC, OUT.distanceQCFile);
    for areaIndex = 1:2
        writetable(area(areaIndex).pairTable, OUT.pairFiles{areaIndex});
    end
    set(fig, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 11.8 6.1], 'PaperSize', [11.8 6.1]);
    print(fig, OUT.figureFile, '-dpng', '-r220');
    print(fig, OUT.figurePdf, '-dpdf', '-painters');
    fprintf('Saved %s\nSaved %s\nSaved %s\n', OUT.resultFile, ...
        OUT.summaryFile, OUT.distanceQCFile);
    fprintf('Saved %s\nSaved %s\nSaved %s\nSaved %s\n', ...
        OUT.pairFiles{1}, OUT.pairFiles{2}, OUT.figureFile, OUT.figurePdf);
end
if ~opt.Visible
    close(fig);
else
    OUT.figure = fig;
end
end

function D = areaDefinition(name, pairName, sourceArea, globalOffset, ...
        geometry, attention)
D = struct('name', name, 'pairName', pairName, ...
    'sourceArea', sourceArea, 'globalOffset', globalOffset, ...
    'geometry', geometry, 'attention', attention);
end

function A = analyzeArea(D, sites, positions, pValue, dprime, centersPx, ...
        pixelsPerDegree, residualAll, stimulus, roles, opt, randomSeed)
nSites = numel(sites);
upper = triu(true(nSites), 1);
distanceDeg = pairDistances(centersPx, pixelsPerDegree);
moments = accumulateMoments(residualAll(positions, :), stimulus, roles);
same = finishMoments(moments.same, opt.MinTrialsPerCondition, upper);
different = finishMoments(moments.different, ...
    opt.MinTrialsPerCondition, upper);
targetTarget = finishMoments(moments.targetTarget, ...
    opt.MinTrialsPerCondition, upper);
distractorDistractor = finishMoments(moments.distractorDistractor, ...
    opt.MinTrialsPerCondition, upper);

sameDifferentEligible = upper & isfinite(same.r) & isfinite(different.r);
threeEligible = upper & isfinite(targetTarget.r) & ...
    isfinite(distractorDistractor.r) & isfinite(different.r);
assert(any(sameDifferentEligible(:)) && any(threeEligible(:)), ...
    '%s lacks matched pairs for the requested comparisons.', D.name);

rSame = same.r;
rDifferent = different.r;
rTT = targetTarget.r;
rDD = distractorDistractor.r;
rng(randomSeed, 'twister');
bootstrapSameDifferent = withinSiteBootstrap( ...
    {rSame, rDifferent}, sameDifferentEligible, ...
    opt.NumBootstraps, sprintf('%s same/different', D.name));
bootstrapThree = withinSiteBootstrap({rTT, rDD, rDifferent}, ...
    threeEligible, opt.NumBootstraps, sprintf('%s TT/DD/different', D.name));

sameDifferentR = categoryMeans({rSame, rDifferent}, ...
    sameDifferentEligible);
sameDifferentMedianR = categoryMedians({rSame, rDifferent}, ...
    sameDifferentEligible);
sameDifferentSEM = bootstrapSEM(bootstrapSameDifferent);
[sameDifferentCILow, sameDifferentCIHigh] = ...
    bootstrapIntervals(bootstrapSameDifferent);
sameMinusDifferentZ = mean(safeAtanh(rSame(sameDifferentEligible)) - ...
    safeAtanh(rDifferent(sameDifferentEligible)));
sameMinusDifferentP = bootstrapWaldP(sameMinusDifferentZ, ...
    bootstrapSameDifferent(:, 1) - bootstrapSameDifferent(:, 2));

threeConditionR = categoryMeans({rTT, rDD, rDifferent}, threeEligible);
threeConditionMedianR = categoryMedians({rTT, rDD, rDifferent}, ...
    threeEligible);
threeConditionSEM = bootstrapSEM(bootstrapThree);
[threeConditionCILow, threeConditionCIHigh] = ...
    bootstrapIntervals(bootstrapThree);
ttMinusDdZ = mean(safeAtanh(rTT(threeEligible)) - ...
    safeAtanh(rDD(threeEligible)));
ttMinusDifferentZ = mean(safeAtanh(rTT(threeEligible)) - ...
    safeAtanh(rDifferent(threeEligible)));
ddMinusDifferentZ = mean(safeAtanh(rDD(threeEligible)) - ...
    safeAtanh(rDifferent(threeEligible)));
threeConditionP = [ ...
    bootstrapWaldP(ttMinusDdZ, ...
    bootstrapThree(:, 1) - bootstrapThree(:, 2)), ...
    bootstrapWaldP(ttMinusDifferentZ, ...
    bootstrapThree(:, 1) - bootstrapThree(:, 3)), ...
    bootstrapWaldP(ddMinusDifferentZ, ...
    bootstrapThree(:, 2) - bootstrapThree(:, 3))];

zDelta = safeAtanh(rSame) - safeAtanh(rDifferent);
[distanceCorrelation, distanceCorrelationP] = correlationWithP( ...
    distanceDeg(sameDifferentEligible), zDelta(sameDifferentEligible));
distanceSlope = linearSlope(distanceDeg(sameDifferentEligible), ...
    zDelta(sameDifferentEligible));
distanceBins = summarizeDistanceBins(D.pairName, distanceDeg, rSame, ...
    rDifferent, sameDifferentEligible, opt.NumDistanceBins);

[row, column] = find(upper);
site1 = sites(row);
site2 = sites(column);
assert(all(site1 < site2) && all(site1 ~= site2), ...
    '%s pair indexing includes a self-pair or reversed pair.', D.name);
unorderedKey = sort([site1 site2], 2);
assert(size(unique(unorderedKey, 'rows'), 1) == numel(row) && ...
    numel(row) == nchoosek(nSites, 2), ...
    '%s pair indexing contains duplicates or omissions.', D.name);
pairTable = table(site1, site2, pValue(row), pValue(column), ...
    dprime(row), dprime(column), distanceDeg(upper), ...
    same.n(upper), different.n(upper), targetTarget.n(upper), ...
    distractorDistractor.n(upper), rSame(upper), rDifferent(upper), ...
    rTT(upper), rDD(upper), zDelta(upper), ...
    sameDifferentEligible(upper), threeEligible(upper), ...
    'VariableNames', {'Site1Global', 'Site2Global', ...
    'Site1AttentionP', 'Site2AttentionP', 'Site1AttentionDprime', ...
    'Site2AttentionDprime', 'RFDistanceDeg', 'NSame', 'NDifferent', ...
    'NTT', 'NDD', 'RSame', 'RDifferent', 'RTT', 'RDD', ...
    'SameMinusDifferentFisherZ', 'EligibleSameDifferent', ...
    'EligibleTTDDDifferent'});

sameValid = upper & isfinite(rSame);
differentValid = upper & isfinite(rDifferent);
A = struct();
A.name = D.name;
A.pairName = D.pairName;
A.sites = sites;
A.positionsInResidual = positions;
A.attentionPValue = pValue;
A.attentionDprime = dprime;
A.rfCentersPx = centersPx;
A.roles = roles;
A.nSelectedSites = nSites;
A.nPossibleUniquePairs = nchoosek(nSites, 2);
A.same = same;
A.different = different;
A.targetTarget = targetTarget;
A.distractorDistractor = distractorDistractor;
A.sameDifferentEligible = sameDifferentEligible;
A.threeConditionEligible = threeEligible;
A.nSameDifferentPairs = nnz(sameDifferentEligible);
A.nThreeConditionPairs = nnz(threeEligible);
A.sameDifferentR = sameDifferentR;
A.sameDifferentMedianR = sameDifferentMedianR;
A.sameDifferentSEM = sameDifferentSEM;
A.sameDifferentCILow = sameDifferentCILow;
A.sameDifferentCIHigh = sameDifferentCIHigh;
A.sameMinusDifferentFisherZ = sameMinusDifferentZ;
A.sameMinusDifferentP = sameMinusDifferentP;
A.threeConditionR = threeConditionR;
A.threeConditionMedianR = threeConditionMedianR;
A.threeConditionSEM = threeConditionSEM;
A.threeConditionCILow = threeConditionCILow;
A.threeConditionCIHigh = threeConditionCIHigh;
A.threeConditionP = threeConditionP;
A.bootstrapSameDifferentFisherZ = bootstrapSameDifferent;
A.bootstrapThreeConditionFisherZ = bootstrapThree;
A.rfDistanceDeg = distanceDeg;
A.distanceSameValid = distanceSummary(distanceDeg(sameValid));
A.distanceDifferentValid = distanceSummary(distanceDeg(differentValid));
A.distanceMatched = distanceSummary(distanceDeg(sameDifferentEligible));
A.distanceDeltaFisherZCorrelation = distanceCorrelation;
A.distanceDeltaFisherZCorrelationP = distanceCorrelationP;
A.distanceDeltaFisherZSlopePerDegree = distanceSlope;
A.distanceBins = distanceBins;
A.pairTable = pairTable;
A.validation = struct('nSelfPairs', nnz(row == column), ...
    'nUniquePairs', size(unique(unorderedKey, 'rows'), 1), ...
    'expectedUniquePairs', nchoosek(nSites, 2), ...
    'rolesSymmetricByDefinition', true, ...
    'upperTriangleOnly', true);

fprintf(['  %s matched same/different: N=%d pairs, Fisher means r ' ...
    '%.4f/%.4f, medians %.4f/%.4f, p=%.4g.\n'], D.name, ...
    A.nSameDifferentPairs, A.sameDifferentR, ...
    A.sameDifferentMedianR, A.sameMinusDifferentP);
fprintf(['  %s matched TT/DD/different: N=%d pairs, Fisher means r ' ...
    '%.4f/%.4f/%.4f, medians %.4f/%.4f/%.4f, p(TT-DD)=%.4g, ' ...
    'p(TT-diff)=%.4g, p(DD-diff)=%.4g.\n'], D.name, ...
    A.nThreeConditionPairs, A.threeConditionR, ...
    A.threeConditionMedianR, A.threeConditionP);
fprintf(['  %s eligible RF distance median %.3f deg [%.3f, %.3f]; ' ...
    'corr(distance, delta z)=%.3f, p=%.4g.\n'], D.name, ...
    A.distanceMatched.median, A.distanceMatched.minimum, ...
    A.distanceMatched.maximum, A.distanceDeltaFisherZCorrelation, ...
    A.distanceDeltaFisherZCorrelationP);
fprintf(['  %s RF distance by independently valid condition: same ' ...
    'N=%d, median %.3f deg [IQR %.3f-%.3f]; different N=%d, ' ...
    'median %.3f deg [IQR %.3f-%.3f].\n'], D.name, ...
    A.distanceSameValid.n, A.distanceSameValid.median, ...
    A.distanceSameValid.quartile1, A.distanceSameValid.quartile3, ...
    A.distanceDifferentValid.n, A.distanceDifferentValid.median, ...
    A.distanceDifferentValid.quartile1, ...
    A.distanceDifferentValid.quartile3);
end

function role = capsuleRoles(Tall, localSites)
nStimuli = numel(Tall);
role = zeros(numel(localSites), nStimuli, 'int8');
for stimulus = 1:nStimuli
    T = Tall(stimulus).T;
    assignment = string(T.assignment(localSites));
    overlap = false(numel(localSites), 1);
    if ismember('overlap', T.Properties.VariableNames)
        overlap = logical(T.overlap(localSites));
    end
    role(assignment == "target" & ~overlap, stimulus) = 1;
    role(assignment == "distractor" & ~overlap, stimulus) = -1;
end
end

function centers = rfCenters(Tall, localSites)
T = Tall(1).T;
centers = [double(T.x_px(localSites)), double(T.y_px(localSites))];
assert(all(isfinite(centers), 'all'), 'A selected RF center is non-finite.');
end

function distance = pairDistances(centers, pixelsPerDegree)
dx = centers(:, 1) - centers(:, 1)';
dy = centers(:, 2) - centers(:, 2)';
distance = hypot(dx, dy) / pixelsPerDegree;
end

function M = accumulateMoments(response, stimulus, role)
nSites = size(response, 1);
upper = triu(true(nSites), 1);
M.same = emptyMoments(nSites);
M.different = emptyMoments(nSites);
M.targetTarget = emptyMoments(nSites);
M.distractorDistractor = emptyMoments(nSites);
for stim = unique(stimulus(:))'
    trials = stimulus == stim;
    if ~any(trials)
        continue;
    end
    responseStim = response(:, trials);
    siteRole = double(role(:, stim));
    roleProduct = siteRole * siteRole';
    sameFull = roleProduct == 1;
    differentFull = roleProduct == -1;
    targetFull = siteRole == 1 & siteRole' == 1;
    distractorFull = siteRole == -1 & siteRole' == -1;
    assert(isequal(sameFull, sameFull') && ...
        isequal(differentFull, differentFull') && ...
        isequal(targetFull, targetFull') && ...
        isequal(distractorFull, distractorFull'), ...
        'Within-area condition masks are not symmetric.');
    M.same = addMoments(M.same, responseStim, sameFull & upper);
    M.different = addMoments(M.different, responseStim, ...
        differentFull & upper);
    M.targetTarget = addMoments(M.targetTarget, responseStim, ...
        targetFull & upper);
    M.distractorDistractor = addMoments(M.distractorDistractor, ...
        responseStim, distractorFull & upper);
end
end

function M = emptyMoments(nSites)
z = zeros(nSites, nSites);
M = struct('n', z, 'sumX', z, 'sumY', z, 'sumXX', z, ...
    'sumYY', z, 'sumXY', z);
end

function M = addMoments(M, response, mask)
if ~any(mask(:))
    return;
end
n = size(response, 2);
sumX = sum(response, 2);
sumXX = sum(response.^2, 2);
M.n = M.n + n * mask;
M.sumX = M.sumX + mask .* sumX;
M.sumY = M.sumY + mask .* sumX';
M.sumXX = M.sumXX + mask .* sumXX;
M.sumYY = M.sumYY + mask .* sumXX';
M.sumXY = M.sumXY + mask .* (response * response');
end

function R = finishMoments(M, minimumN, upper)
n = M.n;
covNumerator = M.sumXY - M.sumX .* M.sumY ./ max(n, 1);
varX = M.sumXX - M.sumX.^2 ./ max(n, 1);
varY = M.sumYY - M.sumY.^2 ./ max(n, 1);
denominator = sqrt(max(varX, 0) .* max(varY, 0));
valid = upper & n >= minimumN & denominator > 0;
r = nan(size(n));
r(valid) = covNumerator(valid) ./ denominator(valid);
r(valid) = min(1, max(-1, r(valid)));
R = struct('n', n, 'r', r);
end

function bootstrap = withinSiteBootstrap(r, eligible, nBootstraps, label)
nSites = size(eligible, 1);
nCategories = numel(r);
z = cellfun(@safeAtanh, r, 'UniformOutput', false);
bootstrap = nan(nBootstraps, nCategories);
timer = tic;
lastReport = 0;
for sample = 1:nBootstraps
    count = accumarray(randi(nSites, nSites, 1), 1, [nSites 1]);
    weight = count * count';
    use = eligible & weight > 0;
    if any(use(:))
        for category = 1:nCategories
            bootstrap(sample, category) = weightedMean( ...
                z{category}(use), weight(use));
        end
    end
    elapsed = toc(timer);
    if elapsed - lastReport >= 20 || sample == nBootstraps
        remaining = elapsed * (nBootstraps / sample - 1);
        fprintf('  %s bootstrap %d/%d; ETA %.0f s.\n', ...
            label, sample, nBootstraps, remaining);
        lastReport = elapsed;
    end
end
assert(all(any(isfinite(bootstrap), 1)), ...
    'A bootstrap category has no finite samples.');
end

function means = categoryMeans(r, eligible)
means = nan(1, numel(r));
for category = 1:numel(r)
    z = safeAtanh(r{category});
    means(category) = tanh(mean(z(eligible)));
end
end

function medians = categoryMedians(r, eligible)
medians = nan(1, numel(r));
for category = 1:numel(r)
    medians(category) = median(r{category}(eligible));
end
end

function [low, high] = bootstrapIntervals(bootstrapZ)
nCategories = size(bootstrapZ, 2);
low = nan(1, nCategories);
high = nan(1, nCategories);
for category = 1:nCategories
    interval = percentile(tanh(bootstrapZ(:, category)), [2.5 97.5]);
    low(category) = interval(1);
    high(category) = interval(2);
end
end

function sem = bootstrapSEM(bootstrapZ)
sem = nan(1, size(bootstrapZ, 2));
for category = 1:size(bootstrapZ, 2)
    values = tanh(bootstrapZ(:, category));
    sem(category) = std(values(isfinite(values)), 0);
end
end

function p = bootstrapWaldP(observedDifference, bootstrapDifferences)
bootstrapDifferences = bootstrapDifferences(isfinite(bootstrapDifferences));
standardError = std(bootstrapDifferences, 0);
if standardError == 0
    p = double(observedDifference == 0);
else
    p = erfc(abs(observedDifference / standardError) / sqrt(2));
end
end

function value = weightedMean(x, weight)
value = sum(x .* weight) / sum(weight);
end

function z = safeAtanh(r)
limit = 1 - 1e-12;
z = atanh(min(limit, max(-limit, r)));
z(~isfinite(r)) = NaN;
end

function summary = distanceSummary(distance)
distance = double(distance(isfinite(distance)));
assert(~isempty(distance), 'Cannot summarize an empty distance set.');
quartile = percentile(distance, [25 75]);
summary = struct('n', numel(distance), 'minimum', min(distance), ...
    'quartile1', quartile(1), 'median', median(distance), ...
    'quartile3', quartile(2), 'maximum', max(distance));
end

function [r, p] = correlationWithP(x, y)
[correlation, pValue] = corrcoef(double(x(:)), double(y(:)));
r = correlation(1, 2);
p = pValue(1, 2);
end

function slope = linearSlope(x, y)
x = double(x(:));
y = double(y(:));
x = x - mean(x);
y = y - mean(y);
slope = sum(x .* y) / sum(x.^2);
end

function T = summarizeDistanceBins(pairName, distance, rSame, rDifferent, ...
        eligible, nBins)
linear = find(eligible & isfinite(distance));
[~, order] = sort(distance(linear));
linear = linear(order);
edge = round(linspace(0, numel(linear), nBins + 1));
areaPair = strings(nBins, 1);
distanceBin = (1:nBins)';
nPairs = zeros(nBins, 1);
distanceLowDeg = nan(nBins, 1);
distanceMedianDeg = nan(nBins, 1);
distanceHighDeg = nan(nBins, 1);
rSameMean = nan(nBins, 1);
rDifferentMean = nan(nBins, 1);
deltaR = nan(nBins, 1);
for bin = 1:nBins
    use = linear((edge(bin) + 1):edge(bin + 1));
    areaPair(bin) = pairName;
    nPairs(bin) = numel(use);
    distanceLowDeg(bin) = min(distance(use));
    distanceMedianDeg(bin) = median(distance(use));
    distanceHighDeg(bin) = max(distance(use));
    rSameMean(bin) = tanh(mean(safeAtanh(rSame(use))));
    rDifferentMean(bin) = tanh(mean(safeAtanh(rDifferent(use))));
    deltaR(bin) = rSameMean(bin) - rDifferentMean(bin);
end
T = table(areaPair, distanceBin, nPairs, distanceLowDeg, ...
    distanceMedianDeg, distanceHighDeg, rSameMean, rDifferentMean, deltaR, ...
    'VariableNames', {'AreaPair', 'DistanceBin', 'NPairs', ...
    'DistanceLowDeg', 'DistanceMedianDeg', 'DistanceHighDeg', ...
    'RSame', 'RDifferent', 'DeltaR'});
end

function T = makeSummary(area)
areaPair = string({area.pairName})';
nSelectedSites = [area.nSelectedSites]';
nPossibleUniquePairs = [area.nPossibleUniquePairs]';
nSameDifferentPairs = [area.nSameDifferentPairs]';
nThreeConditionPairs = [area.nThreeConditionPairs]';
sameDifferentR = vertcat(area.sameDifferentR);
sameDifferentMedianR = vertcat(area.sameDifferentMedianR);
sameDifferentSEM = vertcat(area.sameDifferentSEM);
sameDifferentCILow = vertcat(area.sameDifferentCILow);
sameDifferentCIHigh = vertcat(area.sameDifferentCIHigh);
threeR = vertcat(area.threeConditionR);
threeMedianR = vertcat(area.threeConditionMedianR);
threeSEM = vertcat(area.threeConditionSEM);
threeCILow = vertcat(area.threeConditionCILow);
threeCIHigh = vertcat(area.threeConditionCIHigh);
pSameDifferent = [area.sameMinusDifferentP]';
threeP = vertcat(area.threeConditionP);
nSameValidPairs = arrayfun(@(x) x.distanceSameValid.n, area);
sameValidDistanceMedianDeg = arrayfun( ...
    @(x) x.distanceSameValid.median, area);
sameValidDistanceQ1Deg = arrayfun( ...
    @(x) x.distanceSameValid.quartile1, area);
sameValidDistanceQ3Deg = arrayfun( ...
    @(x) x.distanceSameValid.quartile3, area);
nDifferentValidPairs = arrayfun(@(x) x.distanceDifferentValid.n, area);
differentValidDistanceMedianDeg = arrayfun( ...
    @(x) x.distanceDifferentValid.median, area);
differentValidDistanceQ1Deg = arrayfun( ...
    @(x) x.distanceDifferentValid.quartile1, area);
differentValidDistanceQ3Deg = arrayfun( ...
    @(x) x.distanceDifferentValid.quartile3, area);
distanceMedianDeg = arrayfun(@(x) x.distanceMatched.median, area);
distanceQ1Deg = arrayfun(@(x) x.distanceMatched.quartile1, area);
distanceQ3Deg = arrayfun(@(x) x.distanceMatched.quartile3, area);
distanceDeltaCorrelation = [area.distanceDeltaFisherZCorrelation]';
distanceDeltaCorrelationP = [area.distanceDeltaFisherZCorrelationP]';
distanceDeltaSlopePerDeg = [area.distanceDeltaFisherZSlopePerDegree]';
T = table(areaPair, nSelectedSites, nPossibleUniquePairs, ...
    nSameDifferentPairs, nThreeConditionPairs, ...
    sameDifferentR(:, 1), sameDifferentMedianR(:, 1), ...
    sameDifferentSEM(:, 1), sameDifferentCILow(:, 1), ...
    sameDifferentCIHigh(:, 1), sameDifferentR(:, 2), ...
    sameDifferentMedianR(:, 2), sameDifferentSEM(:, 2), ...
    sameDifferentCILow(:, 2), sameDifferentCIHigh(:, 2), ...
    pSameDifferent, threeR(:, 1), threeMedianR(:, 1), threeSEM(:, 1), ...
    threeCILow(:, 1), threeCIHigh(:, 1), threeR(:, 2), ...
    threeMedianR(:, 2), threeSEM(:, 2), threeCILow(:, 2), ...
    threeCIHigh(:, 2), threeR(:, 3), threeMedianR(:, 3), ...
    threeSEM(:, 3), threeCILow(:, 3), threeCIHigh(:, 3), ...
    threeP(:, 1), threeP(:, 2), threeP(:, 3), nSameValidPairs, ...
    sameValidDistanceMedianDeg, sameValidDistanceQ1Deg, ...
    sameValidDistanceQ3Deg, nDifferentValidPairs, ...
    differentValidDistanceMedianDeg, differentValidDistanceQ1Deg, ...
    differentValidDistanceQ3Deg, distanceMedianDeg, ...
    distanceQ1Deg, distanceQ3Deg, distanceDeltaCorrelation, ...
    distanceDeltaCorrelationP, distanceDeltaSlopePerDeg, ...
    'VariableNames', {'AreaPair', 'NSelectedSites', ...
    'NPossibleUniquePairs', 'NSameDifferentPairs', ...
    'NThreeConditionPairs', 'RSame', 'RSameMedian', 'RSameSEM', ...
    'RSameCILow', 'RSameCIHigh', 'RDifferent', 'RDifferentMedian', ...
    'RDifferentSEM', 'RDifferentCILow', 'RDifferentCIHigh', ...
    'PSameVsDifferent', 'RTT', 'RTTMedian', 'RTTSEM', 'RTTCILow', ...
    'RTTCIHigh', 'RDD', 'RDDMedian', 'RDDSEM', 'RDDCILow', ...
    'RDDCIHigh', 'RDifferentMatched', 'RDifferentMatchedMedian', ...
    'RDifferentMatchedSEM', 'RDifferentMatchedCILow', ...
    'RDifferentMatchedCIHigh', 'PTTVsDD', 'PTTVsDifferent', ...
    'PDDVsDifferent', 'NSameValidPairs', ...
    'SameValidRFDistanceMedianDeg', 'SameValidRFDistanceQ1Deg', ...
    'SameValidRFDistanceQ3Deg', 'NDifferentValidPairs', ...
    'DifferentValidRFDistanceMedianDeg', ...
    'DifferentValidRFDistanceQ1Deg', 'DifferentValidRFDistanceQ3Deg', ...
    'RFDistanceMedianDeg', 'RFDistanceQ1Deg', ...
    'RFDistanceQ3Deg', 'DistanceDeltaFisherZCorrelation', ...
    'DistanceDeltaFisherZCorrelationP', ...
    'DistanceDeltaFisherZSlopePerDeg'});
end

function fig = makeFigure(OUT, opt)
visibility = 'off';
if opt.Visible
    visibility = 'on';
end
fig = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [80 80 1540 760], 'Name', ...
    'Within-area capsule-conditioned noise correlation', ...
    'NumberTitle', 'off');
sameColor = [0.84 0.25 0.12];
ttColor = [0.82 0.18 0.12];
ddColor = [0.93 0.57 0.12];
differentColor = [0.10 0.39 0.70];
groupX = 1:2;

ax1 = axes('Parent', fig, 'Position', [0.075 0.18 0.39 0.68]);
hold(ax1, 'on');
values = vertcat(OUT.area.sameDifferentR);
errors = vertcat(OUT.area.sameDifferentSEM);
offset = [-0.18 0.18];
colors = {sameColor, differentColor};
handles = gobjects(1, 2);
for category = 1:2
    x = groupX + offset(category);
    handles(category) = bar(ax1, x, values(:, category), 0.32, ...
        'FaceColor', colors{category}, 'EdgeColor', 'none');
    errorbar(ax1, x, values(:, category), errors(:, category), ...
        errors(:, category), 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.2, 'CapSize', 7);
end
set(ax1, 'XTick', groupX, ...
    'XTickLabel', pairLabels(OUT.area, false));
ylabel(ax1, 'Mean residual noise correlation, r');
title(ax1, 'Same versus different capsule');
legend(ax1, handles, {'Same capsule', 'Different capsules'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
styleAxes(ax1);

ax2 = axes('Parent', fig, 'Position', [0.555 0.18 0.39 0.68]);
hold(ax2, 'on');
threeValues = vertcat(OUT.area.threeConditionR);
threeErrors = vertcat(OUT.area.threeConditionSEM);
offset3 = [-0.26 0 0.26];
colors3 = {ttColor, ddColor, differentColor};
handles3 = gobjects(1, 3);
for category = 1:3
    x = groupX + offset3(category);
    handles3(category) = bar(ax2, x, threeValues(:, category), 0.24, ...
        'FaceColor', colors3{category}, 'EdgeColor', 'none');
    errorbar(ax2, x, threeValues(:, category), ...
        threeErrors(:, category), threeErrors(:, category), 'k', ...
        'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 7);
end
set(ax2, 'XTick', groupX, 'XTickLabel', pairLabels(OUT.area, true));
ylabel(ax2, 'Mean residual noise correlation, r');
title(ax2, 'Target-target, distractor-distractor, and different');
legend(ax2, handles3, {'TT', 'DD', 'Different'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
styleAxes(ax2);

yMinimum = min([0; values(:) - errors(:); ...
    threeValues(:) - threeErrors(:)]) - 0.012;
yRange = max([values(:) + errors(:); ...
    threeValues(:) + threeErrors(:)]) - yMinimum;
bracketStep = max(0.008, 0.09 * yRange);
leftY = max(values + errors, [], 2) + 0.45 * bracketStep;
rightY = max(threeValues + threeErrors, [], 2) + 0.45 * bracketStep;
yMaximum = max([leftY; rightY + 2 * bracketStep]) + 1.3 * bracketStep;
ylim(ax1, [yMinimum yMaximum]);
ylim(ax2, [yMinimum yMaximum]);
plot(ax1, xlim(ax1), [0 0], ':', 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off');
plot(ax2, xlim(ax2), [0 0], ':', 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off');

for areaIndex = 1:2
    addPBracket(ax1, groupX(areaIndex) + [-0.18 0.18], ...
        leftY(areaIndex), OUT.area(areaIndex).sameMinusDifferentP, ...
        0.18 * bracketStep);
    pValue = OUT.area(areaIndex).threeConditionP;
    addPBracket(ax2, groupX(areaIndex) + [-0.26 0], ...
        rightY(areaIndex), pValue(1), 0.18 * bracketStep);
    addPBracket(ax2, groupX(areaIndex) + [0 0.26], ...
        rightY(areaIndex) + bracketStep, pValue(3), ...
        0.18 * bracketStep);
    addPBracket(ax2, groupX(areaIndex) + [-0.26 0.26], ...
        rightY(areaIndex) + 2 * bracketStep, pValue(2), ...
        0.18 * bracketStep);
end

annotation(fig, 'textbox', [0.02 0.93 0.96 0.055], ...
    'String', ['Nilson within-area capsule-conditioned residual noise ' ...
    'correlations, 300-500 ms'], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none', ...
    'FontSize', 15.5, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.035 0.025 0.93 0.06], ...
    'String', sprintf(['Unique non-self pairs; both sites p<%.2f and T>D. ' ...
    'Error bars: within-area site-bootstrap SEM (%d resamples). ' ...
    'Two-sided bootstrap-Wald p-values, uncorrected. Minimum %d trials ' ...
    'per condition. RF separation is reported separately as QC.'], ...
    opt.AttentionAlpha, opt.NumBootstraps, opt.MinTrialsPerCondition), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 9.2, 'Color', [0.25 0.25 0.25]);
end

function labels = pairLabels(area, useThree)
labels = cell(numel(area), 1);
for index = 1:numel(area)
    if useThree
        n = area(index).nThreeConditionPairs;
    else
        n = area(index).nSameDifferentPairs;
    end
    labels{index} = sprintf('%s (N=%d pairs)', area(index).pairName, n);
end
end

function addPBracket(ax, x, y, p, capHeight)
if p < 0.05
    color = [0.72 0.08 0.08];
    fontWeight = 'bold';
else
    color = [0.35 0.35 0.35];
    fontWeight = 'normal';
end
plot(ax, [x(1) x(1) x(2) x(2)], ...
    [y - capHeight y y y - capHeight], '-', 'Color', color, ...
    'LineWidth', 1, 'HandleVisibility', 'off');
text(ax, mean(x), y + 0.12 * capHeight, formatPValue(p), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Color', color, 'FontSize', 8.6, 'FontWeight', fontWeight);
end

function label = formatPValue(p)
if p < 0.001
    label = 'p<0.001';
else
    label = sprintf('p=%.3f', p);
end
end

function styleAxes(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11.5, ...
    'LineWidth', 1, 'YGrid', 'on', 'GridAlpha', 0.12, ...
    'Layer', 'top', 'XLim', [0.5 2.5]);
end

function q = percentile(x, requested)
x = sort(double(x(isfinite(x))));
assert(~isempty(x), 'Cannot compute a percentile of an empty vector.');
q = nan(size(requested));
for index = 1:numel(requested)
    position = 1 + (numel(x) - 1) * requested(index) / 100;
    lower = floor(position);
    upper = ceil(position);
    fraction = position - lower;
    q(index) = x(lower) * (1 - fraction) + x(upper) * fraction;
end
end
