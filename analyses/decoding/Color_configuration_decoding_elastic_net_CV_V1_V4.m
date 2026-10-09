function OUT = Color_configuration_decoding_elastic_net_CV_V1_V4(varargin)
%COLOR_CONFIGURATION_DECODING_ELASTIC_NET_CV_V1_V4 Tune color decoders.
%
% Uses the same trials, RF-on-curve availability, response window, and
% minimum site threshold as the attention decoder. The two classes are
% opposite yellow-purple configurations, with target color and target side
% balanced within each class. RF color comes from
% Tall(stim).T.center_color, not ALLMAT.

p = inputParser;
p.addParameter('Monkey', 'N', @(x) any(strcmpi(string(x), ["N" "F"])));
p.addParameter('ResponseWindow', [300 500], ...
    @(x) isnumeric(x) && isequal(size(x), [1 2]) && ...
    all(isfinite(x)) && x(1) < x(2));
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
p.addParameter('RandomSeed', 140926, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('CoefficientTolerance', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ChunkTrials', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
opt.AlphaGrid = unique(double(opt.AlphaGrid(:)'));
monkey = upper(char(opt.Monkey));

cfg = config();
if monkey == 'F'
    baseFile = fullfile(cfg.resultsDir, ...
        'Attention_decoder_Figaro_capsules_snr0p5_minSites1_F.mat');
    base = loadFigaroFeatureSource(baseFile);
else
    baseFile = fullfile(cfg.resultsDir, ...
        'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');
    base = loadFeatureSource(baseFile);
end

OUT = struct();
OUT.description = ['Repeated quartet-balanced cross-validation for ' ...
    'elastic-net logistic color-configuration decoders in V1 and V4.'];
OUT.monkey = 'Mr Nilson';
if monkey == 'F'; OUT.monkey = 'Figaro'; end
OUT.classDefinition = ['Within each eight-stimulus block, positions ' ...
    '[1 3 6 8] are canonical color configuration A (+1), and positions ' ...
    '[2 4 5 7] are the opposite configuration B (-1).'];
OUT.targetBalancing = ['Target color and target side are balanced within ' ...
    'both classes; fold-specific quartet means equally weight all four ' ...
    'configuration-by-target-color cells.'];
OUT.colorSource = sprintf('Tall_(region)_lines_%s.mat: T.center_color', monkey);
OUT.featureSource = baseFile;
OUT.responseCentering = 'training-fold quartet means';
OUT.responseWindow = double(opt.ResponseWindow);
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
OUT.singletonPolicy = ['A singleton configuration-by-target-color cell remains ' ...
    'in every training fold; that singleton trial receives no cross-validated score.'];

colorV1 = prepareColorConfigurationFeatures(base.V1, 'V1', cfg, ...
    OUT.responseWindow, opt.ChunkTrials, monkey);
colorV4 = prepareColorConfigurationFeatures(base.V4, 'V4', cfg, ...
    OUT.responseWindow, opt.ChunkTrials, monkey);
OUT.V1 = crossValidateRegion(colorV1, opt, opt.RandomSeed + 1000);
OUT.V4 = crossValidateRegion(colorV4, opt, opt.RandomSeed + 2000);

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

fileStem = resultFileStem(OUT.responseWindow);
if monkey == 'F'
    fileStem = strrep(fileStem, '_minSites20_N', '_FigaroSNR0p5_minSites1_F');
end
resultFile = fullfile(cfg.resultsDir, [fileStem '.mat']);
figureFile = fullfile(cfg.resultsDir, [fileStem '.png']);
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

function base = loadFigaroFeatureSource(baseFile)
assert(isfile(baseFile), 'Missing Figaro attention feature source %s.', baseFile);
data = load(baseFile, 'OUT');
F = data.OUT;
assert(F.snrThreshold == 0.5 && F.minSitesPerQuartet == 1 && ...
    isequal(F.responseWindow, [300 500]), ...
    'Figaro attention feature source does not match the expected selection.');
base = struct('curveMarginDeg', 0, 'minSitesPerQuartet', 1);
for region = {'V1', 'V4'}
    name = region{1};
    A = F.(name);
    selected = find(A.eligibleSite(:));
    D = struct();
    D.region = name;
    D.siteIndexLocal = selected(:)';
    D.siteIndexGlobal = double(A.globalSites(selected(:)))';
    D.trialIndex = double(A.trialIndex(:));
    D.stimulus = double(A.stimulus(:));
    D.quartet = double(A.quartet(:));
    D.baseline = double(A.baseline(:));
    D.responseScale = double(A.responseScale(:));
    role = double(A.roleByStimulus(selected, :));
    D.availability = role(:, D.stimulus)' ~= 0;
    hasFeature = any(D.availability, 2);
    D.trialIndex = D.trialIndex(hasFeature);
    D.stimulus = D.stimulus(hasFeature);
    D.quartet = D.quartet(hasFeature);
    D.availability = D.availability(hasFeature, :);
    fprintf('%s Figaro RF-on-curve trials: %d/%d across %d quartets.\n', ...
        name, nnz(hasFeature), numel(hasFeature), numel(unique(D.quartet)));
    base.(name) = D;
end
end

function base = loadFeatureSource(baseFile)
needsRebuild = true;
if isfile(baseFile)
    data = load(baseFile, 'OUT');
    base = data.OUT;
    needsRebuild = ~isfield(base.V1, 'normalizedResponse') || ...
        ~isfield(base.V1, 'availability') || ...
        ~isfield(base.V4, 'normalizedResponse') || ...
        ~isfield(base.V4, 'availability');
end

if needsRebuild
    fprintf('Rebuilding the shared RF-on-curve feature source.\n');
    base = Attention_decoding_sparse_linear_V1_V4( ...
        'MakeFigure', false, 'SaveOutputs', true);
end
end

function D = prepareColorConfigurationFeatures( ...
        source, region, cfg, responseWindow, chunkTrials, monkey)
if monkey == 'F'
    geometryDir = fullfile(fileparts(cfg.resultsRoot), 'Figaro files');
else
    geometryDir = cfg.matDir;
end
switch region
    case 'V1'
        tallData = load(fullfile(geometryDir, ...
            sprintf('Tall_V1_lines_%s.mat', monkey)), 'Tall_V1');
        Tall = tallData.Tall_V1;
    case 'V4'
        tallData = load(fullfile(geometryDir, ...
            sprintf('Tall_V4_lines_%s.mat', monkey)), 'Tall_V4');
        Tall = tallData.Tall_V4;
    otherwise
        error('Unknown region %s.', region);
end

stimulusNumbers = arrayfun(@(x) x.stimNum, Tall(:));
[sortedStimulusNumbers, order] = sort(stimulusNumbers(:));
assert(isequal(sortedStimulusNumbers(:)', 1:numel(Tall)), ...
    '%s Tall stimulus numbers must cover 1:%d.', region, numel(Tall));
Tall = Tall(order);

stimulus = double(source.stimulus(:));
position = mod(stimulus - 1, 8) + 1;
classLabel = -ones(numel(stimulus), 1);
classLabel(ismember(position, [1 3 6 8])) = 1;
targetColor = -ones(numel(stimulus), 1);
targetColor(position >= 5) = 1;
targetSide = ones(numel(stimulus), 1);
targetSide(mod(stimulus, 2) == 0) = -1;
assert(all(classLabel == -targetColor .* targetSide), ...
    '%s configuration labels do not match the orthogonal 2x2 design.', region);

nLocalSites = numel(source.siteIndexLocal);
nStimuli = numel(Tall);
colorRoleByStimulus = zeros(nLocalSites, nStimuli);
for stim = 1:nStimuli
    T = Tall(stim).T;
    assert(ismember('center_color', ...
        T.Properties.VariableNames), ...
        '%s Tall stimulus %d lacks center_color.', region, stim);
    color = string(T.center_color(source.siteIndexLocal));
    colorRoleByStimulus(color == "yellowArm", stim) = 1;
    colorRoleByStimulus(color == "purple", stim) = -1;
end

colorRole = colorRoleByStimulus(:, stimulus)';
availability = logical(source.availability);
assert(all(colorRole(availability) ~= 0), ...
    '%s has RF-on-curve entries without a yellow/purple label.', region);
canonicalColorSign = classLabel .* colorRole;
canonicalColorSign(~availability) = 0;
validateComplementColors(colorRoleByStimulus, region);
validateConfigurationSigns(colorRoleByStimulus, source, region);
validateTrialBalance(source.quartet, classLabel, targetColor, targetSide);

D = source;
D.normalizedResponse = readNormalizedResponses( ...
    source, cfg, responseWindow, chunkTrials, monkey);
D.responseWindowRequested = responseWindow;
D.classLabel = classLabel;
D.targetColor = targetColor;
D.targetSide = targetSide;
D.colorRole = colorRole;
D.canonicalSign = canonicalColorSign;
D.positiveClassStimulusPosition = 'positions [1 3 6 8] within each 8-stimulus block';
D.negativeClassStimulusPosition = 'positions [2 4 5 7] within each 8-stimulus block';
D.colorRoleDefinition = '+1 yellowArm, -1 purple, 0 background/unavailable';
D.featureDefinition = ['Training-fold quartet-centered normalized response ' ...
    'multiplied by color-configuration class x local RF-center color. ' ...
    'Background and overlap entries are zero.'];
D.colorSource = 'Tall(stim).T.center_color';
D.featureMatrix = makeFoldFeatures(D, true(numel(classLabel), 1));

variablePredictor = std(D.featureMatrix, 0, 1) > 0;
D.featureMatrix = D.featureMatrix(:, variablePredictor);
D.availability = D.availability(:, variablePredictor);
D.normalizedResponse = D.normalizedResponse(:, variablePredictor);
D.colorRole = D.colorRole(:, variablePredictor);
D.canonicalSign = D.canonicalSign(:, variablePredictor);
D.siteIndexLocal = D.siteIndexLocal(variablePredictor);
D.siteIndexGlobal = D.siteIndexGlobal(variablePredictor);
assert(all(isfinite(D.featureMatrix(:))), ...
    '%s color feature matrix contains non-finite values.', region);

fprintf(['%s color features: %d trials, %d quartets, %d candidate sites; ' ...
    'labels from Tall center_color; window %g-%g ms.\n'], ...
    region, numel(classLabel), numel(unique(D.quartet)), ...
    size(D.featureMatrix, 2), responseWindow(1), responseWindow(2));
end

function normalizedResponse = readNormalizedResponses( ...
        source, cfg, responseWindow, chunkTrials, monkey)
dataFolder = 'Mr Nilson';
if monkey == 'F'; dataFolder = 'Figaro'; end
dataDir = fullfile(cfg.dataRoot, dataFolder);
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
tb = double(m2.tb);
tb = tb(:)';
timeMask = tb >= responseWindow(1) & tb <= responseWindow(2);
timeSamples = find(timeMask);
assert(~isempty(timeSamples), ...
    'Response window %g-%g ms does not overlap tb.', ...
    responseWindow(1), responseWindow(2));

trialIndex = double(source.trialIndex(:));
globalSites = double(source.siteIndexGlobal(:));
localSites = double(source.siteIndexLocal(:));
switch source.region
    case 'V1'
        regionSites = 1:512;
    case 'V4'
        if monkey == 'F'
            regionSites = 833:1024;
        else
            regionSites = 513:768;
        end
    otherwise
        error('Unknown region %s.', source.region);
end
assert(isequal(globalSites, regionSites(localSites)'), ...
    '%s local/global site mapping is inconsistent.', source.region);
baseline = double(source.baseline(localSites));
responseScale = double(source.responseScale(localSites));
nTrials = numel(trialIndex);
nSites = numel(globalSites);
nRegionSites = numel(regionSites);
[~, nAllTrials, ~] = size(m1, 'normMUA');
rowByGlobalTrial = zeros(nAllTrials, 1);
rowByGlobalTrial(trialIndex) = 1:nTrials;
normalizedResponse = nan(nTrials, nSites);

for firstGlobalTrial = 1:chunkTrials:nAllTrials
    lastGlobalTrial = min(nAllTrials, firstGlobalTrial + chunkTrials - 1);
    globalRange = firstGlobalTrial:lastGlobalTrial;
    selected = rowByGlobalTrial(globalRange) > 0;
    if ~any(selected)
        continue;
    end
    rows = rowByGlobalTrial(globalRange(selected));
    raw = double(m1.normMUA(regionSites, globalRange, timeSamples));
    trialResponse = reshape(mean(raw, 3, 'omitnan'), ...
        nRegionSites, numel(globalRange));
    trialResponse = trialResponse(localSites, :);
    normalized = (trialResponse - baseline) ./ responseScale;
    normalizedResponse(rows, :) = normalized(:, selected)';
end

assert(all(any(isfinite(normalizedResponse), 2)), ...
    'Some selected trials have no finite responses in the requested window.');
end

function validateComplementColors(colorRoleByStimulus, region)
nMismatch = 0;
nCompared = 0;
for stimulusA = 1:size(colorRoleByStimulus, 2)
    if mod(stimulusA - 1, 8) + 1 <= 4
        stimulusB = stimulusA + 4;
        roleA = colorRoleByStimulus(:, stimulusA);
        roleB = colorRoleByStimulus(:, stimulusB);
        compare = roleA ~= 0 & roleB ~= 0;
        nMismatch = nMismatch + nnz(roleA(compare) ~= -roleB(compare));
        nCompared = nCompared + nnz(compare);
    end
end
assert(nMismatch == 0, ...
    '%s has %d non-complementary RF color labels.', region, nMismatch);
fprintf('%s color-complement validation: %d comparisons, 0 mismatches.\n', ...
    region, nCompared);
end

function validateConfigurationSigns(colorRoleByStimulus, source, region)
nStimuli = size(colorRoleByStimulus, 2);
nSites = size(colorRoleByStimulus, 1);
activeByStimulus = false(nSites, nStimuli);
for stimulus = unique(source.stimulus(:))'
    rows = find(source.stimulus == stimulus);
    reference = logical(source.availability(rows(1), :));
    assert(all(all(source.availability(rows, :) == reference, 2)), ...
        '%s availability varies across trials of stimulus %d.', region, stimulus);
    activeByStimulus(:, stimulus) = reference(:);
end

classByStimulus = -ones(1, nStimuli);
position = mod((1:nStimuli) - 1, 8) + 1;
classByStimulus(ismember(position, [1 3 6 8])) = 1;
canonicalByStimulus = colorRoleByStimulus .* classByStimulus;
nCompared = 0;
nMismatch = 0;
for quartet = unique(source.quartet(:))'
    members = quartetMembers(quartet);
    reference = members(1);
    for member = members(2:end)
        compare = activeByStimulus(:, reference) & ...
            activeByStimulus(:, member);
        nCompared = nCompared + nnz(compare);
        nMismatch = nMismatch + nnz(canonicalByStimulus(compare, reference) ~= ...
            canonicalByStimulus(compare, member));
    end
end
assert(nMismatch == 0, ...
    '%s has %d inconsistent canonical color-configuration signs.', ...
    region, nMismatch);
fprintf('%s canonical configuration-sign validation: %d comparisons, 0 mismatches.\n', ...
    region, nCompared);
end

function validateTrialBalance(quartet, classLabel, targetColor, targetSide)
assert(all(classLabel == -targetColor .* targetSide), ...
    'Configuration labels do not match target color x target side.');
for quartetValue = unique(quartet(:))'
    for classValue = [-1 1]
        for colorValue = [-1 1]
            n = nnz(quartet == quartetValue & classLabel == classValue & ...
                targetColor == colorValue);
            assert(n >= 1, ...
                'Quartet %d, class %d, target color %d has no trials.', ...
                quartetValue, classValue, colorValue);
        end
        assert(all(ismember([-1 1], unique(targetSide( ...
            quartet == quartetValue & classLabel == classValue))')), ...
            'Quartet %d class %d does not balance target side.', ...
            quartetValue, classValue);
    end
end
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

function R = crossValidateRegion(D, opt, seed)
region = D.region;
fprintf('\n%s color elastic-net cross-validation: %d trials, %d sites.\n', ...
    region, numel(D.classLabel), size(D.normalizedResponse, 2));

fullTrainingMask = true(numel(D.classLabel), 1);
fullX = makeFoldFeatures(D, fullTrainingMask);
maxDifference = max(abs(fullX(:) - D.featureMatrix(:)));
assert(maxDifference < 1e-10, ...
    '%s full-data quartet centering differs by %.3g.', region, maxDifference);

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
        D.quartet, D.classLabel, D.targetColor, ...
        opt.NumFolds, seed + repeatIdx);
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
        correct = double(signClass(scores) == D.classLabel);
        correct(~isfinite(scores)) = NaN;
        accuracyByRepeat(alphaIdx, :, repeatIdx) = ...
            mean(correct, 1, 'omitnan');
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
crossValidatedScore = mean(selectedScores, 2, 'omitnan');
validCV = isfinite(crossValidatedScore);
crossValidatedAccuracy = mean(signClass(crossValidatedScore(validCV)) == ...
    D.classLabel(validCV));

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
R.targetColor = D.targetColor;
R.targetSide = D.targetSide;
R.colorRole = D.colorRole;
R.siteIndexLocal = D.siteIndexLocal;
R.siteIndexGlobal = D.siteIndexGlobal;
R.availability = D.availability;
R.foldIndex = foldIndex;
R.alwaysTrainingTrial = all(foldIndex == 0, 2);
R.nCrossValidatedTrials = nnz(validCV);
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
R.crossValidatedConfigurationAlignedScore = ...
    D.classLabel .* crossValidatedScore;
R.crossValidatedAccuracyOfMeanScore = crossValidatedAccuracy;
R.finalCoefficient = finalB;
R.finalIntercept = finalStats.Intercept;
R.finalNonzero = finalNonzero;
R.finalClassScore = finalClassScore;
R.finalConfigurationAlignedScore = D.classLabel .* finalClassScore;
R.finalInSampleAccuracy = mean( ...
    signClass(finalClassScore) == D.classLabel);
R.finalClassScoreNormalized = normalizedScore;
R.finalConfigurationAlignedScoreNormalized = ...
    D.classLabel .* normalizedScore;
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
    cellMean = nan(2, 2, nSites);
    for classIdx = 1:2
        classValue = [-1 1];
        for colorIdx = 1:2
            colorValue = [-1 1];
            use = trainingMask & D.quartet == quartet & ...
                D.classLabel == classValue(classIdx) & ...
                D.targetColor == colorValue(colorIdx);
            assert(any(use), ...
                ['Training fold lacks quartet %d, class %d, ' ...
                'target color %d.'], ...
                quartet, classValue(classIdx), colorValue(colorIdx));
            cellMean(classIdx, colorIdx, :) = ...
                mean(D.normalizedResponse(use, :), 1, 'omitnan');
        end
    end
    trainingMean(quartet, :) = ...
        reshape(mean(cellMean, [1 2]), 1, nSites);
end

centered = D.normalizedResponse - trainingMean(D.quartet, :);
X = centered .* D.canonicalSign;
X(~D.availability | ~isfinite(X)) = 0;
end

function foldIndex = balancedFoldIndex( ...
        quartet, classLabel, targetColor, nFolds, seed)
rng(seed, 'twister');
foldIndex = zeros(numel(classLabel), 1);
for quartetValue = unique(quartet(:))'
    for classValue = [-1 1]
        for colorValue = [-1 1]
            trials = find(quartet == quartetValue & ...
                classLabel == classValue & targetColor == colorValue);
            if numel(trials) == 1
                foldIndex(trials) = 0;
                continue;
            end
            trials = trials(randperm(numel(trials)));
            foldOrder = randperm(nFolds);
            labels = repmat(foldOrder, 1, ceil(numel(trials) / nFolds));
            foldIndex(trials) = labels(1:numel(trials));
        end
    end
end
assert(all(foldIndex >= 0), 'Invalid fold assignment.');
for fold = 1:nFolds
    assert(all(ismember([-1 1], unique(classLabel(foldIndex == fold))')), ...
        'Fold %d does not contain both classes.', fold);
    isTrain = foldIndex ~= fold;
    for quartetValue = unique(quartet(:))'
        for classValue = [-1 1]
            for colorValue = [-1 1]
                assert(any(isTrain & quartet == quartetValue & ...
                    classLabel == classValue & targetColor == colorValue), ...
                    ['Training fold %d lacks quartet %d, class %d, ' ...
                    'target color %d.'], ...
                    fold, quartetValue, classValue, colorValue);
            end
        end
    end
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
predicted = nan(size(score));
predicted(isfinite(score) & score >= 0) = 1;
predicted(isfinite(score) & score < 0) = -1;
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
bar(accuracy, 'FaceColor', [0.48 0.35 0.62]);
hold on;
errorbar(1:2, accuracy, accuracySE, '.k', 'LineWidth', 1.4);
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('Cross-validated accuracy (%)');
ylim([50 100]);
title('One-SE selected models');
formatAxes(gca);

subplot(2, 2, 4);
siteCount = [OUT.V1.finalNonzero, OUT.V4.finalNonzero];
bar(siteCount, 'FaceColor', [0.86 0.68 0.16]);
set(gca, 'XTickLabel', {'V1', 'V4'});
ylabel('Nonzero sites in final fit');
title('Final decoder population');
formatAxes(gca);

sgtitle(sprintf(['%s elastic-net color-configuration decoder, %g-%g ms: %dx%d-fold ' ...
    'quartet-balanced cross-validation'], OUT.monkey, OUT.responseWindow(1), ...
    OUT.responseWindow(2), OUT.numRepeats, OUT.numFolds), ...
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

function stem = resultFileStem(responseWindow)
stem = 'Color_configuration_decoder_elastic_net_CV_V1_V4';
if ~isequal(responseWindow, [300 500])
    windowText = sprintf('_%gto%gms', responseWindow(1), responseWindow(2));
    windowText = strrep(windowText, '.', 'p');
    stem = [stem windowText];
end
stem = [stem '_quartetCentered_minSites20_N'];
end

function printSummary(OUT)
fprintf('\nColor-configuration elastic-net cross-validation summary\n');
regions = {'V1', 'V4'};
for idx = 1:numel(regions)
    R = OUT.(regions{idx});
    fprintf(['  %s: alpha %.3g, lambda %.4g; CV %.2f%% +/- %.2f%%; ' ...
        'mean fold sites %.1f; final sites %d; final in-sample %.2f%%; ' ...
        '%d/%d trials cross-validated.\n'], ...
        regions{idx}, R.selectedAlpha, R.selectedLambda, ...
        100 * R.selectedMeanAccuracy, 100 * R.selectedAccuracySE, ...
        R.selectedMeanNonzero, R.finalNonzero, ...
        100 * R.finalInSampleAccuracy, R.nCrossValidatedTrials, ...
        numel(R.trialIndex));
end
end
