function OUT = Attention_decoding_elastic_net_CV_IT(varargin)
%ATTENTION_DECODING_ELASTIC_NET_CV_IT Combined IT attention decoder.
%
% All valid IT sites contribute on every trial. For each site and quartet,
% a training-fold tuning model combines:
%   1) RF-center proximity to the nearest point on the target versus
%      distractor capsule;
%   2) tuning to the distractor-to-target screen direction.
%
% The fitted combination supplies one continuous alignment value per site
% and quartet. Elastic-net logistic regression then learns one decoder
% coefficient per site from the aligned, quartet-centered trial responses.

p = inputParser;
p.addParameter('Monkey','N',@(x) any(strcmpi(string(x),["N" "F"])));
p.addParameter('Config',struct(),@isstruct);
p.addParameter('CandidateSiteGlobal',[],@(x) isnumeric(x)&&isvector(x)||isempty(x));
p.addParameter('ResponseWindow', [300 500], ...
    @(x) isnumeric(x) && isequal(size(x), [1 2]) && x(1) < x(2));
p.addParameter('Days', [1 2], ...
    @(x) isnumeric(x) && isvector(x) && ~isempty(x));
p.addParameter('OnlyCorrect', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SiteSelection', 'functional', ...
    @(x) any(strcmpi(string(x), ["functional" "visual"])));
p.addParameter('AlignmentMode', 'tuning', ...
    @(x) any(strcmpi(string(x), ["tuning" "rfRole"])));
p.addParameter('ResponseSNRThreshold', 0.7, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ResponseTopKQuartets', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
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
p.addParameter('VarianceFloorFraction', 1e-3, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('CoefficientTolerance', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('ChunkTrials', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
opt.AlphaGrid = unique(double(opt.AlphaGrid(:)'));
siteSelection = lower(string(opt.SiteSelection));

assert(exist('lassoglm', 'file') == 2, ...
    'Attention_decoding_elastic_net_CV_IT requires lassoglm.');

cfg = config();
if ~isempty(fieldnames(opt.Config)); cfg=opt.Config; end
opt.Monkey=upper(char(opt.Monkey));
opt.AlignmentMode=char(opt.AlignmentMode);
monkeyName='Mr Nilson'; if opt.Monkey=='F'; monkeyName='Figaro'; end
D = buildITData(cfg, opt);

fprintf('\nIT combined attention decoder: %d trials, %d sites, %d quartets.\n', ...
    numel(D.classLabel), size(D.normalizedResponse, 2), ...
    numel(unique(D.quartet)));
fprintf('  Valid RF centers for proximity: %d / %d sites.\n', ...
    nnz(D.hasRFcenter), numel(D.hasRFcenter));
if ~isempty(opt.CandidateSiteGlobal)
    fprintf('  Fixed previous candidate pool: %d / %d sites.\n',D.nIncludedSites,D.nAllITSites);
elseif siteSelection == "functional"
    fprintf(['  Functional criterion: visual OR attention OR color; ' ...
        '%d / %d sites retained.\n'], D.nIncludedSites, D.nAllITSites);
else
    fprintf(['  Visual-response criterion: top-%d quartet |SNR| > %.2f; ' ...
        '%d / %d sites retained.\n'], opt.ResponseTopKQuartets, ...
        opt.ResponseSNRThreshold, D.nIncludedSites, D.nAllITSites);
end
fprintf('  Response availability: every included site is present on every trial.\n');

R = crossValidateIT(D, opt);

OUT = struct();
OUT.description = ['Combined nearest-capsule RF proximity and target-side ' ...
    'direction IT attention decoder with one elastic-net coefficient per site.'];
OUT.alignmentMode = opt.AlignmentMode;
OUT.monkey = monkeyName;
OUT.responseWindow = opt.ResponseWindow;
OUT.days = opt.Days;
OUT.onlyCorrect = opt.OnlyCorrect;
OUT.siteSelection = char(siteSelection);
OUT.responseSNRThreshold = opt.ResponseSNRThreshold;
OUT.responseTopKQuartets = opt.ResponseTopKQuartets;
if siteSelection == "functional"
    OUT.siteInclusion = ['Valid spontaneous normalization, complete trial ' ...
        'responses, and visual OR attention OR color functional criterion; ' ...
        'only triple-null sites are excluded.'];
    OUT.siteSelectionCaveat = ['Attention and color inclusion masks were ' ...
        'estimated from the full dataset and are not nested within the ' ...
        'decoder cross-validation.'];
else
    OUT.siteInclusion = ['Valid spontaneous normalization, complete trial ' ...
        'responses, and historical top-K absolute quartet-response SNR criterion.'];
    OUT.siteSelectionCaveat = '';
end
if ~isempty(opt.CandidateSiteGlobal)
    OUT.siteSelection='fixedCandidates';
    OUT.siteInclusion='Explicit previously selected candidate IDs; valid normalization and complete responses.';
    OUT.siteSelectionCaveat='Candidate IDs are inherited from the original Figaro comparison.';
    OUT.candidateSiteGlobal=double(opt.CandidateSiteGlobal(:));
end
OUT.responseCentering = 'training-fold quartet means';
OUT.responseNormalization = 'spontaneous SD';
OUT.geometryModel = ['Per-site additive model of A-minus-B attention response: ' ...
    'cos(D-to-T angle), sin(D-to-T angle), and signed RF proximity to the ' ...
    'nearest target versus distractor capsule points.'];
OUT.proximityDefinition = ['d(RF, nearest distractor capsule point) - ' ...
    'd(RF, nearest target capsule point); positive means target is nearer.'];
OUT.alignmentNormalization = ['Predicted tuning excluding the intercept, divided ' ...
    'by its RMS across quartets.'];
if strcmpi(opt.AlignmentMode,'rfRole')
    OUT.description = ['RF-role-only IT attention decoder with one ' ...
        'elastic-net coefficient per site.'];
    OUT.geometryModel = ['RF on target: +1; RF on distractor: -1; ' ...
        'background or overlap: 0. No fitted direction or proximity terms.'];
    OUT.alignmentNormalization = 'Binary RF role; no RMS normalization.';
end
OUT.classDefinition = ['Within each quartet, odd stimulus numbers are class +1 ' ...
    'and even stimulus numbers are class -1; the two color-swapped stimuli ' ...
    'within each class are pooled.'];
OUT.alphaGrid = opt.AlphaGrid;
OUT.numLambda = opt.NumLambda;
OUT.lambdaRatio = opt.LambdaRatio;
OUT.numFolds = opt.NumFolds;
OUT.numRepeats = opt.NumRepeats;
OUT.randomSeed = opt.RandomSeed;
OUT.selectionRule = ['Sparsest alpha/lambda combination within one standard ' ...
    'error of the highest mean cross-validated accuracy.'];
OUT.IT = R;
if isfield(D, 'functionalSelection')
    OUT.ITFunctionalSelection = D.functionalSelection;
end
OUT.IT.targetDirectionDegByQuartet = D.targetDirectionDegByQuartet;
OUT.IT.proximityAdvantageBySiteQuartet = D.proximityAdvantageBySiteQuartet;
OUT.IT.hasRFcenter = D.hasRFcenter;
if strcmpi(opt.AlignmentMode,'rfRole')
    OUT.IT.rfRoleBySiteQuartet = D.rfRoleBySiteQuartet;
    OUT.IT.spontaneousSD = D.noiseSD;
end

if siteSelection == "functional"
    stem = ['Attention_decoder_elastic_net_CV_IT_combined_' ...
        'nearestCurve_ITfunctional_N'];
else
    thresholdTag = strrep(sprintf('%.2g', opt.ResponseSNRThreshold), '.', 'p');
    stem = sprintf(['Attention_decoder_elastic_net_CV_IT_combined_' ...
        'nearestCurve_respSNR%s_N'], thresholdTag);
end
stem=[stem(1:end-1) opt.Monkey];
if ~isempty(opt.CandidateSiteGlobal); stem=[stem '_fixedCandidates']; end
if strcmpi(opt.AlignmentMode,'rfRole'); stem=[stem '_rfRoleOnly']; end
if opt.Monkey=='F' && isequal(double(opt.Days(:).'),[1 2 3])
    stem = [stem '_days123'];
end
resultFile = fullfile(cfg.resultsDir, [stem '.mat']);
figureFile = fullfile(cfg.resultsDir, [stem '.png']);
OUT.resultFile = resultFile;
OUT.figureFile = figureFile;

fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT);
end

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

function D = buildITData(cfg, opt)
nStimuli = 384; if opt.Monkey=='F'; nStimuli=216; end
nQuartets = nStimuli/4;

geometryFile = fullfile(cfg.matDir, ['Tall_IT_lines_' opt.Monkey '.mat']);
responseSummaryFile = fullfile(cfg.matDir, ['SNR_capsules_' opt.Monkey '_d12.mat']);
if opt.Monkey=='F' && ~isfile(geometryFile)
    geometryFile = fullfile(fileparts(cfg.resultsRoot), ...
        'Figaro files','Tall_IT_lines_F.mat');
end
if opt.Monkey=='F' && ~isfile(responseSummaryFile)
    responseSummaryFile = fullfile(cfg.resultsDir, ...
        'SNR_capsules_F_d12.mat');
end
assert(isfile(geometryFile), 'Missing %s. Run Line_Stimuli_IT first.', geometryFile);
assert(isfile(responseSummaryFile), 'Missing %s.', responseSummaryFile);

if strcmpi(opt.AlignmentMode,'rfRole')
    G = load(geometryFile, 'Tall_IT', 'RFrange');
else
    G = load(geometryFile, 'Tall_IT', 'ALLCOORDS', ...
        'RTAB384', 'RFrange');
end
S = load(responseSummaryFile, 'R');
assert(numel(G.Tall_IT) == nStimuli, 'Geometry stimulus count differs from the selected monkey.');
assert(numel(G.RFrange) > 0, 'Tall_IT_lines_N.mat has an empty RFrange.');

globalSites = double(G.RFrange(:));
nSites = numel(globalSites);
assert(all(size(S.R.meanAct, 1) >= globalSites), ...
    'The response summary does not contain every IT channel.');

timeIdx = find(all(abs(double(S.R.timeWindows) - opt.ResponseWindow) < 1e-9, 2), 1);
if isempty(timeIdx)
    error('Requested response window %s is absent from SNR_capsules_N_d12.mat.', ...
        mat2str(opt.ResponseWindow));
end

snr = compute_snr_per_color_region(S.R, G.Tall_IT, globalSites);
baselineAll = double(snr.muSpont(:));
noiseSDAll = double(snr.sdSpont(:));
validNoise = isfinite(baselineAll) & isfinite(noiseSDAll) & noiseSDAll > 0;
[visualResponseScore, visuallyResponsive] = computeVisualResponseSelection( ...
    S.R, globalSites, baselineAll, noiseSDAll, ...
    opt.ResponseTopKQuartets, opt.ResponseSNRThreshold, nQuartets);

if strcmpi(opt.AlignmentMode,'rfRole')
    rfRole = buildRFRole(G.Tall_IT, nQuartets, nSites);
    targetDirectionDeg = nan(nQuartets,1);
    proximityAdvantage = nan(nSites,nQuartets);
    hasRFcenter = any(rfRole~=0,2);
    geometryQC = struct('nSitesWithRFcenter',nnz(hasRFcenter), ...
        'nSitesWithoutRFcenter',nnz(~hasRFcenter));
else
    [targetDirectionDeg, proximityAdvantage, hasRFcenter, geometryQC] = ...
        buildGeometry(G.Tall_IT, G.ALLCOORDS, G.RTAB384, nQuartets);
end

monkeyName='Mr Nilson'; dayCol=11;
if opt.Monkey=='F'; monkeyName='Figaro'; dayCol=10; end
dataDir = fullfile(cfg.dataRoot, monkeyName);
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:)';
[nChannels, nAllTrials, nTimes] = size(m1, 'normMUA');

assert(nChannels >= max(globalSites), 'normMUA does not contain every IT channel.');
assert(size(ALLMAT, 1) == nAllTrials, 'ALLMAT and normMUA trial counts differ.');
assert(numel(tb) == nTimes, 'tb and normMUA time dimensions differ.');
assert(size(ALLMAT, 2) >= dayCol, 'Missing session column in ALLMAT.');

stimulusAll = double(ALLMAT(:, 1));
quartetAll = stimulusToQuartet(stimulusAll);
classAll = ones(nAllTrials, 1);
classAll(mod(stimulusAll, 2) == 0) = -1;

includeTrial = ismember(ALLMAT(:, dayCol), opt.Days(:));
if opt.OnlyCorrect
    includeTrial = includeTrial & ALLMAT(:, 9) == 1;
end
includeTrial = includeTrial & isfinite(stimulusAll) & ...
    stimulusAll >= 1 & stimulusAll <= nStimuli & ...
    stimulusAll == floor(stimulusAll);

timeMask = tb >= opt.ResponseWindow(1) & tb <= opt.ResponseWindow(2);
assert(any(timeMask), 'Response window %s does not overlap tb.', ...
    mat2str(opt.ResponseWindow));
timeSamples = find(timeMask);

trialIndex = find(includeTrial);
nTrials = numel(trialIndex);
normalizedAll = nan(nAllTrials, nSites);

fprintf('Reading %d IT trials in chunks of %d (%g-%g ms).\n', ...
    nTrials, opt.ChunkTrials, tb(timeSamples(1)), tb(timeSamples(end)));
for firstTrial = 1:opt.ChunkTrials:nAllTrials
    lastTrial = min(nAllTrials, firstTrial + opt.ChunkTrials - 1);
    globalRange = firstTrial:lastTrial;
    localInclude = includeTrial(globalRange);
    if ~any(localInclude)
        continue;
    end

    raw = double(m1.normMUA(globalSites, globalRange, timeSamples));
    trialResponse = reshape(mean(raw, 3, 'omitnan'), nSites, numel(globalRange));
    normalized = (trialResponse - baselineAll) ./ noiseSDAll;
    normalizedAll(globalRange(localInclude), :) = normalized(:, localInclude)';
end

responseComplete = all(isfinite(normalizedAll(trialIndex, :)), 1)';
if ~isempty(opt.CandidateSiteGlobal)
    functionalSelection=[];
    baseSelection=ismember(globalSites,opt.CandidateSiteGlobal);
elseif lower(string(opt.SiteSelection)) == "functional"
    functionalSelection = IT_functional_inclusion_N('Monkey',opt.Monkey, ...
        'Config',cfg,'VisualSNRThreshold',opt.ResponseSNRThreshold);
    [found, position] = ismember(globalSites, ...
        functionalSelection.siteIndexGlobal);
    assert(all(found), 'The decoder and functional IT site lists differ.');
    baseSelection = functionalSelection.included(position);
else
    functionalSelection = [];
    baseSelection = visuallyResponsive;
end
eligibleSite = validNoise & responseComplete & baseSelection;
if strcmpi(opt.AlignmentMode,'rfRole')
    eligibleSite = eligibleSite & any(rfRole ~= 0, 2);
end
if any(validNoise & ~responseComplete)
    fprintf('Excluding %d IT sites with non-finite trial responses.\n', ...
        nnz(validNoise & ~responseComplete));
end
assert(any(eligibleSite), 'No IT sites have valid normalization and complete trial responses.');

normalizedResponse = normalizedAll(trialIndex, eligibleSite);
assert(all(isfinite(normalizedResponse(:))), ...
    'The selected IT response matrix contains non-finite values.');

stimulus = stimulusAll(trialIndex);
quartet = quartetAll(trialIndex);
classLabel = classAll(trialIndex);
validateTrialBalance(quartet, classLabel, opt.NumFolds, nQuartets);

D = struct();
D.region = 'IT';
D.trialIndex = trialIndex;
D.stimulus = stimulus;
D.quartet = quartet;
D.classLabel = classLabel;
D.siteIndexLocal = find(eligibleSite);
D.siteIndexGlobal = globalSites(eligibleSite);
D.normalizedResponse = normalizedResponse;
D.baseline = baselineAll(eligibleSite);
D.noiseSD = noiseSDAll(eligibleSite);
D.visualResponseScore = visualResponseScore(eligibleSite);
D.visuallyResponsiveAllSites = visuallyResponsive;
D.targetDirectionDegByQuartet = targetDirectionDeg;
D.proximityAdvantageBySiteQuartet = proximityAdvantage(eligibleSite, :);
D.proximityZByQuartetSite = zscoreColumns(proximityAdvantage(eligibleSite, :)');
if strcmpi(opt.AlignmentMode,'rfRole')
    D.rfRoleBySiteQuartet = rfRole(eligibleSite, :);
end
D.hasRFcenter = hasRFcenter(eligibleSite);
D.geometryQC = geometryQC;
D.responseWindowSampled = [tb(timeSamples(1)), tb(timeSamples(end))];
D.nAllITSites = nSites;
D.nValidNoiseSites = nnz(validNoise);
D.nCompleteSites = nnz(responseComplete);
D.nVisuallyResponsiveSites = nnz(visuallyResponsive & validNoise);
D.nFunctionallyIncludedSites = nnz(baseSelection & validNoise);
D.nIncludedSites = nnz(eligibleSite);
if ~isempty(functionalSelection)
    D.functionalSelection = functionalSelection;
end
end

function [score, responsive] = computeVisualResponseSelection( ...
        R, globalSites, baseline, noiseSD, topK, threshold, nQuartets)
nSites = numel(globalSites);
quartetEarly = nan(nSites, nQuartets);
quartetLate = nan(nSites, nQuartets);
nTrials = double(R.nTrials);

for site = 1:nSites
    globalSite = globalSites(site);
    if isvector(nTrials)
        nTrialSite = double(nTrials(:));
    else
        nTrialSite = double(nTrials(globalSite, :))';
    end
    early = squeeze(double(R.meanAct(globalSite, :, 2)))';
    late = squeeze(double(R.meanAct(globalSite, :, 3)))';

    for q = 1:nQuartets
        members = quartetMembers(q);
        quartetEarly(site, q) = weightedStimulusMean( ...
            early(members), nTrialSite(members));
        quartetLate(site, q) = weightedStimulusMean( ...
            late(members), nTrialSite(members));
    end
end

quartetEarly = (quartetEarly - baseline) ./ noiseSD;
quartetLate = (quartetLate - baseline) ./ noiseSD;
topEarly = topKAbsoluteMean(quartetEarly, topK);
topLate = topKAbsoluteMean(quartetLate, topK);
score = max([topEarly, topLate], [], 2, 'omitnan');
responsive = isfinite(score) & score > threshold;
end

function mu = weightedStimulusMean(response, nTrials)
valid = isfinite(response) & isfinite(nTrials) & nTrials > 0;
if any(valid)
    mu = sum(response(valid) .* nTrials(valid)) / sum(nTrials(valid));
else
    mu = NaN;
end
end

function value = topKAbsoluteMean(responseByQuartet, topK)
nSites = size(responseByQuartet, 1);
value = nan(nSites, 1);
for site = 1:nSites
    x = abs(responseByQuartet(site, :));
    x = sort(x(isfinite(x)), 'descend');
    if ~isempty(x)
        value(site) = mean(x(1:min(topK, numel(x))));
    end
end
end

function [thetaDeg, proximityAdvantage, hasRFcenter, QC] = ...
        buildGeometry(TallIT, ALLCOORDS, RTAB384, nQuartets)
nSites = height(TallIT(1).T);
nStimuli = numel(TallIT);
assert(nStimuli == 4 * nQuartets, 'Expected four stimuli per quartet.');
assert(size(RTAB384, 1) == nStimuli, ...
    'RTAB384 and Tall_IT must contain the same number of stimuli.');

xPx = double(TallIT(1).T.x_px(:));
yPx = double(TallIT(1).T.y_px(:));
imageWidth = 1024;
imageHeight = 768;
xSample = round(xPx);
ySample = round(yPx);
hasRFcenter = isfinite(xPx) & isfinite(yPx) & ...
    xSample >= 1 & xSample <= imageWidth & ...
    ySample >= 1 & ySample <= imageHeight;
sampleIndex = nan(nSites, 1);
sampleIndex(hasRFcenter) = sub2ind([imageHeight imageWidth], ...
    ySample(hasRFcenter), xSample(hasRFcenter));

targetAdvByStimulus = nan(nSites, nStimuli);
targetDirectionByStimulus = nan(nStimuli, 1);
fprintf('Computing nearest target/distractor capsule distances for %d stimuli.\n', ...
    nStimuli);
for stim = 1:nStimuli
    f = sprintf('stim_%d', stim);
    targetEnd = double(ALLCOORDS.(f).t_fig(:))';
    distractorEnd = double(ALLCOORDS.(f).t_back(:))';

    [~, masks] = render_stim_with_masks2( ...
        ALLCOORDS, RTAB384, stim, 'DrawDots', false);
    distanceTarget = bwdist(masks.figArm);
    distanceDistractor = bwdist(masks.backArm);
    targetAdvByStimulus(hasRFcenter, stim) = ...
        distanceDistractor(sampleIndex(hasRFcenter)) - ...
        distanceTarget(sampleIndex(hasRFcenter));
    targetDirectionByStimulus(stim) = mod(atan2d( ...
        targetEnd(2) - distractorEnd(2), ...
        targetEnd(1) - distractorEnd(1)), 360);
end

thetaDeg = nan(nQuartets, 1);
proximityAdvantage = nan(nSites, nQuartets);
maxPairDirectionError = 0;
maxOppositeDirectionError = 0;
for q = 1:nQuartets
    members = quartetMembers(q);
    pairA = members([1 3]);
    pairB = members([2 4]);
    thetaA = targetDirectionByStimulus(pairA);
    thetaB = targetDirectionByStimulus(pairB);

    pairDirectionError = max(abs(angleDifference(thetaA, thetaA(1))));
    oppositeDirectionError = max(abs(abs(angleDifference(thetaB, thetaA(1))) - 180));
    maxPairDirectionError = max(maxPairDirectionError, pairDirectionError);
    maxOppositeDirectionError = max(maxOppositeDirectionError, oppositeDirectionError);
    assert(pairDirectionError < 1e-6, ...
        'Class +1 target directions differ within quartet %d.', q);
    assert(oppositeDirectionError < 1e-6, ...
        'Target directions are not opposite within quartet %d.', q);

    thetaDeg(q) = thetaA(1);
    proxA = mean(targetAdvByStimulus(:, pairA), 2, 'omitnan');
    proxB = mean(targetAdvByStimulus(:, pairB), 2, 'omitnan');
    proximityAdvantage(:, q) = 0.5 * (proxA - proxB);
end

proximityAdvantage(~hasRFcenter, :) = NaN;
QC = struct();
QC.maxWithinPairDirectionErrorDeg = maxPairDirectionError;
QC.maxOppositeDirectionErrorDeg = maxOppositeDirectionError;
QC.nSitesWithRFcenter = nnz(hasRFcenter);
QC.nSitesWithoutRFcenter = nnz(~hasRFcenter);
QC.proximityDefinition = ['signed difference between distance to the nearest ' ...
    'point on the complete distractor and target capsule masks'];
QC.capsuleMaskDefinition = ['Complete figArm/target and backArm/distractor ' ...
    'capsules, including their geometric overlap.'];
end

function R = crossValidateIT(D, opt)
fullTrainingMask = true(numel(D.classLabel), 1);
[fullX, fullTuning] = makeFoldFeatures(D, fullTrainingMask, opt);
assert(all(isfinite(fullX(:))), 'Full IT feature matrix contains non-finite values.');

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
cvTimer = tic;
completedFolds = 0;
totalFolds = opt.NumRepeats * opt.NumFolds;

for repeatIdx = 1:opt.NumRepeats
    foldIndex(:, repeatIdx) = balancedFoldIndex( ...
        D.quartet, D.classLabel, opt.NumFolds, opt.RandomSeed + repeatIdx);
    fprintf('  repeat %d/%d\n', repeatIdx, opt.NumRepeats);
    for fold = 1:opt.NumFolds
        isTest = foldIndex(:, repeatIdx) == fold;
        isTrain = ~isTest;
        X = makeFoldFeatures(D, isTrain, opt);
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
        completedFolds = completedFolds + 1;
        elapsed = toc(cvTimer);
        eta = elapsed * (totalFolds / completedFolds - 1);
        fprintf(['    fold %d/%d complete (%d/%d total); elapsed %.0f s, ' ...
            'ETA %.0f s\n'], fold, opt.NumFolds, completedFolds, ...
            totalFolds, elapsed, eta);
    end
end

accuracyByRepeat = nan(nAlpha, nLambda, opt.NumRepeats);
accuracyByFold = nan(nAlpha, nLambda, opt.NumRepeats, opt.NumFolds);
for repeatIdx = 1:opt.NumRepeats
    for alphaIdx = 1:nAlpha
        scores = scoreByTrial(:, :, alphaIdx, repeatIdx);
        accuracyByRepeat(alphaIdx, :, repeatIdx) = ...
            mean(signClass(scores) == D.classLabel, 1);
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
crossValidatedScore = mean(selectedScores, 2);

selectedAlpha = opt.AlphaGrid(alphaIdx);
selectedLambda = lambdaGrid(alphaIdx, lambdaIdx);
[finalB, finalStats] = lassoglm(fullX, D.classLabel == 1, 'binomial', ...
    'Alpha', selectedAlpha, 'Lambda', selectedLambda, ...
    'Standardize', true, 'Weights', fullWeights, 'CV', 'resubstitution');
finalClassScore = fullX * finalB + finalStats.Intercept;

assert(numel(finalB) == size(D.normalizedResponse, 2), ...
    'The final decoder must contain exactly one coefficient per included IT site.');

R = struct();
R.region = D.region;
R.trialIndex = D.trialIndex;
R.stimulus = D.stimulus;
R.quartet = D.quartet;
R.classLabel = D.classLabel;
R.siteIndexLocal = D.siteIndexLocal;
R.siteIndexGlobal = D.siteIndexGlobal;
R.nAllITSites = D.nAllITSites;
R.nValidNoiseSites = D.nValidNoiseSites;
R.nCompleteSites = D.nCompleteSites;
R.nVisuallyResponsiveSites = D.nVisuallyResponsiveSites;
R.nFunctionallyIncludedSites = D.nFunctionallyIncludedSites;
R.nIncludedSites = D.nIncludedSites;
R.visualResponseScore = D.visualResponseScore;
R.foldIndex = foldIndex;
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
R.crossValidatedTargetAlignedScore = D.classLabel .* crossValidatedScore;
R.crossValidatedAccuracyOfMeanScore = mean( ...
    signClass(crossValidatedScore) == D.classLabel);
R.finalCoefficient = finalB;
R.finalIntercept = finalStats.Intercept;
R.finalNonzero = countNonzero(finalB, opt.CoefficientTolerance);
R.finalClassScore = finalClassScore;
R.finalTargetAlignedScore = D.classLabel .* finalClassScore;
R.finalInSampleAccuracy = mean(signClass(finalClassScore) == D.classLabel);
R.fullDataFeatureMatrix = fullX;
R.fullDataQuartetMeanResponse = fullTuning.quartetMeanResponse;
R.fullDataAlignmentByQuartetSite = fullTuning.alignmentByQuartetSite;
R.fullDataTuningCoefficient = fullTuning.coefficient;
R.fullDataTuningR2 = fullTuning.r2;
R.fullDataTuningStatus = fullTuning.status;
R.fullDataDirectionAmplitude = hypot( ...
    fullTuning.coefficient(:, 2), fullTuning.coefficient(:, 3));
R.fullDataDirectionPreferenceDeg = mod(atan2d( ...
    fullTuning.coefficient(:, 3), fullTuning.coefficient(:, 2)), 360);
R.fullDataProximityCoefficient = fullTuning.coefficient(:, 4);
R.responseWindowSampled = D.responseWindowSampled;
R.geometryQC = D.geometryQC;
R.oneCoefficientPerSite = true;

fprintf(['IT selected alpha %.3g, lambda %.4g (fraction %.4g): ' ...
    'CV %.2f%% +/- %.2f%%, mean %.1f sites; final %d sites.\n'], ...
    R.selectedAlpha, R.selectedLambda, R.selectedLambdaFraction, ...
    100 * R.selectedMeanAccuracy, 100 * R.selectedAccuracySE, ...
    R.selectedMeanNonzero, R.finalNonzero);
end

function [X, tuning] = makeFoldFeatures(D, trainingMask, opt)
nQuartets = max(D.quartet);
nSites = size(D.normalizedResponse, 2);
quartetMean = nan(nQuartets, nSites);
delta = nan(nQuartets, nSites);
deltaVariance = nan(nQuartets, nSites);

for q = 1:nQuartets
    inQuartet = trainingMask & D.quartet == q;
    positive = inQuartet & D.classLabel == 1;
    negative = inQuartet & D.classLabel == -1;
    assert(any(positive) && any(negative), ...
        'Training fold lacks one class for quartet %d.', q);

    quartetMean(q, :) = mean(D.normalizedResponse(inQuartet, :), 1);
    muPositive = mean(D.normalizedResponse(positive, :), 1);
    muNegative = mean(D.normalizedResponse(negative, :), 1);
    delta(q, :) = muPositive - muNegative;
    deltaVariance(q, :) = ...
        var(D.normalizedResponse(positive, :), 0, 1) ./ nnz(positive) + ...
        var(D.normalizedResponse(negative, :), 0, 1) ./ nnz(negative);
end

if strcmpi(opt.AlignmentMode,'rfRole')
    alignment = double(D.rfRoleBySiteQuartet.');
    centered = D.normalizedResponse - quartetMean(D.quartet, :);
    X = centered .* alignment(D.quartet, :);
    tuning = struct('quartetMeanResponse',quartetMean, ...
        'coefficient',zeros(nSites,4), ...
        'alignmentByQuartetSite',alignment, ...
        'r2',nan(nSites,1), ...
        'status',repmat("rf_role",nSites,1));
    return
end

theta = double(D.targetDirectionDegByQuartet(:));
directionDesign = [cosd(theta), sind(theta)];
coefficient = zeros(nSites, 4);
alignment = zeros(nQuartets, nSites);
r2 = nan(nSites, 1);
status = strings(nSites, 1);

for site = 1:nSites
    proximity = D.proximityZByQuartetSite(:, site);
    hasProximity = any(isfinite(proximity));
    if hasProximity
        proximity(~isfinite(proximity)) = 0;
        design = [ones(nQuartets, 1), directionDesign, proximity];
    else
        design = [ones(nQuartets, 1), directionDesign];
    end

    y = delta(:, site);
    variance = deltaVariance(:, site);
    valid = all(isfinite(design), 2) & isfinite(y) & ...
        isfinite(variance) & variance >= 0;
    if nnz(valid) < size(design, 2) + 1
        status(site) = "too_few_quartets";
        continue;
    end

    positiveVariance = variance(valid & variance > 0);
    if isempty(positiveVariance)
        varianceFloor = 1;
    else
        varianceFloor = max(median(positiveVariance) * ...
            opt.VarianceFloorFraction, eps);
    end
    variance = max(variance, varianceFloor);
    sqrtWeight = sqrt(1 ./ variance(valid));
    weightedDesign = design(valid, :) .* sqrtWeight;
    weightedResponse = y(valid) .* sqrtWeight;
    beta = weightedDesign \ weightedResponse;

    coefficient(site, 1:numel(beta)) = beta(:)';
    predicted = design * beta;
    predictedNoIntercept = predicted - beta(1);
    rmsPrediction = sqrt(mean(predictedNoIntercept(valid).^2));
    if ~(isfinite(rmsPrediction) && rmsPrediction > eps)
        status(site) = "constant_prediction";
        continue;
    end
    alignment(:, site) = predictedNoIntercept ./ rmsPrediction;

    w = 1 ./ variance(valid);
    weightedMean = sum(w .* y(valid)) / sum(w);
    rssNull = sum(w .* (y(valid) - weightedMean).^2);
    rssFull = sum(w .* (y(valid) - predicted(valid)).^2);
    if rssNull > 0
        r2(site) = 1 - rssFull / rssNull;
    end
    status(site) = "ok";
end

centered = D.normalizedResponse - quartetMean(D.quartet, :);
X = centered .* alignment(D.quartet, :);
X(~isfinite(X)) = 0;

tuning = struct();
tuning.quartetMeanResponse = quartetMean;
tuning.deltaResponse = delta;
tuning.deltaVariance = deltaVariance;
tuning.coefficient = coefficient;
tuning.coefficientNames = {'intercept','directionCos','directionSin','proximity'};
tuning.alignmentByQuartetSite = alignment;
tuning.r2 = r2;
tuning.status = status;
end

function rfRole = buildRFRole(Tall, nQuartets, nSites)
rfRole = zeros(nSites,nQuartets);
for q = 1:nQuartets
    members = quartetMembers(q);
    for j = 1:4
        T = Tall(members(j)).T;
        assignment = string(T.assignment(:));
        assert(numel(assignment)==nSites && ismember('overlap', ...
            T.Properties.VariableNames));
        role = zeros(nSites,1);
        role(assignment=="target") = 1;
        role(assignment=="distractor") = -1;
        role(logical(T.overlap(:))) = 0;
        if mod(members(j),2)==0; role = -role; end
        if j==1
            rfRole(:,q) = role;
        else
            assert(isequal(rfRole(:,q),role), ...
                'RF roles differ within quartet %d.',q);
        end
    end
end
end

function Z = zscoreColumns(X)
Z = nan(size(X));
for col = 1:size(X, 2)
    x = double(X(:, col));
    valid = isfinite(x);
    if nnz(valid) < 2
        continue;
    end
    sd = std(x(valid), 0);
    if isfinite(sd) && sd > 0
        Z(valid, col) = (x(valid) - mean(x(valid))) ./ sd;
    end
end
end

function quartet = stimulusToQuartet(stimulus)
quartet = nan(size(stimulus));
valid = isfinite(stimulus) & stimulus >= 1 & stimulus <= 384 & ...
    stimulus == floor(stimulus);
block = floor((stimulus(valid) - 1) / 8);
withinBlock = mod(stimulus(valid) - 1, 8) + 1;
quartet(valid) = 2 * block + 1 + (withinBlock >= 3 & withinBlock <= 4 | ...
    withinBlock >= 7 & withinBlock <= 8);
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

function validateTrialBalance(quartet, classLabel, ~, nQuartets)
assert(all(ismember(classLabel, [-1 1])), 'Class labels must be -1 or +1.');
assert(isequal(unique(quartet(:))', 1:nQuartets), ...
    'Included trials must cover all %d quartets.', nQuartets);
for q = 1:nQuartets
    for classValue = [-1 1]
        n = nnz(quartet == q & classLabel == classValue);
        assert(n >= 2, ...
            'Quartet %d class %d has %d trials; need at least 2.', ...
            q, classValue, n);
    end
end
end

function foldIndex = balancedFoldIndex(quartet, classLabel, nFolds, seed)
rng(seed, 'twister');
foldIndex = zeros(numel(classLabel), 1);
for quartetValue = unique(quartet(:))'
    for classValue = [-1 1]
        trials = find(quartet == quartetValue & classLabel == classValue);
        trials = trials(randperm(numel(trials)));
        foldOrder = randperm(nFolds);
        labels = repmat(foldOrder, 1, ceil(numel(trials) / nFolds));
        foldIndex(trials) = labels(1:numel(trials));
    end
end
assert(all(foldIndex > 0), 'Some trials were not assigned to a fold.');
for fold = 1:nFolds
    isTrain = foldIndex ~= fold;
    for quartetValue = unique(quartet(:))'
        assert(all(ismember([-1 1], ...
            unique(classLabel(isTrain & quartet == quartetValue))')), ...
            'Training fold %d lacks one class for quartet %d.', ...
            fold, quartetValue);
    end
end
end

function weights = classBalancingWeights(classLabel)
n = numel(classLabel);
nPositive = nnz(classLabel == 1);
nNegative = nnz(classLabel == -1);
assert(nPositive > 0 && nNegative > 0, 'Both classes must contain trials.');
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
predicted = ones(size(score));
predicted(score < 0) = -1;
end

function d = angleDifference(a, b)
d = mod(double(a) - double(b) + 180, 360) - 180;
end

function fig = makeFigure(OUT)
R = OUT.IT;
fig = figure('Color', 'w', 'Position', [80 80 1380 800]);
colors = lines(numel(R.alphaGrid));

subplot(2, 2, 1);
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
yline(50, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xlabel('Mean sites with nonzero weight');
ylabel('Cross-validated accuracy (%)');
ylim([50 100]);
legend('Location', 'southeast');
title(sprintf('Selected alpha %.2g, lambda fraction %.3g', ...
    R.selectedAlpha, R.selectedLambdaFraction));
formatAxes(gca);

subplot(2, 2, 2);
score = R.crossValidatedTargetAlignedScore;
histogram(score(isfinite(score)), 35, 'FaceColor', [0.24 0.48 0.58], ...
    'EdgeColor', 'none');
xline(0, ':k');
xlabel('Cross-validated target-aligned score S');
ylabel('Trials');
title(sprintf('Accuracy %.2f%%', 100 * R.crossValidatedAccuracyOfMeanScore));
formatAxes(gca);

subplot(2, 2, 3);
plot(R.siteIndexGlobal, R.finalCoefficient, '.', ...
    'Color', [0.72 0.28 0.16], 'MarkerSize', 11);
yline(0, ':', 'Color', [0.4 0.4 0.4]);
xlabel('Global IT site');
ylabel('Final elastic-net coefficient');
title(sprintf('%d / %d nonzero site weights', ...
    R.finalNonzero, R.nIncludedSites));
formatAxes(gca);

subplot(2, 2, 4);
directionAmplitude = R.fullDataDirectionAmplitude;
proximityCoefficient = abs(R.fullDataProximityCoefficient);
hasRF = R.hasRFcenter;
scatter(directionAmplitude(~hasRF), proximityCoefficient(~hasRF), 22, ...
    [0.65 0.65 0.65], 'filled', 'DisplayName', 'No RF center');
hold on;
scatter(directionAmplitude(hasRF), proximityCoefficient(hasRF), 22, ...
    [0.20 0.55 0.38], 'filled', 'DisplayName', 'RF center available');
xlabel('Direction tuning amplitude');
ylabel('|Proximity coefficient|');
title('Direction and nearest-curve components');
legend('Location', 'best');
formatAxes(gca);

sgtitle(sprintf(['%s IT nearest-curve attention decoder, %g-%g ms: ' ...
    '%dx%d-fold cross-validation'], OUT.monkey, OUT.responseWindow(1), ...
    OUT.responseWindow(2), OUT.numRepeats, OUT.numFolds), ...
    'FontWeight', 'bold');
end

function formatAxes(ax)
grid(ax, 'on');
box(ax, 'off');
set(ax, 'FontSize', 11, 'Layer', 'top');
end

function printSummary(OUT)
R = OUT.IT;
fprintf('\nCombined IT attention-decoder summary\n');
fprintf('  Included sites after all response-quality criteria: %d / %d.\n', ...
    R.nIncludedSites, R.nAllITSites);
fprintf('  Historical visual-response criterion retained: %d / %d sites.\n', ...
    R.nVisuallyResponsiveSites, R.nAllITSites);
if isfield(OUT, 'siteSelection') && strcmpi(OUT.siteSelection, 'functional')
    fprintf(['  Functional inclusion was computed on the full dataset and ' ...
        'is not nested within decoder cross-validation.\n']);
end
fprintf('  Sites with an RF center for proximity: %d / %d.\n', ...
    nnz(R.hasRFcenter), R.nIncludedSites);
fprintf('  Trials: %d; quartets: %d.\n', ...
    numel(R.trialIndex), numel(unique(R.quartet)));
fprintf(['  Selected alpha %.3g, lambda %.4g; CV %.2f%% +/- %.2f%%; ' ...
    'final nonzero weights %d.\n'], ...
    R.selectedAlpha, R.selectedLambda, 100 * R.selectedMeanAccuracy, ...
    100 * R.selectedAccuracySE, R.finalNonzero);
fprintf('  Final in-sample accuracy: %.2f%%.\n', ...
    100 * R.finalInSampleAccuracy);
end
