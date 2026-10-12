function OUT = Attention_decoding_Figaro_capsules(varargin)
%ATTENTION_DECODING_FIGARO_CAPSULES Fixed late-window RF-aligned readouts.
% Weights and scores use the same Figaro trials, so accuracy is in-sample.
% The saved weights, normalization and stimulus roles can be reused later
% with time-resolved responses without refitting the decoder.

p = inputParser;
p.addParameter('SNRThreshold', 0.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('MinSitesPerQuartet', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('ChunkTrials', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('V4QuartetRFOnly', false, ...
    @(x) islogical(x) && isscalar(x));
p.addParameter('ExcludeV4Site926', false, ...
    @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
assert(~opt.V4QuartetRFOnly || opt.SNRThreshold == 0.5, ...
    'The temporary V4 quartet-RF variant uses SNR > 0.5.');
assert(~opt.ExcludeV4Site926 || opt.V4QuartetRFOnly, ...
    'Site 926 exclusion is only defined for the temporary V4 variant.');

cfg = config();
if isempty(opt.OutputDir)
    outputDir = cfg.resultsDir;
else
    outputDir = char(opt.OutputDir);
end
geometryDir = fullfile(fileparts(cfg.resultsRoot), 'Figaro files');
responseFile = firstExisting({ ...
    fullfile(cfg.matDir, 'SNR_capsules_F_d12.mat'), ...
    fullfile(cfg.resultsDir, 'SNR_capsules_F_d12.mat')});
R = load(responseFile, 'R');
R = R.R;
assert(isequal(size(R.meanAct), [1024 216 3]), ...
    'Expected Figaro response summaries for 1024 sites, 216 stimuli, 3 windows.');
timeIdx = find(all(double(R.timeWindows) == [300 500], 2), 1);
assert(~isempty(timeIdx), 'Figaro response summary lacks 300-500 ms.');

dataDir = fullfile(cfg.dataRoot, 'Figaro');
m1 = matfile(fullfile(dataDir, 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(dataDir, 'ObjAtt_lines_MUA_trials.mat'));
ALLMAT = double(m2.ALLMAT);
tb = double(m2.tb);
tb = tb(:);
[nChannels, nTrials, nTimes] = size(m1, 'normMUA');
assert(nChannels == 1024 && size(ALLMAT, 1) == nTrials && ...
    size(ALLMAT, 2) == 10 && numel(tb) == nTimes, ...
    'Unexpected Figaro trial matrix dimensions.');
stimulus = ALLMAT(:, 1);
includeTrial = ALLMAT(:, 9) == 1 & isfinite(stimulus) & ...
    stimulus >= 1 & stimulus <= 216 & stimulus == floor(stimulus);
assert(isequal(accumarray(stimulus(includeTrial), 1, [216 1]), ...
    double(R.nTrials(:))), ...
    'Saved Figaro response summaries do not match the included trials.');
timeSamples = find(tb >= 300 & tb <= 500);
assert(~isempty(timeSamples), 'Figaro trial time base lacks 300-500 ms.');
trialQuartet = nan(nTrials, 1);
trialQuartet(includeTrial) = stimulusToQuartet(stimulus(includeTrial));

OUT = struct();
OUT.monkey = 'Figaro';
OUT.responseWindow = [300 500];
OUT.sampledWindow = [tb(timeSamples(1)) tb(timeSamples(end))];
OUT.responseCentering = 'quartet';
OUT.colorBalance = 'Matched target/distractor trial counts per color';
OUT.siteSelection = 'Best color/window visual SNR > threshold; valid attention d-prime';
OUT.v4QuartetRFOnly = opt.V4QuartetRFOnly;
OUT.excludedV4Site = [];
if opt.ExcludeV4Site926
    OUT.excludedV4Site = 926;
end
if opt.V4QuartetRFOnly
    OUT.siteSelection = ['V4: best color/window OR best quartet/window ' ...
        'visual SNR > 0.5; valid attention d-prime and response scale. ' ...
        'V1 and IT use their original selection.'];
end
OUT.scoreDefinition = ['Sum over selected sites of d-prime * RF role * ' ...
    'quartet-centered normalized response, divided by sum(abs(d-prime)).'];
OUT.accuracyDefinition = 'Fraction of included, scored trials with score > 0 (in-sample)';
OUT.latencyUse = ['Keep d-prime and RF roles fixed; recompute each ' ...
    'time bin''s site-by-quartet response mean before scoring.'];
OUT.snrThreshold = opt.SNRThreshold;
OUT.minSitesPerQuartet = opt.MinSitesPerQuartet;
OUT.nIncludedTrials = nnz(includeTrial);
OUT.daysIncluded = unique(ALLMAT(includeTrial, 10))';
OUT.onlyCorrect = true;
OUT.nStimuli = 216;
OUT.nQuartets = 54;
OUT.source = struct('geometryDir', geometryDir, 'responseFile', responseFile, ...
    'trialDataDir', dataDir);

regions = {'V1', 'V4', 'IT'};
for iRegion = 1:numel(regions)
    region = regions{iRegion};
    file = firstExisting({ ...
        fullfile(cfg.matDir, sprintf('Tall_%s_lines_F.mat', region)), ...
        fullfile(geometryDir, sprintf('Tall_%s_lines_F.mat', region))});
    geometry = load(file);
    Tall = geometry.(sprintf('Tall_%s', region));
    assert(numel(Tall) == 216 && ...
        isequal(sort([Tall.stimNum]), 1:216), ...
        '%s geometry does not cover Figaro stimuli 1:216.', region);
    [~, order] = sort([Tall.stimNum]);
    Tall = Tall(order);
    if strcmp(region, 'V1')
        globalSites = (1:512)';
    else
        globalSites = double(geometry.RFrange(:));
    end
    assert(height(Tall(1).T) == numel(globalSites), ...
        '%s geometry rows and site indices differ.', region);
    fprintf('Figaro %s: %d sites, %d included trials.\n', ...
        region, numel(globalSites), OUT.nIncludedTrials);
    OUT.(region) = decodeRegion(region, Tall, globalSites, R, ...
        m1, stimulus, trialQuartet, includeTrial, timeIdx, timeSamples, opt);
    OUT.(region).geometryFile = file;
end

baseName = 'Attention_decoder_Figaro_capsules';
if opt.SNRThreshold ~= 0.7 || opt.MinSitesPerQuartet ~= 1
    thresholdLabel = strrep(sprintf('%g', opt.SNRThreshold), '.', 'p');
    baseName = sprintf('%s_snr%s_minSites%d', baseName, ...
        thresholdLabel, opt.MinSitesPerQuartet);
end
if opt.V4QuartetRFOnly
    baseName = [baseName '_V4quartetRF_temp'];
end
if opt.ExcludeV4Site926
    baseName = [baseName '_minusSite926'];
end
OUT.resultFile = fullfile(outputDir, [baseName '_F.mat']);
OUT.figureFile = fullfile(outputDir, [baseName '_F.png']);
fig = [];
if opt.MakeFigure
    fig = makeFigure(OUT, regions);
end
if opt.SaveOutputs
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end
    save(OUT.resultFile, 'OUT', '-v7.3');
    if ~isempty(fig)
        print(fig, OUT.figureFile, '-dpng', '-r180');
    end
end
if ~isempty(fig)
    close(fig);
end
end

function D = decodeRegion(region, Tall, globalSites, R, m1, stimulus, ...
        trialQuartet, includeTrial, timeIdx, timeSamples, opt)
nSites = numel(globalSites);
nTrials = numel(stimulus);
localR = R;
localR.meanAct = R.meanAct(globalSites, :, :);
localR.meanSqAct = R.meanSqAct(globalSites, :, :);
SNR = compute_snr_per_color_sites(localR, Tall, (1:nSites)', ...
    'Verbose', false);
bestSNR = max([SNR.yellowEarly(:), SNR.yellowLate(:), ...
    SNR.purpleEarly(:), SNR.purpleLate(:)], [], 2, 'omitnan');
baseline = double(SNR.muSpont(:));
topResponse = max([SNR.muYellowEarly(:), SNR.muYellowLate(:), ...
    SNR.muPurpleEarly(:), SNR.muPurpleLate(:)], [], 2, 'omitnan');
responseScale = topResponse - baseline;
responseScale(~isfinite(responseScale) | responseScale <= 0) = NaN;

attOpts = struct('v1Sites', 1:nSites, 'timeIdx', timeIdx, ...
    'excludeOverlap', true, 'verbose', false);
attention = attention_modulation_V1_3bin(localR, Tall, SNR, attOpts);
dprime = double(attention.dprime(:));
bestQuartetSNR = nan(nSites,1);
passesVisual = bestSNR > opt.SNRThreshold;
if strcmp(region,'V4') && opt.V4QuartetRFOnly
    nTrialsPerStimulus = double(R.nTrials(:));
    quartetScores = nan(nSites,54,2);
    for q = 1:54
        members = quartetMembers(q);
        for w = 1:2
            values = double(localR.meanAct(:,members,w+1));
            weights = repmat(nTrialsPerStimulus(members).',nSites,1);
            validResponse = isfinite(values) & weights>0;
            weights(~validResponse) = 0;
            values(~validResponse) = 0;
            mu = sum(weights.*values,2)./sum(weights,2);
            quartetScores(:,q,w) = (mu-SNR.muSpont)./SNR.sdSpont;
        end
    end
    flat = reshape(quartetScores,nSites,[]);
    bestQuartetSNR = max(flat,[],2,'omitnan');
    bestQuartetSNR(all(~isfinite(flat),2)) = NaN;
    passesVisual = passesVisual | ...
        bestQuartetSNR > opt.SNRThreshold;
end
eligibleSite = passesVisual & ...
    attention.validSite(:) & isfinite(dprime) & ...
    isfinite(responseScale);
if strcmp(region,'V4') && opt.V4QuartetRFOnly
    assert(nnz(eligibleSite)==84, ...
        'The expected 84-site V4 capsule selection has changed.');
end
if strcmp(region,'V4') && opt.ExcludeV4Site926
    excluded = globalSites==926;
    assert(nnz(excluded)==1 && eligibleSite(excluded), ...
        'V4 site 926 is not uniquely selected.');
    eligibleSite(excluded) = false;
end

roleByStimulus = zeros(nSites, 216, 'int8');
for stim = 1:216
    T = Tall(stim).T;
    assignment = string(T.assignment);
    role = int8(assignment == "target") - int8(assignment == "distractor");
    role(logical(T.overlap)) = 0;
    roleByStimulus(:, stim) = role;
end
% The four quartet members must use the same RF sites with reversed roles.
nMismatchedSiteQuartets = 0;
for q = 1:54
    members = quartetMembers(q);
    r = double(roleByStimulus(:, members));
    valid = r(:, 1) ~= 0 & r(:, 1) == -r(:, 2) & ...
        r(:, 1) == r(:, 3) & r(:, 2) == r(:, 4);
    nMismatchedSiteQuartets = nMismatchedSiteQuartets + ...
        nnz(eligibleSite & any(r ~= 0, 2) & ~valid);
    r(~valid, :) = 0;
    roleByStimulus(:, members) = int8(r);
end
quartetSiteCount = zeros(54, 1);
for q = 1:54
    members = quartetMembers(q);
    quartetSiteCount(q) = nnz(eligibleSite & ...
        roleByStimulus(:, members(1)) ~= 0);
end

normalized = nan(nSites, nTrials);
for firstTrial = 1:opt.ChunkTrials:nTrials
    lastTrial = min(nTrials, firstTrial + opt.ChunkTrials - 1);
    trialRange = firstTrial:lastTrial;
    raw = double(m1.normMUA(globalSites, trialRange, timeSamples));
    response = reshape(mean(raw, 3, 'omitnan'), nSites, numel(trialRange));
    normalized(:, trialRange) = (response - baseline) ./ responseScale;
end
[centered, quartetMean, quartetTrialCount] = ...
    subtract_quartet_mean_response(normalized, trialQuartet, includeTrial, 54);

score = nan(nTrials, 1);
nSitesUsed = zeros(nTrials, 1);
for stim = 1:216
    rows = find(includeTrial & stimulus == stim);
    if isempty(rows)
        continue;
    end
    siteMask = eligibleSite & roleByStimulus(:, stim) ~= 0;
    weight = dprime(siteMask) .* double(roleByStimulus(siteMask, stim));
    response = centered(siteMask, rows);
    available = isfinite(response);
    response(~available) = 0;
    denominator = sum(abs(weight) .* available, 1);
    numerator = sum(weight .* response, 1);
    valid = denominator > 0 & sum(available, 1) >= opt.MinSitesPerQuartet;
    score(rows(valid)) = numerator(valid) ./ denominator(valid);
    nSitesUsed(rows) = sum(available, 1);
end
coverage = false(nTrials, 1);
coverage(includeTrial) = ...
    quartetSiteCount(trialQuartet(includeTrial)) >= opt.MinSitesPerQuartet;
scored = includeTrial & isfinite(score) & coverage;

D = struct();
D.region = region;
D.globalSites = globalSites;
D.baseline = baseline;
D.responseScale = responseScale;
D.bestSNR = bestSNR;
D.bestQuartetSNR = bestQuartetSNR;
D.eligibleSite = eligibleSite;
D.attentionDprime = dprime;
D.attentionPValue = attention.pValueTD(:);
D.roleByStimulus = roleByStimulus;
D.quartetSiteCount = quartetSiteCount;
D.quartetMeanResponse = quartetMean;
D.quartetTrialCount = quartetTrialCount;
D.nMismatchedSiteQuartets = nMismatchedSiteQuartets;
D.trialIndex = find(scored);
D.stimulus = stimulus(scored);
D.quartet = trialQuartet(scored);
D.score = score(scored);
D.nSitesUsed = nSitesUsed(scored);
D.nVisualSites = nnz(passesVisual);
D.nAddedEligibleSites = nnz(eligibleSite & ...
    ~(bestSNR > opt.SNRThreshold));
D.nEligibleSites = nnz(eligibleSite);
D.nScored = nnz(scored);
D.nUnscored = nnz(includeTrial & ~scored);
D.accuracyInSample = mean(score(scored) > 0);
fprintf(['  %s: visual %d, eligible %d, scored %d/%d, ' ...
    'median sites %.1f, accuracy %.1f%% (in-sample).\n'], region, ...
    D.nVisualSites, D.nEligibleSites, D.nScored, nnz(includeTrial), ...
    median(D.nSitesUsed), 100 * D.accuracyInSample);
end

function q = stimulusToQuartet(stimulus)
block = floor((stimulus - 1) / 8);
position = mod(stimulus - 1, 8) + 1;
q = 2 * block + 1 + ismember(position, [3 4 7 8]);
end

function members = quartetMembers(q)
block = floor((q - 1) / 2);
if mod(q, 2) == 1
    offsets = [1 2 5 6];
else
    offsets = [3 4 7 8];
end
members = 8 * block + offsets;
end

function path = firstExisting(candidates)
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        path = candidates{i};
        return;
    end
end
error('Missing input file. Checked: %s', strjoin(candidates, ', '));
end

function fig = makeFigure(OUT, regions)
fig = figure('Color', 'w', 'Position', [100 100 1350 420]);
colors = [0.16 0.43 0.66; 0.20 0.59 0.43; 0.70 0.33 0.25];
for i = 1:numel(regions)
    ax = subplot(1, 3, i, 'Parent', fig);
    D = OUT.(regions{i});
    if isempty(D.score)
        title(ax, sprintf('%s: no scored trials', regions{i}));
        continue;
    end
    histogram(ax, D.score, 'BinMethod', 'fd', ...
        'FaceColor', colors(i, :), 'EdgeColor', 'none');
    hold(ax, 'on');
    xline(ax, 0, '--k');
    xlabel(ax, 'Attention score');
    ylabel(ax, 'Trials');
    title(ax, sprintf('%s: %.1f%% (%d trials; %d sites)', ...
        regions{i}, 100 * D.accuracyInSample, ...
        D.nScored, D.nEligibleSites));
    grid(ax, 'on');
    box(ax, 'off');
end
sgtitle(fig, sprintf(['Figaro capsule attention decoders, 300-500 ms ' ...
    '(SNR > %.2f; in-sample)'], OUT.snrThreshold));
end
