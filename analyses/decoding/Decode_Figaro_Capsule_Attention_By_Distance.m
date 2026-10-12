function OUT = Decode_Figaro_Capsule_Attention_By_Distance(varargin)
%DECODE_FIGARO_CAPSULE_ATTENTION_BY_DISTANCE Apply fixed decoder weights in hextiles.
% RF distances are pooled over eligible site-quartets before binning.
% Hextile contributions retain the full decoder denominator and add up to
% the full decoder contrast; they are not divided by hextile site counts.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'core', 'rf'));
cfg = config();
p = inputParser;
p.addParameter('OutputDir', cfg.resultsDir, @(x) ischar(x) || isstring(x));
p.addParameter('BootstrapIterations', 300, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x == floor(x));
p.addParameter('BootstrapUnit', 'quartet', @(x) ischar(x) || isstring(x));
p.addParameter('RandomSeed', 270926, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('ReplotFrom', '', @(x) ischar(x) || isstring(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('V4QuartetRFOnly', false, ...
    @(x) islogical(x) && isscalar(x));
p.addParameter('ExcludeV4Site926', false, ...
    @(x) islogical(x) && isscalar(x));
p.addParameter('ElasticNet', false, ...
    @(x) islogical(x) && isscalar(x));
p.addParameter('ITRFRoleOnly', false, ...
    @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
if opt.ElasticNet && ismember('OutputDir',p.UsingDefaults)
    opt.OutputDir = fullfile(cfg.resultsRoot,'FigaroDecoderCrossover');
end
opt.BootstrapUnit = validatestring(char(opt.BootstrapUnit), ...
    {'quartet', 'site'});
assert(~opt.ExcludeV4Site926 || opt.V4QuartetRFOnly, ...
    'Site 926 exclusion requires the temporary V4 decoder.');
assert(~opt.ElasticNet || ...
    ~(opt.V4QuartetRFOnly || opt.ExcludeV4Site926), ...
    'The elastic-net comparison uses the original V1/V4 site sets.');
assert(~opt.ITRFRoleOnly || opt.ElasticNet, ...
    'RF-role-only IT comparison requires elastic-net V1/V4 weights.');
if ~isempty(opt.ReplotFrom)
    A = load(opt.ReplotFrom, 'OUT');
    OUT = A.OUT;
    if exist(opt.OutputDir, 'dir') ~= 7
        mkdir(opt.OutputDir);
    end
    savedUnit = 'quartet';
    if isfield(OUT.options, 'BootstrapUnit')
        savedUnit = OUT.options.BootstrapUnit;
    end
    variant = isfield(OUT,'v4QuartetRFOnly') && OUT.v4QuartetRFOnly;
    excludeSite = isfield(OUT,'excludedV4Site') && ...
        isequal(OUT.excludedV4Site,926);
    OUT.figureFile = fullfile(opt.OutputDir, ...
        [outputStem(savedUnit,variant,excludeSite,opt.ElasticNet, ...
        opt.ITRFRoleOnly) '.png']);
    makeFigure(OUT.result, OUT.timeMs, OUT.figureFile, ...
        savedUnit,variant,excludeSite,opt.ElasticNet,opt.ITRFRoleOnly);
    return;
end

decoderStem = 'Attention_decoder_Figaro_capsules_snr0p5_minSites1';
if opt.V4QuartetRFOnly
    decoderStem = [decoderStem '_V4quartetRF_temp'];
end
if opt.ExcludeV4Site926
    decoderStem = [decoderStem '_minusSite926'];
end
decoderFile = fullfile(cfg.resultsDir,[decoderStem '_F.mat']);
A = load(decoderFile, 'OUT');
decoder = A.OUT;
elasticFile = '';
rfItFile = '';
if opt.ElasticNet
    elasticFile = fullfile(cfg.resultsRoot,'FigaroDecoderCrossover', ...
        'figaro_capsules_elastic_net_best_quartet_IT_snr0p5.mat');
    A = load(elasticFile,'OUT');
    elastic = A.OUT;
    assert(elastic.snrThreshold==0.5 && ...
        elastic.nIncludedTrials==decoder.nIncludedTrials, ...
        'The elastic-net decoder is not aligned to this capsule session.');
end
if opt.ITRFRoleOnly
    rfItFile = fullfile(cfg.resultsDir, ...
        ['Attention_decoder_elastic_net_CV_IT_combined_' ...
        'nearestCurve_ITfunctional_F_rfRoleOnly_days123.mat']);
    A = load(rfItFile,'OUT');
    rfIt = A.OUT.IT;
    assert(numel(rfIt.trialIndex)==decoder.nIncludedTrials && ...
        isequal(double(rfIt.trialIndex(:)), ...
        double(elastic.IT.trialIndex(:))), ...
        'The RF-only IT decoder uses different Figaro trials.');
end
assert(decoder.snrThreshold == 0.5 && decoder.nStimuli == 216 && ...
    decoder.nQuartets == 54 && isequal(decoder.responseWindow, [300 500]) && ...
    (~opt.V4QuartetRFOnly || ...
    (decoder.v4QuartetRFOnly && ...
    decoder.V4.nEligibleSites==84-double(opt.ExcludeV4Site926))) && ...
    (~opt.ExcludeV4Site926 || isequal(decoder.excludedV4Site,926)), ...
    'Expected the fixed 0.5-SNR Figaro 300-500 ms decoder.');

m = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
trialFile = fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat');
M = load(trialFile, 'tb');
sampleTime = double(M.tb(:).');
assert(isequal(sampleTime, -200:499), ...
    'Expected Figaro trial samples at -200:499 ms.');
binTime = sampleTime(1:10:end) + 4.5;
baselineBins = binTime >= -150 & binTime < 0;
lateBins = binTime >= 300 & binTime < 500;

if opt.V4QuartetRFOnly
    names = {'V4'};
else
    names = {'V1', 'V4', 'IT'};
end
result = cell(numel(names), 1);
rows = cell(numel(names), 1);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng));
rng(opt.RandomSeed, 'twister');
for a = 1:numel(names)
    name = names{a};
    D = decoder.(name);
    geometry = load(D.geometryFile, ['Tall_' name]);
    Tall = geometry.(['Tall_' name]);
    if opt.ElasticNet
        E = elastic.(name);
        if opt.ITRFRoleOnly && strcmp(name,'IT')
            E = struct('sites',rfIt.siteIndexGlobal, ...
                'spontaneousSD',rfIt.spontaneousSD, ...
                'weights',rfIt.finalCoefficient, ...
                'signByQuartet',rfIt.rfRoleBySiteQuartet, ...
                'trialIndex',rfIt.trialIndex, ...
                'trialQuartet',rfIt.quartet, ...
                'trialLabel',rfIt.classLabel==1, ...
                'trialScore',rfIt.finalClassScore);
        end
        globalSite = double(E.sites(:));
        [known,selected] = ismember(globalSite,double(D.globalSites(:)));
        assert(all(known),'%s elastic sites lack capsule geometry.',name);
        responseScale = double(E.spontaneousSD(:));
        weight = double(E.weights(:));
        roleByStimulus = double(D.roleByStimulus(selected,:));
        trial = double(E.trialIndex(:));
        quartet = double(E.trialQuartet(:));
        label = 2*double(E.trialLabel(:))-1;
    else
        selected = find(D.eligibleSite(:));
        globalSite = double(D.globalSites(selected));
        responseScale = double(D.responseScale(selected));
        weight = double(D.attentionDprime(selected));
        roleByStimulus = double(D.roleByStimulus(selected, :));
        trial = double(D.trialIndex(:));
        stimulus = double(D.stimulus(:));
        quartet = double(D.quartet(:));
        label = ones(numel(trial), 1);
        label(mod(stimulus, 2) == 0) = -1;
    end
    eligible = unique(quartet);
    nQ = numel(eligible);
    nSites = numel(selected);
    nBins = numel(binTime);
    assert(numel(Tall) == 216 && all(isfinite(weight)) && ...
        all(isfinite(responseScale) & responseScale > 0));

    distance = nan(nSites, nQ);
    siteSign = zeros(nSites, nQ);
    fullDenominator = nan(nQ, 1);
    for k = 1:nQ
        q = eligible(k);
        members = quartetMembers(q);
        first = Tall(members(1)).T;
        x = double(first.along_GC(selected));
        role = roleByStimulus(:, members(1));
        active = role ~= 0;
        valid = active & isfinite(x) & ~logical(first.overlap(selected));
        assert(all(valid == active), ...
            '%s quartet %d has an active RF without a valid distance.', name, q);
        for j = 2:4
            other = Tall(members(j)).T;
            xOther = double(other.along_GC(selected));
            assert(all(isfinite(xOther(valid)) & ...
                ~logical(other.overlap(selected(valid))) & ...
                abs(xOther(valid) - x(valid)) < 1e-8), ...
                '%s quartet %d changes RF distance.', name, q);
        end
        distance(valid, k) = x(valid);
        if opt.ElasticNet
            siteSign(:,k) = double(E.signByQuartet(:,q));
            fullDenominator(k) = 1;
            assert(all(siteSign(active,k)==role(active)), ...
                '%s elastic RF signs differ from the saved geometry.',name);
        else
            siteSign(:, k) = role;
            fullDenominator(k) = sum(abs(weight(active)));
        end
        assert(fullDenominator(k) > 0, ...
            '%s quartet %d has no nonzero decoder weight.', name, q);
    end

    valid = isfinite(distance);
    nValid = nnz(valid);
    [~, order] = sortrows([distance(valid) (1:nValid).'], [1 2]);
    group = zeros(nValid, 1, 'uint8');
    group(order) = uint8(floor((0:nValid-1) * 6 / nValid) + 1);
    groupMap = zeros(nSites, nQ, 'uint8');
    groupMap(valid) = group;
    groupCounts = accumarray(double(group), 1, [6 1]);
    medianDistance = accumarray(double(group), distance(valid), ...
        [6 1], @median);
    assert(all(diff(medianDistance) > 0), ...
        '%s hextile distances are not ordered.', name);

    positiveRows = cell(nQ, 1);
    negativeRows = cell(nQ, 1);
    for k = 1:nQ
        positiveRows{k} = find(quartet == eligible(k) & label == 1);
        negativeRows{k} = find(quartet == eligible(k) & label == -1);
        assert(~isempty(positiveRows{k}) && ~isempty(negativeRows{k}), ...
            '%s quartet %d lacks a complementary trial set.', name, eligible(k));
    end

    contribution = zeros(nSites, nQ, nBins, 'single');
    for s = 1:nSites
        if abs(weight(s)) <= 1e-8
            continue;
        end
        raw = double(squeeze(m.normMUA(globalSite(s), :, :)));
        raw = raw(trial, :);
        assert(isequal(size(raw), [numel(trial) 700]));
        binned = squeeze(mean(reshape(raw, numel(trial), 10, nBins), ...
            2, 'omitnan')) / responseScale(s);
        for k = 1:nQ
            if siteSign(s, k) == 0
                continue;
            end
            delta = mean(binned(positiveRows{k}, :), 1, 'omitnan') - ...
                mean(binned(negativeRows{k}, :), 1, 'omitnan');
            contribution(s, k, :) = single( ...
                weight(s) * siteSign(s, k) / fullDenominator(k) * delta);
        end
    end

    savedContrast = nan(nQ, 1);
    reconstructedContrast = nan(nQ, 1);
    if opt.ElasticNet
        classScore = double(E.trialScore(:));
    else
        classScore = label .* double(D.score(:));
    end
    for k = 1:nQ
        savedContrast(k) = mean(classScore(positiveRows{k})) - ...
            mean(classScore(negativeRows{k}));
        x = squeeze(sum(double(contribution(:, k, lateBins)), 1));
        reconstructedContrast(k) = mean(x(:));
    end
    decoderCheck = max(abs(savedContrast - reconstructedContrast));
    assert(decoderCheck < 1e-4, ...
        '%s binned time courses do not reproduce the saved decoder.', name);
    fprintf('%s full-decoder contrast check: %.3g\n', name, decoderCheck);

    qGroup = nan(nQ, 6, nBins);
    weightedSites = zeros(nQ, 6);
    unassigned = zeros(nQ,nBins);
    unassignedWeightedSites = zeros(nQ,1);
    for k = 1:nQ
        offObject = groupMap(:,k)==0 & abs(weight)>1e-8 & ...
            siteSign(:,k)~=0;
        unassignedWeightedSites(k) = nnz(offObject);
        unassigned(k,:) = reshape(sum(double( ...
            contribution(offObject,k,:)),1),1,nBins);
        for g = 1:6
            keep = groupMap(:, k) == g;
            weightedSites(k, g) = nnz(keep & abs(weight) > 1e-8);
            if any(keep)
                qGroup(k, g, :) = sum(double(contribution(keep, k, :)), 1);
            end
        end
    end
    if opt.ElasticNet
        assigned = reshape(sum(qGroup,2,'omitnan'),nQ,nBins);
        full = reshape(sum(double(contribution),1),nQ,nBins);
        assert(max(abs(assigned(:)+unassigned(:)-full(:)))<1e-5, ...
            '%s RF groups and off-object sites do not reconstruct the decoder.',name);
        assignedLateContrast = mean(assigned(:,lateBins),'all','omitnan');
        unassignedLateContrast = mean(unassigned(:,lateBins),'all','omitnan');
        fprintf('%s late contrast: RF groups %.3f, no-RF %.3f\n', ...
            name,assignedLateContrast,unassignedLateContrast);
    else
        assignedLateContrast = NaN;
        unassignedLateContrast = NaN;
    end
    curve = reshape(mean(qGroup, 1, 'omitnan'), 6, nBins);
    curve = curve - mean(curve(:, baselineBins), 2);
    fit = fitGroups(curve, binTime, medianDistance);
    assert(fit.exitflag > 0 && all(isfinite(fit.summary.t50)), ...
        '%s sigmoid fit did not yield six finite t50 values.', name);
    line = fitSlope(medianDistance, fit.summary.t50);

    bootT50 = nan(opt.BootstrapIterations, 6);
    bootSlope = nan(opt.BootstrapIterations, 1);
    for b = 1:opt.BootstrapIterations
        if strcmp(opt.BootstrapUnit, 'quartet')
            draw = randi(nQ, nQ, 1);
            bootGroup = qGroup(draw, :, :);
        else
            siteDraw = randi(nSites, nSites, 1);
            multiplicity = accumarray(siteDraw, 1, [nSites 1]);
            bootGroup = nan(nQ, 6, nBins);
            for k = 1:nQ
                for g = 1:6
                    keep = groupMap(:, k) == g & multiplicity > 0;
                    if any(keep)
                        weightedContribution = double(contribution(keep, k, :)) .* ...
                            reshape(multiplicity(keep), [], 1, 1);
                        bootGroup(k, g, :) = sum(weightedContribution, 1);
                    end
                end
            end
        end
        bootCurve = reshape(mean(bootGroup, 1, 'omitnan'), 6, nBins);
        bootCurve = bootCurve - mean(bootCurve(:, baselineBins), 2);
        try
            bootFit = fitGroups(bootCurve, binTime, medianDistance);
            if bootFit.exitflag > 0 && all(isfinite(bootFit.summary.t50))
                bootT50(b, :) = bootFit.summary.t50(:).';
                bootLine = fitSlope(medianDistance, bootFit.summary.t50);
                bootSlope(b) = bootLine.slope;
            end
        catch
            % Failed resamples have no timing estimate.
        end
    end
    t50CI = nan(6, 2);
    for g = 1:6
        values = bootT50(isfinite(bootT50(:, g)), g);
        if ~isempty(values)
            t50CI(g, :) = prctile(values, [2.5 97.5]);
        end
    end
    validSlope = bootSlope(isfinite(bootSlope));
    slopeCI = [NaN NaN];
    if ~isempty(validSlope)
        slopeCI = prctile(validSlope, [2.5 97.5]);
    end
    result{a} = struct('area', string(name), ...
        'eligibleQuartets', eligible, 'medianDistance', medianDistance, ...
        'groupSiteQuartetCounts', groupCounts, ...
        'groupBySiteQuartet', groupMap, ...
        'selectedSiteLocal', selected, 'globalSite', globalSite, ...
        'fullDenominator', fullDenominator, ...
        'weightedSitesByQuartet', weightedSites, ...
        'quartetDifference', qGroup, 'curve', curve, 'fit', fit, ...
        'unassignedQuartetDifference',unassigned, ...
        'unassignedWeightedSitesByQuartet',unassignedWeightedSites, ...
        'assignedLateContrast',assignedLateContrast, ...
        'unassignedLateContrast',unassignedLateContrast, ...
        'line', line, 't50CI', t50CI, 'slopeCI', slopeCI, ...
        'bootstrapSlope', bootSlope, ...
        'decoderReconstructionError', decoderCheck);
    tauAtLowerBound = fit.sharedTau <= fit.options.minTau + 1e-6;
    rows{a} = table(string(name), nQ, nSites, nValid, ...
        line.slope, slopeCI(1), slopeCI(2), line.rSquared, ...
        fit.sharedTau, tauAtLowerBound, nnz(isfinite(bootSlope)), ...
        'VariableNames', {'area', 'nQuartets', 'nDecoderSites', ...
        'nEligibleSiteQuartets', 'slopeMsPerWidth', 'slopeCILow', ...
        'slopeCIHigh', 'timingRegressionRSquared', 'sharedTauMs', ...
        'tauAtLowerBound', 'nValidBootstrap'});
    fprintf(['%s: %d quartets, %d site-quartets, median %d ' ...
        'unassigned weighted sites; slope %.2f [%.2f %.2f] ms/width\n'], ...
        name, nQ, nValid, median(unassignedWeightedSites), ...
        line.slope, slopeCI(1), slopeCI(2));
end

OUT = struct();
OUT.result = result;
OUT.summary = vertcat(rows{:});
OUT.timeMs = binTime;
OUT.options = opt;
OUT.sourceDecoderFile = decoderFile;
if opt.ElasticNet
    OUT.sourceDecoderFile = elasticFile;
end
if opt.ITRFRoleOnly
    OUT.sourceITDecoderFile = rfItFile;
end
OUT.v4QuartetRFOnly = opt.V4QuartetRFOnly;
OUT.excludedV4Site = [];
if isfield(decoder,'excludedV4Site')
    OUT.excludedV4Site = decoder.excludedV4Site;
end
OUT.description = ['Fixed SNR>0.5 Figaro attention-decoder weights in six ' ...
    'pooled site-quartet RF distance hextiles. Distance is Tall.along_GC, ' ...
    'the cue-to-RF projection divided by capsule width, not an exact ' ...
    'growth-cone path. Non-overlap RFs only. Site contributions retain ' ...
    'the full quartet decoder denominator; hextiles are not renormalized. ' ...
    'Quartets are averaged equally. Sigmoids fit 80-450 ms; observed ' ...
    'curves receive five-bin smoothing for display only. ' ...
    opt.BootstrapUnit ' bootstrap keeps full-data decoder weights fixed; ' ...
    'site bootstrap carries all quartets of each sampled site together.'];
if opt.ExcludeV4Site926
    OUT.description = [OUT.description ...
        ' This sensitivity refits the separate V4 decoder with site 926 ' ...
        'excluded; individual trials are not removed.'];
end
if opt.ElasticNet
    OUT.description = ['Fixed Figaro capsule elastic-net weights in the ' ...
        'same six pooled RF-site/quartet distance hextiles and timing ' ...
        'fit as the d-prime decoder. Distance is Tall.along_GC, a ' ...
        'cue-to-RF projection divided by capsule width. Training ' ...
        'uses the full data and IT best-quartet SNR > 0.5. IT sites ' ...
        'without a localized RF contribute to the full decoder but ' ...
        'not to a distance group; the unassigned contribution is ' ...
        'saved separately. Class contrasts cancel the elastic-net ' ...
        'intercept and quartet centering. Site bootstrap keeps the ' ...
        'fixed full-data weights and all quartets of each sampled ' ...
        'site together.'];
end
if opt.ITRFRoleOnly
    OUT.description = ['Fixed Figaro capsule elastic-net weights in six ' ...
        'pooled RF site/quartet distance hextiles. IT was retrained ' ...
        'using only binary target/distractor RF roles, with background ' ...
        'and overlap set to zero; 0.5 functional inclusion was used. ' ...
        'V1/V4 retain the existing elastic-net decoders. The ' ...
        opt.BootstrapUnit ' bootstrap keeps full-data weights fixed.'];
end
stem = outputStem(opt.BootstrapUnit,opt.V4QuartetRFOnly, ...
    opt.ExcludeV4Site926,opt.ElasticNet,opt.ITRFRoleOnly);
OUT.resultFile = fullfile(opt.OutputDir, [stem '.mat']);
OUT.summaryFile = fullfile(opt.OutputDir, [stem '.csv']);
OUT.figureFile = fullfile(opt.OutputDir, [stem '.png']);
if opt.SaveOutputs || opt.MakeFigure
    if exist(opt.OutputDir, 'dir') ~= 7
        mkdir(opt.OutputDir);
    end
end
if opt.ElasticNet && (opt.SaveOutputs || opt.MakeFigure)
    assert(~isfile(OUT.resultFile) && ~isfile(OUT.summaryFile) && ...
        ~isfile(OUT.figureFile), ...
        'Elastic-net distance outputs already exist; refusing overwrite.');
end
if opt.MakeFigure
    makeFigure(result, binTime, OUT.figureFile, ...
        opt.BootstrapUnit,opt.V4QuartetRFOnly, ...
        opt.ExcludeV4Site926,opt.ElasticNet,opt.ITRFRoleOnly);
end
if opt.SaveOutputs
    writetable(OUT.summary, OUT.summaryFile);
    save(OUT.resultFile, 'OUT', '-v7.3');
end
disp(OUT.summary);
end

function stem = outputStem(bootstrapUnit,v4QuartetRFOnly,excludeSite926,elasticNet,itRfRoleOnly)
stem = 'Figaro_capsule_attention_decoder_six_distance_bins_snr0p5_F';
if strcmp(bootstrapUnit, 'site')
    stem = 'Figaro_capsule_attention_decoder_six_distance_bins_siteBootstrap_snr0p5_F';
end
if v4QuartetRFOnly
    stem = [stem '_V4quartetRF_temp'];
end
if excludeSite926
    stem = [stem '_minusSite926'];
end
if elasticNet
    stem = [stem '_elasticNet_ITbestQuartet'];
end
if itRfRoleOnly
    stem = [stem '_IT_rfRoleOnly'];
end
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

function fit = fitGroups(Y, t, x)
G = struct('groupMeanSigned', Y, ...
    'timeWindows', [t(:) t(:)], ...
    'groupSummary', table((1:6).', x(:), x(:), ...
    'VariableNames', {'groupIdx', 'alongMin', 'alongMax'}));
fit = fit_grouped_alonggc_sigmoid_sharedslope(G, ...
    'fitStartMs', 80, 'fitEndMs', 450, ...
    'smoothW', 1, 'plotMode', 'none', ...
    'useT50SEWeightedRegression', false, 'verbose', false);
end

function line = fitSlope(x, y)
x = double(x(:));
y = double(y(:));
X = [ones(numel(x), 1) x];
beta = X \ y;
residual = y - X * beta;
total = sum((y - mean(y)).^2);
line = struct('intercept', beta(1), 'slope', beta(2), ...
    'rSquared', 1 - sum(residual.^2) / max(total, eps));
end

function makeFigure(result, t, figureFile, bootstrapUnit, ...
        v4QuartetRFOnly, excludeSite926,elasticNet,itRfRoleOnly)
h = figure('Color', 'w', 'Visible', 'off', ...
    'Position', [50 50 1900 460*numel(result)]);
layout = tiledlayout(h, numel(result), 2, 'Padding', 'compact', ...
    'TileSpacing', 'compact');
colors = [0.07 0.35 0.82; ...
          0.00 0.60 0.73; ...
          0.04 0.57 0.31; ...
          0.48 0.62 0.07; ...
          0.67 0.48 0.02; ...
          0.79 0.12 0.10];
visibleTime = t >= 0 & t <= 500;
low = -0.12;
high = 1.2;
for a = 1:numel(result)
    Z = result{a};
    fitMax = max(Z.fit.yFit(:, visibleTime), [], 2);
    observed = movmean(Z.curve ./ fitMax, 5, 2, 'omitnan');
    values = observed(:, visibleTime);
    low = min(low, min(values(:)));
    high = max(high, max(values(:)));
end
sharedYLim = [floor(10 * low) / 10, ceil(10 * high) / 10];
for a = 1:numel(result)
    Z = result{a};
    ax = nexttile(layout, 2*a-1);
    hold(ax, 'on');
    fitMax = max(Z.fit.yFit(:, visibleTime), [], 2);
    for g = 1:6
        observed = movmean(Z.curve(g, :) / fitMax(g), ...
            5, 2, 'omitnan');
        plot(ax, t, observed, 'Color', ...
            0.78 * colors(g, :) + 0.22 * [1 1 1], ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax, t, Z.fit.yFit(g, :) / fitMax(g), ...
            'Color', colors(g, :), 'LineWidth', 2, ...
            'DisplayName', sprintf('Group %d', g));
    end
    xline(ax, 0, ':', 'HandleVisibility', 'off');
    xlim(ax, [0 500]);
    ylim(ax, sharedYLim);
    grid(ax, 'on');
    box(ax, 'off');
    ylabel(ax, 'Fit normalized to max = 1');
    weightedQuartets = sum(Z.weightedSitesByQuartet > 0, 1);
    title(ax, sprintf('%s | %d quartets | weighted %s', ...
        Z.area, numel(Z.eligibleQuartets), mat2str(weightedQuartets)), ...
        'FontWeight', 'bold');
    if a == 1
        legend(ax, 'Location', 'southeast');
    end
    ax = nexttile(layout, 2*a);
    hold(ax, 'on');
    for g = 1:6
        x = Z.medianDistance(g);
        y = Z.fit.summary.t50(g);
        ci = Z.t50CI(g, :);
        if all(isfinite(ci))
            errorbar(ax, x, y, max(y-ci(1), 0), max(ci(2)-y, 0), ...
                'o', 'Color', colors(g, :), ...
                'MarkerFaceColor', colors(g, :), ...
                'MarkerSize', 7, 'CapSize', 0);
        else
            plot(ax, x, y, 'o', 'Color', colors(g, :), ...
                'MarkerFaceColor', colors(g, :), 'MarkerSize', 7);
        end
    end
    xx = linspace(min(Z.medianDistance), max(Z.medianDistance), 100);
    plot(ax, xx, Z.line.intercept + Z.line.slope * xx, ...
        'k--', 'LineWidth', 1.5);
    grid(ax, 'on');
    box(ax, 'off');
    ylabel(ax, 'Fitted t_{50} (ms)');
    if Z.fit.sharedTau <= Z.fit.options.minTau + 1e-6
        tauNote = sprintf('shared tau %.1f ms (lower bound)', Z.fit.sharedTau);
    else
        tauNote = sprintf('shared tau %.1f ms', Z.fit.sharedTau);
    end
    title(ax, sprintf('%s | slope %+.2f [%.2f, %.2f] ms/width\n%s', ...
        Z.area, Z.line.slope, Z.slopeCI(1), Z.slopeCI(2), tauNote), ...
        'FontWeight', 'bold');
    if a == numel(result)
        xlabel(ax, 'Cue-to-RF along-capsule distance / width');
    end
end
if v4QuartetRFOnly
    variantLabel = 'temporary V4 best-color OR best-quartet RF decoder';
else
    variantLabel = 'fixed attention-decoder weights';
end
if excludeSite926
    variantLabel = [variantLabel ' excluding site 926'];
end
if elasticNet
    variantLabel = 'fixed elastic-net weights | IT best-quartet SNR';
end
if itRfRoleOnly
    variantLabel = ['fixed elastic-net weights | IT functional ' ...
        'RF-role-only; V1/V4 SNR > 0.5'];
end
if itRfRoleOnly
    inclusionLabel = '';
else
    inclusionLabel = 'SNR > 0.5 | ';
end
title(layout, ['Figaro capsules | ' inclusionLabel variantLabel ...
    ' | six pooled distance groups | overlap excluded | ' ...
    bootstrapUnit ' bootstrap'], ...
    'FontWeight', 'bold');
print(h, figureFile, '-dpng', '-r180');
close(h);
end
