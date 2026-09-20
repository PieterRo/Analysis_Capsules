function SNR = compute_snr_per_color_region(R, Tall, globalSites)
% COMPUTE_SNR_PER_COLOR_REGION
% Compute the established per-color SNR normalization for an arbitrary
% contiguous or non-contiguous set of channels. Rows in Tall(stim).T must
% correspond one-to-one with globalSites.

WIN_SPONT = 1; % -200..0 ms
WIN_EARLY = 2; % 40..240 ms
WIN_LATE = 3;  % 300..500 ms

COL_YELLOW = "yellowArm";
COL_PURPLE = "purple";
minTrialsPerStim = 1;
minTotalTrialsPerColor = 1;
useBesselCorrection = true;

globalSites = double(globalSites(:));
nSites = numel(globalSites);
nTrials = double(R.nTrials(:));
nStim = numel(nTrials);

[nSitesTotal, nStimData, nWin] = size(R.meanAct);
assert(nStimData == nStim, ...
    'R.meanAct has %d stimuli but R.nTrials has %d.', nStimData, nStim);
assert(all(size(R.meanSqAct) == size(R.meanAct)), ...
    'R.meanSqAct must match R.meanAct.');
assert(nWin >= WIN_LATE, 'Expected at least three response windows.');
assert(max(globalSites) <= nSitesTotal, ...
    'Requested channel %d exceeds the %d channels in R.', ...
    max(globalSites), nSitesTotal);
assert(numel(Tall) == nStim, ...
    'Tall has %d stimuli but R has %d.', numel(Tall), nStim);

if isfield(Tall, 'stimNum')
    tallStimNums = arrayfun(@(x) x.stimNum, Tall(:));
    [sortedStimNums, order] = sort(tallStimNums(:));
    assert(all(sortedStimNums(:).' == 1:nStim), ...
        'Tall.stimNum must cover 1:%d exactly.', nStim);
    Tall = Tall(order);
end

T0 = Tall(1).T;
assert(istable(T0) && height(T0) == nSites, ...
    'Tall(stim).T must have %d rows matching globalSites.', nSites);
varNames = string(T0.Properties.VariableNames);
colorColumn = find(varNames == "center_color", 1);
if isempty(colorColumn)
    colorColumn = find(contains(lower(varNames), "center") & ...
        contains(lower(varNames), "color"), 1);
end
assert(~isempty(colorColumn), 'Tall tables must contain center_color.');

stimOk = nTrials >= minTrialsPerStim;
muSpont = nan(nSites,1);
sdSpont = nan(nSites,1);
nSpontTrials = nan(nSites,1);
isYellow = false(nSites,nStim);
isPurple = false(nSites,nStim);

for stimIdx = 1:nStim
    T = Tall(stimIdx).T;
    assert(istable(T) && height(T) == nSites, ...
        'Tall(%d).T must have %d rows.', stimIdx, nSites);
    colorLabels = strtrim(string(T{:,colorColumn}));
    isYellow(:,stimIdx) = colorLabels == COL_YELLOW;
    isPurple(:,stimIdx) = colorLabels == COL_PURPLE;
end

muYellowEarly = nan(nSites,1);
muYellowLate = nan(nSites,1);
muPurpleEarly = nan(nSites,1);
muPurpleLate = nan(nSites,1);
nYellowTrials = zeros(nSites,1);
nPurpleTrials = zeros(nSites,1);

for localSite = 1:nSites
    globalSite = globalSites(localSite);

    muBase = squeeze(R.meanAct(globalSite,:,WIN_SPONT)).';
    meanSqBase = squeeze(R.meanSqAct(globalSite,:,WIN_SPONT)).';
    validBase = stimOk & isfinite(muBase) & isfinite(meanSqBase);
    nBase = sum(nTrials(validBase));
    nSpontTrials(localSite) = nBase;
    if nBase > 1
        mu = sum(nTrials(validBase).*muBase(validBase))/nBase;
        ex2 = sum(nTrials(validBase).*meanSqBase(validBase))/nBase;
        varPop = max(0, ex2 - mu^2);
        if useBesselCorrection
            varPop = varPop*(nBase/(nBase-1));
        end
        muSpont(localSite) = mu;
        sdSpont(localSite) = sqrt(varPop);
    end

    responseEarly = squeeze(R.meanAct(globalSite,:,WIN_EARLY)).';
    responseLate = squeeze(R.meanAct(globalSite,:,WIN_LATE)).';

    useYellow = stimOk & isYellow(localSite,:).';
    nYellow = sum(nTrials(useYellow));
    nYellowTrials(localSite) = nYellow;
    if nYellow >= minTotalTrialsPerColor
        muYellowEarly(localSite) = ...
            sum(nTrials(useYellow).*responseEarly(useYellow))/nYellow;
        muYellowLate(localSite) = ...
            sum(nTrials(useYellow).*responseLate(useYellow))/nYellow;
    end

    usePurple = stimOk & isPurple(localSite,:).';
    nPurple = sum(nTrials(usePurple));
    nPurpleTrials(localSite) = nPurple;
    if nPurple >= minTotalTrialsPerColor
        muPurpleEarly(localSite) = ...
            sum(nTrials(usePurple).*responseEarly(usePurple))/nPurple;
        muPurpleLate(localSite) = ...
            sum(nTrials(usePurple).*responseLate(usePurple))/nPurple;
    end
end

SNR = struct();
SNR.globalSites = globalSites;
SNR.muSpont = muSpont;
SNR.sdSpont = sdSpont;
SNR.nSpontTrials = nSpontTrials;
SNR.muYellowEarly = muYellowEarly;
SNR.muYellowLate = muYellowLate;
SNR.muPurpleEarly = muPurpleEarly;
SNR.muPurpleLate = muPurpleLate;
SNR.nYellowTrials = nYellowTrials;
SNR.nPurpleTrials = nPurpleTrials;

SNR.yellowEarly = (muYellowEarly - muSpont)./sdSpont;
SNR.yellowLate = (muYellowLate - muSpont)./sdSpont;
SNR.purpleEarly = (muPurpleEarly - muSpont)./sdSpont;
SNR.purpleLate = (muPurpleLate - muSpont)./sdSpont;

badNoise = ~isfinite(sdSpont) | sdSpont <= 0;
colorFields = {'yellowEarly','yellowLate','purpleEarly','purpleLate'};
for i = 1:numel(colorFields)
    SNR.(colorFields{i})(badNoise) = NaN;
end
SNR.yellowEarly(nYellowTrials < minTotalTrialsPerColor) = NaN;
SNR.yellowLate(nYellowTrials < minTotalTrialsPerColor) = NaN;
SNR.purpleEarly(nPurpleTrials < minTotalTrialsPerColor) = NaN;
SNR.purpleLate(nPurpleTrials < minTotalTrialsPerColor) = NaN;

end
