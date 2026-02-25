function deltaQ = compute_deltaQ_quartet(tb, Rdata, Tall_V1, SNRn, siteRange, quartetMembers, excludeOverlap)
% COMPUTE_DELTAQ_QUARTET
% Compute per-site, per-quartet T-D deltas for a given time bin.

nSites = numel(siteRange);
nQuartets = size(quartetMembers, 1);

bAll = SNRn.muSpont(:);
topMat = [SNRn.muYellowEarly(:), SNRn.muYellowLate(:), SNRn.muPurpleEarly(:), SNRn.muPurpleLate(:)];
muTopAll = max(topMat, [], 2);
b = bAll(siteRange);
scale = muTopAll(siteRange) - b;
scale(~isfinite(scale) | scale <= 0) = NaN;

nTrials = Rdata.nTrials;
if isvector(nTrials)
    nTrials = nTrials(:)';  % 1x384
    perSiteTrials = false;
else
    perSiteTrials = true;
end

requiredVarsQ = ["assignment","overlap"];
deltaQ = nan(nSites, nQuartets);
for qIdx = 1:nQuartets
    stimsQ = quartetMembers(qIdx,:);
    sumT = zeros(nSites,1); NT = zeros(nSites,1);
    sumD = zeros(nSites,1); ND = zeros(nSites,1);

    for stimNum = stimsQ
        if ~isfield(Tall_V1(stimNum),'T') || ~istable(Tall_V1(stimNum).T)
            continue;
        end
        Tq = Tall_V1(stimNum).T;
        if height(Tq) < max(siteRange) || any(~ismember(requiredVarsQ, string(Tq.Properties.VariableNames)))
            continue;
        end

        asg = string(Tq.assignment(siteRange));
        isT = (asg == "target");
        isD = (asg == "distractor");
        isBG = (asg == "background");
        isOV = Tq.overlap(siteRange) ~= 0;

        EX = squeeze(double(Rdata.meanAct(siteRange, stimNum, tb)));
        EY = (EX - b) ./ scale;

        if ~perSiteTrials
            wAll = repmat(double(nTrials(stimNum)), nSites, 1);
        else
            wAll = double(nTrials(siteRange, stimNum));
        end

        if excludeOverlap
            good = isfinite(EY) & isfinite(wAll) & (wAll > 0) & ~isBG & ~isOV;
        else
            good = isfinite(EY) & isfinite(wAll) & (wAll > 0) & ~isBG;
        end

        mT = good & isT;
        if any(mT)
            w = wAll(mT);
            sumT(mT) = sumT(mT) + w .* EY(mT);
            NT(mT) = NT(mT) + w;
        end

        mD = good & isD;
        if any(mD)
            w = wAll(mD);
            sumD(mD) = sumD(mD) + w .* EY(mD);
            ND(mD) = ND(mD) + w;
        end
    end

    hasBoth = (NT > 0) & (ND > 0);
    muTq = nan(nSites,1); muDq = nan(nSites,1);
    muTq(hasBoth) = sumT(hasBoth) ./ NT(hasBoth);
    muDq(hasBoth) = sumD(hasBoth) ./ ND(hasBoth);
    deltaQ(:, qIdx) = muTq - muDq;
end
end
