function h = plot_projected_attentiondiff_on_example_stim(stimID_example, OUT, Tall_V1, ALLCOORDS, optsPlot)
% PLOT_PROJECTED_ATTENTIONDIFF_ON_EXAMPLE_STIM
% Single still frame: project attention effect onto target/distractor curve regions.
%
% Required significance rule:
%   significant sites are selected by OUT.pValueTD < 0.05 (default threshold).
%
% Required routing rule:
%   delta = T - D per site
%   - on target-assigned RFs:     value = max(0,  delta)
%   - on distractor-assigned RFs: value = max(0, -delta)
%
% Inputs:
%   stimID_example : stimulus ID in 1..384
%   OUT            : output struct from attention_modulation_V1_3bin
%   Tall_V1        : struct array, Tall_V1(stim).T table contains RF assignment/overlap/GC coords
%   ALLCOORDS      : geometry struct
%   optsPlot       : struct with plotting/config options:
%       .RTAB384   (required)
%       .pThresh   (default 0.05)
%       .siteRange (default 1:512)
%       .markerSize (default 18)
%       .alpha      (default 0.75)   % constant alpha
%       .robustPct  (default 95)     % color max percentile
%
% Safety clause:
%   If OUT lacks clear per-site TD magnitude fields (muT/muD), prints fieldnames/sizes and errors.

if nargin < 5 || isempty(optsPlot), optsPlot = struct(); end
if ~isfield(optsPlot,'RTAB384')
    error('optsPlot.RTAB384 is required.');
end
if ~isfield(optsPlot,'pThresh'), optsPlot.pThresh = 0.05; end
if ~isfield(optsPlot,'siteRange'), optsPlot.siteRange = 1:512; end
if ~isfield(optsPlot,'markerSize'), optsPlot.markerSize = 18; end
if ~isfield(optsPlot,'alpha'), optsPlot.alpha = 0.75; end
if ~isfield(optsPlot,'robustPct'), optsPlot.robustPct = 95; end
if ~isfield(optsPlot,'stimIdx'), optsPlot.stimIdx = []; end
if ~isfield(optsPlot,'cMaxFixed'), optsPlot.cMaxFixed = []; end
if ~isfield(optsPlot,'bgColor'), optsPlot.bgColor = [0.5 0.5 0.5]; end
if ~isfield(optsPlot,'cLow'), optsPlot.cLow = [0.62 0.62 0.62]; end
if ~isfield(optsPlot,'cHigh'), optsPlot.cHigh = [0.85 0.05 0.05]; end
if ~isfield(optsPlot,'projectionMode'), optsPlot.projectionMode = 'quartet_pooled'; end

RTAB384 = optsPlot.RTAB384;
siteRange = optsPlot.siteRange(:)';

if stimID_example < 1 || stimID_example > numel(Tall_V1)
    error('stimID_example out of range.');
end

% ---- Safety check for effect magnitude fields ----
if ~(isfield(OUT,'muT') && isfield(OUT,'muD'))
    fprintf('OUT is missing muT/muD. Available fields:\n');
    fn = fieldnames(OUT);
    disp(fn);
    fprintf('Candidate field sizes:\n');
    for k = 1:numel(fn)
        v = OUT.(fn{k});
        if isnumeric(v) || islogical(v)
            sz = size(v);
            fprintf('  %s: [%s]\n', fn{k}, num2str(sz));
        end
    end
    error(['Cannot compute TD magnitude safely for plotting. ' ...
           'Expected per-site fields OUT.muT and OUT.muD for window 3.']);
end

if ~isfield(OUT,'pValueTD')
    error('OUT.pValueTD is required for significance masking.');
end

% Optional sanity check on expected time window
if isfield(OUT,'timeWindowUsed') && numel(OUT.timeWindowUsed)==2
    tw = OUT.timeWindowUsed(:)';
    if ~(abs(tw(1)-300) < 1e-9 && abs(tw(2)-500) < 1e-9)
        fprintf('Warning: OUT.timeWindowUsed is [%g %g], not [300 500].\n', tw(1), tw(2));
    end
end

% ---- Site-level masks ----
sigMaskAll = isfinite(OUT.pValueTD(:)) & (OUT.pValueTD(:) < optsPlot.pThresh);
nSigSites = sum(sigMaskAll);
fprintf('nSigSites = %d\n', nSigSites);

% ---- Per-site effect ----
deltaAll = OUT.muT(:) - OUT.muD(:);  % T - D

% ---- Canonical geometry (same conventions as activity projection) ----
W = 1024; H = 768;
toPx = @(q) [q(1) + W/2, H/2 - q(2)];

fieldName = sprintf('stim_%d', stimID_example);
s    = double(ALLCOORDS.(fieldName).s(:))';
tFig = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack= double(ALLCOORDS.(fieldName).t_back(:))';

s_px = toPx(s);
tT   = toPx(tFig);
tD   = toPx(tBack);

uT = (tT - s_px) / norm(tT - s_px);
uD = (tD - s_px) / norm(tD - s_px);
nT = perpTowardOther(s_px, tT, tD);
nD = perpTowardOther(s_px, tD, tT);

widthEx = double(RTAB384(stimID_example,7));

% ---- Project coordinates ----
X = [];
Y = [];
V = [];
nTargetContrib = 0;
nDistrContrib  = 0;
valTargetAll = [];
valDistrAll  = [];

siteGlobal = siteRange(:);
sigMask = sigMaskAll(siteGlobal);
delta = deltaAll(siteGlobal);

requiredVars = ["assignment","overlap","along_GC","perp_signed_GC"];
if ~isfield(Tall_V1(stimID_example),'T') || ~istable(Tall_V1(stimID_example).T)
    error('Tall_V1(%d).T missing/invalid.', stimID_example);
end
TtabEx = Tall_V1(stimID_example).T;
if height(TtabEx) < max(siteRange) || any(~ismember(requiredVars, string(TtabEx.Properties.VariableNames)))
    error('Tall_V1(%d).T lacks required rows/vars.', stimID_example);
end

mode = lower(string(optsPlot.projectionMode));
if mode == "quartet_pooled"
    % Quartet-conditioned pooling:
    % - Hard-coded quartets per block of 8 stimuli:
    %   [1 2 5 6] and [3 4 7 8] (+ block offset)
    % - For each site and quartet, compute pooled T and D from selected time bin
    %   using Tall_V1 assignment/overlap and normalized response from Rdata/SNRnorm.
    % - Use sign of quartet delta to route each site-stim contribution to target/distractor.
    if ~isfield(optsPlot,'Rdata') || ~isfield(optsPlot,'SNRnorm') || ~isfield(optsPlot,'timeIdx')
        error(['projectionMode=quartet_pooled requires optsPlot.Rdata, ' ...
               'optsPlot.SNRnorm, and optsPlot.timeIdx.']);
    end
    Rdata = optsPlot.Rdata;
    SNRn = optsPlot.SNRnorm;
    tb = optsPlot.timeIdx;
    assert(isfield(Rdata,'meanAct') && isfield(Rdata,'nTrials'), ...
        'optsPlot.Rdata must contain meanAct and nTrials.');
    assert(size(Rdata.meanAct,2) == 384 && size(Rdata.meanAct,3) >= tb, ...
        'Rdata.meanAct must be [nSites x 384 x nTime] with nTime >= optsPlot.timeIdx.');

    % Build normalization constants
    bAll = SNRn.muSpont(:);
    topMat = [SNRn.muYellowEarly(:), SNRn.muYellowLate(:), SNRn.muPurpleEarly(:), SNRn.muPurpleLate(:)];
    muTopAll = max(topMat, [], 2);
    b = bAll(siteRange);
    scale = muTopAll(siteRange) - b;
    scale(~isfinite(scale) | scale <= 0) = NaN;

    % Trials shape
    nTrials = Rdata.nTrials;
    if isvector(nTrials)
        nTrials = nTrials(:)';  % 1x384
        perSiteTrials = false;
    else
        perSiteTrials = true;
    end

    % Quartet map
    nStim = 384;
    nBlocks = nStim / 8;
    nQuartets = nBlocks * 2;
    stimToQuartet = zeros(1, nStim);
    quartetMembers = zeros(nQuartets, 4);
    q = 0;
    for bIdx = 0:(nBlocks-1)
        base = 8*bIdx;
        q = q + 1;
        quartetMembers(q,:) = base + [1 2 5 6];
        stimToQuartet(base + [1 2 5 6]) = q;
        q = q + 1;
        quartetMembers(q,:) = base + [3 4 7 8];
        stimToQuartet(base + [3 4 7 8]) = q;
    end

    % Per-site, per-quartet pooled delta
    nSites = numel(siteRange);
    deltaQ = nan(nSites, nQuartets);
    requiredVarsQ = ["assignment","overlap"];
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

            good = isfinite(EY) & isfinite(wAll) & (wAll > 0) & ~isBG & ~isOV;

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

    % Accumulate geometry across selected stimuli
    if isempty(optsPlot.stimIdx)
        stimList = 1:numel(Tall_V1);
    else
        stimList = optsPlot.stimIdx(:)';
    end
    for stimNum = stimList
        if stimNum < 1 || stimNum > numel(Tall_V1)
            continue;
        end
        qIdx = stimToQuartet(stimNum);
        if qIdx < 1
            continue;
        end
        if ~isfield(Tall_V1(stimNum),'T') || ~istable(Tall_V1(stimNum).T)
            continue;
        end
        Ttab = Tall_V1(stimNum).T;
        if height(Ttab) < max(siteRange) || any(~ismember(requiredVars, string(Ttab.Properties.VariableNames)))
            continue;
        end

        assign = string(Ttab.assignment(siteRange));
        isTargetAssign = (assign == "target");
        isDistrAssign  = (assign == "distractor");
        isOverlap = Ttab.overlap(siteRange) ~= 0;
        dHere = deltaQ(:, qIdx);

        validBase = sigMask & ~isOverlap & (isTargetAssign | isDistrAssign) & isfinite(dHere);

        valTarget = zeros(size(siteGlobal));
        valDistr  = zeros(size(siteGlobal));
        idxT = validBase & isTargetAssign;
        idxD = validBase & isDistrAssign;
        valTarget(idxT) = max(0,  dHere(idxT));
        valDistr(idxD)  = max(0, -dHere(idxD));

        plotT = idxT & (valTarget > 0);
        plotD = idxD & (valDistr > 0);

        nTargetContrib = nTargetContrib + nnz(plotT);
        nDistrContrib  = nDistrContrib + nnz(plotD);
        valTargetAll = [valTargetAll; valTarget(plotT)]; %#ok<AGROW>
        valDistrAll  = [valDistrAll;  valDistr(plotD)]; %#ok<AGROW>

        if any(plotT)
            along = double(Ttab.along_GC(siteRange(plotT))) * widthEx;
            perp  = double(Ttab.perp_signed_GC(siteRange(plotT))) * widthEx;
            p = s_px + along.*uT + perp.*nT;
            X = [X; p(:,1)]; %#ok<AGROW>
            Y = [Y; p(:,2)]; %#ok<AGROW>
            V = [V; valTarget(plotT)]; %#ok<AGROW>
        end

        if any(plotD)
            along = double(Ttab.along_GC(siteRange(plotD))) * widthEx;
            perp  = double(Ttab.perp_signed_GC(siteRange(plotD))) * widthEx;
            p = s_px + along.*uD + perp.*nD;
            X = [X; p(:,1)]; %#ok<AGROW>
            Y = [Y; p(:,2)]; %#ok<AGROW>
            V = [V; valDistr(plotD)]; %#ok<AGROW>
        end
    end

elseif mode == "site_pooled_template"
    % New default:
    % - one pooled T-D effect per site (from OUT)
    % - one template location per site (stimID_example geometry)
    % - sign routes to target arm (T>D) or mirrored distractor arm (D>T)
    assignEx = string(TtabEx.assignment(siteRange));
    isCurveEx = (assignEx == "target") | (assignEx == "distractor");
    isOverlapEx = TtabEx.overlap(siteRange) ~= 0;
    validBase = sigMask & isCurveEx & ~isOverlapEx & isfinite(delta);

    valTarget = zeros(size(siteGlobal));
    valDistr  = zeros(size(siteGlobal));
    valTarget(validBase) = max(0,  delta(validBase));
    valDistr(validBase)  = max(0, -delta(validBase));

    plotT = validBase & (valTarget > 0);
    plotD = validBase & (valDistr > 0);

    nTargetContrib = nnz(plotT);
    nDistrContrib  = nnz(plotD);
    valTargetAll = valTarget(plotT);
    valDistrAll  = valDistr(plotD);

    if any(plotT)
        along = double(TtabEx.along_GC(siteRange(plotT))) * widthEx;
        perp  = double(TtabEx.perp_signed_GC(siteRange(plotT))) * widthEx;
        p = s_px + along.*uT + perp.*nT;
        X = [X; p(:,1)]; %#ok<AGROW>
        Y = [Y; p(:,2)]; %#ok<AGROW>
        V = [V; valTarget(plotT)]; %#ok<AGROW>
    end

    if any(plotD)
        along = double(TtabEx.along_GC(siteRange(plotD))) * widthEx;
        perp  = double(TtabEx.perp_signed_GC(siteRange(plotD))) * widthEx;
        p = s_px + along.*uD + perp.*nD;
        X = [X; p(:,1)]; %#ok<AGROW>
        Y = [Y; p(:,2)]; %#ok<AGROW>
        V = [V; valDistr(plotD)]; %#ok<AGROW>
    end

elseif mode == "all_stimuli"
    % Legacy behavior: accumulate geometry across selected stimuli.
    if isempty(optsPlot.stimIdx)
        stimList = 1:numel(Tall_V1);
    else
        stimList = optsPlot.stimIdx(:)';
    end

    for stimNum = stimList
        if stimNum < 1 || stimNum > numel(Tall_V1)
            continue;
        end
        if ~isfield(Tall_V1(stimNum),'T') || ~istable(Tall_V1(stimNum).T)
            continue;
        end
        Ttab = Tall_V1(stimNum).T;
        if height(Ttab) < max(siteRange)
            continue;
        end
        if any(~ismember(requiredVars, string(Ttab.Properties.VariableNames)))
            continue;
        end

        assign = string(Ttab.assignment(siteRange));
        isTargetAssign = (assign == "target");
        isDistrAssign  = (assign == "distractor");
        isOverlap = Ttab.overlap(siteRange) ~= 0;

        validBase = sigMask & ~isOverlap & (isTargetAssign | isDistrAssign) & isfinite(delta);

        valTarget = zeros(size(siteGlobal));
        valDistr  = zeros(size(siteGlobal));
        idxT = validBase & isTargetAssign;
        idxD = validBase & isDistrAssign;
        valTarget(idxT) = max(0,  delta(idxT));
        valDistr(idxD)  = max(0, -delta(idxD));

        plotT = idxT & (valTarget > 0);
        plotD = idxD & (valDistr > 0);

        nTargetContrib = nTargetContrib + nnz(plotT);
        nDistrContrib  = nDistrContrib + nnz(plotD);
        valTargetAll = [valTargetAll; valTarget(plotT)]; %#ok<AGROW>
        valDistrAll  = [valDistrAll;  valDistr(plotD)]; %#ok<AGROW>

        if any(plotT)
            along = double(Ttab.along_GC(siteRange(plotT))) * widthEx;
            perp  = double(Ttab.perp_signed_GC(siteRange(plotT))) * widthEx;
            p = s_px + along.*uT + perp.*nT;
            X = [X; p(:,1)]; %#ok<AGROW>
            Y = [Y; p(:,2)]; %#ok<AGROW>
            V = [V; valTarget(plotT)]; %#ok<AGROW>
        end

        if any(plotD)
            along = double(Ttab.along_GC(siteRange(plotD))) * widthEx;
            perp  = double(Ttab.perp_signed_GC(siteRange(plotD))) * widthEx;
            p = s_px + along.*uD + perp.*nD;
            X = [X; p(:,1)]; %#ok<AGROW>
            Y = [Y; p(:,2)]; %#ok<AGROW>
            V = [V; valDistr(plotD)]; %#ok<AGROW>
        end
    end
else
    error('Unknown optsPlot.projectionMode: %s', optsPlot.projectionMode);
end

fprintf('stim %d contributing assignments (target curve): %d\n', stimID_example, nTargetContrib);
fprintf('stim %d contributing assignments (distractor curve): %d\n', stimID_example, nDistrContrib);
printStats('target', valTargetAll);
printStats('distractor', valDistrAll);

% ---- Robust color scaling (95th percentile of plotted positives) ----
if ~isempty(optsPlot.cMaxFixed) && isfinite(optsPlot.cMaxFixed) && optsPlot.cMaxFixed > 0
    cMax = optsPlot.cMaxFixed;
elseif isempty(V)
    cMax = 1;
else
    cMax = prctile(V, optsPlot.robustPct);
    if ~isfinite(cMax) || cMax <= 0
        cMax = max(V);
    end
    if ~isfinite(cMax) || cMax <= 0
        cMax = 1;
    end
end

t = min(max(V ./ cMax, 0), 1);
C = (1-t).*optsPlot.cLow + t.*optsPlot.cHigh;

% ---- Plot single panel ----
fig = figure('Color',optsPlot.bgColor);
ax = axes('Position',[0 0 1 1]); hold(ax,'on');

img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, stimID_example);
imshow(img, 'Parent', ax, 'InitialMagnification','fit');
set(ax,'Position',[0 0 1 1], 'Color', optsPlot.bgColor);
axis(ax,'ij');

if ~isempty(X)
    hSc = scatter(ax, X, Y, optsPlot.markerSize, C, 'filled');
    hSc.MarkerFaceAlpha = optsPlot.alpha;
    if isprop(hSc,'MarkerEdgeAlpha')
        hSc.MarkerEdgeAlpha = optsPlot.alpha;
    end
end

xlim(ax, [1 W]);
ylim(ax, [1 H]);
axis(ax,'equal');
set(ax,'YDir','reverse');

hFrame = rectangle(ax,'Position',[0.5 0.5 W H], 'EdgeColor',[0.85 0.85 0.85], 'LineWidth',1);
uistack(hFrame,'top');

title(ax, sprintf('Attention effect on stim %d (pTD < %.3f)', stimID_example, optsPlot.pThresh), ...
    'Color','w','FontSize',14);

h = struct();
h.fig = fig;
h.ax = ax;
h.nSigSites = nSigSites;
h.nTargetContrib = nTargetContrib;
h.nDistrContrib = nDistrContrib;
h.cMax = cMax;
h.stimID = stimID_example;
h.V = V;
end

function printStats(label, v)
if isempty(v)
    fprintf('%s plotted positive values: min=NaN, median=NaN, max=NaN\n', label);
else
    fprintf('%s plotted positive values: min=%.6g, median=%.6g, max=%.6g\n', ...
        label, min(v), median(v), max(v));
end
end

function n = perpTowardOther(s_px, t_arm, t_other)
v = t_arm - s_px;
u = v / norm(v);
w = t_other - s_px;
w_perp = w - dot(w,u)*u;
n = w_perp / norm(w_perp);
end
