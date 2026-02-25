function make_attention_movie_wrapper_safe(outFile, ...
    Tall_V1, ALLCOORDS, RTAB384, stimID_example, R_resp, SNRnorm, OUT3, varargin)
% MAKE_ATTENTION_MOVIE_WRAPPER_SAFE
% Render a time-resolved movie of attention effect (T-D) in 10 ms bins.
%
% Hard significance rule:
%   site is included iff OUT3.pValueTD < pThresh (default 0.05)
%
% Inputs:
%   outFile         : output .mp4 path
%   Tall_V1         : stimulus/site assignment table struct (1x384)
%   ALLCOORDS       : stimulus geometry struct
%   RTAB384         : stimulus metadata table
%   stimID_example  : stimulus ID for rendering frame (1..384)
%   R_resp          : response struct with many time windows (Resp_capsules_*)
%   SNRnorm         : normalization struct (muSpont/mu* fields)
%   OUT3            : baseline OUT struct providing pValueTD mask
%
% Optional name/value pairs:
%   'pThresh'         (0.05)
%   'frameRate'       (10)
%   'excludeOverlap'  (true)
%   'siteRange'       (1:512)
%   'alpha'           (0.50)
%   'markerSize'      (5)
%   'robustPct'       (95)
%   'bgColor'         ([0.5 0.5 0.5])
%   'cLow'            ([0.55 0.55 0.55])
%   'cHigh'           ([0.85 0.05 0.05])
%   'stimIdx'         ([])
%   'hideStimPreZero' (true)
%   'verbose'         (true)
%   'useWaitbar'      (true)
%   'neighborN'       (30)   % RF-center nearest significant neighbors
%   'pixelNeighborN'  (0)    % image-space KNN averaging after affine projection
%   'pixelSmoothSigma' (0)   % image-space smoothing after affine projection (pixels)
%   'alphaByValue'    (true) % fixed per-point fading by value (same for all frames)
%   'alphaValueGamma' (2.5)  % nonlinearity for value-based alpha
%   'alphaValueMinScale' (0.03) % minimum alpha scale for weakest values
%   'alphaFloorSoft'  (0.20) % softness around pre-stim alpha floor threshold
%   'hotScale'         (true) % red->yellow->white for strongest values
%   'colorHotMaxFactor' (3.0) % allow values > cMax to become hotter
%   'preStimCalibratedAlpha' (true) % set alpha floor from pre-stim distribution
%   'preStimPercentile' (99.7) % percentile used to set alpha floor
%   'cMaxFixed'        ([])   % if set, use this shared color reference for all frames
%   'alphaValueFloorFixed' ([]) % if set, use this fixed alpha floor for all frames

% ---- Defaults ----
optsMovie = struct();
optsMovie.pThresh = 0.05;
optsMovie.frameRate = 10;
optsMovie.excludeOverlap = true;
optsMovie.siteRange = 1:512;
optsMovie.alpha = 0.50;
optsMovie.markerSize = 5;
optsMovie.robustPct = 95;
optsMovie.bgColor = [0.5 0.5 0.5];
optsMovie.cLow = [0.55 0.55 0.55];
optsMovie.cHigh = [0.85 0.05 0.05];
optsMovie.stimIdx = [];
optsMovie.hideStimPreZero = true;
optsMovie.verbose = true;
optsMovie.useWaitbar = true;
optsMovie.neighborN = 30;
optsMovie.pixelNeighborN = 0;
optsMovie.pixelSmoothSigma = 0;
optsMovie.alphaByValue = true;
optsMovie.alphaValueGamma = 2.5;
optsMovie.alphaValueMinScale = 0.03;
optsMovie.alphaFloorSoft = 0.20;
optsMovie.hotScale = true;
optsMovie.colorHotMaxFactor = 3.0;
optsMovie.preStimCalibratedAlpha = true;
optsMovie.preStimPercentile = 99.7;
optsMovie.cMaxFixed = [];
optsMovie.alphaValueFloorFixed = [];

% ---- Parse name/value ----
if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k+1};
    if ~ischar(name)
        error('Optional argument names must be char.');
    end
    if ~isfield(optsMovie, name)
        error('Unknown option: %s', name);
    end
    optsMovie.(name) = value;
end

% ---- Validate inputs ----
if ~isstruct(R_resp) || ~isfield(R_resp, 'timeWindows') || ~isfield(R_resp, 'meanAct')
    error(['R_resp must be a response struct with fields timeWindows and meanAct. ' ...
           'Use the struct loaded from Resp_capsules_* (many windows).']);
end
nFrames = size(R_resp.timeWindows,1);
if nFrames <= 3
    error(['R_resp.timeWindows has %d bins. ' ...
           'For this movie use the response struct from Resp_capsules_* (e.g., ~70 bins).'], nFrames);
end
if size(R_resp.meanAct,3) ~= nFrames
    error('Mismatch: size(R_resp.meanAct,3)=%d but size(R_resp.timeWindows,1)=%d.', ...
        size(R_resp.meanAct,3), nFrames);
end

if ~isstruct(OUT3) || ~isfield(OUT3,'pValueTD')
    error('OUT3 must contain pValueTD for fixed significance masking.');
end
if numel(OUT3.pValueTD) < max(optsMovie.siteRange)
    error('OUT3.pValueTD is smaller than max(siteRange).');
end

siteRange = optsMovie.siteRange(:)';
sigMask = isfinite(OUT3.pValueTD(siteRange)) & (OUT3.pValueTD(siteRange) < optsMovie.pThresh);
nSigSites = nnz(sigMask);

% RF-center neighbor map (visual space) for smoothing
Tref = Tall_V1(stimID_example).T;
[rfX, rfY] = get_rf_centers(Tref, siteRange);
neighborIdx = compute_neighbor_idx(rfX, rfY, sigMask, optsMovie.neighborN);

if optsMovie.verbose
    fprintf('\n============================================\n');
    fprintf('Starting attention-effect movie rendering (%d frames)\n', nFrames);
    fprintf('Output file: %s\n', outFile);
    fprintf('Example stimulus: %d\n', stimID_example);
    fprintf('Significant sites (fixed mask): %d / %d\n', nSigSites, numel(siteRange));
    fprintf('Neighbor smoothing: N=%d (RF-center distance)\n', optsMovie.neighborN);
    fprintf('Pixel neighbor averaging: N=%d (post-affine)\n', optsMovie.pixelNeighborN);
    fprintf('Pixel smoothing sigma: %.3g px (post-affine)\n', optsMovie.pixelSmoothSigma);
    fprintf('Value-based alpha: %d (gamma=%g, minScale=%g)\n', ...
        optsMovie.alphaByValue, optsMovie.alphaValueGamma, optsMovie.alphaValueMinScale);
    fprintf('============================================\n');
end

wb = [];
if optsMovie.useWaitbar
    try
        wb = waitbar(0, 'Precomputing attention bins...');
    catch
        wb = [];
    end
end

% ---- Precompute quartet-pooled deltas per frame + shared color scale ----
% Quartets per block of 8 stimuli: [1 2 5 6] and [3 4 7 8]
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

deltaQ_byFrame = cell(nFrames,1);
deltaAbs = [];
frameStrength = nan(nFrames,1);
preStimVals = [];
tp = tic;
for tb = 1:nFrames
    deltaQ = compute_deltaQ_quartet(tb, R_resp, Tall_V1, SNRnorm, siteRange, quartetMembers, optsMovie.excludeOverlap);
    deltaQ_byFrame{tb} = deltaQ;

    d = abs(deltaQ(sigMask, :));
    d = d(:);
    d = d(isfinite(d));
    if ~isempty(d)
        deltaAbs = [deltaAbs; d]; %#ok<AGROW>
        frameStrength(tb) = prctile(d, optsMovie.robustPct);
        if R_resp.timeWindows(tb,2) <= 0
            preStimVals = [preStimVals; d]; %#ok<AGROW>
        end
    end

    if optsMovie.verbose
        tElapsed = toc(tp);
        pct = 100 * tb / nFrames;
        tPer = tElapsed / tb;
        eta = tPer * (nFrames - tb);
        fprintf('[Precompute] %2d/%2d | %5.1f%% | elapsed %.1fs | ETA %.1fs\n', ...
            tb, nFrames, pct, tElapsed, eta);
    end
    if ~isempty(wb)
        waitbar(0.45 * tb / nFrames, wb, sprintf('Precomputing bins... %d/%d', tb, nFrames));
    end
end

if ~isempty(optsMovie.cMaxFixed) && isfinite(optsMovie.cMaxFixed) && optsMovie.cMaxFixed > 0
    cShared = optsMovie.cMaxFixed;
else
    if isempty(deltaAbs)
        cShared = 1;
    else
        cShared = prctile(deltaAbs, optsMovie.robustPct);
        if ~isfinite(cShared) || cShared <= 0
            cShared = max(deltaAbs);
        end
        if ~isfinite(cShared) || cShared <= 0
            cShared = 1;
        end
    end
end

if optsMovie.verbose
    fprintf('Shared color scale cMax: %.6g\n', cShared);
end

if ~isempty(optsMovie.alphaValueFloorFixed) && isfinite(optsMovie.alphaValueFloorFixed) && optsMovie.alphaValueFloorFixed >= 0
    alphaFloor = optsMovie.alphaValueFloorFixed;
    if optsMovie.verbose
        fprintf('Using fixed alpha floor: %.6g\n', alphaFloor);
    end
else
    alphaFloor = 0;
    if optsMovie.preStimCalibratedAlpha
        if isempty(preStimVals)
            alphaFloor = 0;
            if optsMovie.verbose
                fprintf('Pre-stim alpha floor disabled: no pre-stim values found.\n');
            end
        else
            alphaFloor = prctile(preStimVals, optsMovie.preStimPercentile);
            if ~isfinite(alphaFloor) || alphaFloor < 0
                alphaFloor = 0;
            end
            if optsMovie.verbose
                fprintf('Pre-stim alpha floor: p%d = %.6g\n', optsMovie.preStimPercentile, alphaFloor);
            end
        end
    end
end

% ---- Video setup ----
v = VideoWriter(outFile, 'MPEG-4');
v.FrameRate = optsMovie.frameRate;
open(v);

tStart = tic;

for tb = 1:nFrames
    if optsMovie.verbose
        tElapsed = toc(tStart);
        pct = 100 * tb / nFrames;
        tPer = tElapsed / max(tb-1,1);
        if tb == 1
            eta = NaN;
        else
            eta = tPer * (nFrames - tb + 1);
        end
        fprintf('Frame %2d / %2d  |  %5.1f%%  |  %4d-%4d ms  |  elapsed: %.1fs\n', ...
            tb, nFrames, pct, ...
            R_resp.timeWindows(tb,1), R_resp.timeWindows(tb,2), tElapsed);
        if isfinite(eta)
            fprintf('                 ETA %.1fs\n', eta);
        end
    end

    optsPlot = struct();
    optsPlot.RTAB384 = RTAB384;
    optsPlot.pThresh = optsMovie.pThresh;
    optsPlot.siteRange = siteRange;
    optsPlot.markerSize = optsMovie.markerSize;
    optsPlot.alpha = optsMovie.alpha;
    optsPlot.robustPct = optsMovie.robustPct;
    optsPlot.stimIdx = optsMovie.stimIdx;
    optsPlot.cMaxFixed = cShared;
    optsPlot.bgColor = optsMovie.bgColor;
    optsPlot.cLow = optsMovie.cLow;
    optsPlot.cHigh = optsMovie.cHigh;
    optsPlot.projectionMode = 'quartet_pooled';
    optsPlot.deltaQ = deltaQ_byFrame{tb};
    optsPlot.stimToQuartet = stimToQuartet;
    optsPlot.neighborN = optsMovie.neighborN;
    optsPlot.pixelNeighborN = optsMovie.pixelNeighborN;
    optsPlot.pixelSmoothSigma = optsMovie.pixelSmoothSigma;
    optsPlot.neighborIdx = neighborIdx;
    optsPlot.alphaByValue = optsMovie.alphaByValue;
    optsPlot.alphaValueGamma = optsMovie.alphaValueGamma;
    optsPlot.alphaValueMinScale = optsMovie.alphaValueMinScale;
    optsPlot.alphaFloorSoft = optsMovie.alphaFloorSoft;
    optsPlot.alphaValueFloor = alphaFloor;
    optsPlot.hotScale = optsMovie.hotScale;
    optsPlot.colorHotMaxFactor = optsMovie.colorHotMaxFactor;

    h = plot_projected_attentiondiff_on_example_stim( ...
        stimID_example, OUT3, Tall_V1, ALLCOORDS, optsPlot);

    % Optional: hide stimulus outlines in pre-zero frames
    isPreZero = (R_resp.timeWindows(tb,2) <= 0);
    if optsMovie.hideStimPreZero && isPreZero
        ax = h.ax;
        hImg = findobj(ax, 'Type', 'image');
        if ~isempty(hImg)
            hImg = hImg(1);
            C = hImg.CData;
            [hh, ww, ~] = size(C);
            g = uint8(round(255 * optsMovie.bgColor(1)));
            hImg.CData = repmat(g, [hh, ww, 3]);
        end
        set(findobj(ax, 'Type', 'line'),  'Visible', 'off');
        set(findobj(ax, 'Type', 'patch'), 'Visible', 'off');
    end

    figure(h.fig);
    set(h.fig, 'Name', 'Attention TD Movie (10 ms bins)', 'NumberTitle', 'off');
    title(h.ax, sprintf('Attention effect (T-D) | frame %d/%d | %d-%d ms | stim %d', ...
        tb, nFrames, R_resp.timeWindows(tb,1), R_resp.timeWindows(tb,2), stimID_example), ...
        'Color','w','FontSize',14);

    text(h.ax, 0.02, 0.98, sprintf('t = %.0f ms', mean(R_resp.timeWindows(tb,:))), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize', 14, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'BackgroundColor','k', ...
        'Margin', 4);

    drawnow;
    frame = getframe(h.fig);
    writeVideo(v, frame);
    close(h.fig);
    if ~isempty(wb)
        waitbar(0.45 + 0.55 * tb / nFrames, wb, sprintf('Rendering frames... %d/%d', tb, nFrames));
    end
end

close(v);
if ~isempty(wb)
    try
        close(wb);
    catch
    end
end

if optsMovie.verbose
    fprintf('============================================\n');
    fprintf('Attention movie finished in %.1f seconds\n', toc(tStart));
    fprintf('Saved to: %s\n', outFile);
    fprintf('============================================\n\n');
end
end

function [x, y] = get_rf_centers(Ttab, siteRange)
vn = string(Ttab.Properties.VariableNames);

if ismember("x_deg", vn) && ismember("y_deg", vn)
    x = double(Ttab.x_deg(siteRange));
    y = double(Ttab.y_deg(siteRange));
    return;
end
if ismember("x", vn) && ismember("y", vn)
    x = double(Ttab.x(siteRange));
    y = double(Ttab.y(siteRange));
    return;
end
if ismember("x_px", vn) && ismember("y_px", vn)
    x = double(Ttab.x_px(siteRange));
    y = double(Ttab.y_px(siteRange));
    return;
end

error(['Could not find RF center coordinates in Tall_V1(stim).T. ' ...
       'Expected (x_deg,y_deg), (x,y), or (x_px,y_px).']);
end

function nbrIdx = compute_neighbor_idx(x, y, sigMask, neighborN)
n = numel(x);
nbrIdx = nan(n, neighborN);

sigIdx = find(sigMask(:));
if isempty(sigIdx)
    return;
end

xs = x(sigIdx); ys = y(sigIdx);
for i = 1:n
    if ~isfinite(x(i)) || ~isfinite(y(i))
        continue;
    end
    d2 = (xs - x(i)).^2 + (ys - y(i)).^2;
    [~, ord] = sort(d2, 'ascend');
    k = min(neighborN, numel(ord));
    nbrIdx(i,1:k) = sigIdx(ord(1:k));
end
end

function deltaQ = compute_deltaQ_quartet(tb, Rdata, Tall_V1, SNRn, siteRange, quartetMembers, excludeOverlap)
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
