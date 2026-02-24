%% Weighted pair-wise color decoding over time (no fitting, complementary balanced)
% Each complementary pair is a sample.
% For each pair and each time bin:
%   score = weighted mean over sites of (preferredColorResponse - nonpreferredColorResponse)
% Sites on gray are ignored (missing dimensions).
% Weights w(site) = abs(dprime(site)) from ColorTune.(refWin).dprime.

% =========================
% REQUIRED INPUTS
% =========================
% R.meanAct    : [1024 x 384 x nBins]  (10 ms bins)
% R.meanSqAct  : [1024 x 384 x nBins]
% R.nTrials    : [384 x 1]
% R.timeWindows: [nBins x 2]  e.g. -200..500 in 10 ms
% Tall_V1      : [384 x 1] struct with fields stimNum and T (table 512 rows) with center_color
% ColorTune.early (and/or .late) with fields colorIndex, dprime
% Optional: SNR.sdSpont (512x1) from your previous SNR script

% =========================
% SETTINGS
% =========================
refWin = "early";      % "early" or "late" -> which ColorTune.* fields define preference + dprime weights
capWeightPct = 95;     % cap |d'| weights at this percentile to avoid a few huge weights dominating
useSNRsdSpontIfPresent = true;  % if SNR.sdSpont exists, use it
computeSpontSDIfMissing = true; % if no SNR.sdSpont, compute from R pre-0 bins
normalizeBySpontSD = true;      % normalize per-site contributions by spont SD

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

% =========================
% BASIC DIMS
% =========================
nStim = numel(R.nTrials);
assert(nStim == 384, 'Expected 384 stimuli.');
nTrials = double(R.nTrials(:));

[nCh, nStim2, nBins] = size(R.meanAct);
assert(nStim2 == nStim, 'R.meanAct stimulus dim mismatch.');
assert(all(size(R.meanSqAct) == size(R.meanAct)), 'R.meanSqAct must match R.meanAct size.');
assert(size(R.timeWindows,1) == nBins, 'R.timeWindows rows must equal #bins.');

v1Sites = 1:512;      % mapping: site i -> channel i in R
tCenters = mean(R.timeWindows, 2);  % nBins x 1

% =========================
% Choose preference + weights from ColorTune
% =========================
switch lower(refWin)
    case "early"
        ci = ColorTune.early.colorIndex(:);   % 512x1
        dp = ColorTune.early.dprime(:);       % 512x1
    case "late"
        ci = ColorTune.late.colorIndex(:);
        dp = ColorTune.late.dprime(:);
    otherwise
        error('refWin must be "early" or "late".');
end

% Preference direction: CI>0 => prefers yellow, CI<0 => prefers purple
prefIsYellow = (ci > 0);
prefIsYellow(~isfinite(ci)) = false; % arbitrary, but those will get ~0 weight anyway

% Weights: abs(d')
w = abs(dp);
w(~isfinite(w)) = 0;

% Cap extreme weights robustly
wPos = w(w>0);
if ~isempty(wPos)
    wCap = prctile(wPos, capWeightPct);
    w = min(w, wCap);
end

fprintf('Weight summary |d''|: median=%.3f, 95%%=%.3f, max(after cap)=%.3f\n', ...
    median(w(w>0),'omitnan'), prctile(w(w>0),95), max(w));

% If a site has 0 weight, it contributes nothing (no need to threshold by CI)
useSites = find(w > 0);
fprintf('Using %d sites with nonzero weight\n', numel(useSites));

% =========================
% Sort Tall_V1 by stimNum and extract center_color per site x stim
% =========================
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..384.');
TallSorted = Tall_V1(order);

T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
colIdx = find(vn=="center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column in Tall_V1(stim).T.');

CC = strings(512, nStim);
for stim = 1:nStim
    CC(:,stim) = strtrim(string(TallSorted(stim).T{:, colIdx})); % 512x1
end

% =========================
% Build complementary pairs (1<->5),(2<->6),(3<->7),(4<->8), repeating
% =========================
pairsA = []; pairsB = [];
for a = 1:nStim
    pos = mod(a-1,8) + 1;
    if pos <= 4
        pairsA(end+1,1) = a; %#ok<AGROW>
        pairsB(end+1,1) = a + 4; %#ok<AGROW>
    end
end
nPairs = numel(pairsA);
assert(nPairs == 192, 'Expected 192 complementary pairs, got %d.', nPairs);
fprintf('Using %d complementary pairs\n', nPairs);

% =========================
% Spontaneous SD per site (for normalization)
% =========================
sdSpont = ones(512,1);

if normalizeBySpontSD
    if useSNRsdSpontIfPresent && exist('SNR','var') && isfield(SNR,'sdSpont') && numel(SNR.sdSpont) >= 512
        sdSpont = SNR.sdSpont(1:512);
        sdSpont(~isfinite(sdSpont) | sdSpont<=0) = 1;
        fprintf('Using sdSpont from SNR.sdSpont\n');

    elseif computeSpontSDIfMissing
        % Compute from pre-0 bins in the 10 ms representation
        isSpontBin = (R.timeWindows(:,2) <= 0);
        assert(any(isSpontBin), 'No spontaneous bins found (timeWindows end<=0).');

        bIdx = find(isSpontBin);

        % We pool over stimuli using trial weights (no need to balance here)
        wStim = nTrials;
        wStim = wStim / sum(wStim);

        fprintf('Computing sdSpont from %d pre-0 bins using R.meanAct/meanSqAct...\n', numel(bIdx));

        for site = 1:512
            ch = v1Sites(site);

            varBins = nan(numel(bIdx),1);
            for ib = 1:numel(bIdx)
                tb = bIdx(ib);

                m  = squeeze(R.meanAct(ch,:,tb)).';    % 384x1
                m2 = squeeze(R.meanSqAct(ch,:,tb)).';

                good = isfinite(m) & isfinite(m2) & (nTrials>0);
                if ~any(good), continue; end

                ww = wStim(good); ww = ww / sum(ww);
                mu = sum(ww .* m(good));
                ex2 = sum(ww .* m2(good));
                varBins(ib) = max(0, ex2 - mu^2);
            end

            sd = sqrt(mean(varBins, 'omitnan'));
            if isfinite(sd) && sd>0
                sdSpont(site) = sd;
            else
                sdSpont(site) = 1;
            end
        end
        fprintf('Done computing sdSpont\n');
    else
        fprintf('normalizeBySpontSD=true but no sdSpont available; using sdSpont=1\n');
    end
end

% =========================
% MAIN: score per pair per time bin (weighted)
% =========================
scorePerPair = nan(nPairs, nBins);
nUsedSitesPerPair = nan(nPairs,1);

for ip = 1:nPairs
    a = pairsA(ip);
    b = pairsB(ip);

    % Pull responses for all V1 channels for these stimuli
    Ra = squeeze(R.meanAct(v1Sites, a, :)); % 512 x nBins
    Rb = squeeze(R.meanAct(v1Sites, b, :)); % 512 x nBins

    sumScore = zeros(1, nBins);
    sumW     = zeros(1, nBins);

    nUsed = 0;

    for site = useSites(:).'
        ca = CC(site,a);
        cb = CC(site,b);

        % Ignore if either member is gray (missing) or unknown
        if ca==COL_G || cb==COL_G || ca=="" || cb==""
            continue;
        end

        % Determine which stim has Yellow and which has Purple for this site
        if ca==COL_Y && cb==COL_P
            rY = Ra(site,:); rP = Rb(site,:);
        elseif ca==COL_P && cb==COL_Y
            rY = Rb(site,:); rP = Ra(site,:);
        else
            % unexpected label combo (e.g. not a clean Y/P swap), skip
            continue;
        end

        wi = w(site);
        if wi <= 0
            continue;
        end

        % Signed difference according to preference
        if prefIsYellow(site)
            d = (rY - rP);
        else
            d = (rP - rY);
        end

        % Normalize by spont SD (optional)
        d = d ./ sdSpont(site);

        sumScore = sumScore + wi * d;
        sumW     = sumW + wi;

        nUsed = nUsed + 1;
    end

    nUsedSitesPerPair(ip) = nUsed;

    ok = sumW > 0;
    tmp = nan(1, nBins);
    tmp(ok) = sumScore(ok) ./ sumW(ok);
    scorePerPair(ip,:) = tmp;
end

fprintf('Median #sites contributing per pair: %.1f\n', median(nUsedSitesPerPair, 'omitnan'));

% =========================
% Average across pairs
% =========================
mScore = mean(scorePerPair, 1, 'omitnan');
semScore = std(scorePerPair, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(scorePerPair),1));

figure; hold on;
plot(tCenters, mScore, 'LineWidth', 2);
fill([tCenters; flipud(tCenters)], [(mScore-semScore)'; flipud((mScore+semScore)')], ...
     'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
xline(0,'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Weighted signed score (pref - nonpref), a.u.');
title(sprintf('Weighted pair-wise color decoding (%s weights |d''|), Npairs=%d', refWin, nPairs));
grid on;

% Optional: quick pre/post sanity
tbPre  = find(tCenters < 0, 1, 'last');
tbPost = find(tCenters >= 100, 1, 'first');
if ~isempty(tbPre) && ~isempty(tbPost)
    fprintf('Mean score pre (%.0f ms): %.4f\n', tCenters(tbPre), mean(scorePerPair(:,tbPre),'omitnan'));
    fprintf('Mean score post (%.0f ms): %.4f\n', tCenters(tbPost), mean(scorePerPair(:,tbPost),'omitnan'));
end


