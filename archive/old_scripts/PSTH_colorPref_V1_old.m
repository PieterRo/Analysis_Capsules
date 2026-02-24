%% Balanced time course: preferred vs non-preferred color (10 ms bins)
% Requires:
%   R.meanAct    [1024 x 384 x nBins]  (you have 70 bins)
%   R.nTrials    [384 x 1]
%   R.timeWindows [nBins x 2]          (e.g. -200..500 in 10 ms steps)
%   Tall_V1(384x1).T is a table with center_color per site (row) (512 rows)
%   ColorTune.early.colorIndex or ColorTune.late.colorIndex (512x1) from previous analysis

% -----------------------
% SETTINGS
% -----------------------
useColorIndexFrom = "early";  % "early" or "late"
snrThr = SNRthr;                % optional extra gate: require SNR.bestSNR > snrThr
useSNRgate = true;

COL_Y = "yellowArm";
COL_P = "purple";
COL_G = "gray";

% -----------------------
% BASIC DIMS
% -----------------------
nStim = numel(R.nTrials);
nTrials = double(R.nTrials(:));
assert(nStim == 384, 'Expected 384 stimuli.');

[nCh, nStim2, nBins] = size(R.meanAct);
assert(nStim2 == nStim, 'R.meanAct stimulus dim mismatch.');
assert(size(R.timeWindows,1) == nBins, 'R.timeWindows rows must equal #bins.');

% V1 sites: rows 1..512 in R, and rows 1..512 in Tall_V1(stim).T
v1Sites = 1:512;

% time axis (bin centers)
tCenters = mean(R.timeWindows, 2);   % nBins x 1

% -----------------------
% Choose which colorIndex to use for tuning selection + preference
% -----------------------
switch lower(useColorIndexFrom)
    case "early"
        ci = ColorTune.early.colorIndex(:);  % 512x1
    case "late"
        ci = ColorTune.late.colorIndex(:);   % 512x1
    otherwise
        error('useColorIndexFrom must be "early" or "late".');
end

% Optional: also require SNR threshold
if useSNRgate
    keep = isfinite(ci) & (abs(ci) > ciThr) & isfinite(bestSNR(:)) & (bestSNR(:) > snrThr);
else
    keep = isfinite(ci) & (abs(ci) > ciThr);
end

keepSites = find(keep);   % indices 1..512
fprintf('Selected %d / 512 sites (abs(CI)>%.2f)%s\n', numel(keepSites), ciThr, ...
    ternary(useSNRgate, sprintf(' & bestSNR>%.2f', snrThr), ''));

% Preferred color per site
prefIsYellow = ci > 0;    % CI = (Y-P)/... ; positive -> prefers yellow

% -----------------------
% Sort Tall_V1 by stimNum (safety)
% -----------------------
TallStimNums = arrayfun(@(x) x.stimNum, Tall_V1(:));
[sortedStimNums, order] = sort(TallStimNums(:));
assert(all(sortedStimNums(:).' == 1:nStim), 'Tall_V1.stimNum should cover 1..384.');
TallSorted = Tall_V1(order);

% Identify center_color column in the table
T0 = TallSorted(1).T;
assert(istable(T0), 'Tall_V1(stim).T must be a table.');
vn = string(T0.Properties.VariableNames);
colIdx = find(vn == "center_color", 1);
if isempty(colIdx)
    colIdx = find(contains(lower(vn), "center") & contains(lower(vn), "color"), 1);
end
assert(~isempty(colIdx), 'Could not find center_color column.');

% Build CC(site,stim): 512 x 384 string matrix
CC = strings(512, nStim);
for stim = 1:nStim
    labs = string(TallSorted(stim).T{:, colIdx});   % 512x1
    CC(:,stim) = strtrim(labs);
end

% -----------------------
% Complement pairs within each block of 8: (1<->5),(2<->6),(3<->7),(4<->8), repeat
% We'll store only "first half" (pos 1..4) and pair with +4
% -----------------------
pairsA = zeros(nStim/2,1);
pairsB = zeros(nStim/2,1);
k = 0;
for a = 1:nStim
    pos = mod(a-1, 8) + 1;
    if pos <= 4
        k = k + 1;
        pairsA(k) = a;
        pairsB(k) = a + 4;
    end
end
nPairs = k;  % should be 192

% -----------------------
% Compute balanced time courses per site: muYellow(t), muPurple(t)
% then map to preferred/nonpreferred per site
% -----------------------
muPref = nan(numel(keepSites), nBins);
muNon  = nan(numel(keepSites), nBins);

for ii = 1:numel(keepSites)
    site = keepSites(ii);     % 1..512 index into Tall / V1 list
    ch   = v1Sites(site);     % row in R.meanAct (here identity)

    % accumulators over time bins
    sumY  = zeros(nBins,1);  sumY2 = zeros(nBins,1);  NY = 0;
    sumP  = zeros(nBins,1);  sumP2 = zeros(nBins,1);  NP = 0;

    for ip = 1:nPairs
        a = pairsA(ip);
        b = pairsB(ip);

        na = nTrials(a); nb = nTrials(b);
        nEff = min(na, nb);
        if nEff <= 0
            continue;
        end

        ca = CC(site,a);
        cb = CC(site,b);

        % skip if RF is background (gray) in this pair
        if (ca == COL_G) && (cb == COL_G)
            continue;
        end

        % time courses for these stimuli (nBins x 1)
        ma  = squeeze(R.meanAct(ch,a,:));
        mb  = squeeze(R.meanAct(ch,b,:));
        m2a = squeeze(R.meanSqAct(ch,a,:));
        m2b = squeeze(R.meanSqAct(ch,b,:));

        % Add balanced contributions for stimulus a
        if ca == COL_Y
            sumY  = sumY  + nEff * ma;
            sumY2 = sumY2 + nEff * m2a;
            NY = NY + nEff;
        elseif ca == COL_P
            sumP  = sumP  + nEff * ma;
            sumP2 = sumP2 + nEff * m2a;
            NP = NP + nEff;
        end

        % Add balanced contributions for stimulus b
        if cb == COL_Y
            sumY  = sumY  + nEff * mb;
            sumY2 = sumY2 + nEff * m2b;
            NY = NY + nEff;
        elseif cb == COL_P
            sumP  = sumP  + nEff * mb;
            sumP2 = sumP2 + nEff * m2b;
            NP = NP + nEff;
        end
    end

    if NY < 1 || NP < 1
        continue;
    end

    muY = (sumY / NY).';   % 1 x nBins
    muP = (sumP / NP).';   % 1 x nBins

    if prefIsYellow(site)
        muPref(ii,:) = muY;
        muNon(ii,:)  = muP;
    else
        muPref(ii,:) = muP;
        muNon(ii,:)  = muY;
    end
end

% Drop sites that ended up NaN everywhere (e.g., no usable trials)
good = any(isfinite(muPref), 2) & any(isfinite(muNon), 2);
muPref = muPref(good,:);
muNon  = muNon(good,:);
fprintf('After requiring usable balanced trials: N=%d sites\n', size(muPref,1));

% -----------------------
% Plot mean ± SEM over sites
% -----------------------
mPref = mean(muPref, 1, 'omitnan');
mNon  = mean(muNon,  1, 'omitnan');

semPref = std(muPref, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(muPref),1));
semNon  = std(muNon,  0, 1, 'omitnan') ./ sqrt(sum(isfinite(muNon),1));

figure; hold on;

% mean lines
plot(tCenters, mPref, 'LineWidth', 2);
plot(tCenters, mNon,  'LineWidth', 2);

% shaded error (simple patch)
fill([tCenters; flipud(tCenters)], [ (mPref-semPref)'; flipud((mPref+semPref)') ], ...
     'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');  % color auto will be black-ish
fill([tCenters; flipud(tCenters)], [ (mNon-semNon)';  flipud((mNon+semNon)')  ], ...
     'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none');

xline(0,'k-');
xlabel('Time from stimulus onset (ms)');
ylabel('Mean response (a.u.)');
title(sprintf('Balanced time course: preferred vs nonpreferred (abs(CI)>%.2f, N=%d)', ciThr, size(muPref,1)));
legend('Preferred','Nonpreferred','Location','best');
grid on;

%% Helper (inline ternary)
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end