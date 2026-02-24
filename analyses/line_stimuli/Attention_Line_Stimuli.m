% Attention_Line_Stimuli
% Attention modulation summary, still-frame diagnostics, and optional movie wrapper.

%% Run switches
RUN_HISTOGRAM = true;
RUN_QC = true;
RUN_STILLS = true;
RUN_MOVIE = false;

%% Central parameters
P = struct();
P.timeIdx3 = 3;                      % 3-bin dataset: 300-500 ms
P.targetWin10 = [300 310];           % short window target (nearest bin)
P.targetWinPre = [-100 -90];         % pre-stim control window (nearest bin)
P.pThresh = 0.05;                    % hard significance rule
P.siteRange = 1:512;
P.stimIDExample = 38;
P.idxClip = 2;                       % histogram clipping
P.plotMarkerSize = 5;
P.plotAlpha = 0.50;                  % constant across windows
P.plotRobustPct = 95;
P.plotBgColor = [0.5 0.5 0.5];
P.plotCLow = [0.55 0.55 0.55];       % weak effects blend into background
P.plotCHigh = [0.85 0.05 0.05];

%% Required context
cfg = config();
assert(exist('RTAB384','var') == 1, 'RTAB384 must exist in workspace.');
assert(exist('ALLCOORDS','var') == 1, 'ALLCOORDS must exist in workspace.');

if ~exist('Tall_V1', 'var') || ~isstruct(Tall_V1)
    S = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));  % loads Tall_V1
    Tall_V1 = S.Tall_V1;
end

%% Load 3-bin response + normalization and get baseline OUT3
S = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'));   % loads R
assert(isfield(S, 'R') && isstruct(S.R), 'SNR_capsules_N_d12.mat must contain struct R.');
R3 = S.R;
assert(isfield(R3, 'timeWindows') && size(R3.timeWindows,1) == 3, ...
    'Expected 3 time windows in SNR_capsules_N_d12.mat.');

S = load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat')); % loads SNR
assert(isfield(S, 'SNR') && isstruct(S.SNR) && isfield(S.SNR, 'muSpont'), ...
    'SNR_V1_byColor_byWindow.mat must contain normalization struct SNR.');
SNRnorm = S.SNR;

saveTag = 'OUT_attention_modulation_3bin_timeIdx3';
outFile = fullfile(cfg.resultsDir, [saveTag '.mat']);

if exist(outFile, 'file')
    S = load(outFile, 'OUT');
    OUT3 = S.OUT;
    fprintf('Loaded baseline OUT from: %s\n', outFile);
else
    opts3 = struct('timeIdx', P.timeIdx3, 'excludeOverlap', true, 'verbose', true);
    OUT3 = attention_modulation_V1_3bin(R3, Tall_V1, SNRnorm, opts3);
    meta = struct();
    meta.created = datestr(now, 30);
    meta.script = mfilename;
    meta.timeIdx = opts3.timeIdx;
    meta.excludeOverlap = opts3.excludeOverlap;
    meta.note = 'Use OUT.pValueTD from this file as fixed significance mask.';
    save(outFile, 'OUT', 'meta', '-v7.3');
    fprintf('Computed and saved baseline OUT to: %s\n', outFile);
end

fprintf('Median d'': %.3f\n', median(OUT3.dprime, 'omitnan'));
fprintf('Median index: %.3f\n', median(OUT3.index, 'omitnan'));

%% Histogram: attention index distribution
if RUN_HISTOGRAM
    idxAll = OUT3.index;
    isSig = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh);

    vAllRaw = idxAll(isfinite(idxAll));
    vSigRaw = idxAll(isSig & isfinite(idxAll));
    vAll = max(min(vAllRaw, P.idxClip), -P.idxClip);
    vSig = max(min(vSigRaw, P.idxClip), -P.idxClip);

    fprintf('Clipped attention index for display at +/-%.1f (all: %d, sig: %d values clipped)\n', ...
        P.idxClip, nnz(abs(vAllRaw) > P.idxClip), nnz(abs(vSigRaw) > P.idxClip));

    figure('Color','w');
    hold on;
    histogram(vAll, 30, 'FaceColor', [0.80 0.80 0.80], 'EdgeColor', 'none');
    histogram(vSig, 30, 'FaceColor', [0.85 0.20 0.20], 'EdgeColor', 'none');
    xlabel('Attention index');
    ylabel('Number of sites');
    title(sprintf('Attention index (all vs significant, pTD < %.3f)', P.pThresh));
    legend(sprintf('All sites (N=%d)', numel(vAll)), ...
           sprintf('Significant sites (N=%d)', numel(vSig)), ...
           'Location','best');
    grid on;
end

%% Load high-resolution response windows
if RUN_QC || RUN_STILLS || RUN_MOVIE
    if ~exist('R_resp', 'var') || ~isstruct(R_resp)
        S = load(fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat'));  % loads R
        R_resp = S.R;
    end

    assert(isfield(R_resp, 'timeWindows') && size(R_resp.timeWindows,2) == 2, ...
        'Resp_capsules_N_d12.mat must contain R.timeWindows as [nWindows x 2].');
    assert(size(R_resp.timeWindows,1) > 3, ...
        ['Resp_capsules_N_d12.mat appears to have only %d windows. ' ...
         'Expected many short windows for 10 ms plotting/movie.'], size(R_resp.timeWindows,1));

    [~, timeIdx10] = min(sum(abs(R_resp.timeWindows - P.targetWin10),2));
    win10 = R_resp.timeWindows(timeIdx10,:);
    assert((win10(2)-win10(1)) <= 25, ...
        'Selected 10 ms bin is [%d %d] ms (duration %.1f ms).', win10(1), win10(2), win10(2)-win10(1));

    [~, timeIdxPre] = min(sum(abs(R_resp.timeWindows - P.targetWinPre),2));
    winPre = R_resp.timeWindows(timeIdxPre,:);
    assert((winPre(2)-winPre(1)) <= 25, ...
        'Selected pre-stim bin is [%d %d] ms (duration %.1f ms).', winPre(1), winPre(2), winPre(2)-winPre(1));

    opts10 = struct('timeIdx',timeIdx10,'excludeOverlap',true,'verbose',false);
    optsPre = struct('timeIdx',timeIdxPre,'excludeOverlap',true,'verbose',false);
    OUT10 = attention_modulation_V1_3bin(R_resp, Tall_V1, SNRnorm, opts10);
    OUTpre = attention_modulation_V1_3bin(R_resp, Tall_V1, SNRnorm, optsPre);

    OUT10plot = OUT10;
    OUT10plot.pValueTD = OUT3.pValueTD;   % hard significance rule
    OUTprePlot = OUTpre;
    OUTprePlot.pValueTD = OUT3.pValueTD;  % hard significance rule
end

%% QC diagnostics
if RUN_QC
    d3 = OUT3.muT - OUT3.muD;
    d10 = OUT10.muT - OUT10.muD;
    dpre = OUTpre.muT - OUTpre.muD;

    ok10 = isfinite(d3) & isfinite(d10);
    if any(ok10)
        r10 = corr(d3(ok10), d10(ok10));
        dd10 = abs(d3(ok10) - d10(ok10));
        fprintf('Delta(T-D) 3-bin vs 10-ms: corr=%.4f, median|diff|=%.6g, max|diff|=%.6g\n', ...
            r10, median(dd10), max(dd10));
    else
        fprintf('Delta(T-D) 3-bin vs 10-ms unavailable (no finite overlap).\n');
    end

    okPre = isfinite(d3) & isfinite(dpre);
    if any(okPre)
        rPre = corr(d3(okPre), dpre(okPre));
        ddPre = abs(d3(okPre) - dpre(okPre));
        fprintf('Delta(T-D) 3-bin vs pre-stim: corr=%.4f, median|diff|=%.6g, max|diff|=%.6g\n', ...
            rPre, median(ddPre), max(ddPre));
    end

    fprintf('nSig pre-stim (pTD < %.2f): %d / %d sites\n', P.pThresh, ...
        nnz(isfinite(OUTpre.pValueTD) & OUTpre.pValueTD < P.pThresh), numel(OUTpre.pValueTD));

    sigPre = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh) & isfinite(dpre);
    nPrePos = nnz(sigPre & (dpre > 0));   % T > D
    nPreNeg = nnz(sigPre & (dpre < 0));   % D > T
    nPreTot = nPrePos + nPreNeg;
    if nPreTot > 0
        fracPos = nPrePos / nPreTot;
        fracNeg = nPreNeg / nPreTot;
        fprintf(['Pre-stim sign balance (sig sites): T>D: %d, D>T: %d, ' ...
                 'frac(T>D)=%.3f, frac(D>T)=%.3f\n'], nPrePos, nPreNeg, fracPos, fracNeg);
        if abs(fracPos - 0.5) > 0.15
            warning(['Pre-stim sign balance is noticeably asymmetric. ' ...
                     'Check normalization/trial selection if this persists.']);
        end
    else
        fprintf('Pre-stim sign balance (sig sites): no non-zero d(T-D) values.\n');
    end
end

%% Shared visualization scale (across 10-ms and pre-stim still frames)
if RUN_STILLS
    d10 = OUT10.muT - OUT10.muD;
    dpre = OUTpre.muT - OUTpre.muD;
    sigMask = isfinite(OUT3.pValueTD) & (OUT3.pValueTD < P.pThresh);
    poolVals = [abs(d10(sigMask)); abs(dpre(sigMask))];
    poolVals = poolVals(isfinite(poolVals) & poolVals > 0);

    if isempty(poolVals)
        cShared = 1;
    else
        cShared = prctile(poolVals, P.plotRobustPct);
        if ~isfinite(cShared) || cShared <= 0
            cShared = max(poolVals);
        end
        if ~isfinite(cShared) || cShared <= 0
            cShared = 1;
        end
    end
    fprintf('Using shared color scale cMax=%.6g for still frames.\n', cShared);

    optsPlot = struct();
    optsPlot.RTAB384 = RTAB384;
    optsPlot.pThresh = P.pThresh;
    optsPlot.siteRange = P.siteRange;
    optsPlot.markerSize = P.plotMarkerSize;
    optsPlot.alpha = P.plotAlpha;
    optsPlot.robustPct = P.plotRobustPct;
    optsPlot.cMaxFixed = cShared;
    optsPlot.bgColor = P.plotBgColor;
    optsPlot.cLow = P.plotCLow;
    optsPlot.cHigh = P.plotCHigh;

    % ===== SMALL WINDOW FRAME PLOT (~10 ms) =====
    hSmall = plot_projected_attentiondiff_on_example_stim( ...
        P.stimIDExample, OUT10plot, Tall_V1, ALLCOORDS, optsPlot);
    figure(hSmall.fig);
    set(hSmall.fig, 'Name', 'Attention TD | Small window (~10 ms)', 'NumberTitle', 'off');
    title(hSmall.ax, sprintf('Small-window attention map | [%d %d] ms | stim %d', ...
        win10(1), win10(2), P.stimIDExample), 'Color','w');
    fprintf('Using 10 ms bin %d: [%d %d] ms\n', timeIdx10, win10(1), win10(2));

    % ===== PRE-STIM (SPONTANEOUS) FRAME PLOT =====
    hPre = plot_projected_attentiondiff_on_example_stim( ...
        P.stimIDExample, OUTprePlot, Tall_V1, ALLCOORDS, optsPlot);
    figure(hPre.fig);
    set(hPre.fig, 'Name', 'Attention TD | Pre-stim control', 'NumberTitle', 'off');
    title(hPre.ax, sprintf('Pre-stim control map | [%d %d] ms | stim %d', ...
        winPre(1), winPre(2), P.stimIDExample), 'Color','w');
    fprintf('Using pre-stim bin %d: [%d %d] ms\n', timeIdxPre, winPre(1), winPre(2));
    fprintf('Pre-stim plotted contributions: target=%d, distractor=%d\n', ...
        hPre.nTargetContrib, hPre.nDistrContrib);
end

%% Optional movie rendering
if RUN_MOVIE
    outMovie = fullfile(cfg.resultsDir, 'V1_attentiondiff_movie_10ms.mp4');
    make_attention_movie_wrapper_safe( ...
        outMovie, Tall_V1, ALLCOORDS, RTAB384, P.stimIDExample, R_resp, SNRnorm, OUT3);
end
