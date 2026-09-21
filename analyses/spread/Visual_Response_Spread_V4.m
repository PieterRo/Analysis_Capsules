% Visual_Response_Spread_V4
% Plot the spatial distribution of the mean V4 visual response over a
% broad time window by averaging the corresponding short response bins.

%% Paths and parameters
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'analyses')));
addpath(genpath(fullfile(repoRoot, 'core')));
addpath(genpath(fullfile(repoRoot, 'utils')));

cfg = config();

P = struct();
P.windowMs = [300 500];
P.exampleStimulus = 38;
P.globalSites = (513:768).';
P.snrThreshold = 0.7;
P.pTDThreshold = 0.05;
P.minMatchedTrials = 20;
P.pColorThreshold = 0.05;
P.markerSize = 5;
P.onlyOnObjects = false;
P.compareBroadBin = true;
P.saveFigure = false;
P.plotMode = 'smoothed';
P.smoothSigmaPx = 15;
P.smoothSupportFraction = 0.01;
P.smoothActivityAlphaGain = 3;
P.smoothFramePaddingPx = 100;
P.contourLineWidth = 4;
P.coverageSiteIdx = (1:256).';

%% Load geometry and response data
load(fullfile(cfg.logsDir, 'ObjAtt_lines_monkeyN_20220201_B1.mat')); % ALLCOORDS
load(fullfile(cfg.logsDir, 'RTAB384.mat'));                         % RTAB384

S = load(fullfile(cfg.matDir, 'Tall_V4_lines_N.mat'));
Tall_V4 = S.Tall_V4;

S = load(fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat'));
R_resp = S.R;

S = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'));
R3 = S.R;

S = load(fullfile(cfg.matDir, 'ColorTune_balanced_V4_N.mat'));
ColorTune = S.ColorTune;

assert(numel(P.globalSites) == 256, 'Expected 256 Nilson V4 channels.');
assert(height(Tall_V4(1).T) == numel(P.globalSites), ...
    'Tall_V4 rows must correspond to global channels 513:768.');
assert(isequal(double(ColorTune.RFrange(:)), P.globalSites), ...
    'ColorTune V4 channel mapping does not match global channels 513:768.');

%% Compute V4 normalization using the established V1 definition
SNR_V4 = compute_snr_per_color_region(R3, Tall_V4, P.globalSites);
SNRmat = [SNR_V4.yellowEarly, SNR_V4.yellowLate, ...
          SNR_V4.purpleEarly, SNR_V4.purpleLate];
bestSNR = max(SNRmat, [], 2, 'omitnan');

validStoredSNR = isfinite(bestSNR) & isfinite(ColorTune.bestSNR(:));
if any(validStoredSNR)
    snrDifference = abs(bestSNR(validStoredSNR) - ...
        double(ColorTune.bestSNR(validStoredSNR)));
    fprintf(['Recomputed versus stored V4 bestSNR: median |difference| = %.6g, ' ...
             'max |difference| = %.6g.\n'], ...
        median(snrDifference), max(snrDifference));
end

%% Apply the same site-inclusion rule as the V1 spread analysis
R3_V4 = R3;
R3_V4.meanAct = R3.meanAct(P.globalSites,:,:);
R3_V4.meanSqAct = R3.meanSqAct(P.globalSites,:,:);
if ~isvector(R3.nTrials)
    R3_V4.nTrials = R3.nTrials(P.globalSites,:);
end

optsTD = struct('v1Sites', 1:256, 'timeIdx', 3, ...
    'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R3_V4, Tall_V4, SNR_V4, optsTD);
matchedN = OUTtd.wY + OUTtd.wP;

isMain = isfinite(bestSNR) & (bestSNR > P.snrThreshold);
isRescueTD = isfinite(OUTtd.pValueTD) & ...
    (OUTtd.pValueTD < P.pTDThreshold) & ...
    (matchedN >= P.minMatchedTrials);
colorP = ColorTune.early.p(:);
isRescueColor = isfinite(colorP) & (colorP < P.pColorThreshold);
keepSites = find(isMain | isRescueTD | isRescueColor);

fprintf('SNR-selected V4 sites: %d / 256\n', nnz(isMain));
fprintf('TD rescue-only V4 sites: %d\n', nnz(isRescueTD & ~isMain));
fprintf('Color rescue-only V4 sites: %d\n', ...
    nnz(isRescueColor & ~isMain & ~isRescueTD));
fprintf('Total included V4 sites: %d / 256\n', numel(keepSites));

%% Average the short movie bins within the requested window
binMask = R_resp.timeWindows(:,1) >= P.windowMs(1) & ...
          R_resp.timeWindows(:,2) <= P.windowMs(2);
assert(any(binMask), ...
    'No response bins found within [%g %g] ms.', P.windowMs);

selectedWindows = double(R_resp.timeWindows(binMask,:));
binWidths = selectedWindows(:,2) - selectedWindows(:,1);
assert(all(abs(binWidths - binWidths(1)) < 1e-9), ...
    'Selected response bins have unequal widths; use a weighted mean.');
assert(abs(selectedWindows(1,1) - P.windowMs(1)) < 1e-9 && ...
       abs(selectedWindows(end,2) - P.windowMs(2)) < 1e-9, ...
    'Selected bins do not fully cover the requested time window.');

R_spread = struct();
R_spread.meanAct = mean( ...
    R_resp.meanAct(P.globalSites,:,binMask), 3, 'omitnan');
R_spread.timeWindows = P.windowMs;

fprintf('Averaging %d bins covering %.0f-%.0f ms.\n', ...
    nnz(binMask), selectedWindows(1,1), selectedWindows(end,2));

%% Optional check against the independently computed broad response bin
if P.compareBroadBin
    broadIdx = find(all(abs(double(R3.timeWindows) - P.windowMs) < 1e-9, 2), ...
        1, 'first');
    if isempty(broadIdx)
        warning('No matching broad [%g %g] ms bin found in R3.', P.windowMs);
    else
        dBroad = double(R_spread.meanAct) - ...
            double(R3.meanAct(P.globalSites,:,broadIdx));
        fprintf(['Short-bin mean versus broad bin: median |difference| = %.6g, ' ...
                 'max |difference| = %.6g.\n'], ...
            median(abs(dBroad(:)), 'omitnan'), ...
            max(abs(dBroad(:)), [], 'omitnan'));
    end
end

%% Plot using local V4 rows mapped to global response channels above
hSpread = plot_projected_activity_on_example_stim( ...
    Tall_V4, ALLCOORDS, RTAB384, P.exampleStimulus, R_spread, SNR_V4, ...
    'TimeBin', 1, ...
    'UseOnlyV1', false, ...
    'OnlyOnObjects', P.onlyOnObjects, ...
    'MarkerSize', P.markerSize, ...
    'PlotMode', P.plotMode, ...
    'SmoothSigmaPx', P.smoothSigmaPx, ...
    'SmoothSupportFraction', P.smoothSupportFraction, ...
    'SmoothActivityAlphaGain', P.smoothActivityAlphaGain, ...
    'SmoothFramePaddingPx', P.smoothFramePaddingPx, ...
    'ContourLineWidth', P.contourLineWidth, ...
    'CoverageSiteIdx', P.coverageSiteIdx, ...
    'SiteIdx', keepSites);

windowLabel = sprintf('V4 | %.0f-%.0f ms', P.windowMs(1), P.windowMs(2));
set(hSpread.fig, 'Name', ...
    sprintf('V4 mean visual response | %.0f-%.0f ms', P.windowMs), ...
    'NumberTitle', 'on', ...
    'InvertHardcopy', 'off');
text(hSpread.ax, hSpread.xLimits(1)+14, hSpread.yLimits(1)+14, windowLabel, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'Color', [0.1 0.1 0.1]);
text(hSpread.ax, hSpread.xLimits(2)-14, hSpread.yLimits(2)-14, ...
    sprintf('N = %d projected points', hSpread.nRasterPoints), ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'FontName', 'Helvetica', ...
    'FontSize', 12, ...
    'Color', [0.25 0.25 0.25]);

%% Optional figure export to the Dropbox results tree
if P.saveFigure
    spreadResultsDir = fullfile(cfg.resultsDir, 'spread');
    if exist(spreadResultsDir, 'dir') ~= 7
        mkdir(spreadResultsDir);
    end
    fileStem = sprintf('V4_visual_response_spread_%g_%gms', P.windowMs);
    savefig(hSpread.fig, fullfile(spreadResultsDir, [fileStem '.fig']));
    print(hSpread.fig, fullfile(spreadResultsDir, [fileStem '.png']), ...
        '-dpng', '-r300');
    fprintf('Saved spread figure to: %s\n', spreadResultsDir);
end
