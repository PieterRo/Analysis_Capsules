% Visual_Response_Spread_V1
% Plot the spatial distribution of the mean V1 visual response over a
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
P.snrThreshold = 0.7;
P.pTDThreshold = 0.05;
P.minMatchedTrials = 20;
P.pColorThreshold = 0.05;
P.markerSize = 5;
P.onlyOnObjects = false;
P.compareBroadBin = true;
P.saveFigure = false;

%% Load geometry and response data
load(fullfile(cfg.logsDir, 'ObjAtt_lines_monkeyN_20220201_B1.mat')); % ALLCOORDS
load(fullfile(cfg.logsDir, 'RTAB384.mat'));                         % RTAB384

S = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));
Tall_V1 = S.Tall_V1;

S = load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat'));
SNR = S.SNR;

S = load(fullfile(cfg.matDir, 'Resp_capsules_N_d12.mat'));
R_resp = S.R;

S = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'));
R3 = S.R;

S = load(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'));
ColorTune = S.ColorTune;

%% Apply the established V1 site-inclusion rule
SNRmat = [SNR.yellowEarly, SNR.yellowLate, ...
          SNR.purpleEarly, SNR.purpleLate];
bestSNR = max(SNRmat, [], 2, 'omitnan');

optsTD = struct('timeIdx', 3, 'excludeOverlap', true, 'verbose', false);
OUTtd = attention_modulation_V1_3bin(R3, Tall_V1, SNR, optsTD);
matchedN = OUTtd.wY + OUTtd.wP;

isMain = isfinite(bestSNR) & (bestSNR > P.snrThreshold);
isRescueTD = isfinite(OUTtd.pValueTD) & ...
    (OUTtd.pValueTD < P.pTDThreshold) & ...
    (matchedN >= P.minMatchedTrials);
colorP = ColorTune.early.p(1:512);
isRescueColor = isfinite(colorP(:)) & (colorP(:) < P.pColorThreshold);

keepSites = find(isMain | isRescueTD | isRescueColor);

fprintf('SNR-selected sites: %d / 512\n', nnz(isMain));
fprintf('TD rescue-only sites: %d\n', nnz(isRescueTD & ~isMain));
fprintf('Color rescue-only sites: %d\n', ...
    nnz(isRescueColor & ~isMain & ~isRescueTD));
fprintf('Total included V1 sites: %d / 512\n', numel(keepSites));

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
R_spread.meanAct = mean(R_resp.meanAct(:,:,binMask), 3, 'omitnan');
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
        dBroad = double(R_spread.meanAct(1:512,:)) - ...
            double(R3.meanAct(1:512,:,broadIdx));
        fprintf(['Short-bin mean versus broad bin: median |difference| = %.6g, ' ...
                 'max |difference| = %.6g.\n'], ...
            median(abs(dBroad(:)), 'omitnan'), ...
            max(abs(dBroad(:)), [], 'omitnan'));
    end
end

%% Plot the spatial spread of the mean visual response
hSpread = plot_projected_activity_on_example_stim( ...
    Tall_V1, ALLCOORDS, RTAB384, P.exampleStimulus, R_spread, SNR, ...
    'TimeBin', 1, ...
    'UseOnlyV1', true, ...
    'OnlyOnObjects', P.onlyOnObjects, ...
    'MarkerSize', P.markerSize, ...
    'SiteIdx', keepSites);

windowLabel = sprintf('%.0f-%.0f ms', P.windowMs(1), P.windowMs(2));
set(hSpread.fig, 'Name', ...
    sprintf('V1 mean visual response | %s', windowLabel), ...
    'NumberTitle', 'on', ...
    'InvertHardcopy', 'off');
text(hSpread.ax, 14, 14, windowLabel, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'Color', 'w');

%% Optional figure export to the Dropbox results tree
if P.saveFigure
    spreadResultsDir = fullfile(cfg.resultsDir, 'spread');
    if exist(spreadResultsDir, 'dir') ~= 7
        mkdir(spreadResultsDir);
    end
    fileStem = sprintf('V1_visual_response_spread_%g_%gms', P.windowMs);
    savefig(hSpread.fig, fullfile(spreadResultsDir, [fileStem '.fig']));
    print(hSpread.fig, fullfile(spreadResultsDir, [fileStem '.png']), ...
        '-dpng', '-r300');
    fprintf('Saved spread figure to: %s\n', spreadResultsDir);
end
