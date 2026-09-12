% Attention_Modulation_Spread_Activity_Inclusion_V1_V4
% V1/V4 attention d-prime using the visual-activity site-inclusion mask.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'analyses')));
addpath(genpath(fullfile(repoRoot, 'core')));
addpath(genpath(fullfile(repoRoot, 'utils')));
cfg = config();

P = struct();
P.timeIdx = 3;
P.windowMs = [300 500];
P.exampleStimulus = 38;
P.snrThreshold = 0.7;
P.pTDThreshold = 0.05;
P.minMatchedTrials = 20;
P.pColorThreshold = 0.05;
P.smoothSigmaPx = 15;
P.smoothSupportFraction = 0.01;
P.framePaddingPx = 100;
P.clipRange = [-0.5 0.5];
P.alphaGain = 2.5;
P.contourLineWidth = 4;

load(fullfile(cfg.logsDir, 'ObjAtt_lines_monkeyN_20220201_B1.mat')); % ALLCOORDS
load(fullfile(cfg.logsDir, 'RTAB384.mat'));                         % RTAB384
S = load(fullfile(cfg.matDir, 'SNR_capsules_N_d12.mat'));
R3 = S.R;
assert(isequal(double(R3.timeWindows(P.timeIdx,:)), P.windowMs), ...
    'Expected the selected response bin to be 300-500 ms.');

%% V1: reproduce the visual-response inclusion mask
S = load(fullfile(cfg.matDir, 'Tall_V1_lines_N.mat'));
TallV1 = S.Tall_V1;
S = load(fullfile(cfg.matDir, 'SNR_V1_byColor_byWindow.mat'));
SNRv1 = S.SNR;
S = load(fullfile(cfg.matDir, 'ColorTune_balanced_V1.mat'));
ColorTuneV1 = S.ColorTune;

R3V1 = subset_response(R3, 1:512);
optsTD = struct('v1Sites', 1:512, 'timeIdx', P.timeIdx, ...
    'excludeOverlap', true, 'verbose', false);
OUTv1 = attention_modulation_V1_3bin(R3V1, TallV1, SNRv1, optsTD);
[keepV1, maskStatsV1] = activity_inclusion_mask( ...
    SNRv1, OUTv1, ColorTuneV1, P);
print_selection('V1', keepV1, maskStatsV1, 512);

hV1 = plot_projected_attention_dprime_on_example_stim( ...
    TallV1, ALLCOORDS, RTAB384, P.exampleStimulus, R3V1, SNRv1, ...
    'TimeBin', P.timeIdx, 'SiteIdx', keepV1, ...
    'CoverageSiteIdx', (1:512).', 'ExcludeOverlap', true, ...
    'SmoothSigmaPx', P.smoothSigmaPx, ...
    'SmoothSupportFraction', P.smoothSupportFraction, ...
    'FramePaddingPx', P.framePaddingPx, 'ClipRange', P.clipRange, ...
    'AlphaGain', P.alphaGain, 'ContourLineWidth', P.contourLineWidth);
add_panel_labels(hV1, 'V1', P.windowMs);

%% V4: reproduce the visual-response inclusion mask
S = load(fullfile(cfg.matDir, 'Tall_V4_lines_N.mat'));
TallV4 = S.Tall_V4;
S = load(fullfile(cfg.matDir, 'ColorTune_balanced_V4_N.mat'));
ColorTuneV4 = S.ColorTune;
globalV4 = (513:768).';
assert(isequal(double(ColorTuneV4.RFrange(:)), globalV4), ...
    'ColorTune V4 channel mapping does not match global channels 513:768.');

SNRv4 = compute_snr_per_color_region(R3, TallV4, globalV4);
R3V4 = subset_response(R3, globalV4);
optsTD = struct('v1Sites', 1:256, 'timeIdx', P.timeIdx, ...
    'excludeOverlap', true, 'verbose', false);
OUTv4 = attention_modulation_V1_3bin(R3V4, TallV4, SNRv4, optsTD);
[keepV4, maskStatsV4] = activity_inclusion_mask( ...
    SNRv4, OUTv4, ColorTuneV4, P);
print_selection('V4', keepV4, maskStatsV4, 256);

hV4 = plot_projected_attention_dprime_on_example_stim( ...
    TallV4, ALLCOORDS, RTAB384, P.exampleStimulus, R3V4, SNRv4, ...
    'TimeBin', P.timeIdx, 'SiteIdx', keepV4, ...
    'CoverageSiteIdx', (1:256).', 'ExcludeOverlap', true, ...
    'SmoothSigmaPx', P.smoothSigmaPx, ...
    'SmoothSupportFraction', P.smoothSupportFraction, ...
    'FramePaddingPx', P.framePaddingPx, 'ClipRange', P.clipRange, ...
    'AlphaGain', P.alphaGain, 'ContourLineWidth', P.contourLineWidth);
add_panel_labels(hV4, 'V4', P.windowMs);

%% Assemble and save a distinct comparison figure
figCombined = figure('Color','w', ...
    'Name','Nilson attention d-prime | activity-selected sites', ...
    'NumberTitle','off','InvertHardcopy','off', ...
    'Position',[50 100 1600 650]);
axV1 = axes('Parent',figCombined,'Position',[0.01 0.02 0.485 0.96]);
copy_panel(hV1.ax,axV1);
axV4 = axes('Parent',figCombined,'Position',[0.505 0.02 0.485 0.96]);
copy_panel(hV4.ax,axV4);

scaleAx = axes('Parent',figCombined,'Position',[0.31 0.905 0.17 0.025]);
image(scaleAx,[hV1.colorScaleValues(1) hV1.colorScaleValues(end)], ...
    [0 1],hV1.colorScaleRGB);
set(scaleAx,'Box','on','Color','w','FontName','Helvetica','FontSize',9, ...
    'Layer','top','TickDir','out', ...
    'XLim',[hV1.colorScaleValues(1) hV1.colorScaleValues(end)], ...
    'XTick',[-0.5 -0.25 0 0.25 0.5], ...
    'YDir','normal','YLim',[0 1],'YTick',[]);
xlabel(scaleAx,'Attentional modulation d'' (T-D)', ...
    'FontName','Helvetica','FontSize',9);

close(hV1.fig);
close(hV4.fig);
drawnow;

hAttentionActivityInclusion = struct('fig',figCombined,'axV1',axV1, ...
    'axV4',axV4,'scaleAx',scaleAx,'keepV1',keepV1,'keepV4',keepV4, ...
    'maskStatsV1',maskStatsV1,'maskStatsV4',maskStatsV4, ...
    'nPointsV1',hV1.nRasterPoints,'nPointsV4',hV4.nRasterPoints, ...
    'nBackgroundV1',hV1.nBackgroundSamples, ...
    'nBackgroundV4',hV4.nBackgroundSamples, ...
    'fieldRangeV1',hV1.fieldRange,'fieldRangeV4',hV4.fieldRange, ...
    'fractionClippedV1',hV1.fractionClipped, ...
    'fractionClippedV4',hV4.fractionClipped);

outDir = fullfile(cfg.resultsDir,'spread');
if exist(outDir,'dir')~=7, mkdir(outDir); end
stem = 'V1_V4_attentional_modulation_dprime_activity_inclusion_300_500ms';
figPath = fullfile(outDir,[stem '.fig']);
pngPath = fullfile(outDir,[stem '.png']);
savefig(figCombined,figPath);
print(figCombined,pngPath,'-dpng','-r300');
fprintf('Saved activity-inclusion attention figure:\n  %s\n  %s\n', ...
    figPath,pngPath);

function Rlocal = subset_response(R,siteIdx)
Rlocal = R;
Rlocal.meanAct = R.meanAct(siteIdx,:,:);
Rlocal.meanSqAct = R.meanSqAct(siteIdx,:,:);
if ~isvector(R.nTrials)
    Rlocal.nTrials = R.nTrials(siteIdx,:);
end
end

function [keepSites,stats] = activity_inclusion_mask(SNR,OUT,ColorTune,P)
snrMatrix = [SNR.yellowEarly,SNR.yellowLate, ...
    SNR.purpleEarly,SNR.purpleLate];
bestSNR = max(snrMatrix,[],2,'omitnan');
matchedN = OUT.wY+OUT.wP;
isMain = isfinite(bestSNR) & bestSNR>P.snrThreshold;
isRescueTD = isfinite(OUT.pValueTD) & OUT.pValueTD<P.pTDThreshold & ...
    matchedN>=P.minMatchedTrials;
colorP = double(ColorTune.early.p(:));
isRescueColor = isfinite(colorP) & colorP<P.pColorThreshold;
keepSites = find(isMain | isRescueTD | isRescueColor);
stats = struct('nMain',nnz(isMain), ...
    'nTDRescueOnly',nnz(isRescueTD & ~isMain), ...
    'nColorRescueOnly',nnz(isRescueColor & ~isMain & ~isRescueTD));
end

function print_selection(region,keepSites,stats,nTotal)
fprintf(['Activity-inclusion %s sites: %d / %d ' ...
    '(SNR=%d, TD rescue-only=%d, color rescue-only=%d).\n'], ...
    region,numel(keepSites),nTotal,stats.nMain, ...
    stats.nTDRescueOnly,stats.nColorRescueOnly);
end

function add_panel_labels(h,region,windowMs)
label = sprintf('%s | T-D d'' | activity-selected | %.0f-%.0f ms', ...
    region,windowMs);
text(h.ax,h.xLimits(1)+14,h.yLimits(1)+14,label, ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Helvetica','FontSize',14,'FontWeight','bold', ...
    'Color',[0.1 0.1 0.1]);
text(h.ax,h.xLimits(2)-14,h.yLimits(2)-14, ...
    sprintf('N = %d projected points',h.nRasterPoints), ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontName','Helvetica','FontSize',12,'Color',[0.25 0.25 0.25]);
end

function copy_panel(sourceAx,targetAx)
copyobj(allchild(sourceAx),targetAx);
set(targetAx,'Color','w','XLim',sourceAx.XLim,'YLim',sourceAx.YLim, ...
    'YDir',sourceAx.YDir,'DataAspectRatio',sourceAx.DataAspectRatio, ...
    'Visible','off');
axis(targetAx,'image');
axis(targetAx,'off');
end
