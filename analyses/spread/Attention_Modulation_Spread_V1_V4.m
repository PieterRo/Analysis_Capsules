% Attention_Modulation_Spread_V1_V4
% Create one V1/V4 figure of spatial attentional modulation d-prime.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'analyses')));
addpath(genpath(fullfile(repoRoot, 'core')));
addpath(genpath(fullfile(repoRoot, 'utils')));
cfg = config();

run(fullfile(scriptDir, 'Attention_Modulation_Spread_V1.m'));
hV1 = hAttentionSpread;
v1Stats = attention_stats(hV1);

run(fullfile(scriptDir, 'Attention_Modulation_Spread_V4.m'));
hV4 = hAttentionSpread;
v4Stats = attention_stats(hV4);

figCombined = figure('Color','w', ...
    'Name','Nilson attentional-modulation spread | V1 and V4', ...
    'NumberTitle','off','InvertHardcopy','off', ...
    'Position',[50 100 1600 650]);
axV1 = axes('Parent',figCombined,'Position',[0.01 0.02 0.485 0.96]);
copy_panel(hV1.ax,axV1);
axV4 = axes('Parent',figCombined,'Position',[0.505 0.02 0.485 0.96]);
copy_panel(hV4.ax,axV4);

scaleAx = axes('Parent',figCombined,'Position',[0.31 0.905 0.17 0.025]);
image(scaleAx,[v1Stats.colorScaleValues(1) v1Stats.colorScaleValues(end)], ...
    [0 1],v1Stats.colorScaleRGB);
set(scaleAx,'Box','on','Color','w','FontName','Helvetica','FontSize',9, ...
    'Layer','top','TickDir','out', ...
    'XLim',[v1Stats.colorScaleValues(1) v1Stats.colorScaleValues(end)], ...
    'XTick',[-0.5 -0.25 0 0.25 0.5], ...
    'YDir','normal','YLim',[0 1],'YTick',[]);
xlabel(scaleAx,'Attentional modulation d'' (T-D)', ...
    'FontName','Helvetica','FontSize',9);

close(hV1.fig);
close(hV4.fig);
drawnow;

hAttentionSpreadCombined = struct('fig',figCombined,'axV1',axV1, ...
    'axV4',axV4,'scaleAx',scaleAx,'v1',v1Stats,'v4',v4Stats);

outDir = fullfile(cfg.resultsDir,'spread');
if exist(outDir,'dir')~=7, mkdir(outDir); end
stem = 'V1_V4_attentional_modulation_dprime_300_500ms';
figPath = fullfile(outDir,[stem '.fig']);
pngPath = fullfile(outDir,[stem '.png']);
savefig(figCombined,figPath);
print(figCombined,pngPath,'-dpng','-r300');
fprintf('Saved combined attention-spread figure:\n  %s\n  %s\n', ...
    figPath,pngPath);

function stats = attention_stats(h)
stats = struct('nSites',h.nSites,'nSamples',h.nSamples, ...
    'nTargetSamples',h.nTargetSamples, ...
    'nDistractorSamples',h.nDistractorSamples, ...
    'nBackgroundSamples',h.nBackgroundSamples, ...
    'nRasterPoints',h.nRasterPoints, ...
    'nCoveragePoints',h.nCoveragePoints, ...
    'fieldRange',h.fieldRange,'fractionClipped',h.fractionClipped, ...
    'clipRange',h.clipRange,'colorScaleValues',h.colorScaleValues, ...
    'colorScaleRGB',h.colorScaleRGB,'timeWindow',h.timeWindow);
end

function copy_panel(sourceAx,targetAx)
copyobj(allchild(sourceAx),targetAx);
set(targetAx,'Color','w','XLim',sourceAx.XLim,'YLim',sourceAx.YLim, ...
    'YDir',sourceAx.YDir,'DataAspectRatio',sourceAx.DataAspectRatio, ...
    'Visible','off');
axis(targetAx,'image');
axis(targetAx,'off');
end
