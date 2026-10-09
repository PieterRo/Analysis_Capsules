function build_supp_fig1_panel_a_peaknorm(resultsDir,figureDir)
cfg = config();
if nargin<1 || isempty(resultsDir); resultsDir = cfg.resultsDir; end
if nargin<2 || isempty(figureDir)
    figureDir = fullfile(fileparts(cfg.resultsRoot),'Figures');
end
if ~exist(figureDir,'dir'); mkdir(figureDir); end
monkeys = {'N','F'};
areas = {'V1','V4','IT'};
areaColor = [0.00 0.61 0.88; 0.98 0.56 0.00; 0.93 0.06 0.12];
data = cell(2,3);
latency = nan(2,3);
normalizationPeak = nan(2,3);
for m=1:2
    for a=1:3
        s = load(fullfile(resultsDir,sprintf( ...
            'SuppFig1_color_preferred_sites_%s_%s.mat',monkeys{m},areas{a})),'result');
        z = s.result;
        smoothPreferred = movmean(double(z.muPreferred(:)'),3);
        postStimulus = double(z.timeMs(:)') >= 0 & ...
            double(z.timeMs(:)') <= 500;
        normalizationPeak(m,a) = max(smoothPreferred(postStimulus));
        assert(isfinite(normalizationPeak(m,a)) && ...
            normalizationPeak(m,a)>0);
        z.muPreferred = double(z.muPreferred)/normalizationPeak(m,a);
        z.muOther = double(z.muOther)/normalizationPeak(m,a);
        data{m,a} = z;
        x = load(fullfile(resultsDir,sprintf( ...
            'Color_decoder_latency_fit_%s_%s_2000_-100to300.mat', ...
            monkeys{m},areas{a})),'R');
        latency(m,a) = x.R.t50;
        assert(isequal(s.result.timeMs,x.R.timeMs));
    end
end

values = [];
for m=1:2
    for a=1:3
        z = data{m,a};
        values = [values z.muPreferred(:)' z.muOther(:)']; %#ok<AGROW>
    end
end
low = min(values); high = max(values);
span = max(high-low,0.1);
yLim = repmat([min(low-0.12*span,-0.10*high), ...
    high+0.13*span],3,1);

fig = figure('Visible','off','Color','white', ...
    'Units','pixels','Position',[100 100 1190 1684]);
set(fig,'PaperUnits','points','PaperSize',[595 842], ...
    'PaperPosition',[0 0 595 842], ...
    'InvertHardcopy','off','Renderer','painters');
axesX = [0.12 0.41 0.70];
axesY = [0.82 0.65];
axesW = 0.225; axesH = 0.105;
for m=1:2
    for a=1:3
        z = data{m,a};
        ax = axes('Parent',fig,'Position',[axesX(a) axesY(m) axesW axesH]);
        hold(ax,'on');
        t = double(z.timeMs(:)');
        preferred = movmean(double(z.muPreferred(:)'),3);
        other = movmean(double(z.muOther(:)'),3);
        plot(ax,t,preferred,'-','Color',areaColor(a,:), ...
            'LineWidth',1.4);
        plot(ax,t,other,'k--','LineWidth',1.1);
        plot(ax,[-200 500],[0 0],'k-','LineWidth',0.5);
        plot(ax,[0 0],yLim(a,:),'--','Color',[0.45 0.45 0.45], ...
            'LineWidth',0.55);
        xlim(ax,[-200 500]); ylim(ax,yLim(a,:));
        set(ax,'XTick',[],'YTick',[],'Box','off', ...
            'XColor','none','YColor','none');
        text(ax,0.02,0.86,sprintf('N = %d',z.nSignificantSites), ...
            'Units','normalized','HorizontalAlignment','left', ...
            'FontSize',7,'Color',[0.1 0.1 0.1]);
        if a==1
            text(ax,-0.28,0.5,sprintf('Monkey %s',monkeys{m}), ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontSize',7,'Rotation',90,'Clipping','off');
        end
        arrowX = axesX(a)+axesW*(latency(m,a)+200)/700;
        arrowY = axesY(m);
        annotation(fig,'arrow',[arrowX arrowX], ...
            [arrowY-0.023 arrowY-0.006], ...
            'Color',areaColor(a,:),'LineWidth',1.0, ...
            'HeadLength',6,'HeadWidth',5);
        annotation(fig,'textbox', ...
            [arrowX-0.047 arrowY-0.046 0.094 0.019], ...
            'String',sprintf('%.0f ms',latency(m,a)), ...
            'LineStyle','none','HorizontalAlignment','center', ...
            'FontSize',7,'Color',areaColor(a,:));
    end
end
for a=1:3
    annotation(fig,'textbox',[axesX(a) 0.938 axesW 0.020], ...
        'String',areas{a},'LineStyle','none', ...
        'HorizontalAlignment','center','FontSize',10, ...
        'Color',areaColor(a,:));
end
annotation(fig,'textbox',[0.037 0.930 0.025 0.028], ...
    'String','a','LineStyle','none','FontSize',10,'FontWeight','bold');

% Small vector swatches keep the color convention explicit without a box.
annotation(fig,'line',[0.12 0.15],[0.968 0.968], ...
    'Color',areaColor(1,:),'LineWidth',1.3);
annotation(fig,'textbox',[0.152 0.957 0.11 0.020], ...
    'String','Preferred color','LineStyle','none','FontSize',7);
annotation(fig,'line',[0.285 0.315],[0.968 0.968], ...
    'Color',[0 0 0],'LineWidth',1.1,'LineStyle','--');
annotation(fig,'textbox',[0.317 0.957 0.10 0.020], ...
    'String','Other color','LineStyle','none','FontSize',7);

svg = fullfile(figureDir,'Supp Fig 1_panel_a_peaknorm.svg');
pdf = fullfile(figureDir,'Supp Fig 1_panel_a_peaknorm.pdf');
png = fullfile(figureDir,'Supp Fig 1_panel_a_peaknorm.png');
print(fig,svg,'-dsvg','-painters');
print(fig,pdf,'-dpdf','-painters');
print(fig,png,'-dpng','-r150');
close(fig);
fprintf('Saved %s\nSaved %s\nSaved %s\n',svg,pdf,png);
for m=1:2
    for a=1:3
        z = data{m,a};
        fprintf('%s %s: N=%d/%d, window %s, decoder t50=%.1f ms, smoothed peak=%.3f spontaneous SD\n', ...
            monkeys{m},areas{a},z.nSignificantSites,z.nCandidateSites, ...
            mat2str(z.selectionWindowMs),latency(m,a),normalizationPeak(m,a));
    end
end
end
