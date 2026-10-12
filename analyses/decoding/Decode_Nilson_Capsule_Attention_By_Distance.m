function OUT = Decode_Nilson_Capsule_Attention_By_Distance(varargin)
% Apply fixed capsule attention-decoder weights in six cue-distance groups.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
addpath(fullfile(repoRoot,'core','rf'));
addpath(fullfile(repoRoot,'core','metrics'));
addpath(fullfile(repoRoot,'core','snr'));
cfg = config();
p = inputParser;
addParameter(p,'OutputDir',cfg.resultsDir);
addParameter(p,'BootstrapIterations',300);
addParameter(p,'RandomSeed',270926);
addParameter(p,'ReplotFrom','');
addParameter(p,'ITRFRoleOnly',false,@(x) islogical(x) && isscalar(x));
parse(p,varargin{:});
opt = p.Results;
stem = 'Nilson_capsule_attention_decoder_six_distance_bins_N';
if opt.ITRFRoleOnly; stem = [stem '_itRfRoleOnly']; end
assert(opt.BootstrapIterations>=0 && ...
    opt.BootstrapIterations==floor(opt.BootstrapIterations));
if ~isempty(opt.ReplotFrom)
    A = load(opt.ReplotFrom,'OUT'); OUT = A.OUT;
    if ~isfolder(opt.OutputDir); mkdir(opt.OutputDir); end
    OUT.figureFile = fullfile(opt.OutputDir,[stem '.png']);
    makeFigure(OUT.result,OUT.timeMs,OUT.figureFile,opt.ITRFRoleOnly);
    return
end

lowerFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat');
lowerBaseFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');
itFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_elastic_net_CV_IT_combined_nearestCurve_ITfunctional_N.mat');
if opt.ITRFRoleOnly
    itFile = fullfile(cfg.resultsDir, ...
        ['Attention_decoder_elastic_net_CV_IT_combined_' ...
        'nearestCurve_ITfunctional_N_rfRoleOnly.mat']);
end
A = load(lowerFile,'OUT'); lower = A.OUT;
A = load(lowerBaseFile,'OUT'); lowerBase = A.OUT;
A = load(itFile,'OUT'); it = A.OUT;
m = matfile(fullfile(cfg.dataRoot,'Mr Nilson', ...
    'ObjAtt_lines_normMUA.mat'));
trialFile = fullfile(cfg.dataRoot,'Mr Nilson', ...
    'ObjAtt_lines_MUA_trials.mat');
M = load(trialFile,'tb');
sampleTime = double(M.tb(:).');
assert(isequal(sampleTime,-200:499));
binTime = sampleTime(1:10:end)+4.5;
baselineBins = binTime>=-150 & binTime<0;
lateBins = binTime>=300 & binTime<500;
names = {'V1','V4','IT'};
result = cell(3,1);
rows = cell(3,1);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng));
rng(opt.RandomSeed,'twister');
for area = 1:3
    name = names{area};
    if area==3
        R = it.IT;
        D = [];
    else
        R = lower.(name);
        D = lowerBase.(name);
        assert(isequal(R.siteIndexGlobal,D.siteIndexGlobal) && ...
            isequal(R.trialIndex,D.trialIndex));
    end
    T = load(fullfile(cfg.matDir,['Tall_' name '_lines_N.mat']), ...
        ['Tall_' name]);
    Tall = T.(['Tall_' name]);
    nSites = numel(R.siteIndexGlobal);
    localSite = double(R.siteIndexLocal(:));
    globalSite = double(R.siteIndexGlobal(:));
    weight = double(R.finalCoefficient(:));
    trial = double(R.trialIndex(:));
    quartet = double(R.quartet(:));
    label = double(R.classLabel(:));
    eligible = unique(quartet);
    nQ = numel(eligible);
    assert(numel(Tall)==384 && numel(weight)==nSites && ...
        numel(localSite)==nSites && numel(trial)==numel(quartet) && ...
        all(ismember(label,[-1 1])));
    if area==3
        S = load(fullfile(cfg.matDir,'SNR_capsules_N_d12.mat'),'R');
        snr = compute_snr_per_color_region( ...
            S.R,Tall,(769:1024).');
        responseScale = double(snr.sdSpont(localSite));
        assert(isequal(size(R.fullDataAlignmentByQuartetSite), ...
            [96 nSites]));
    else
        responseScale = double(D.responseScale(localSite));
    end
    assert(all(isfinite(responseScale) & responseScale>0));

    distance = nan(nSites,nQ);
    siteSign = zeros(nSites,nQ);
    for k = 1:nQ
        q = eligible(k);
        members = quartetMembers(q);
        first = Tall(members(1)).T;
        x = double(first.along_GC(localSite));
        overlap = logical(first.overlap(localSite));
        valid = isfinite(x) & ~overlap;
        for j = 2:4
            other = Tall(members(j)).T;
            xOther = double(other.along_GC(localSite));
            assert(all(isfinite(xOther(valid)) & ...
                ~other.overlap(localSite(valid)) & ...
                abs(xOther(valid)-x(valid))<1e-8), ...
                '%s quartet %d changed non-overlap RF distance.', ...
                name,q);
        end
        distance(valid,k) = x(valid);
        if area==3
            siteSign(:,k) = double( ...
                R.fullDataAlignmentByQuartetSite(q,:)).';
        else
            inQ = find(quartet==q);
            signRows = double(D.canonicalSign(inQ,:));
            assert(all(all(signRows==signRows(1,:))) && ...
                all(signRows(1,valid)~=0), ...
                '%s quartet %d has inconsistent decoder signs.',name,q);
            siteSign(:,k) = signRows(1,:).';
        end
    end
    valid = isfinite(distance);
    nValid = nnz(valid);
    [~,order] = sortrows([distance(valid) (1:nValid).'],[1 2]);
    group = zeros(nValid,1,'uint8');
    group(order) = uint8(floor((0:nValid-1)*6/nValid)+1);
    groupMap = zeros(nSites,nQ,'uint8');
    groupMap(valid) = group;
    groupCounts = accumarray(double(group),1,[6 1]);
    medianDistance = accumarray(double(group), ...
        distance(valid),[6 1],@median);
    assert(all(diff(medianDistance)>0));

    nTrials = numel(trial);
    nBins = numel(binTime);
    contribution = zeros(nSites,nQ,nBins,'single');
    positiveRows = cell(nQ,1);
    negativeRows = cell(nQ,1);
    for k = 1:nQ
        positiveRows{k} = find(quartet==eligible(k) & label==1);
        negativeRows{k} = find(quartet==eligible(k) & label==-1);
        assert(~isempty(positiveRows{k}) && ...
            ~isempty(negativeRows{k}));
    end
    for s = 1:nSites
        if abs(weight(s))<=1e-8; continue; end
        raw = double(squeeze(m.normMUA(globalSite(s),:,:)));
        raw = raw(trial,:);
        assert(isequal(size(raw),[nTrials 700]));
        binned = squeeze(mean(reshape(raw,nTrials,10,nBins), ...
            2,'omitnan'))/responseScale(s);
        for k = 1:nQ
            delta = mean(binned(positiveRows{k},:),1)- ...
                mean(binned(negativeRows{k},:),1);
            contribution(s,k,:) = single( ...
                weight(s)*siteSign(s,k)*delta);
        end
    end

    savedContrast = nan(nQ,1);
    reconstructedContrast = nan(nQ,1);
    for k = 1:nQ
        pos = positiveRows{k}; neg = negativeRows{k};
        savedContrast(k) = mean(R.finalClassScore(pos))- ...
            mean(R.finalClassScore(neg));
        x = squeeze(sum(double(contribution(:,k,lateBins)),1));
        reconstructedContrast(k) = mean(x(:));
    end
    decoderCheck = max(abs(savedContrast-reconstructedContrast));
    assert(decoderCheck<1e-4, ...
        '%s time courses do not reproduce saved decoder contrasts.',name);
    fprintf('%s decoder contrast check: %.3g\n',name,decoderCheck);

    qGroup = nan(nQ,6,nBins);
    candidateSites = zeros(nQ,6);
    weightedSites = zeros(nQ,6);
    for k = 1:nQ
        for g = 1:6
            keep = groupMap(:,k)==g;
            candidateSites(k,g) = nnz(keep);
            weightedSites(k,g) = nnz(keep & ...
                abs(weight.*siteSign(:,k))>1e-8);
            if ~any(keep); continue; end
            qGroup(k,g,:) = sum(double( ...
                contribution(keep,k,:)),1);
        end
    end
    curve = reshape(mean(qGroup,1,'omitnan'),6,nBins);
    curve = curve-mean(curve(:,baselineBins),2);
    fit = fitGroups(curve,binTime,medianDistance);
    assert(fit.exitflag>0 && ...
        all(isfinite(fit.summary.t50)));
    line = fitSlope(medianDistance,fit.summary.t50);

    bootT50 = nan(opt.BootstrapIterations,6);
    bootSlope = nan(opt.BootstrapIterations,1);
    for b = 1:opt.BootstrapIterations
        draw = randi(nQ,nQ,1);
        bootCurve = reshape(mean(qGroup(draw,:,:), ...
            1,'omitnan'),6,nBins);
        bootCurve = bootCurve-mean( ...
            bootCurve(:,baselineBins),2);
        try
            bootFit = fitGroups(bootCurve,binTime,medianDistance);
            if bootFit.exitflag>0 && ...
                    all(isfinite(bootFit.summary.t50))
                bootT50(b,:) = bootFit.summary.t50(:).';
                bootLine = fitSlope(medianDistance, ...
                    bootFit.summary.t50);
                bootSlope(b) = bootLine.slope;
            end
        catch
            % A nonconvergent resample contributes no timing estimate.
        end
    end
    t50CI = nan(6,2);
    for g = 1:6
        values = bootT50(isfinite(bootT50(:,g)),g);
        if ~isempty(values)
            t50CI(g,:) = prctile(values,[2.5 97.5]);
        end
    end
    validSlope = bootSlope(isfinite(bootSlope));
    slopeCI = [NaN NaN];
    if ~isempty(validSlope)
        slopeCI = prctile(validSlope,[2.5 97.5]);
    end
    result{area} = struct('area',string(name), ...
        'eligibleQuartets',eligible,'medianDistance',medianDistance, ...
        'groupSiteQuartetCounts',groupCounts, ...
        'groupBySiteQuartet',groupMap, ...
        'candidateSitesByQuartet',candidateSites, ...
        'weightedSitesByQuartet',weightedSites, ...
        'quartetDifference',qGroup,'curve',curve,'fit',fit, ...
        'line',line,'t50CI',t50CI,'slopeCI',slopeCI, ...
        'bootstrapSlope',bootSlope, ...
        'decoderReconstructionError',decoderCheck);
    rows{area} = table(string(name),nQ,nnz(abs(weight)>1e-8), ...
        nValid,line.slope,slopeCI(1),slopeCI(2), ...
        line.rSquared,nnz(isfinite(bootSlope)), ...
        'VariableNames',{'area','nQuartets', ...
        'nNonzeroDecoderWeights','nEligibleSiteQuartets', ...
        'slopeMsPerWidth','slopeCILow','slopeCIHigh', ...
        'timingRegressionRSquared','nValidBootstrap'});
    fprintf('%s: %d quartets, %d eligible site-quartets; slope %.2f [%.2f %.2f] ms/width\n', ...
        name,nQ,nValid,line.slope,slopeCI(1),slopeCI(2));
end

summary = vertcat(rows{:});
if ~isfolder(opt.OutputDir); mkdir(opt.OutputDir); end
figureFile = fullfile(opt.OutputDir,[stem '.png']);
summaryFile = fullfile(opt.OutputDir,[stem '.csv']);
resultFile = fullfile(opt.OutputDir,[stem '.mat']);
makeFigure(result,binTime,figureFile,opt.ITRFRoleOnly);
OUT = struct('result',{result},'summary',summary, ...
    'timeMs',binTime,'options',opt, ...
    'sourceDecoderFiles',{{lowerFile,itFile}}, ...
    'sourceBaseFile',lowerBaseFile, ...
    'resultFile',resultFile,'summaryFile',summaryFile, ...
    'figureFile',figureFile, ...
    'description',['Fixed full-data V1, V4, and IT attention-decoder ' ...
    'weights applied within six pooled RF distance groups. ' ...
    'Distance is Tall.along_GC, the cue-to-RF projection along ' ...
    'a capsule divided by stimulus-specific capsule width; it ' ...
    'is not an exact growth-cone shortest path. Groups use all ' ...
    'eligible non-overlap RF site-quartets before applying ' ...
    'decoder weights. Scores sum weighted site contributions ' ...
    'without dividing by site count, then average quartets ' ...
    'equally. IT sites without RF centers are not distance ' ...
    'eligible. Sigmoids fit 80-450 ms; observed traces are ' ...
    'smoothed by five bins for display only. Decoder weights ' ...
    'were trained on the full dataset; quartet bootstrap keeps ' ...
    'those weights fixed.']);
if opt.ITRFRoleOnly
    OUT.description = [OUT.description ' IT alignment and training use ' ...
        'only binary RF target/distractor roles; background and ' ...
        'overlap have zero alignment.'];
end
writetable(summary,summaryFile);
save(resultFile,'OUT','-v7.3');
disp(summary);
end

function members = quartetMembers(q)
block = floor((q-1)/2);
if mod(q,2)==1
    offset = [1 2 5 6];
else
    offset = [3 4 7 8];
end
members = 8*block+offset;
end

function fit = fitGroups(Y,t,x)
G = struct('groupMeanSigned',Y, ...
    'timeWindows',[t(:) t(:)], ...
    'groupSummary',table((1:6).',x(:),x(:), ...
    'VariableNames',{'groupIdx','alongMin','alongMax'}));
fit = fit_grouped_alonggc_sigmoid_sharedslope(G, ...
    'fitStartMs',80,'fitEndMs',450, ...
    'smoothW',1,'plotMode','none', ...
    'useT50SEWeightedRegression',false,'verbose',false);
end

function line = fitSlope(x,y)
x = double(x(:)); y = double(y(:));
X = [ones(numel(x),1) x];
beta = X\y;
residual = y-X*beta;
total = sum((y-mean(y)).^2);
line = struct('intercept',beta(1),'slope',beta(2), ...
    'rSquared',1-sum(residual.^2)/max(total,eps));
end

function makeFigure(result,t,figureFile,itRfRoleOnly)
h = figure('Color','w','Visible','off', ...
    'Position',[50 50 1900 1360]);
layout = tiledlayout(h,3,2,'Padding','compact', ...
    'TileSpacing','compact');
colors = [0.07 0.35 0.82; ...
          0.00 0.60 0.73; ...
          0.04 0.57 0.31; ...
          0.48 0.62 0.07; ...
          0.67 0.48 0.02; ...
          0.79 0.12 0.10];
for a = 1:3
    Z = result{a};
    ax = nexttile(layout,2*a-1);
    hold(ax,'on');
    fitMax = max(Z.fit.yFit(:,t>=0 & t<=500),[],2);
    for g = 1:6
        observed = movmean(Z.curve(g,:)/fitMax(g), ...
            5,2,'omitnan');
        plot(ax,t,observed,'Color', ...
            0.78*colors(g,:)+0.22*[1 1 1], ...
            'LineWidth',1.2,'HandleVisibility','off');
        plot(ax,t,Z.fit.yFit(g,:)/fitMax(g), ...
            'Color',colors(g,:),'LineWidth',2, ...
            'DisplayName',sprintf('Group %d',g));
    end
    xline(ax,0,':','HandleVisibility','off');
    xlim(ax,[0 500]); ylim(ax,[-0.12 1.2]);
    grid(ax,'on'); box(ax,'off');
    ylabel(ax,'Fit normalized to max = 1');
    weightedQuartets = sum(Z.weightedSitesByQuartet>0,1);
    title(ax,sprintf('%s | %d total quartets | quartets with weighted sites %s', ...
        Z.area,numel(Z.eligibleQuartets), ...
        mat2str(weightedQuartets)), ...
        'FontWeight','bold');
    if a==1; legend(ax,'Location','southeast'); end
    ax = nexttile(layout,2*a);
    hold(ax,'on');
    for g = 1:6
        x = Z.medianDistance(g);
        y = Z.fit.summary.t50(g);
        ci = Z.t50CI(g,:);
        if all(isfinite(ci))
            errorbar(ax,x,y,max(y-ci(1),0),max(ci(2)-y,0), ...
                'o','Color',colors(g,:), ...
                'MarkerFaceColor',colors(g,:), ...
                'MarkerSize',7,'CapSize',0);
        else
            plot(ax,x,y,'o','Color',colors(g,:), ...
                'MarkerFaceColor',colors(g,:),'MarkerSize',7);
        end
    end
    xx = linspace(min(Z.medianDistance), ...
        max(Z.medianDistance),100);
    plot(ax,xx,Z.line.intercept+Z.line.slope*xx, ...
        'k--','LineWidth',1.5);
    grid(ax,'on'); box(ax,'off');
    ylabel(ax,'Fitted t_{50} (ms)');
    title(ax,sprintf('%s | slope %+.2f [%.2f, %.2f] ms/width', ...
        Z.area,Z.line.slope,Z.slopeCI(1),Z.slopeCI(2)), ...
        'FontWeight','bold');
    if a==3
        xlabel(ax,'Cue-to-RF along-capsule distance / width');
    end
end
figureTitle = ['Nilson capsules | fixed attention-decoder weights | ' ...
    'six pooled distance groups | overlap excluded'];
if itRfRoleOnly
    figureTitle = [figureTitle ' | RF-role-only IT'];
end
title(layout,figureTitle, ...
    'FontWeight','bold');
print(h,figureFile,'-dpng','-r180');
close(h);
end
