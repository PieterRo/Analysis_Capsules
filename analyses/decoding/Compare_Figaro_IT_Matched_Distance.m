function OUT = Compare_Figaro_IT_Matched_Distance(decoderFile,varargin)
% Apply Nilson's distance grouping, timing fit and fixed-weight site bootstrap.
repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(repoRoot);addpath(genpath(fullfile(repoRoot,'core')));
cfg=config();p=inputParser;
p.addParameter('Config',struct(),@isstruct);p.addParameter('OutputDir',fullfile(repoRoot,'tmp','figaro_nilson_matched'),@(x)ischar(x)||isstring(x));
p.addParameter('BootstrapIterations',10000,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x==floor(x));p.addParameter('RandomSeed',270926,@isscalar);
p.addParameter('Monkey','F',@(x)any(strcmpi(string(x),["N" "F"])));
p.addParameter('Area','IT',@(x)any(strcmpi(string(x),["V1" "V4" "IT"])));
p.parse(varargin{:});opt=p.Results;if ~isempty(fieldnames(opt.Config));cfg=opt.Config;end
suffix=upper(char(opt.Monkey));monkeyName='Figaro';nStimuli=216;
if suffix=='N';monkeyName='Mr Nilson';nStimuli=384;end
name=upper(char(opt.Area));A=load(decoderFile,'OUT');R=A.OUT.(name);
Geo=load(fullfile(cfg.matDir,['Tall_' name '_lines_' suffix '.mat']));Tall=Geo.(['Tall_' name]);
m=matfile(fullfile(cfg.dataRoot,monkeyName,'ObjAtt_lines_normMUA.mat'));M=load(fullfile(cfg.dataRoot,monkeyName,'ObjAtt_lines_MUA_trials.mat'),'tb');
sampleTime=double(M.tb(:)');assert(isequal(sampleTime,-200:499));binTime=sampleTime(1:10:end)+4.5;
baselineBins=binTime>=-150&binTime<0;lateBins=binTime>=300&binTime<500;
oldRng=rng;restoreRng=onCleanup(@()rng(oldRng));rng(opt.RandomSeed,'twister');
    nSites = numel(R.siteIndexGlobal);
    localSite = double(R.siteIndexLocal(:));
    globalSite = double(R.siteIndexGlobal(:));
    weight = double(R.finalCoefficient(:));
    trial = double(R.trialIndex(:));
    quartet = double(R.quartet(:));
    label = double(R.classLabel(:));
    eligible = unique(quartet);
    nQ = numel(eligible);
    assert(numel(Tall)==nStimuli && numel(weight)==nSites && ...
        numel(localSite)==nSites && numel(trial)==numel(quartet) && ...
        all(ismember(label,[-1 1])));
    S = load(fullfile(cfg.matDir,['SNR_capsules_' suffix '_d12.mat']),'R');
    if isfield(R,'distanceResponseScale')
        responseScale=double(R.distanceResponseScale(:));
    else
        snr = compute_snr_per_color_region(S.R,Tall,Geo.RFrange(:));
        responseScale = double(snr.sdSpont(localSite));
    end
    assert(isequal(size(R.fullDataAlignmentByQuartetSite),[nStimuli/4 nSites]));
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
        siteSign(:,k)=double(R.fullDataAlignmentByQuartetSite(q,:)).';
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
    nAllTrials=size(m,'normMUA'); nAllTrials=nAllTrials(2);
    trialRow=zeros(nAllTrials,1);trialRow(trial)=1:nTrials;
    sums=zeros(nSites,nQ,2,nBins);counts=zeros(nQ,2);
    for first=1:200:nAllTrials
        last=min(first+199,nAllTrials);ix=first:last;use=trialRow(ix)>0;
        if ~any(use);continue;end
        channelSpan=min(globalSite):max(globalSite);
        raw=double(m.normMUA(channelSpan,first:last,:));
        raw=raw(globalSite-channelSpan(1)+1,:,:);
        binned=mean(reshape(raw,nSites,numel(ix),10,nBins),3,'omitnan');
        binned=reshape(binned,nSites,numel(ix),nBins);rows=trialRow(ix(use));
        for kk=1:nQ
            for cls=1:2
                lab=3-2*cls;choose=use;choose(use)=quartet(rows)==eligible(kk)&label(rows)==lab;
                if any(choose)
                    sums(:,kk,cls,:)=sums(:,kk,cls,:)+reshape(sum(binned(:,choose,:),2),nSites,1,1,nBins);
                    counts(kk,cls)=counts(kk,cls)+nnz(choose);
                end
            end
        end
        fprintf('IT time-course raw trials %d/%d\n',last,nAllTrials);
    end
    for ss=1:nSites
        for kk=1:nQ
            delta=reshape(sums(ss,kk,1,:)/counts(kk,1)-sums(ss,kk,2,:)/counts(kk,2),1,nBins)/responseScale(ss);
            contribution(ss,kk,:)=single(weight(ss)*siteSign(ss,kk)*delta);
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
    % Algebraically identical site bootstrap, with sums cached per site/group.
    assert(all(isfinite(contribution(:))));
    siteGroupSum=zeros(nSites,6,nBins);availability=false(nSites,nQ,6);
    for gg=1:6
        availability(:,:,gg)=groupMap==gg;
        siteGroupSum(:,gg,:)=sum(double(contribution).*reshape(availability(:,:,gg),nSites,nQ,1),2);
    end
    for b = 1:opt.BootstrapIterations
        if mod(b,100)==0;fprintf('IT site bootstrap %d/%d\n',b,opt.BootstrapIterations);end
        draw=randi(nSites,nSites,1);multiplicity=accumarray(draw,1,[nSites 1]);
        denom=sum(reshape((multiplicity'>0)*double(reshape(availability,nSites,[]))>0,nQ,6),1)';
        bootCurve=reshape(multiplicity'*reshape(siteGroupSum,nSites,[]),6,nBins)./denom;
        if b<=20
            reference=nan(nQ,6,nBins);
            for kk=1:nQ
                for gg=1:6
                    keep=groupMap(:,kk)==gg & multiplicity>0;
                    if any(keep);reference(kk,gg,:)=sum(double(contribution(keep,kk,:)).*reshape(multiplicity(keep),[],1,1),1);end
                end
            end
            refCurve=reshape(mean(reference,1,'omitnan'),6,nBins);
            assert(isequal(isnan(refCurve),isnan(bootCurve))&&max(abs(refCurve(:)-bootCurve(:)),[],'omitnan')<1e-10);
        end
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
    Z = struct('area',string(name), ...
        'eligibleQuartets',eligible,'medianDistance',medianDistance, ...
        'groupSiteQuartetCounts',groupCounts, ...
        'groupBySiteQuartet',groupMap,'globalSite',globalSite,'weights',weight, ...
        'candidateSitesByQuartet',candidateSites, ...
        'weightedSitesByQuartet',weightedSites, ...
        'quartetDifference',qGroup,'curve',curve,'fit',fit, ...
        'line',line,'t50CI',t50CI,'slopeCI',slopeCI, ...
        'bootstrapSlope',bootSlope, ...
        'decoderReconstructionError',decoderCheck);
    summary = table(string(name),nQ,nnz(abs(weight)>1e-8), ...
        nValid,line.slope,slopeCI(1),slopeCI(2), ...
        line.rSquared,nnz(isfinite(bootSlope)), ...
        'VariableNames',{'area','nQuartets', ...
        'nNonzeroDecoderWeights','nEligibleSiteQuartets', ...
        'slopeMsPerWidth','slopeCILow','slopeCIHigh', ...
        'timingRegressionRSquared','nValidBootstrap'});
    fprintf('%s: %d quartets, %d eligible site-quartets; slope %.2f [%.2f %.2f] ms/width\n', ...
        name,nQ,nValid,line.slope,slopeCI(1),slopeCI(2));
valid=Z.bootstrapSlope(isfinite(Z.bootstrapSlope));
Z.nValidBootstrap=numel(valid);
Z.nOppositeSign=min(nnz(valid<=0),nnz(valid>=0));
Z.slopePValue=min(1,2*Z.nOppositeSign/max(numel(valid),1));
if isempty(valid);Z.slopePValue=NaN;end
Z.pResolution=2/max(numel(valid),1);
OUT=struct('result',Z,'summary',summary,'timeMs',binTime,'options',opt,'sourceDecoderFile',decoderFile,'bootstrapUnit','site');
[~,stem]=fileparts(decoderFile);
if ~strcmp(name,'IT');stem=[stem '_' name];end
OUT.resultFile=fullfile(char(opt.OutputDir),[stem '_distance_sitebootstrap.mat']);save(OUT.resultFile,'OUT','-v7.3');
fprintf('%s %s slope %.4f, 95%% CI [%.4f %.4f], bootstrap p %.8g, %d valid\n',monkeyName,name,Z.line.slope,Z.slopeCI(1),Z.slopeCI(2),Z.slopePValue,numel(valid));
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
