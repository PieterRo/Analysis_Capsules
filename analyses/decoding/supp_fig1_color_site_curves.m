function supp_fig1_color_site_curves(monkey,area,outDir)
% Matched complementary-pair color responses in attention-decoder sites.
cfg = config();
if nargin<3 || isempty(outDir); outDir = cfg.resultsDir; end
monkey = upper(char(monkey)); area = upper(char(area));
assert(ismember(monkey,{'N','F'}) && ismember(area,{'V1','V4','IT'}));
if strcmp(monkey,'F')
    dataDir = 'Figaro';
    geometryDir = fullfile(fileparts(cfg.resultsRoot),'Figaro files');
    lowerName = 'Color_configuration_decoder_elastic_net_CV_V1_V4_quartetCentered_FigaroSNR0p5_minSites1_F.mat';
    attentionName = 'Attention_decoder_Figaro_capsules_snr0p5_minSites1_F.mat';
    itName = 'Color_configuration_decoder_elastic_net_CV_IT_combined_ITfunctional_F.mat';
    snrFile = fullfile(cfg.resultsDir,'SNR_capsules_F_d12.mat');
else
    dataDir = 'Mr Nilson';
    geometryDir = cfg.matDir;
    lowerName = 'Color_configuration_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat';
    attentionName = 'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat';
    itName = 'Color_configuration_decoder_elastic_net_CV_IT_combined_ITfunctional_N.mat';
    snrFile = fullfile(cfg.matDir,'SNR_capsules_N_d12.mat');
end
geometryPath = fullfile(geometryDir,sprintf('Tall_%s_lines_%s.mat',area,monkey));
if strcmp(area,'IT')
    geometry = load(geometryPath,['Tall_' area],'RFrange');
else
    geometry = load(geometryPath,['Tall_' area]);
end
Tall = geometry.(['Tall_' area]);
[~,order] = sort([Tall.stimNum]); Tall = Tall(order);
assert(isequal([Tall.stimNum],1:numel(Tall)));
if strcmp(area,'IT')
    C = load(fullfile(cfg.resultsDir,itName),'OUT');
    R = C.OUT.IT;
    S = load(snrFile,'R');
    snr = compute_snr_per_color_region(S.R,Tall,double(geometry.RFrange(:)));
    scale = double(snr.sdSpont(double(R.siteIndexLocal(:))));
else
    C = load(fullfile(cfg.resultsDir,lowerName),'OUT');
    R = C.OUT.(area);
    A = load(fullfile(cfg.resultsDir,attentionName),'OUT');
    attention = A.OUT.(area);
    if strcmp(monkey,'F')
        expected = double(attention.globalSites(attention.eligibleSite(:)));
        assert(isequal(sort(double(R.siteIndexGlobal(:))),sort(expected)));
        [found,index] = ismember(double(R.siteIndexGlobal(:)), ...
            double(attention.globalSites(:)));
        assert(all(found));
    else
        assert(isequal(double(R.siteIndexGlobal(:)), ...
            double(attention.siteIndexGlobal(:))));
        index = double(R.siteIndexLocal(:));
    end
    scale = double(attention.responseScale(index));
end
globalSites = double(R.siteIndexGlobal(:));
localSites = double(R.siteIndexLocal(:));
nSites = numel(globalSites); nStim = numel(Tall);
assert(all(isfinite(scale) & scale>0));
roleByStim = zeros(nSites,nStim,'int8');
for stim=1:nStim
    color = string(Tall(stim).T.center_color(localSites));
    roleByStim(color=="yellowArm",stim)=1;
    roleByStim(color=="purple",stim)=-1;
end
if isfield(R,'colorRole')
    assert(isequal(double(roleByStim(:,double(R.stimulus(:))))', ...
        double(R.colorRole)));
end
M = load(fullfile(cfg.dataRoot,dataDir,'ObjAtt_lines_MUA_trials.mat'), ...
    'tb','ALLMAT');
tb = double(M.tb(:)');
assert(isequal(tb,-200:499));
timeMs = mean(reshape(tb,10,[]),1);
trial = double(R.trialIndex(:));
stimulus = double(R.stimulus(:));
day = double(M.ALLMAT(trial,10));
days = unique(day)';
assert(all(M.ALLMAT(trial,9)==1));
pre = timeMs>=-150 & timeMs<0;
early = timeMs>=40 & timeMs<240;
late = timeMs>=300 & timeMs<500;
if strcmp(area,'IT'); selectedWindow = late; else; selectedWindow = early; end
m = matfile(fullfile(cfg.dataRoot,dataDir,'ObjAtt_lines_normMUA.mat'));
preferred = nan(nSites,numel(timeMs));
other = nan(nSites,numel(timeMs));
pValue = nan(nSites,1);
meanDifference = nan(nSites,1);
nPairDays = zeros(nSites,1);
preferredColor = zeros(nSites,1);
for s=1:nSites
    raw = double(squeeze(m.normMUA(globalSites(s),:,:)));
    binned = squeeze(mean(reshape(raw(trial,:),numel(trial),10,[]),2,'omitnan'));
    assert(all(isfinite(binned(:))));
    for d=days
        rows = day==d;
        z = binned(rows,pre);
        dayBaseline = mean(z(:));
        binned(rows,:) = (binned(rows,:)-dayBaseline)/scale(s);
    end
    sumY = zeros(1,numel(timeMs)); sumP = sumY;
    weightTotal = 0; weightedDifference = 0;
    pairDifference = [];
    for block=0:(nStim/8-1)
        for offset=1:4
            stimA = 8*block+offset; stimB = stimA+4;
            roleA = double(roleByStim(s,stimA));
            roleB = double(roleByStim(s,stimB));
            if roleA==0 || roleB~=-roleA; continue; end
            for d=days
                ia = find(stimulus==stimA & day==d);
                ib = find(stimulus==stimB & day==d);
                nEff = min(numel(ia),numel(ib));
                if nEff==0; continue; end
                meanA = mean(binned(ia,:),1);
                meanB = mean(binned(ib,:),1);
                if roleA==1
                    y = meanA; p = meanB;
                else
                    y = meanB; p = meanA;
                end
                sumY = sumY+nEff*y;
                sumP = sumP+nEff*p;
                weightTotal = weightTotal+nEff;
                diffValue = mean(y(selectedWindow)-p(selectedWindow));
                pairDifference(end+1,1) = diffValue; %#ok<AGROW>
                weightedDifference = weightedDifference+nEff*diffValue;
            end
        end
    end
    nPairDays(s) = numel(pairDifference);
    if weightTotal==0 || numel(pairDifference)<3; continue; end
    [~,pValue(s)] = ttest(pairDifference);
    meanDifference(s) = weightedDifference/weightTotal;
    preferredColor(s) = sign(meanDifference(s));
    avgY = sumY/weightTotal; avgP = sumP/weightTotal;
    if preferredColor(s)>=0
        preferred(s,:) = avgY; other(s,:) = avgP;
    else
        preferred(s,:) = avgP; other(s,:) = avgY;
    end
    if mod(s,20)==0 || s==nSites
        fprintf('%s %s site %d/%d\n',monkey,area,s,nSites);
    end
end
significant = isfinite(pValue) & pValue<0.05;
assert(nnz(significant)>0,'%s %s has no significant sites',monkey,area);
muPreferred = mean(preferred(significant,:),1);
muOther = mean(other(significant,:),1);
semPreferred = std(preferred(significant,:),0,1)/sqrt(nnz(significant));
semOther = std(other(significant,:),0,1)/sqrt(nnz(significant));
result = struct('monkey',monkey,'area',area,'globalSites',globalSites, ...
    'timeMs',timeMs,'pValue',pValue,'significant',significant, ...
    'preferredColor',preferredColor,'meanDifference',meanDifference, ...
    'nPairDays',nPairDays,'nCandidateSites',nSites, ...
    'nSignificantSites',nnz(significant), ...
    'nYellowPreferred',nnz(significant & preferredColor>0), ...
    'nPurplePreferred',nnz(significant & preferredColor<0), ...
    'selectionWindowMs',ternary(strcmp(area,'IT'),[300 500],[40 240]), ...
    'muPreferred',muPreferred,'muOther',muOther, ...
    'semPreferred',semPreferred,'semOther',semOther, ...
    'selectionRule','Two-sided paired t test across complementary stimulus-pair by recording-day differences, p<0.05 uncorrected.', ...
    'normalizationRule','Recording-day spontaneous baseline subtraction, division by site spontaneous SD; pair-day means weighted by min trial count, sites averaged equally.');
outFile = fullfile(outDir,sprintf('SuppFig1_color_preferred_sites_%s_%s.mat',monkey,area));
save(outFile,'result','-v7.3');
fprintf('%s %s: significant %d/%d (%d yellow, %d purple), median pair-days %.0f; saved %s\n', ...
    monkey,area,nnz(significant),nSites,result.nYellowPreferred, ...
    result.nPurplePreferred,median(nPairDays(significant)),outFile);
end

function value = ternary(condition,yes,no)
if condition; value = yes; else; value = no; end
end
