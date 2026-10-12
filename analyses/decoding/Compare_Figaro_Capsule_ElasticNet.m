function OUT = Compare_Figaro_Capsule_ElasticNet(varargin)
% CompHalf-style full-data elastic-net attention decoder for Figaro capsules.
% This is a separate comparison; it does not replace the d-prime decoder.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
cfg = config();
p = inputParser;
addParameter(p,'OutputDir',fullfile(cfg.resultsRoot,'FigaroDecoderCrossover'));
addParameter(p,'SNRThreshold',0.5);
addParameter(p,'Alpha',0.25);
addParameter(p,'NumLambda',12);
addParameter(p,'LambdaRatio',0.05);
addParameter(p,'SaveOutputs',true);
parse(p,varargin{:});
opt = p.Results;
assert(opt.SNRThreshold==0.5 && exist('lassoglm','file')==2);

oldFile = fullfile(cfg.resultsDir, ...
    'Attention_decoder_Figaro_capsules_snr0p5_minSites1_F.mat');
responseFile = firstExisting({fullfile(cfg.matDir,'SNR_capsules_F_d12.mat'), ...
    fullfile(cfg.resultsDir,'SNR_capsules_F_d12.mat')});
assert(isfile(oldFile) && isfile(responseFile));
A = load(oldFile,'OUT'); old = A.OUT;
A = load(responseFile,'R'); R = A.R;
m = matfile(fullfile(cfg.dataRoot,'Figaro','ObjAtt_lines_normMUA.mat'));
t = matfile(fullfile(cfg.dataRoot,'Figaro','ObjAtt_lines_MUA_trials.mat'));
info = double(t.ALLMAT);
tb = double(t.tb);
tb = tb(:);
[nAllSites,nAllTrials,nTimes] = size(m,'normMUA');
assert(nAllSites==1024 && nAllTrials==size(info,1) && ...
    nTimes==numel(tb) && isequal(size(R.meanAct),[1024 216 3]));
use = info(:,9)==1 & info(:,1)>=1 & info(:,1)<=216;
trialIndex = find(use);
stimulus = info(use,1);
quartet = stimulusToQuartet(stimulus);
label = mod(stimulus,2)==1;
lateSamples = find(tb>=300 & tb<=500);
assert(~isempty(lateSamples) && ...
    isequal(accumarray(stimulus,1,[216 1]),double(R.nTrials(:))));

regions = {'V1','V4','IT'};
OUT = struct('description', ...
    ['CompHalf-style elastic-net logistic decoder, fit and scored on ' ...
     'the same correct Figaro capsule trials. Features are quartet-' ...
     'centered late responses in spontaneous-SD units, multiplied by ' ...
     'the RF target/distractor sign. V1/V4 off-curve RFs are masked. ' ...
     'As in the CompHalf decoder, unlocalized IT sites have sign +1.'], ...
    'siteSelection', ...
    ['V1/V4: original best-color early/late SNR > 0.5 selection; ' ...
     'IT: best-quartet early/late SNR > 0.5, no RF required.'], ...
    'sourceDecoderFile',oldFile,'sourceResponseFile',responseFile, ...
    'snrThreshold',opt.SNRThreshold,'alpha',opt.Alpha, ...
    'numLambda',opt.NumLambda,'lambdaRatio',opt.LambdaRatio, ...
    'responseWindow',[300 500],'nIncludedTrials',numel(trialIndex));
for a = 1:3
    name = regions{a};
    D = old.(name);
    sites = double(D.globalSites(:));
    localR = R;
    localR.meanAct = R.meanAct(sites,:,:);
    localR.meanSqAct = R.meanSqAct(sites,:,:);
    geomFile = D.geometryFile;
    G = load(geomFile);
    Tall = G.(sprintf('Tall_%s',name));
    [~,order] = sort([Tall.stimNum]);
    Tall = Tall(order);
    SNR = compute_snr_per_color_sites(localR,Tall, ...
        (1:numel(sites))','Verbose',false);
    if strcmp(name,'IT')
        bestQuartetSNR = quartetSNR(localR,SNR,R.nTrials);
        eligible = bestQuartetSNR>opt.SNRThreshold;
    else
        bestQuartetSNR = nan(numel(sites),1);
        eligible = D.eligibleSite(:);
    end
    eligible = eligible & isfinite(SNR.muSpont) & ...
        isfinite(SNR.sdSpont) & SNR.sdSpont>0;
    siteRows = find(eligible);
    assert(~isempty(siteRows),'No selected %s sites.',name);
    selectedSites = sites(siteRows);
    role = double(D.roleByStimulus(siteRows,:));
    signByQuartet = zeros(numel(siteRows),54);
    for q = 1:54
        members = quartetMembers(q);
        signByQuartet(:,q) = role(:,members(1));
    end
    if strcmp(name,'IT')
        signByQuartet(signByQuartet==0) = 1;
    end
    late = zeros(numel(trialIndex),numel(siteRows),'single');
    for i = 1:numel(siteRows)
        raw = double(squeeze(m.normMUA(selectedSites(i), ...
            trialIndex,lateSamples)));
        assert(size(raw,1)==numel(trialIndex));
        late(:,i) = single((mean(raw,2,'omitnan')- ...
            SNR.muSpont(siteRows(i)))/SNR.sdSpont(siteRows(i)));
        if i==1 || mod(i,40)==0 || i==numel(siteRows)
            fprintf('Capsule elastic %s site %d/%d\n', ...
                name,i,numel(siteRows));
        end
    end
    means = zeros(54,numel(siteRows));
    for q = 1:54
        means(q,:) = mean(late(quartet==q,:),1);
    end
    X = (double(late)-means(quartet,:)).*signByQuartet(:,quartet).';
    assert(all(isfinite(X(:))),'Nonfinite %s features.',name);
    variable = std(X,0,1)>1e-8;
    assert(nnz(variable)>=2,'Too few variable %s sites.',name);
    weights = zeros(numel(siteRows),1);
    [B,stats] = lassoglm(X(:,variable),label,'binomial', ...
        'Alpha',opt.Alpha,'NumLambda',opt.NumLambda, ...
        'LambdaRatio',opt.LambdaRatio,'CV','resubstitution');
    weights(variable) = B(:,1);
    intercept = stats.Intercept(1);
    score = X*weights+intercept;
    Dnew = struct('region',name,'sites',selectedSites, ...
        'eligibleSiteLocal',eligible,'bestQuartetSNR',bestQuartetSNR, ...
        'baseline',SNR.muSpont(siteRows), ...
        'spontaneousSD',SNR.sdSpont(siteRows), ...
        'signByQuartet',signByQuartet,'weights',weights, ...
        'intercept',intercept,'lambda',stats.Lambda(1), ...
        'nSelected',numel(siteRows),'nNonzero',nnz(abs(weights)>1e-8), ...
        'nNoRFSelected',nnz(all(role==0,2)), ...
        'accuracyInSample',mean((score>0)==label), ...
        'nTrials',numel(label),'trialIndex',trialIndex, ...
        'trialQuartet',quartet,'trialStimulus',stimulus, ...
        'trialLabel',label,'trialScore',score);
    OUT.(name) = Dnew;
    fprintf('%s: selected %d, nonzero %d, no RF %d, in-sample %.1f%%\n', ...
        name,Dnew.nSelected,Dnew.nNonzero,Dnew.nNoRFSelected, ...
        100*Dnew.accuracyInSample);
end
OUT.resultFile = fullfile(char(opt.OutputDir), ...
    'figaro_capsules_elastic_net_best_quartet_IT_snr0p5.mat');
OUT.figureFile = fullfile(char(opt.OutputDir), ...
    'figaro_capsules_elastic_net_best_quartet_IT_snr0p5.png');
if opt.SaveOutputs
    if ~isfolder(char(opt.OutputDir)); mkdir(char(opt.OutputDir)); end
    assert(~isfile(OUT.resultFile) && ~isfile(OUT.figureFile), ...
        'Crossover output exists; refusing to overwrite it.');
    save(OUT.resultFile,'OUT','-v7.3');
    plotSummary(OUT,regions);
end
end

function best = quartetSNR(R,SNR,nTrials)
nSites = size(R.meanAct,1);
scores = nan(nSites,54,2);
nTrials = double(nTrials(:));
for q = 1:54
    members = quartetMembers(q);
    for w = 1:2
        values = double(R.meanAct(:,members,w+1));
        weights = repmat(nTrials(members).',nSites,1);
        good = isfinite(values) & weights>0;
        values(~good)=0; weights(~good)=0;
        mu = sum(values.*weights,2)./sum(weights,2);
        scores(:,q,w) = (mu-SNR.muSpont)./SNR.sdSpont;
    end
end
flat = reshape(scores,nSites,[]);
best = max(flat,[],2,'omitnan');
best(all(~isfinite(flat),2))=NaN;
end

function q = stimulusToQuartet(stimulus)
block = floor((stimulus-1)/8);
position = mod(stimulus-1,8)+1;
q = 2*block+1+ismember(position,[3 4 7 8]);
end

function members = quartetMembers(q)
block = floor((q-1)/2);
if mod(q,2)==1; offsets=[1 2 5 6];
else; offsets=[3 4 7 8]; end
members = 8*block+offsets;
end

function path = firstExisting(candidates)
for i=1:numel(candidates)
    if isfile(candidates{i}); path=candidates{i}; return; end
end
error('Required response file is missing.');
end

function plotSummary(OUT,regions)
h = figure('Color','w','Visible','off','Position',[100 100 1150 420]);
t = tiledlayout(h,1,3,'Padding','compact');
for i=1:3
    D = OUT.(regions{i});
    ax=nexttile(t,i);
    histogram(ax,D.trialScore(~D.trialLabel),30,'FaceColor',[0.85 .35 .30]);
    hold(ax,'on');
    histogram(ax,D.trialScore(D.trialLabel),30,'FaceColor',[.20 .48 .75]);
    xline(ax,0,'k:'); grid(ax,'on');
    title(ax,sprintf('%s | %d sites (%d nonzero) | %.1f%%', ...
        regions{i},D.nSelected,D.nNonzero,100*D.accuracyInSample));
    xlabel(ax,'Elastic-net decision score'); ylabel(ax,'Trials');
    legend(ax,{'D','T'},'Location','best');
end
title(t,'Figaro capsules | full-data elastic net | in-sample');
print(h,OUT.figureFile,'-dpng','-r180'); close(h);
end
