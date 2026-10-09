function color_decoder_latency_prepare(monkey,area,outDir)
% Reconstruct fixed-weight full color-decoder time courses per recording site.
cfg = config();
if nargin<3 || isempty(outDir); outDir = cfg.resultsDir; end
monkey = upper(char(monkey));
area = upper(char(area));
assert(ismember(monkey,{'N','F'}) && ismember(area,{'V1','V4','IT'}));
if strcmp(monkey,'F')
    dataFolder = 'Figaro';
    lowerName = 'Color_configuration_decoder_elastic_net_CV_V1_V4_quartetCentered_FigaroSNR0p5_minSites1_F.mat';
    baseName = 'Attention_decoder_Figaro_capsules_snr0p5_minSites1_F.mat';
    snrFile = fullfile(cfg.resultsDir,'SNR_capsules_F_d12.mat');
    geometryDir = fullfile(fileparts(cfg.resultsRoot),'Figaro files');
else
    dataFolder = 'Mr Nilson';
    lowerName = 'Color_configuration_decoder_elastic_net_CV_V1_V4_quartetCentered_minSites20_N.mat';
    baseName = 'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat';
    snrFile = fullfile(cfg.matDir,'SNR_capsules_N_d12.mat');
    geometryDir = cfg.matDir;
end
if strcmp(area,'IT')
    sourceFile = fullfile(cfg.resultsDir, ...
        sprintf('Color_configuration_decoder_elastic_net_CV_IT_combined_ITfunctional_%s.mat',monkey));
    D = load(sourceFile,'OUT');
    R = D.OUT.IT;
    S = load(snrFile,'R');
    H = load(fullfile(geometryDir,sprintf('Tall_IT_lines_%s.mat',monkey)), ...
        'Tall_IT','RFrange');
    snr = compute_snr_per_color_region(S.R,H.Tall_IT,double(H.RFrange(:)));
    baseline = double(snr.muSpont(R.siteIndexLocal(:)));
    scale = double(snr.sdSpont(R.siteIndexLocal(:)));
else
    sourceFile = fullfile(cfg.resultsDir,lowerName);
    D = load(sourceFile,'OUT');
    R = D.OUT.(area);
    B = load(fullfile(cfg.resultsDir,baseName),'OUT');
    B = B.OUT.(area);
    if strcmp(monkey,'F')
        [found,index] = ismember(double(R.siteIndexGlobal(:)),double(B.globalSites(:)));
        assert(all(found));
    else
        assert(isequal(double(R.siteIndexGlobal(:)),double(B.siteIndexGlobal(:))));
        index = double(R.siteIndexLocal(:));
    end
    baseline = double(B.baseline(index));
    scale = double(B.responseScale(index));
end
sites = double(R.siteIndexGlobal(:));
trial = double(R.trialIndex(:));
quartet = double(R.quartet(:));
label = double(R.classLabel(:));
weight = double(R.finalCoefficient(:));
qValues = unique(quartet);
nSites = numel(sites); nQ = numel(qValues);
assert(numel(weight)==nSites && numel(scale)==nSites && ...
    all(isfinite(scale) & scale>0));
T = load(fullfile(cfg.dataRoot,dataFolder,'ObjAtt_lines_MUA_trials.mat'),'tb');
sampleTime = double(T.tb(:)');
assert(isequal(sampleTime,-200:499));
binTime = mean(reshape(sampleTime,10,[]),1);
nBins = numel(binTime);
late = binTime>=300 & binTime<500;
m = matfile(fullfile(cfg.dataRoot,dataFolder,'ObjAtt_lines_normMUA.mat'));
siteContrib = zeros(nSites,nQ,nBins,'single');
for s=1:nSites
    if abs(weight(s))<=1e-8
        continue;
    end
    raw = double(squeeze(m.normMUA(sites(s),:,:)));
    binned = squeeze(mean(reshape(raw(trial,:),numel(trial),10,nBins),2,'omitnan'));
    binned = (binned-baseline(s))/scale(s);
    assert(all(isfinite(binned(:))));
    for k=1:nQ
        rows = find(quartet==qValues(k));
        pos = label(rows)==1;
        neg = label(rows)==-1;
        assert(any(pos) && any(neg));
        if strcmp(area,'IT')
            a = double(R.fullDataAlignmentByQuartetSite(qValues(k),s));
            feature = binned(rows,:)*a;
        else
            targetColor = double(R.targetColor(rows));
            cellMean = zeros(2,2);
            for c=1:2
                for t=1:2
                    use = label(rows)==2*c-3 & targetColor==2*t-3;
                    assert(any(use));
                    z = binned(rows(use),late);
                    cellMean(c,t) = mean(z(:));
                end
            end
            center = mean(cellMean(:));
            sign = label(rows).*double(R.colorRole(rows,s));
            sign(~R.availability(rows,s)) = 0;
            feature = (binned(rows,:)-center).*sign;
        end
        delta = weight(s)*(mean(feature(pos,:),1)-mean(feature(neg,:),1));
        siteContrib(s,k,:) = single(delta);
    end
    if mod(s,20)==0 || s==nSites
        fprintf('%s %s %d/%d sites\n',monkey,area,s,nSites);
    end
end
saved = zeros(nQ,1);
for k=1:nQ
    rows = find(quartet==qValues(k));
    score = double(R.finalClassScore(rows));
    labs = label(rows);
    saved(k) = mean(score(labs==1))-mean(score(labs==-1));
end
reconstructed = reshape(mean(sum(double(siteContrib(:,:,late)),1),3),nQ,1);
maxReconstructionError = max(abs(saved-reconstructed));
assert(maxReconstructionError<1e-4, ...
    '%s %s decoder score mismatch: %g',monkey,area,maxReconstructionError);
outFile = fullfile(outDir,sprintf('Color_decoder_timecourse_sites_%s_%s.mat',monkey,area));
save(outFile,'siteContrib','binTime','qValues','sites','weight', ...
    'maxReconstructionError','sourceFile','-v7.3');
fprintf('%s %s: %d sites (%d nonzero), %d quartets, error %.3g; saved %s\n', ...
    monkey,area,nSites,nnz(abs(weight)>1e-8),nQ,maxReconstructionError,outFile);
end
