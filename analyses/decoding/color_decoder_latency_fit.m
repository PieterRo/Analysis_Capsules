function color_decoder_latency_fit(monkey,area,nBoot,fitStartMs,fitEndMs,outDir)
% Same sigmoid form as capsule attention, with an explicit fitting window.
if nargin<4; fitStartMs = 0; end
if nargin<5; fitEndMs = 200; end
cfg = config();
if nargin<6 || isempty(outDir); outDir = cfg.resultsDir; end
monkey = upper(char(monkey)); area = upper(char(area));
assert(ismember(monkey,{'N','F'}) && ismember(area,{'V1','V4','IT'}));
assert(nBoot>=0 && nBoot==floor(nBoot));
source = fullfile(outDir,sprintf('Color_decoder_timecourse_sites_%s_%s.mat',monkey,area));
D = load(source);
t = double(D.binTime(:)');
nSites = numel(D.sites);
nBins = numel(t);
siteCurve = reshape(mean(double(D.siteContrib),2),nSites,nBins);
baselineMask = t>=-150 & t<0;
siteCurve = siteCurve-mean(siteCurve(:,baselineMask),2);
curve = sum(siteCurve,1);
G = struct('groupMeanSigned',curve,'timeWindows',[t(:) t(:)]);
fit = fit_grouped_alonggc_sigmoid_sharedslope(G, ...
    'fitStartMs',fitStartMs,'fitEndMs',fitEndMs,'smoothW',1,'plotMode','none', ...
    'useT50SEWeightedRegression',false,'verbose',false);
assert(fit.exitflag>0 && isfinite(fit.summary.t50));
fitMask = t>=fitStartMs & t<=fitEndMs;
observed = curve(fitMask);
predicted = fit.yFit(1,fitMask);
rSquared = 1-sum((observed-predicted).^2)/sum((observed-mean(observed)).^2);
bootT50 = nan(nBoot,1);
bootTau = nan(nBoot,1);
oldRng = rng;
restoreRng = onCleanup(@()rng(oldRng));
seed = 91026 + double(monkey) + sum(double(area));
rng(seed,'twister');
for b=1:nBoot
    drawn = randi(nSites,nSites,1);
    multiplicity = accumarray(drawn,1,[nSites 1]);
    sampleCurve = multiplicity'*siteCurve;
    G.groupMeanSigned = sampleCurve;
    try
        sampleFit = fit_grouped_alonggc_sigmoid_sharedslope(G, ...
            'fitStartMs',fitStartMs,'fitEndMs',fitEndMs,'smoothW',1,'plotMode','none', ...
            'useT50SEWeightedRegression',false,'verbose',false);
        if sampleFit.exitflag>0 && isfinite(sampleFit.summary.t50)
            bootT50(b) = sampleFit.summary.t50;
            bootTau(b) = sampleFit.summary.tauShared;
        end
    catch
        % An unconverged resample contributes no latency estimate.
    end
    if mod(b,200)==0 || b==nBoot
        fprintf('%s %s bootstrap %d/%d (%d valid)\n', ...
            monkey,area,b,nBoot,nnz(isfinite(bootT50)));
    end
end
valid = bootT50(isfinite(bootT50));
ci = [NaN NaN];
if ~isempty(valid); ci = prctile(valid,[2.5 97.5]); end
R = struct('monkey',monkey,'area',area,'nSites',nSites, ...
    'nNonzeroWeights',nnz(abs(D.weight(:))>1e-8), ...
    'nQuartets',numel(D.qValues),'timeMs',t,'curve',curve, ...
    'fit',fit,'fitStartMs',fitStartMs,'fitEndMs',fitEndMs,'rSquared',rSquared, ...
    't50',fit.summary.t50,'tau',fit.summary.tauShared, ...
    't50CI',ci,'bootT50',bootT50,'bootTau',bootTau, ...
    'nValidBootstrap',numel(valid),'nBootstrap',nBoot, ...
    'bootstrapUnit','site','weightsRefitted',false, ...
    'sourceFile',source,'reconstructionError',D.maxReconstructionError);
if fitStartMs==0 && fitEndMs==200
    suffix = '';
else
    suffix = sprintf('_%dto%d',fitStartMs,fitEndMs);
end
outFile = fullfile(outDir,sprintf('Color_decoder_latency_fit_%s_%s_%d%s.mat', ...
    monkey,area,nBoot,suffix));
save(outFile,'R','-v7.3');
fprintf('%s %s: t50 %.3f [%.3f %.3f] ms, tau %.3f, R2 %.3f, %d/%d valid; saved %s\n', ...
    monkey,area,R.t50,ci(1),ci(2),R.tau,R.rSquared,numel(valid),nBoot,outFile);
end
