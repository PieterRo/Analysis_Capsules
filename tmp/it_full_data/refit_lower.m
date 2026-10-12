repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(repoRoot);addpath(genpath('core'));cfg=config();root=fullfile(pwd,'tmp','it_full_data');if ~isfolder(root);mkdir(root);end
sourceFile=fullfile(cfg.resultsDir,'Attention_decoder_sparse_linear_V1_V4_quartetCentered_minSites20_N.mat');A=load(sourceFile,'OUT');base=A.OUT;
OUT=struct('description','Full-data Nilson V1/V4 elastic net with the fixed Figaro V1/V4 fitting rule','sourceDecoderFile',sourceFile,'alpha',0.25,'numLambda',12,'lambdaRatio',0.05,'accuracyDefinition','In-sample');
for area={'V1','V4'}
 name=area{1};D=base.(name);X=double(D.featureMatrix);label=double(D.classLabel(:));
 assert(all(isfinite(X(:)))&&size(X,2)==numel(D.siteIndexGlobal));variable=std(X,0,1)>1e-8;
 [B,stats]=lassoglm(X(:,variable),label==1,'binomial','Alpha',0.25,'NumLambda',12,'LambdaRatio',0.05,'CV','resubstitution');assert(stats.Lambda(1)==min(stats.Lambda));
 coefficient=zeros(size(X,2),1);coefficient(variable)=B(:,1);score=X*coefficient+stats.Intercept(1);
 R=struct();for field={'region','siteIndexGlobal','siteIndexLocal','trialIndex','stimulus','quartet','classLabel','responseWindowSampled'};R.(field{1})=D.(field{1});end
 R.finalCoefficient=coefficient;R.finalIntercept=stats.Intercept(1);R.finalClassScore=score;R.finalNonzero=nnz(abs(coefficient)>1e-8);R.finalInSampleAccuracy=mean((score>0)==(label==1));
 R.fullDataFeatureMatrix=X;R.fullDataAlignmentByQuartetSite=zeros(96,numel(coefficient));
 for q=unique(D.quartet(:))';ix=find(D.quartet==q);sg=double(D.canonicalSign(ix,:));assert(all(all(sg==sg(1,:))));R.fullDataAlignmentByQuartetSite(q,:)=sg(1,:);end
 R.distanceResponseScale=double(D.responseScale(D.siteIndexLocal));R.selectedAlpha=0.25;R.selectedLambda=stats.Lambda(1);R.lambdaPath=stats.Lambda;R.variableFeature=variable;
 assert(isequal(R.siteIndexGlobal,D.siteIndexGlobal)&&isequal(R.trialIndex,D.trialIndex));
 OUT.(name)=R;fprintf('Nilson %s full data: %d candidates, %d nonzero, %d trials, in-sample %.2f%%, lambda %.6g\n',name,numel(coefficient),R.finalNonzero,numel(label),100*R.finalInSampleAccuracy,R.selectedLambda);
end
OUT.resultFile=fullfile(root,'Attention_decoder_elastic_net_FullData_V1_V4_N.mat');assert(~isfile(OUT.resultFile));save(OUT.resultFile,'OUT','-v7.3');
