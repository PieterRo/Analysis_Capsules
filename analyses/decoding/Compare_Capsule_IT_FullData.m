function OUT = Compare_Capsule_IT_FullData(sourceFile,varargin)
% Refit saved full-data IT features using Figaro V1/V4's fixed fitting rule.
% Candidate IDs, trial IDs, quartet centering and IT alignment are inherited.
p=inputParser;
p.addParameter('OutputDir',fullfile(pwd,'tmp','it_full_data'),@(x)ischar(x)||isstring(x));
p.addParameter('Alpha',0.25,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
p.addParameter('NumLambda',12,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('LambdaRatio',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.parse(varargin{:});opt=p.Results;
A=load(sourceFile,'OUT');previous=A.OUT;D=previous.IT;
X=double(D.fullDataFeatureMatrix);label=double(D.classLabel(:));
assert(size(X,1)==numel(label)&&size(X,2)==numel(D.siteIndexGlobal));
assert(all(isfinite(X(:)))&&all(ismember(label,[-1 1])));
variable=std(X,0,1)>1e-8;assert(nnz(variable)>=2);
[B,stats]=lassoglm(X(:,variable),label==1,'binomial', ...
    'Alpha',opt.Alpha,'NumLambda',opt.NumLambda, ...
    'LambdaRatio',opt.LambdaRatio,'CV','resubstitution');
assert(stats.Lambda(1)==min(stats.Lambda));
coefficient=zeros(size(X,2),1);coefficient(variable)=B(:,1);
score=X*coefficient+stats.Intercept(1);
fields={'region','trialIndex','stimulus','quartet','classLabel', ...
    'siteIndexLocal','siteIndexGlobal','nAllITSites','nValidNoiseSites', ...
    'nCompleteSites','nVisuallyResponsiveSites','nFunctionallyIncludedSites', ...
    'nIncludedSites','visualResponseScore','fullDataFeatureMatrix', ...
    'fullDataQuartetMeanResponse','fullDataAlignmentByQuartetSite', ...
    'fullDataTuningCoefficient','fullDataTuningR2','fullDataTuningStatus', ...
    'fullDataDirectionAmplitude','fullDataDirectionPreferenceDeg', ...
    'fullDataProximityCoefficient','responseWindowSampled','geometryQC', ...
    'targetDirectionDegByQuartet','proximityAdvantageBySiteQuartet','hasRFcenter'};
R=struct();for i=1:numel(fields);if isfield(D,fields{i});R.(fields{i})=D.(fields{i});end;end
R.finalCoefficient=coefficient;R.finalIntercept=stats.Intercept(1);
R.finalClassScore=score;R.finalTargetAlignedScore=label.*score;
R.finalNonzero=nnz(abs(coefficient)>1e-8);
R.finalInSampleAccuracy=mean((score>0)==(label==1));
R.selectedAlpha=opt.Alpha;R.selectedLambda=stats.Lambda(1);
R.lambdaPath=stats.Lambda;R.coefficientPath=B;R.variableFeature=variable;
R.selectionNote='Fixed alpha; first returned (least-penalized) lambda; no cross-validation or class-balancing weights, matching Figaro V1/V4.';
assert(isequal(R.siteIndexGlobal,D.siteIndexGlobal)&&isequal(R.trialIndex,D.trialIndex));
assert(isequal(R.fullDataAlignmentByQuartetSite,D.fullDataAlignmentByQuartetSite));
OUT=struct('description','Full-data IT refit with the fixed Figaro V1/V4 elastic-net fitting rule', ...
    'monkey',previous.monkey,'responseWindow',previous.responseWindow, ...
    'responseCentering','full-data quartet means','responseNormalization','spontaneous SD', ...
    'sourceDecoderFile',sourceFile,'options',opt,'IT',R, ...
    'selectionRule',R.selectionNote,'accuracyDefinition','In-sample; no held-out performance estimate');
for field={'siteInclusion','siteSelection','responseSNRThreshold','siteSelectionCaveat','geometryModel','alignmentNormalization','ITFunctionalSelection'}
 if isfield(previous,field{1});OUT.(field{1})=previous.(field{1});end
end
OUT.siteSelectionCaveat='Inherited full-data functional site selection; this refit uses no cross-validation.';
[~,stem]=fileparts(sourceFile);stem=strrep(stem,'elastic_net_CV','elastic_net_FullData');
OUT.resultFile=fullfile(char(opt.OutputDir),[stem '.mat']);
assert(~isfile(OUT.resultFile),'Refusing to overwrite an existing comparison.');
if ~isfolder(char(opt.OutputDir));mkdir(char(opt.OutputDir));end
save(OUT.resultFile,'OUT','-v7.3');
fprintf('%s IT full data: %d candidates, %d nonzero, %d trials, in-sample %.2f%%, alpha %.2f, lambda %.6g\n', ...
 OUT.monkey,numel(coefficient),R.finalNonzero,numel(label),100*R.finalInSampleAccuracy,R.selectedAlpha,R.selectedLambda);
end
