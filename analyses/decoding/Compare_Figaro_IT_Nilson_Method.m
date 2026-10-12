function OUT = Compare_Figaro_IT_Nilson_Method(varargin)
% Compare Figaro IT using the same functional masks and decoder as Nilson.
% Outputs are isolated; original Figaro decoders and figure remain intact.
p=inputParser;
p.addParameter('OutputDir',fullfile(pwd,'tmp','figaro_nilson_matched'),@(x)ischar(x)||isstring(x));
p.addParameter('RunFixedPool',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});opt=p.Results;
repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(repoRoot);addpath(genpath(fullfile(repoRoot,'core')));addpath(fullfile(repoRoot,'analyses','line_stimuli'));addpath(fullfile(repoRoot,'analyses','color_tuning'));
cfg=config();sourceCfg=cfg;root=char(opt.OutputDir);if ~isfolder(root);mkdir(root);end
cfg.matDir=fullfile(root,'data_mat');cfg.resultsDir=fullfile(root,'results');if ~isfolder(cfg.matDir);mkdir(cfg.matDir);end;if ~isfolder(cfg.resultsDir);mkdir(cfg.resultsDir);end
geomFile=fullfile(fileparts(sourceCfg.resultsRoot),'Figaro files','Tall_IT_lines_F.mat');
G=load(geomFile);assert(size(G.ALLTAB,1)==216&&numel(G.Tall_IT)==216);
% ALLTAB is Figaro's native stimulus table. The alias changes no values.
G.RTAB384=G.ALLTAB;save(fullfile(cfg.matDir,'Tall_IT_lines_F.mat'),'-struct','G','-v7.3');
responseFile=fullfile(sourceCfg.resultsDir,'SNR_capsules_F_d12.mat');copyfile(responseFile,fullfile(cfg.matDir,'SNR_capsules_F_d12.mat'));
I=load(fullfile(cfg.dataRoot,'Figaro','ObjAtt_lines_MUA_trials.mat'),'ALLMAT');days=unique(I.ALLMAT(:,10))';
A=load(fullfile(sourceCfg.resultsRoot,'FigaroDecoderCrossover','figaro_capsules_elastic_net_best_quartet_IT_snr0p5.mat'),'OUT');previous=A.OUT.IT;
assert(numel(previous.sites)==94&&previous.nNonzero==32);
 assert(isequal(find(I.ALLMAT(:,9)==1),double(previous.trialIndex(:))));
BASE=Attention_TargetSide_Tuning_IT(struct('Monkey',2,'Config',cfg,'makeSummaryFigures',false,'makeExampleFigures',false,'makeAngleReferenceFigure',false,'saveResult',true,'forceRefit',true)); %#ok<NASGU>
ATT=Attention_TargetSide_DistanceControl_IT(struct('Monkey',2,'Config',cfg,'makeSummaryFigure',false,'saveResult',true,'forceRefit',true)); %#ok<NASGU>
% Nilson requires 80/96 usable complementary pairs for color direction.
% Preserve this fraction for Figaro's 54 pairs: ceil(80*54/96)=45.
colorMinPairs=ceil(80*54/96);
CD=Color_direction_tuning_IT('Monkey','F','Config',cfg,'ResponseSNRThreshold',0,'MinPairs',colorMinPairs,'MakeFigure',false);
% A constant spatial rescaling does not alter these regression p-values.
% Use native pixels for selection; no Figaro degrees calibration is assumed.
CP=Color_tuning_IT_nearest_curve('Monkey','F','Config',cfg,'ResponseSNRThreshold',0,'PixelsPerDegree',1,'MakeFigure',false);
CP.selectionDistanceUnits='pixels';
% Save under the shared analysis contract.
OUT=CP;save(CP.resultFile,'OUT');
selection=IT_functional_inclusion_N('Monkey','F','Config',cfg,'VisualSNRThreshold',0.5);
matched=Attention_decoding_elastic_net_CV_IT('Monkey','F','Config',cfg,'Days',days,'ResponseSNRThreshold',0.5,'MakeFigure',false);
fixed=[];
if opt.RunFixedPool
 fixed=Attention_decoding_elastic_net_CV_IT('Monkey','F','Config',cfg,'Days',days,'ResponseSNRThreshold',0.5,'CandidateSiteGlobal',double(previous.sites(:)),'MakeFigure',false);
 assert(isequal(sort(double(fixed.IT.siteIndexGlobal(:))),sort(double(previous.sites(:)))));
 fixed.siteSelection='fixedCandidates';
 fixed.siteInclusion='Original Figaro 94-site candidate IDs; same Nilson alignment and CV fitting.';
 fixed.siteSelectionCaveat='Candidate IDs inherited from the original Figaro selection.';
 fixed.candidateSiteGlobal=double(previous.sites(:));
 OUT=fixed;save(fixed.resultFile,'OUT','-v7.3');
end
OUT=struct('description','Figaro IT matched to Nilson direction/proximity alignment and cross-validated elastic net','matched',matched,'fixedPool',fixed,'functionalSelection',selection,'previous',previous,'geometrySource',geomFile,'responseSource',responseFile,'colorDirectionMinPairs',colorMinPairs,'colorDirection',CD,'colorProximity',CP);
OUT.resultFile=fullfile(root,'Figaro_IT_Nilson_method_comparison.mat');save(OUT.resultFile,'OUT','-v7.3');
fprintf('Matched Figaro IT: %d candidates, %d nonzero, CV accuracy %.2f%%\n',matched.IT.nIncludedSites,matched.IT.finalNonzero,100*matched.IT.crossValidatedAccuracyOfMeanScore);
end
