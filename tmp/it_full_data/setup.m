repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(repoRoot);addpath('analyses/decoding');addpath(genpath('core'));cfg=config();root=fullfile(pwd,'tmp','it_full_data');
N=Compare_Capsule_IT_FullData(fullfile(cfg.resultsDir,'Attention_decoder_elastic_net_CV_IT_combined_nearestCurve_ITfunctional_N.mat'),'OutputDir',root);
F=Compare_Capsule_IT_FullData(fullfile(pwd,'tmp','figaro_nilson_matched','results','Attention_decoder_elastic_net_CV_IT_combined_nearestCurve_ITfunctional_F.mat'),'OutputDir',root);
checkcode('analyses/decoding/Compare_Capsule_IT_FullData.m');
