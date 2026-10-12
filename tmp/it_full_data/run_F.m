repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(repoRoot);addpath('analyses/decoding');addpath(genpath('core'));cfg=config();root=fullfile(pwd,'tmp','it_full_data');cfg.matDir=fullfile(pwd,'tmp','figaro_nilson_matched','data_mat');
file=fullfile(root,'Attention_decoder_elastic_net_FullData_IT_combined_nearestCurve_ITfunctional_F.mat');
Z=Compare_Figaro_IT_Matched_Distance(file,'Monkey','F','Config',cfg,'OutputDir',root,'BootstrapIterations',10000);
