repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(repoRoot);addpath('analyses/decoding');addpath(genpath('core'));cfg=config();root=fullfile(pwd,'tmp','it_full_data');
file=fullfile(root,'Attention_decoder_elastic_net_FullData_IT_combined_nearestCurve_ITfunctional_N.mat');
Z=Compare_Figaro_IT_Matched_Distance(file,'Monkey','N','Config',cfg,'OutputDir',root,'BootstrapIterations',10000);
