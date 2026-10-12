repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));cd(repoRoot);addpath('analyses/decoding');addpath(genpath('core'));cfg=config();root=fullfile(pwd,'tmp','it_full_data');
file=fullfile(root,'Attention_decoder_elastic_net_FullData_V1_V4_N.mat');
for area={'V1','V4'}
 Z=Compare_Figaro_IT_Matched_Distance(file,'Monkey','N','Area',area{1},'Config',cfg,'OutputDir',root,'BootstrapIterations',10000);
end
