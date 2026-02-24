% RUN_PIPELINE
% Entry point to initialize paths/config and manually run selected analyses.

clearvars;
close all;
clc;

% Repository root (this file is expected at the root).
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

% Load configuration from config.m
cfg = config(); %#ok<NASGU>

% -------------------------------------------------------------------------
% Run options by analysis category (uncomment what you want to execute)
% -------------------------------------------------------------------------

% LINE STIMULI
run(fullfile(rootDir, 'analyses', 'line_stimuli', 'Line_Stimuli.m'));
run(fullfile(rootDir, 'analyses', 'line_stimuli', 'Analyse_Line_Stimuli.m'));
run(fullfile(rootDir, 'analyses', 'line_stimuli', 'Attention_Line_Stimuli.m'));

% LATENCY
% run(fullfile(rootDir, 'analyses', 'latency', 'Latency_analysis_SVM.m'));
% run(fullfile(rootDir, 'analyses', 'latency', 'Latency_analysis_color.m'));
% run(fullfile(rootDir, 'analyses', 'latency', 'Latency_analysis_color_integrated.m'));
% run(fullfile(rootDir, 'analyses', 'latency', 'Latency_analysis_execSVM_fast.m'));

% DECODING
% run(fullfile(rootDir, 'analyses', 'decoding', 'Color_SVM.m'));
% run(fullfile(rootDir, 'analyses', 'decoding', 'Color_decoding_V1.m'));
% run(fullfile(rootDir, 'analyses', 'decoding', 'time_resolved_SVM.m'));

% PSTH
% run(fullfile(rootDir, 'analyses', 'psth', 'PSTH_colorPref_V1.m'));
% run(fullfile(rootDir, 'analyses', 'psth', 'PSTH_colorPref_V1_grayDist.m'));

% RF OVERLAY
% run(fullfile(rootDir, 'analyses', 'rf_overlay', 'RF_on_stim_V1.m'));
% run(fullfile(rootDir, 'analyses', 'rf_overlay', 'RF_on_stim.m'));
% run(fullfile(rootDir, 'analyses', 'rf_overlay', 'AlignBitmapsRFs.m'));
% run(fullfile(rootDir, 'analyses', 'rf_overlay', 'AllRFs_AlignBitmaps.m'));

% COLOR TUNING
% run(fullfile(rootDir, 'analyses', 'color_tuning', 'ColorTuning_Capsules.m'));
% run(fullfile(rootDir, 'analyses', 'color_tuning', 'color_tuning_cone.m'));

% STIMULUS GENERATION
% run(fullfile(rootDir, 'analyses', 'stimulus_generation', 'Eval_RANDTAB.m'));
% run(fullfile(rootDir, 'analyses', 'stimulus_generation', 'makle_randatab_shapes_compl.m'));

% EXPLORATORY
% run(fullfile(rootDir, 'analyses', 'exploratory', 'Find_stimulus_pairs.m'));
% run(fullfile(rootDir, 'analyses', 'exploratory', 'Compare_Paolos_SNR_Color_SNR.m'));
% run(fullfile(rootDir, 'analyses', 'exploratory', 'test_incl_criteria.m'));

% MOVIE
% run(fullfile(rootDir, 'analyses', 'movie', 'make_activity_movie.m'));
% run(fullfile(rootDir, 'analyses', 'movie', 'make_activity_movie_wrapper_safe.m'));
