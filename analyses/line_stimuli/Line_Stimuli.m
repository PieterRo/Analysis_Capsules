% 11 Feb 2026 The script generates the stimuli of the line task
% - allows analysis of V1 RFs
% - computes the position of the RFs in pixels and GCs, for every stimulus (N=384)
% - saves the result to a structure Tall_V1, to a file

% the #trials is the same in ALLMAT and in normMUA, so the most important stuff is the meaning of the 10 columns, which is the following:
% 1. #of trial (ignore) 
% 2. #array (the number of the target array of this trial) 
% 3. #cones(= the distance from the first cue in RFs) 
% 4. #attention(2=attention ON the RF) 
% 5. #dir (ignore) 
% 6. #angle (ignore) 
% 7. #array_RF_sz (= the size of the RF) 
% 8. #color (= green or purple, I don’t remember rn which was which)
% 9. #correct (=1 is a correct trial)
% 10. #session (=1 is the first day, and so on)

Monkey = 1; % 1 for Nilson, 2 for Figaro
TabFile = "ObjAtt_lines_monkeyN_20220201_B1";  % made on that date
cfg = config();
Dir = cfg.repoRoot;
LogDir = cfg.logsDir;

cd(Dir);
load(fullfile(LogDir, TabFile));   % load ALLCOORDS - they are from Feb 1st
load(fullfile(LogDir, 'RTAB384.mat')); % load the RANDTAB information for each of 384 stimuli on day 1 and 2 (Jan 2026) 
% Thsi file only contains information about Jan 31st and Feb 1st - see EvalRANDTAB.m

% here we need to evaluate RANDTAB to know 
% - the curve thickness RTAB384(stimNum, 7) and 
% - the color of all curves RTAB384(stimNum, 8): If it is odd (1): colFig = green and colBack = purple.
% now let us have a look at ALLMAT
BmpDir = cfg.repoRoot;
StimDir = cfg.stimDir;

if Monkey == 1
    m1 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    m2 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    m2 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
end

% Look at the relevant arrays in this monkey. Try out array 1
%% ALLMAT(5,:) % stimulus 17 is for array 1 and has 4 GC, it should be on the distractor

RFs
% PlotRFs  % Figure 2 removed (superfluous)

% All V1 RFS ---
RFrange = 1:512;            % V1 RFs 
x_rf = x(RFrange); y_rf = y(RFrange);

% RF analysis
stimNum = 139;
T = rf_table_target_distractor(ALLCOORDS, RTAB384, stimNum, x_rf, y_rf);
idxNaN = isnan(T.x_px) | isnan(T.y_px);
T.assignment(idxNaN) = "NaN";
T = add_arm_projection_metrics(T, ALLCOORDS, stimNum);  % now T has along/proj/perp_signed
T = add_polar_about_s(T, ALLCOORDS, stimNum, 'OnlyBackground', true);
T = add_arc_about_s_edges(T, ALLCOORDS, RTAB384, stimNum);
[T, widthPx] = add_GC_normalization(T, RTAB384, stimNum);

% Example: inspect only object RFs
% Tobj = T(T.assignment=="target" | T.assignment=="distractor", :);
% disp(Tobj(10, {'RF','assignment','along','perp_signed','proj_x','proj_y','along_GC','perp_signed_GC'}));
% plot_stim_with_RFs(ALLCOORDS, RTAB384, stimNum, x_rf(74), y_rf(74));


sites_bg = T.RF(T.assignment == "background");
plot_stim_with_RFs(ALLCOORDS, RTAB384, stimNum, x_rf(sites_bg), y_rf(sites_bg));

% build the table with all RF placements across 384 stimuli, only based on
% ALLCOORDS and RTAB384, right now only for day 1 and 2

% Tall_V1 = build_all_stim_tables(ALLCOORDS, RTAB384, x_rf, y_rf);
% save('Tall_V1_lines.mat', 'Tall_V1');  % I later added _N for Nilson

% after this: find all trials 
