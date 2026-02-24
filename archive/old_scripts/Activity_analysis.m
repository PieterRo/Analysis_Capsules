% In both "ObjAtt_lines_normMUA.mat" files you’ll find 3 variables: “SNR", “normMUA" and “lats". You can ignore “lats”. “SNR” is the signal-to-noise, 
% you can use it to include channels. “normMUA” is the normalized MUA, and it’s a 3-dimensional matrix with the following structure: 
% #channels (= 1024) x #trials (depends on the experiment) x #time-course (=700ms, it’s 1Hz downsampled signal).
% 
% In both "ObjAtt_lines_MUA_trials.mat" files you’ll find 3 variables: “tb", “ALLMUA" and “ALLMAT". You should ignore “ALLMUA”, ‘cause it’s the non-normalized 
% MUA, you want to work with the normalized one (normMUA from before). “ tb” tells you, in milliseconds, to what each of the 700 points in the 3rd dimension 
% of “normMUA” correspond, with respect to the onset of the stimulus. Finally “ ALLMAT” contains all the information about each trial, 
% with the following structure: #trials (depends on the experiment) x #columns
% 
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
% You can select trials and channels in normMUA by considering the values in each column. I attach down here a snippet from one of my scripts to give you an example of loop used to average the response of each individual channel, for each distance, across all trials in which the array containing that channel was the target array (second column in ALLMAT).


% --- Setup: working directory + housekeeping -----------------------------

cd '/Users/pieter/Library/CloudStorage/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Analysis'
% ^ Sets the current working folder. Consider avoiding hard-coded paths
%   (see improvements below).

close all
clearvars
% ^ 'clear all' also clears functions from memory (slow). Prefer 'clearvars'.

% --- Load data -----------------------------------------------------------

load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_normMUA.mat')
% Expects variables like:
%   normMUA: [nChans x nTrials x nTime] normalized multi-unit activity
%   SNR:     [nChans x nTrials] (or similar) signal-to-noise measure

load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_MUA_trials.mat')
% Expects variables like:
%   ALLMAT:  [nTrials x nCols] trial metadata
%   tb:      [1 x nTime] timebase in ms (often -200..500 etc.)

% --- Choose condition / subset ------------------------------------------

RelArray = 2;      % which array (RF location / array ID) you want to analyze
RelCones = 6;      % which cone-distance condition you want to analyze - this is 2, 4 or 6
ChansPerArray = 64;

% --- Select trials by condition (attention on RF vs off RF) --------------

f_att = ALLMAT(:,3) == RelCones & ALLMAT(:,2) == RelArray & ALLMAT(:,4) == 2; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
f_unatt = ALLMAT(:,3) == RelCones & ALLMAT(:,2) == RelArray & ALLMAT(:,4) == 1; % unattended trials: same, but attention==1 ("not on RF")


% --- Select channels belonging to the chosen array -----------------------

first_chan = (RelArray-1)*ChansPerArray + 1;  % first channel index in that array
last_chan  = (RelArray)*ChansPerArray;        % last channel index in that array

SNR_Avg = mean(SNR,2);
% Mean SNR across trials per channel.
% NOTE: if SNR contains NaNs, consider mean(SNR,2,'omitnan').

SNR_th = 1.2;                 % threshold for good enough channels
RelChans = SNR_Avg >= SNR_th; % boolean mask for channels that pass threshold

RelChans(1:first_chan-1) = 0;     % exclude channels in arrays before RelArray
RelChans(last_chan+1:end) = 0;    % exclude channels in arrays after RelArray


% --- Compute grand means -------------------------------------------------
% normMUA(RelChans, f_att, :) gives: [nGoodChans x nAttTrials x nTime]
% nanmean(...,[1,2]) averages across channels and trials -> [1 x 1 x nTime]
% squeeze(...) makes it [nTime x 1] or [1 x nTime] depending on dimensions

mean_att   = squeeze(nanmean(normMUA(RelChans,f_att,:),[1,2]));
mean_unatt = squeeze(nanmean(normMUA(RelChans,f_unatt,:),[1,2]));
tb2 = tb(:);
assert(numel(tb) == numel(mean_att),'Length mismatch: tb and mean_att do not match');

% --- Plot (current version: two lines, no labels/time axis) --------------
figure('Color','w'); clf
hold on;

plot(tb, mean_att,   'LineWidth', 2)
plot(tb, mean_unatt, 'LineWidth', 2)

xline(0,'--k','LineWidth',1)    % stimulus/event onset at t = 0 ms
grid on;
box off;

xlabel('Time (ms)');
ylabel('Normalized MUA');
legend({'Attend RF','Unattend RF'}, 'Location','best');

xlim([tb(1) tb(end)]);
title(sprintf('Array %d | Cones %d | SNR %.1f', RelArray, RelCones, SNR_th), ...
      'Interpreter','none');

  
%% === ONE shared attention decoder across the 3 timing conditions (RelCones = 2,4,6) ===
% IMPORTANT:
%  - timing condition is ALLMAT(:,3)  (RelCones) with values {2,4,6}
%  - attention condition is ALLMAT(:,4): 2 = attend RF, 1 = unattend RF
%  - proj will be [nTrials_kept x nTime]  (rows=trials, cols=time)

% this version does diagonal LDA and does not take the correlations into
% account. 


cones_set = [2 4 6];            % the 3 timing conditions you want to include
use_correct_only = true;        % set false if you want to include incorrect trials too

tb = tb(:);

% ---- trial selection: all relevant trials for this array + any cones in {2,4,6} ----
f_ok = true(size(ALLMAT,1),1);
if use_correct_only
    f_ok = (ALLMAT(:,9) == 1);  % column 9 = correctness
end

f_rel = (ALLMAT(:,2) == RelArray) & ismember(ALLMAT(:,3), cones_set) & f_ok;

% attention label (pooled across the three timings)
y_all = nan(size(ALLMAT,1),1);
y_all(f_rel & (ALLMAT(:,4) == 2)) = 1;   % attend RF
y_all(f_rel & (ALLMAT(:,4) == 1)) = 0;   % unattend RF

keep_trials = f_rel & ~isnan(y_all);

% now restrict everything to the kept trials
y = y_all(keep_trials);                   % [nTr x 1]
cones_raw = ALLMAT(keep_trials,3);        % [nTr x 1], values 2/4/6

% map cones values to timing index 1/2/3 for plotting
timing_k = nan(size(cones_raw));
timing_k(cones_raw == 2) = 1;
timing_k(cones_raw == 4) = 2;
timing_k(cones_raw == 6) = 3;

assert(all(ismember(unique(timing_k), [1 2 3])), 'Unexpected cones values in selected trials.');

% ---- extract data (channels already in RelChans) ----
M = normMUA(RelChans, keep_trials, :);    % [nCh x nTr x nTime]
nCh   = size(M,1);
nTr   = size(M,2);
nTime = size(M,3);

assert(nTime == numel(tb), 'Mismatch: size(normMUA,3) (%d) must equal numel(tb) (%d).', nTime, numel(tb));

% ---- baseline normalize per channel/trial (reliability weighting) ----
base_idx = (tb >= -200) & (tb <= 0);      % adjust baseline if needed
fit_idx  = (tb >=  200) & (tb <= 700);    % window used to fit the shared decoder

mu_base = mean(M(:,:,base_idx), 3, 'omitnan');            % [nCh x nTr]
sd_base = std (M(:,:,base_idx), 0, 3, 'omitnan');         % [nCh x nTr]
sd_base(sd_base==0) = NaN;

Mz = (M - mu_base) ./ sd_base;                            % [nCh x nTr x nTime]

% ---- build ONE shared decoder weight vector w (diagonal-LDA / SNR weighting) ----
Xfit = squeeze(mean(Mz(:,:,fit_idx), 3, 'omitnan'))';      % [nTr x nCh]

mu_att   = mean(Xfit(y==1,:), 1, 'omitnan');
mu_unatt = mean(Xfit(y==0,:), 1, 'omitnan');
v        = var (Xfit, 0, 1, 'omitnan');

lambda = 1e-2;
w = ((mu_att - mu_unatt) ./ (v + lambda))';                % [nCh x 1]
w = w / norm(w);

% ---- apply shared decoder over time ----
% proj(trial, time) = w' * Mz(:,trial,time)
proj = squeeze(sum(Mz .* reshape(w,[nCh 1 1]), 1, 'omitnan'));  % [nTr x nTime]

% sanity checks (helpful while debugging)
assert(size(proj,1) == nTr, 'proj rows must equal nTr (kept trials).');
assert(size(proj,2) == nTime, 'proj columns must equal nTime.');

% ---- average within each timing condition and attention condition ----
att_curve  = nan(3, nTime);
un_curve   = nan(3, nTime);
diff_curve = nan(3, nTime);

for k = 1:3
    trk = (timing_k == k);          % [nTr x 1]

    idx_att   = trk & (y == 1);
    idx_unatt = trk & (y == 0);

    att_curve(k,:)  = mean(proj(idx_att,  :), 1, 'omitnan');
    un_curve(k,:)   = mean(proj(idx_unatt,:), 1, 'omitnan');
    diff_curve(k,:) = att_curve(k,:) - un_curve(k,:);
end

%% --- Temporal smoothing of decoder time courses -------------------------

% choose smoothing width (in ms)
smooth_sigma_ms = 30;                 % 15–30 ms is typical for MUA
dt = mean(diff(tb));                  % time step in ms
smooth_sigma = smooth_sigma_ms / dt;  % sigma in samples

% build Gaussian kernel
half_width = ceil(4 * smooth_sigma);
x = -half_width : half_width;
g = exp(-(x.^2) / (2 * smooth_sigma^2));
g = g / sum(g);                       % normalize area to 1

% apply zero-phase smoothing (conv is symmetric here)
diff_curve_s = nan(size(diff_curve));
for k = 1:3
    diff_curve_s(k,:) = conv(diff_curve(k,:), g, 'same');
end

figure('Color','w'); clf
plot(tb, diff_curve_s(1,:), 'LineWidth', 2); hold on
plot(tb, diff_curve_s(2,:), 'LineWidth', 2);
plot(tb, diff_curve_s(3,:), 'LineWidth', 2);
xline(0,'--k','LineWidth',1);
grid on; box off
xlabel('Time (ms)');
ylabel('Shared-decoder attention index (smoothed)');
legend({'RelCones=2','RelCones=4','RelCones=6'}, 'Location','best');
xlim([tb(1) tb(end)]);
title('Shared attention decoder (Gaussian smoothed)', 'Interpreter','none');

  