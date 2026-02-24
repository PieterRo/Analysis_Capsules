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

%% --- Average attended vs unattended across arrays 2, 5, and 7 (no decoding) ----
% The three attention-timing conditions are set by RelCones (ALLMAT(:,3)).
% Set RelCones = 2, 4, or 6 at the top of the script to view the three timings.

RelCones = 4;                 % timing condition: 2, 4, or 6
arrays_set = [2 5 7];         % arrays to include
ChansPerArray = 64;

% Channel inclusion threshold
SNR_th = 2;

% Optional: keep only correct trials (ALLMAT col 9)
use_correct_only = true;

% --- Time base ---
tb = tb(:);
nTime = numel(tb);
assert(size(normMUA,3) == nTime, 'Mismatch: size(normMUA,3) (%d) must equal numel(tb) (%d).', size(normMUA,3), nTime);

% --- SNR per channel (averaged across trials) ---
SNR_Avg = mean(SNR, 2, 'omitnan');   % [nCh x 1]

% --- Collect per-channel mean time courses (equal weight per included channel) ---
% We first average across trials within each channel, then average across channels.
% This avoids giving extra weight to arrays with more channels/trials.

att_allCh   = [];   % will become [nChTotal_att x nTime]
unatt_allCh = [];   % will become [nChTotal_unatt x nTime]

% For reporting
nCh_used   = nan(numel(arrays_set),1);
nAttTrials = nan(numel(arrays_set),1);
nUnTrials  = nan(numel(arrays_set),1);

for ai = 1:numel(arrays_set)
    RelArray = arrays_set(ai);

    % --- Channel mask for this array ---
    first_chan = (RelArray-1)*ChansPerArray + 1;
    last_chan  = (RelArray)*ChansPerArray;

    RelChans = (SNR_Avg >= SNR_th);
    RelChans(1:first_chan-1) = 0;
    RelChans(last_chan+1:end) = 0;

    nCh_used(ai) = sum(RelChans);

    % --- Trial selection for this array and timing condition ---
    f_ok = true(size(ALLMAT,1),1);
    if use_correct_only
        f_ok = (ALLMAT(:,9) == 1);
    end

    f_att   = (ALLMAT(:,3) == RelCones) & (ALLMAT(:,2) == RelArray) & (ALLMAT(:,4) == 2) & f_ok;
    f_unatt = (ALLMAT(:,3) == RelCones) & (ALLMAT(:,2) == RelArray) & (ALLMAT(:,4) == 1) & f_ok;

    nAttTrials(ai) = sum(f_att);
    nUnTrials(ai)  = sum(f_unatt);

    % --- Per-channel mean time courses (average across trials) ---
    % Result: [nChUsed x nTime]
    if nCh_used(ai) > 0 && nAttTrials(ai) > 0
        att_ch_tc = squeeze(mean(normMUA(RelChans, f_att, :), 2, 'omitnan'));
        att_allCh = [att_allCh; att_ch_tc];
    end

    if nCh_used(ai) > 0 && nUnTrials(ai) > 0
        unatt_ch_tc = squeeze(mean(normMUA(RelChans, f_unatt, :), 2, 'omitnan'));
        unatt_allCh = [unatt_allCh; unatt_ch_tc];
    end
end

% --- Combine across channels (equal weight per included channel) ---
% Average across channels, keep time
mean_att   = mean(att_allCh,   1, 'omitnan')';
mean_unatt = mean(unatt_allCh, 1, 'omitnan')';

% SEM across channels (note: attended and unattended can have different Nch)
nCh_att_eff   = size(att_allCh,   1);
nCh_unatt_eff = size(unatt_allCh, 1);

sem_att   = std(att_allCh,   0, 1, 'omitnan')' ./ sqrt(max(nCh_att_eff,1));
sem_unatt = std(unatt_allCh, 0, 1, 'omitnan')' ./ sqrt(max(nCh_unatt_eff,1));

% --- Report trial/channel counts per array ---
fprintf('\nAveraging across arrays [%s] at RelCones=%d (timing)\n', num2str(arrays_set), RelCones);
for ai = 1:numel(arrays_set)
    fprintf('  Array %d: Nch=%d (SNR>=%.1f), Natt=%d, Nunatt=%d\n', ...
        arrays_set(ai), nCh_used(ai), SNR_th, nAttTrials(ai), nUnTrials(ai));
end

% --- Plot ---
figure('Color','w'); clf
hold on

% Shaded bands (SEM across arrays)
fill([tb; flipud(tb)], [mean_att-sem_att; flipud(mean_att+sem_att)], ...
     0.8*[1 1 1], 'EdgeColor','none', 'FaceAlpha', 0.25);
plot(tb, mean_att, 'LineWidth', 2);

fill([tb; flipud(tb)], [mean_unatt-sem_unatt; flipud(mean_unatt+sem_unatt)], ...
     0.8*[1 1 1], 'EdgeColor','none', 'FaceAlpha', 0.25);
plot(tb, mean_unatt, 'LineWidth', 2);

xline(0,'--k','LineWidth',1);
grid on; box off
xlabel('Time (ms)')
ylabel('Normalized MUA (mean across arrays)')
legend({'Attend ± SEM','Attend','Unattend ± SEM','Unattend'}, 'Location','best')
xlim([tb(1) tb(end)])
title(sprintf('Arrays [%s] | Cones %d | SNR >= %.1f', num2str(arrays_set), RelCones, SNR_th), ...
      'Interpreter','none');
