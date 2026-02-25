%% 23 dec 2025, try to get a good latency estimate, first focusing on Figaro

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

% clearvars
% ^ 'clear all' also clears functions from memory (slow). Prefer 'clearvars'.

Monkey = 1; % 1 for Nilson, 2 for Figaro
SNR_th = 1;                 % activity threshold for good enough SNR 
% 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.075;
% 
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR

cfg = config();
cd(cfg.repoRoot)
if Monkey == 1
    m1 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    SNR_Paolo = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    SNR_Paolo = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
end


% --- Choose condition / subset ------------------------------------------

% Look at the relevant arrays in this monkey
[U, ~, ic] = unique(ALLMAT(:,2));
counts = accumarray(ic, 1);

if Monkey == 1
    RelArray = [2];      % focus on array 2
%    RelArray = [1,2,3,4,5,8];      % which array (RF location / array ID) you want to analyze 
else
    RelArray = [2,5,7];      % which array (RF location / array ID) you want to analyze 
end

NArrays = numel(RelArray);
RelCones = [2,4,6];      % cone-distance conditions
NCones = numel(RelCones);
close all

SmoothW = 10;
ChansPerArray = 64;
AttentionWindow=[300,500];      % was previously 250:500
idxAtt = tb >= AttentionWindow(1) & tb <= AttentionWindow(2);
SNR_Spont=[-200,-1];         % window for spontaneous activity
SpontIndx= tb >= SNR_Spont(1) & tb <= SNR_Spont(2);

SNR_Windows=[[40,200]',[200,500]']; % windows for checking if channel is sufficiently active in one of these windows
SNR_Windows_tb = SNR_Windows - tb(1); % windows in the tb
N_SNR_Win = numel(SNR_Windows(1,:));

if p_val_local
    rel_chans_per_cone = cell(NCones,1);
    rel_arr_per_cone = cell(NCones,1);
end

% look for SNR in the data if required for the analysis
if SelectLocalSNR && ~exist('SNR_local','var')
    SNR_local = nan(1, NChansGlob);
    fprintf('Computing the SNRs of all channels\n');

    for ArLoop = 1:NArrays
        Arr = RelArray(ArLoop);
        fprintf('Processing array %d\n', Arr);
        f_relevant = (ALLMAT(:,2) == Arr);

        first_chan = (Arr-1)*ChansPerArray + 1;
        last_chan  = Arr*ChansPerArray;
        
        if ReadPerTrial
            chIdx = first_chan:last_chan;
            fIdx  = find(f_relevant);      % ensure row vector
            fIdx = fIdx(:)';
            tIdx  = find(SpontIndx);       % time indices (dim 3)
            tIdx = tIdx(:)';

            TrialSpont = NaN(numel(chIdx), numel(fIdx));

            for j = 1:numel(fIdx)
                X = m1.normMUA(chIdx, fIdx(j), tIdx);    % loads: channels × 1 × time
                TrialSpont(:,j) = squeeze(mean(X, 3, 'omitnan'));
            end
        else
            TrialSpont = squeeze(nanmean(normMUA(first_chan:last_chan, f_relevant, SpontIndx), 3));
        end

        sd_per_channel = nanstd(TrialSpont, 0, 2);   % SD across trials, per channel: [nCh × 1]
        sd_per_channel(sd_per_channel == 0) = NaN; % Avoid divide-by-zero
        SNR_loc_win = -inf(1, ChansPerArray); % We'll store the *best* SNR across windows, per channel
        if ReadPerTrial
            chIdx = first_chan:last_chan;
            fIdx  = find(f_relevant);      % ensure row vector
            fIdx = fIdx(:)';

            for WinLoop = 1:N_SNR_Win
                tidx = SNR_Windows_tb(WinLoop,1):SNR_Windows_tb(WinLoop,2);
                sumCh = zeros(numel(chIdx),1,'double');
                cntCh = zeros(numel(chIdx),1,'double');
                for j = 1:numel(fIdx)
                    X = m1.normMUA(chIdx, fIdx(j), tidx);          % [nCh × 1 × nT], tidx is contiguous
                    mask = ~isnan(X);
                    sumCh = sumCh + squeeze(sum(double(X).*double(mask), 3));  % sum over time
                    cntCh = cntCh + squeeze(sum(double(mask), 3));             % count over time
                end

                WinMean = sumCh ./ cntCh;          % mean over (trial,time)
                WinMean(cntCh==0) = NaN;

                WinAct = WinMean - nanmean(TrialSpont, 2);  % subtract spontaneous activity
                SNR_thiswin = WinAct ./ sd_per_channel;              % SNR for this window: [nCh × 1]
                SNR_loc_win = max(SNR_loc_win, SNR_thiswin.');              % Keep the max across windows (row vector 1×nCh)
            end
        else
            for WinLoop = 1:N_SNR_Win
                tidx = SNR_Windows_tb(WinLoop,1):SNR_Windows_tb(WinLoop,2);
                WinAct = squeeze(nanmean(normMUA(first_chan:last_chan, f_relevant, tidx), [2,3]))...
                  -nanmean(TrialSpont, 2);   % subtract spontaneous activity
                SNR_thiswin = WinAct ./ sd_per_channel;              % SNR for this window: [nCh × 1]
                SNR_loc_win = max(SNR_loc_win, SNR_thiswin.');              % Keep the max across windows (row vector 1×nCh)
            end
        end
        SNR_local(first_chan:last_chan) = SNR_loc_win;
    end
end

TMP_AttentionPSTH; % makes the PSTHs, GC=2, channels 99,100
% Xatt: 51 values, channel 100 - trials listed in attIdx
% Xunatt: 50 values, channel 100 - trials listed in unattIdx
ALLMAT(attIdx,1);  % these are the stimuli for the attended conditions
ALLMAT(unattIdx,1); % these are the stimuli for the unattended condition



% Stimulus numbers per trial
stimAtt   = ALLMAT(attIdx,1);
stimUnatt = ALLMAT(unattIdx,1);

uAtt   = unique(stimAtt);
uUnatt = unique(stimUnatt);

fprintf('----- ATTENDED -----\n');
attStimMeans = zeros(numel(uAtt),1);

for i = 1:numel(uAtt)
    s = uAtt(i);
    idx = (stimAtt == s);
    attStimMeans(i) = mean(Xatt(idx));
    fprintf('Stim %3d: mean=%.4f  nTrials=%d\n', ...
        s, attStimMeans(i), sum(idx));
end

fprintf('\n----- UNATTENDED -----\n');
unattStimMeans = zeros(numel(uUnatt),1);

for i = 1:numel(uUnatt)
    s = uUnatt(i);
    idx = (stimUnatt == s);
    unattStimMeans(i) = mean(Xunatt(idx));
    fprintf('Stim %3d: mean=%.4f  nTrials=%d\n', ...
        s, unattStimMeans(i), sum(idx));
end

% Grand means across stimulus means (stimulus-weighted)
grandAtt   = mean(attStimMeans);
grandUnatt = mean(unattStimMeans);

fprintf('\n===== SUMMARY =====\n');
fprintf('Mean across attended stimuli:   %.4f\n', grandAtt);
fprintf('Mean across unattended stimuli: %.4f\n', grandUnatt);
fprintf('Difference (Att - Unatt): %.4f\n', grandAtt - grandUnatt);




site = 100;
timeIdx = 3;
stims = 49:64;

fprintf('Stim  Tall.assignment   meanAct(300-500ms)\n');

for stim = stims
    Ttab = Tall_V1(stim).T;
    a = Ttab.assignment(site);

    if iscell(a), a=a{1}; end
    if isstring(a), a=char(a); end
    if iscategorical(a), a=char(a); end

    val = R.meanAct(site, stim, timeIdx);

    fprintf('%3d   %12s   %.4f\n', stim, a, val);
end



