%% 23 dec 2025, try to get a good latency estimate, focusing on Figaro

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

Monkey = 2; % 1 for Nilson, 2 for Figaro


% work
% cd '/Users/pieter/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Analysis'   % werk Mac
% if Monkey == 1
%     load('/Users/pieter/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_normMUA.mat');
%     load('/Users/pieter/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_MUA_trials.mat');
% else
%     load('/Users/pieter/Dropbox/Pieter/data/Figaro/ObjAtt_lines_normMUA.mat');
%     load('/Users/pieter/Dropbox/Pieter/data/Figaro/ObjAtt_lines_MUA_trials.mat');
% end

% home
% cd '/Users/pieter/Library/CloudStorage/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Analysis'
% if Monkey == 1
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_normMUA.mat')
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_MUA_trials.mat')
% else
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_normMUA.mat')
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_MUA_trials.mat')
% end


[NChansGlob,NTrialsGlob,NTimesGlob]=size(normMUA);


% --- Choose condition / subset ------------------------------------------

% Look at the relevant arrays in this monkey
[U, ~, ic] = unique(ALLMAT(:,2));
counts = accumarray(ic, 1);

% Arrays: Figaro: 2, 5 or 7
% Nilson: 1, 2, 3, 4, 5, 8:  9, 10 are in V4 and 101 are special V4; % arrays 4, 5, 8 may miss the cone=6 condition.

if Monkey == 1
    RelArray = [1,2,3,4,5,8];      % which array (RF location / array ID) you want to analyze 
else
    RelArray = [2,5,7];      % which array (RF location / array ID) you want to analyze 
end

NArrays = numel(RelArray);
RelCones = [2,4,6];      % which cone-distance condition you want to analyze - this is 2, 4 or 6
close all

SNR_th = 0.85;                 % threshold for good enough SNR
p_th = 0.05;
SmoothW = 10;
ChansPerArray = 64;

SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
SNR_Spont=[-200,-1];         % window for spontaneous activity
SNR_Windows=[[40,200]',[200,500]']; % windows for checking if channel is sufficiently active in one of these windows
SNR_Windows_tb = SNR_Windows - tb(1); % windows in the tb
N_SNR_Win = numel(SNR_Windows(1,:));

SpontIndx= tb >= SNR_Spont(1) & tb <= SNR_Spont(2);

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

        TrialSpont = squeeze(nanmean(normMUA(first_chan:last_chan, f_relevant, SpontIndx), 3));
        sd_per_channel = nanstd(TrialSpont, 0, 2);   % SD across trials, per channel: [nCh × 1]
        sd_per_channel(sd_per_channel == 0) = NaN; % Avoid divide-by-zero
        SNR_loc_win = -inf(1, ChansPerArray); % We'll store the *best* SNR across windows, per channel

        for WinLoop = 1:N_SNR_Win
            tidx = SNR_Windows_tb(WinLoop,1):SNR_Windows_tb(WinLoop,2);
            WinAct = squeeze(nanmean(normMUA(first_chan:last_chan, f_relevant, tidx), [2,3]))...
              -nanmean(TrialSpont, 2);   % subtract spontaneous activity
            SNR_thiswin = WinAct ./ sd_per_channel;              % SNR for this window: [nCh × 1]
            SNR_loc_win = max(SNR_loc_win, SNR_thiswin.');              % Keep the max across windows (row vector 1×nCh)
        end
        SNR_local(first_chan:last_chan) = SNR_loc_win;
    end
end


for ConeCond = 1:3           % 1, 2 or 3
    RelCone = RelCones(ConeCond); 

    % select the ones with a significant attention effect:
    AttentionWindow=[250,500];
    idxAtt = tb >= AttentionWindow(1) & tb <= AttentionWindow(2);

    % Loop through the arrays
    selArray   = [];   % array number (Arr)
    selChanG   = [];   % global channel index (1..NChansTotal)
    p_vals    = [];    % array with p-values
    t_vals   = [];     % array with t-values

    for ArLoop=1 : NArrays
        Arr = RelArray(ArLoop);
        f_att = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 2; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
        f_unatt = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 1; % unattended trials: same, but attention==1 ("not on RF")


        % --- Select channels belonging to the chosen array -----------------------

        first_chan = (Arr-1)*ChansPerArray + 1;  % first channel index in that array
        last_chan  = (Arr)*ChansPerArray;        % last channel index in that array

        if SelectLocalSNR         % if set to 1, it will use the response in the object attention task, not the global SNR
            SNR_Avg = SNR_local;
        else
        %    SNR_Avg = mean(SNR,2);         % global SNR value
            SNR_Avg = SNR(:,ConeCond);      % select channels based on this cone condition
        end
        RelChans = SNR_Avg >= SNR_th; % boolean mask for channels that pass threshold
        RelChans(1:first_chan-1) = 0;     % exclude channels in arrays before RelArray
        RelChans(last_chan+1:end) = 0;    % exclude channels in arrays after RelArray
        goodGlobal = find(RelChans);      % global channel indices selected this iteration
        nGoodChans = numel(goodGlobal);

        % append to your running list
        selArray   = [selArray;   repmat(Arr, numel(goodGlobal), 1)];
        selChanG   = [selChanG;   goodGlobal(:)];

        % test significance of the attention effect of the included channels
        for c = 1:nGoodChans
            chIdx = goodGlobal(c);
            Xatt   = squeeze(nanmean(normMUA(chIdx, f_att,  idxAtt), 3));   % [1 × nAttTrials] or [nAttTrials × 1]
            Xunatt = squeeze(nanmean(normMUA(chIdx, f_unatt, idxAtt), 3));  % [1 × nUnattTrials] or [nUnattTrials × 1]
            [h, p, ~, st] = ttest2(Xatt(:), Xunatt(:), 'Vartype','unequal');
            p_vals   = [p_vals;p];
            t_vals   = [t_vals;st.tstat];
        end

        sum(RelChans)
        % --- Compute grand means -------------------------------------------------
        % normMUA(RelChans, f_att, :) gives: [nGoodChans x nAttTrials x nTime]
        % nanmean(...,[1,2]) averages across channels and trials -> [1 x 1 x nTime]
        % squeeze(...) makes it [nTime x 1] or [1 x nTime] depending on dimensions

        if ArLoop == 1
            mean_att_ch   = squeeze(nanmean(normMUA(RelChans,f_att,:),2));
            mean_unatt_ch = squeeze(nanmean(normMUA(RelChans,f_unatt,:),2));
        else
            mean_att_ch   = [mean_att_ch; squeeze(nanmean(normMUA(RelChans,f_att,:),2))];
            mean_unatt_ch = [mean_unatt_ch; squeeze(nanmean(normMUA(RelChans,f_unatt,:),2))];
        end

    end

    mean_att = squeeze(nanmean(mean_att_ch,1));
    mean_unatt = squeeze(nanmean(mean_unatt_ch,1));

    [NChans,NTimes]=size(mean_att_ch);

    % --- Plot (current version: two lines, no labels/time axis) --------------
    % --- One figure with 3 panels (All chans | Sig chans | Fit) --------------
    fig = figure('Color','w'); clf
    fig.Position = [60,583+320*ConeCond,1200,320];
    tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    nexttile(tl, 1);
    hold on

    plot(tb, smoothdata(mean_att, 'movmean', SmoothW), 'LineWidth', 2);
    plot(tb, smoothdata(mean_unatt, 'movmean', SmoothW), 'LineWidth', 2);
    plot(tb, smoothdata(mean_att-mean_unatt, 'movmean', SmoothW), 'k','LineWidth', 2);

    xline(0,'--k','LineWidth',1);    % stimulus/event onset at t = 0 ms
    grid on;
    box off;

    xlabel('Time (ms)');
    ylabel('Normalized MUA');
    legend({'Attend RF','Unattend RF',sprintf('N = %d',NChans)}, 'Location','best');


    xlim([tb(1) tb(end)]);
    title(sprintf('NChans %d | Cones %d | SNR %.1f', NChans, RelCones(ConeCond), SNR_th), ...
          'Interpreter','none');


    % --- Plot average significant channels --------------

    sigChan = p_vals<p_th;
    NSig = sum(sigChan);
    mean_att_sig = squeeze(nanmean(mean_att_ch(sigChan,:),1));
    mean_unatt_sig = squeeze(nanmean(mean_unatt_ch(sigChan,:),1));

    nexttile(tl, 2);
    cla
    hold on

    plot(tb, smoothdata(mean_att_sig, 'movmean', SmoothW),   'LineWidth', 2);
    plot(tb, smoothdata(mean_unatt_sig, 'movmean', SmoothW),   'LineWidth', 2);
    diff_resp = mean_att_sig-mean_unatt_sig;
    plot(tb, smoothdata(diff_resp, 'movmean', SmoothW), 'k','LineWidth', 2);

    xline(0,'--k','LineWidth',1);    % stimulus/event onset at t = 0 ms
    grid on;
    box off;

    xlabel('Time (ms)');
    ylabel('Normalized MUA SIG only');
    legend({'Attend RF','Unattend RF'}, 'Location','best');

    xlim([tb(1) tb(end)]);
    title(sprintf('NChans %d | Cones %d | SNR %.1f', NSig, RelCones(ConeCond), SNR_th), ...
          'Interpreter','none');

    % fit a curve
    % tb: 700-point time base (ms), e.g. 0:699 or -100:599 etc.
    t = tb(:);
    diff_resp = diff_resp(:);

    % ----- Model: dissipating + non-dissipating cumulative Gaussians + baseline -----
    modfun = @(p,t) exGauss_mod(p,t);

    % Parameters p = [mu, sigma, alpha, c, d]
    p0 = [200, 20, 0.01,  0.1,  0.1];                 % initial guess (edit to your data)
    lb = [min(t), 1e-3, 0,   0, 0];        % sigma>0, alpha>=0
    ub = [max(t), 300,  10,    1,  1];

    opts = optimoptions('lsqcurvefit','Display','iter','MaxFunctionEvaluations',2e4);
    pHat = lsqcurvefit(modfun, p0, t, diff_resp, lb, ub, opts);

    % ----- Plot -----
    yFit = modfun(pHat,t);

    nexttile(tl, 3);
    cla
    hold on

    plot(t, smoothdata(diff_resp, 'movmean', SmoothW),    'LineWidth', 2);
    plot(t, yFit, 'LineWidth', 2);
    grid on; box off
    title(sprintf('\\mu=%.1f ms, \\sigma=%.1f ms, \\alpha=%.4f', pHat(1), pHat(2), pHat(3)));
end
  

% SVM
% steps: 
% 1) prepare data matrix, by making speudotrials combining data sets. 
% do SVM





% plot individual channels, array, channel and p-value for the attention effect
 
plotsPerPage = 20;
plotCount=0;
nRows = 5;
nCols = 4;
SmoothW = 15;

% always good to check the underlying data. Here we plot only the included channels, with significant modulation

for ChLoop = 1:NChans
    if sigChan(ChLoop)
        % Start a new "page" every 20 channels
        if mod(plotCount, plotsPerPage) == 0
            figure('Color','w');
            t = tiledlayout(nRows, nCols, ...
                'TileSpacing','compact', ...
                'Padding','compact');
        end
        % increment only for plotted channels
        plotCount = plotCount + 1;

        nexttile
        hold on

        plot(tb, smoothdata(mean_att_ch(ChLoop,:), 'movmean', SmoothW), 'LineWidth', 2);
        plot(tb, smoothdata(mean_unatt_ch(ChLoop,:), 'movmean', SmoothW), 'LineWidth', 2);

        xline(0,'--k','LineWidth',1);
        grid on
        box off
        title(sprintf('Ar %d Ch %d p=%.3f', selArray(ChLoop), selChanG(ChLoop), p_vals(ChLoop)));
    end
end
  

