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

cd(cfg.rootDir)
if Monkey == 1
    m1 = matfile(fullfile(cfg.dataDir, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    SNR_Paolo = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataDir, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataDir, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    SNR_Paolo = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataDir, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 

%     ReadPerTrial = 0; % I can load normMUA at once.   % use these two lines if data from Figaro is loaded at once
%     [NChansGlob,NTrialsGlob,NTimesGlob] = size(normMUA);
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_normMUA.mat')
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_MUA_trials.mat')

end


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
RelCones = [2,4,6];      % cone-distance conditions
NCones = numel(RelCones);
close all

SmoothW = 10;
ChansPerArray = 64;
AttentionWindow=[250,500];
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

% compute significances in the three cone conditions and store in p_vals_global
if ~exist('p_vals_global','var') && ~p_val_local                % if true, looks at the p-value in the actual distance condition
    p_vals_global = nan(NChansGlob,3);
    idxAtt2 = find(idxAtt);
    for ConeCond = 1:3           % 1, 2 or 3
        RelCone = RelCones(ConeCond); 

        for ArLoop=1 : NArrays
            Arr = RelArray(ArLoop);
            f_att = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 2; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
            f_unatt = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 1; % unattended trials: same, but attention==1 ("not on RF")
            attIdx   = find(f_att);
            unattIdx = find(f_unatt);
            
            % --- Select channels belonging to the chosen array -----------------------

            first_chan = (Arr-1)*ChansPerArray + 1;  % first channel index in that array
            last_chan  = (Arr)*ChansPerArray;        % last channel index in that array

            RelChans = ones(1, 1024);        % start with all
            RelChans(1:first_chan-1) = 0;     % exclude channels in arrays before RelArray
            RelChans(last_chan+1:end) = 0;    % exclude channels in arrays after RelArray
            goodGlobal = find(RelChans);      % global channel indices selected this iteration
            nGoodChans = numel(goodGlobal);

            % compute significance of the attention effect of the included channels
            for c = 1:nGoodChans
                chIdx = goodGlobal(c);
                if ReadPerTrial

                    % Compute per-trial mean over time without loading full normMUA
                    Xatt   = NaN(numel(attIdx), 1);
                    for j = 1:numel(attIdx)
                        x = m1.normMUA(chIdx, attIdx(j), idxAtt2);      % [1×1×T]
                        Xatt(j) = mean(x, 3, 'omitnan');
                    end

                    Xunatt = NaN(numel(unattIdx), 1);
                    for j = 1:numel(unattIdx)
                        x = m1.normMUA(chIdx, unattIdx(j), idxAtt2);    % [1×1×T]
                        Xunatt(j) = mean(x, 3, 'omitnan');
                    end
                else 
                    Xatt   = squeeze(nanmean(normMUA(chIdx, f_att,  idxAtt), 3));   % [1 × nAttTrials] or [nAttTrials × 1]
                    Xunatt = squeeze(nanmean(normMUA(chIdx, f_unatt, idxAtt), 3));  % [1 × nUnattTrials] or [nUnattTrials × 1]
                end
                [h, p, ~, st] = ttest2(Xatt(:), Xunatt(:), 'Vartype','unequal');
                p_vals_global(chIdx,ConeCond) = p;
            end
            sum(RelChans)
        end
    end

    for i = 1:size(p_vals_global,1)
        if any(~isnan(p_vals_global(i,:)))
            fprintf('%4d : %8.4f  %8.4f  %8.4f\n', ...
                i, p_vals_global(i,1), p_vals_global(i,2), p_vals_global(i,3));
        end
    end
end


for ConeCond = 1:3           % 1, 2 or 3
    RelCone = RelCones(ConeCond); 

    % Loop through the arrays
    selArray   = [];   % array number (Arr)
    selChanG   = [];   % global channel index (1..NChansTotal)
    p_vals    = [];    % array with p-values
    t_vals   = [];     % array with t-values
    allRelChn = [];    % global indices of the included channels with enough activity

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
        attIdx   = find(f_att);
        unattIdx = find(f_unatt);
        tRange   = find(idxAtt);   % contiguous time range

        for c = 1:nGoodChans
            chIdx = goodGlobal(c);
            if ReadPerTrial
                Xatt = NaN(numel(attIdx),1);
                for j = 1:numel(attIdx)
                    x = m1.normMUA(chIdx, attIdx(j), tRange);   % [1×1×T]
                    Xatt(j) = mean(x, 3, 'omitnan');
                end
                Xunatt = NaN(numel(unattIdx),1);
                for j = 1:numel(unattIdx)
                    x = m1.normMUA(chIdx, unattIdx(j), tRange); % [1×1×T]
                    Xunatt(j) = mean(x, 3, 'omitnan');
                end
            else
                Xatt   = squeeze(nanmean(normMUA(chIdx, f_att,  idxAtt), 3));   % [1 × nAttTrials] or [nAttTrials × 1]
                Xunatt = squeeze(nanmean(normMUA(chIdx, f_unatt, idxAtt), 3));  % [1 × nUnattTrials] or [nUnattTrials × 1]
            end
            [h, p, ~, st] = ttest2(Xatt(:), Xunatt(:), 'Vartype','unequal');
            p_vals   = [p_vals;p];
            t_vals   = [t_vals;st.tstat];
        end
 
        
        sum(RelChans)
        % --- Compute grand means -------------------------------------------------
        % normMUA(RelChans, f_att, :) gives: [nGoodChans x nAttTrials x nTime]
        % nanmean(...,[1,2]) averages across channels and trials -> [1 x 1 x nTime]
        % squeeze(...) makes it [nTime x 1] or [1 x nTime] depending on dimensions

        if nGoodChans>0
            if ReadPerTrial 
                sumA = zeros(nGoodChans, numel(idxAtt), 'double');
                cntA = zeros(nGoodChans, numel(idxAtt), 'double');
                chSorted = goodGlobal(:);
                cuts = [1; find(diff(chSorted)~=1)+1; numel(chSorted)+1]; % Split into contiguous runs so each run is a valid MatFile range
                row0 = 0; 
                for r = 1:numel(cuts)-1
                    chRun = chSorted(cuts(r):cuts(r+1)-1);   % contiguous range values
                    rr    = (row0+1):(row0+numel(chRun));    % rows in sum/cnt for this run
                    row0  = row0 + numel(chRun);
                    for j = 1:numel(attIdx)
                        X = m1.normMUA(chRun(1):chRun(end), attIdx(j), :);  % [nRun×1×nTime]
                        X2 = reshape(X, [numel(chRun) size(X,3)]);         % force [nRun × nTime]
                        sumA(rr,:) = sumA(rr,:) + double(X2);
                        cntA(rr,:) = cntA(rr,:) + double(~isnan(X2));
                    end
                end
                mean_att_tmp = sumA ./ cntA;
                mean_att_tmp(cntA==0) = NaN;

                % --- unattended ---
                sumU = zeros(nGoodChans, numel(idxAtt), 'double');
                cntU = zeros(nGoodChans, numel(idxAtt), 'double');
                row0 = 0; 
                for r = 1:numel(cuts)-1
                    chRun = chSorted(cuts(r):cuts(r+1)-1);   % contiguous range values
                    rr    = (row0+1):(row0+numel(chRun));    % rows in sum/cnt for this run
                    row0  = row0 + numel(chRun);
                    for j = 1:numel(unattIdx)
                        X = m1.normMUA(chRun(1):chRun(end), unattIdx(j), :);  % [nRun×1×nTime]
                        X2 = reshape(X, [numel(chRun) size(X,3)]);         % force [nRun × nTime]
                        sumU(rr,:) = sumU(rr,:) + double(X2);
                        cntU(rr,:) = cntU(rr,:) + double(~isnan(X2));
                    end
                end
                mean_unatt_tmp = sumU ./ cntU;
                mean_unatt_tmp(cntU==0) = NaN;

                % --- concatenate across arrays ---
                if ArLoop == 1
                    mean_att_ch   = mean_att_tmp;
                    mean_unatt_ch = mean_unatt_tmp;
                else
                    mean_att_ch   = [mean_att_ch;   mean_att_tmp];
                    mean_unatt_ch = [mean_unatt_ch; mean_unatt_tmp];
                end
            else
                if ArLoop == 1
                    mean_att_ch   = squeeze(nanmean(normMUA(RelChans,f_att,:),2));
                    mean_unatt_ch = squeeze(nanmean(normMUA(RelChans,f_unatt,:),2));
                else
                    mean_att_ch   = [mean_att_ch; squeeze(nanmean(normMUA(RelChans,f_att,:),2))];
                    mean_unatt_ch = [mean_unatt_ch; squeeze(nanmean(normMUA(RelChans,f_unatt,:),2))];
                end
            end
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

    fig = gcf;   % or better: fig = figure('Color','w'); then tl = tiledlayout(fig,...)

    
    if p_val_local
        infoStr = sprintf([ ...
            'p_{th} = %.3f\n' ...
            'Activity SNR = %.2f\n' ...
            '\\sigma range = [%.1f  %.1f] ms'], ...
            p_th,SNR_th,MinSigma, MaxSigma);
    else
        infoStr = sprintf([ ...
            'p_{th} = [%.3f  %.3f  %.3f]\n' ...
            'Activity SNR = %.2f\n' ...
            '\\sigma range = [%.1f  %.1f] ms'], ...
            p_th_arr(1), p_th_arr(2), p_th_arr(3), ...
            SNR_th,MinSigma, MaxSigma);
        end

    annotation(fig, 'textbox', [0.62 0.65 0.35 0.25], ...
        'String', infoStr, ...
        'FitBoxToText','on', ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.7 0.7 0.7]);
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

    if p_val_local                    % if true, looks at the p-value in the actual distance condition
       sigRelMask = p_vals<p_th;         % these are local, within distance condition 
     else
      sigRelMask = p_vals_global(selChanG,1) < p_th_arr(1) & ... 
          p_vals_global(selChanG,2) < p_th_arr(2) & ...    
          p_vals_global(selChanG,3) < p_th_arr(3);
    end
   
    
    sigChan = find(sigRelMask);
    NSig = numel(sigChan);
    mean_att_sig = squeeze(nanmean(mean_att_ch(sigChan,:),1));
    mean_unatt_sig = squeeze(nanmean(mean_unatt_ch(sigChan,:),1));

    if p_val_local
        rel_chans_per_cone{ConeCond} = selChanG(sigChan);
        rel_arr_per_cone{ConeCond} = selArray(sigChan);
    end
    
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
    lb = [0, MinSigma, 0,   0, 0];        % sigma>0, alpha>=0
    ub = [max(t), MaxSigma,  10,    1,  1];

    opts = optimoptions('lsqcurvefit','Display','iter','MaxFunctionEvaluations',2e4);
    pHat = lsqcurvefit(modfun, p0, t, diff_resp, lb, ub, opts);

    % ----- Plot -----
    yFit = modfun(pHat,t);

    thr = 0.33 * max(yFit);
    idx = find(yFit >= thr, 1, 'first');
    t_33 = t(idx);

    
    nexttile(tl, 3);
    cla
    hold on

    plot(t, smoothdata(diff_resp, 'movmean', SmoothW),    'LineWidth', 2);
    plot(t, yFit, 'LineWidth', 2);
    grid on; box off
    title(sprintf('\\mu=%.1f ms, \\sigma=%.1f ms, 33%% point = %.1f ms', pHat(1), pHat(2), t_33));
end

if ~p_val_local
    SVM_Chans = selChanG(sigChan);
    SVM_Arrs = selArray(sigChan);
end   

% now run the SVM that determines the latency. It creates pseudotrials.  
Latency_analysis_execSVM_fast

% make the pdf file 
make_pdf_clean_v2
