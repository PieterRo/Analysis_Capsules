%% 29 dec 2025 looking at the effect of colour

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

% --- Setup: working directory + housekeeping -----------------------------

% clearvars
% ^ 'clear all' also clears functions from memory (slow). Prefer 'clearvars'.

Monkey = 2; % 1 for Nilson, 2 for Figaro
cfg = config();


% work
cd(cfg.repoRoot)   % work Mac
if Monkey == 1
    load(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    load(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
else
    load(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    load(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
end

% home

cd(cfg.repoRoot)
if Monkey == 1
    m1 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataRoot, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 

%     ReadPerTrial = 0; % I can load normMUA at once.   % use these two lines if data from Figaro is loaded at once
%     [NChansGlob,NTrialsGlob,NTimesGlob] = size(normMUA);
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_normMUA.mat')
%     load('/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data/Figaro/ObjAtt_lines_MUA_trials.mat')

end

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
SelectLocalSNR = 1;    % look at the SNR in the  curve-tracing task
SNR_th = 0.65;                 % activity threshold for good enough SNR 
col_sig_th = 0.1;     % significance threshold of color tuning

SmoothW = 10;
ChansPerArray = 64;
ResponseWindow=[100,500];
idxResp = tb >= ResponseWindow(1) & tb <= ResponseWindow(2);
idxResp2 = find(idxResp);

SNR_Spont=[-200,-1];         % window for spontaneous activity
SpontIndx= tb >= SNR_Spont(1) & tb <= SNR_Spont(2);

SNR_Windows=[[40,200]',[200,500]']; % windows for checking if channel is sufficiently active in one of these windows
SNR_Windows_tb = SNR_Windows - tb(1); % windows in the tb
N_SNR_Win = numel(SNR_Windows(1,:));

if Monkey == 1
    load(fullfile(cfg.matDir, 'Nilson_SNR_local.mat'));
else
    load(fullfile(cfg.matDir, 'Figaro_SNR_local.mat'));
end

%Computes the local SNR, i.e. the reliability of the visual response if necessary
ComputeLocalSNR

Relly = 6;

% compute the color response collapsed across the three cone conditions
p_vals_col = nan(NChansGlob,1);
for ArLoop=1 : NArrays
    Arr = RelArray(ArLoop);
    Arr
    f_col1 = ALLMAT(:,2) == Arr & ALLMAT(:,8) == 2; % trials with color 1
    f_col2 = ALLMAT(:,2) == Arr & ALLMAT(:,8) == 1; % trials with color 2
    Idxcol1 = find(f_col1);
    Idxcol2 = find(f_col2);

    % --- Select channels belonging to the chosen array -----------------------
    if SelectLocalSNR         % if set to 1, it will use the response in the object attention task, not the global SNR
        SNR_Avg = SNR_local;
    else
    %    SNR_Avg = mean(SNR,2);         % global SNR value
        SNR_Avg = SNR(:,ConeCond);      % select channels based on this cone condition
    end

    RelChans = SNR_Avg >= SNR_th; % boolean mask for channels that pass threshold
    first_chan = (Arr-1)*ChansPerArray + 1;  % first channel index in that array
    last_chan  = (Arr)*ChansPerArray;        % last channel index in that array
    RelChans(1:first_chan-1) = 0;     % exclude channels in arrays before RelArray
    RelChans(last_chan+1:end) = 0;    % exclude channels in arrays after RelArray
    goodGlobal = find(RelChans);      % global channel indices selected this iteration
    nGoodChans = numel(goodGlobal);

    % compute significance of the color effect of the included channels
    for c = 1:nGoodChans
        chIdx = goodGlobal(c);
        if ReadPerTrial

            % Compute per-trial mean over time without loading full normMUA
            Xcol1   = NaN(numel(Idxcol1), 1);
            for j = 1:numel(Idxcol1)
                x = m1.normMUA(chIdx, Idxcol1(j), idxResp2);      % [1×1×T]
                Xcol1(j) = mean(x, 3, 'omitnan');
            end

            Xcol2 = NaN(numel(Idxcol2), 1);
            for j = 1:numel(Idxcol2)
                x = m1.normMUA(chIdx, Idxcol2(j), idxResp2);    % [1×1×T]
                Xcol2(j) = mean(x, 3, 'omitnan');
            end
        else 
            Xcol1   = squeeze(nanmean(normMUA(chIdx, f_col1,  idxResp), 3));   % [1 × nAttTrials] or [nAttTrials × 1]
            Xcol2 = squeeze(nanmean(normMUA(chIdx, f_col2, idxResp), 3));  % [1 × nUnattTrials] or [nUnattTrials × 1]
        end
        [h, p, ~, st] = ttest2(Xcol1(:), Xcol2(:), 'Vartype','unequal');
        p_vals_col(chIdx) = p;
    end
    sum(RelChans)
end

col_tune_chans = find(p_vals_col<col_sig_th);

%% ---------- Plot time courses for significant color-tuned channels ----------
% Plots per-channel mean normMUA vs time for the two color conditions (ALLMAT(:,8))
% from -200 to 500 ms, in pages of 4x3 subplots.

% --- time window to plot ---
PlotWindow = [-200, 500]; % ms
idxPlot = find(tb >= PlotWindow(1) & tb <= PlotWindow(2));
tPlot   = tb(idxPlot);

% --- labels (edit if you remember which is green/purple) ---
lab1 = 'ALLMAT(:,8)==2';  % e.g., "green"
lab2 = 'ALLMAT(:,8)==1';  % e.g., "purple"

% --- paging / layout ---
nPerPage = 15; % 5 rows x 3 cols
nCh      = numel(col_tune_chans);
nPages   = ceil(nCh / nPerPage);

for pg = 1:nPages
    fig = figure('Color','w');
    tl  = tiledlayout(5,3,'Padding','compact','TileSpacing','compact');

    % channels on this page
    i1 = (pg-1)*nPerPage + 1;
    i2 = min(pg*nPerPage, nCh);
    chansThis = col_tune_chans(i1:i2);

    for ii = 1:numel(chansThis)
        chIdx = chansThis(ii);

        % infer array from channel index (64 chans per array as used above)
        Arr = ceil(chIdx / ChansPerArray);

        % trial indices for this array + the two color conditions
        f_col1 = (ALLMAT(:,2) == Arr) & (ALLMAT(:,8) == 2);
        f_col2 = (ALLMAT(:,2) == Arr) & (ALLMAT(:,8) == 1);
        Idxcol1 = find(f_col1);
        Idxcol2 = find(f_col2);

        % compute mean time course for each color
        if ReadPerTrial
            % --- streaming mean that avoids loading huge arrays ---
            sum1 = zeros(numel(idxPlot),1);  n1 = zeros(numel(idxPlot),1);
            for j = 1:numel(Idxcol1)
                x = squeeze(m1.normMUA(chIdx, Idxcol1(j), idxPlot)); % [T x 1] after squeeze
                x = x(:);
                ok = ~isnan(x);
                sum1(ok) = sum1(ok) + x(ok);
                n1(ok)   = n1(ok)   + 1;
            end
            mu1 = sum1 ./ max(n1,1);
            mu1(n1==0) = NaN;

            sum2 = zeros(numel(idxPlot),1);  n2 = zeros(numel(idxPlot),1);
            for j = 1:numel(Idxcol2)
                x = squeeze(m1.normMUA(chIdx, Idxcol2(j), idxPlot));
                x = x(:);
                ok = ~isnan(x);
                sum2(ok) = sum2(ok) + x(ok);
                n2(ok)   = n2(ok)   + 1;
            end
            mu2 = sum2 ./ max(n2,1);
            mu2(n2==0) = NaN;
        else
            % --- can load directly if normMUA is in memory ---
            mu1 = squeeze(nanmean(normMUA(chIdx, f_col1, idxPlot), 2)); % [T x 1]
            mu2 = squeeze(nanmean(normMUA(chIdx, f_col2, idxPlot), 2)); % [T x 1]
        end

        % --- plot ---
        nexttile;
        hold on;
        plot(tPlot, mu1, 'LineWidth', 1.2);
        plot(tPlot, mu2, 'LineWidth', 1.2);
        xline(0, '--'); % stimulus onset
        yline(0, ':');  % baseline reference (optional)

        xlim(PlotWindow);
        box off;

        % light annotation
        ttl = sprintf('Ch %d | Arr %d | p=%.2g', chIdx, Arr, p_vals_col(chIdx));
        title(ttl, 'Interpreter','none', 'FontSize', 9);

        if ii == 1
            legend({lab1, lab2}, 'Location','best', 'FontSize', 8);
        end
        xlabel('Time (ms)');
        ylabel('normMUA');
    end

    % page title
    title(tl, sprintf('Color-tuned channels: page %d/%d (n=%d)', pg, nPages, nCh), ...
        'FontWeight','bold');
end


% Color_SVM

