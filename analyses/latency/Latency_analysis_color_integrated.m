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
% You can select trials and channels in normMUA by considering the values in each column. I attach down here a snippet from one of my scripts to give you an example of loop used to average the response of each individual channel, for each distance, across all trials in which the array containing that channel was the target array (second column in ALLMAT).


% --- Setup: working directory + housekeeping -----------------------------

% clearvars
% ^ 'clear all' also clears functions from memory (slow). Prefer 'clearvars'.

Monkey = 2; % 1 for Nilson, 2 for Figaro
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
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    m2 = matfile(fullfile(cfg.dataDir, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataDir, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    SNR = m1.SNR;
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
SNR_th = 1.0;                 % activity threshold for good enough SNR 
col_sig_th = 0.05;     % significance threshold of color tuning

SmoothW = 10;
ChansPerArray = 64;
ResponseWindow=[0,500];
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
% One plot per significant channel (col_tune_chans), two curves for the two colors.
% Pages of 4 rows x 3 columns.

PlotWindow = [-200, 500]; % ms
idxPlot = find(tb >= PlotWindow(1) & tb <= PlotWindow(2));
tPlot   = tb(idxPlot);

lab1 = 'ALLMAT(:,8)==2';  % e.g. "green"
lab2 = 'ALLMAT(:,8)==1';  % e.g. "purple"

nPerPage = 12; % 4 rows x 3 cols
nCh      = numel(col_tune_chans);
nPages   = ceil(nCh / nPerPage);

for pg = 1:nPages
    figure('Color','w');
    tl  = tiledlayout(4,3,'Padding','compact','TileSpacing','compact');

    i1 = (pg-1)*nPerPage + 1;
    i2 = min(pg*nPerPage, nCh);
    chansThis = col_tune_chans(i1:i2);

    for ii = 1:numel(chansThis)
        chIdx = chansThis(ii);
        Arr = ceil(chIdx / ChansPerArray);

        f_col1 = (ALLMAT(:,2) == Arr) & (ALLMAT(:,8) == 2);
        f_col2 = (ALLMAT(:,2) == Arr) & (ALLMAT(:,8) == 1);
        Idxcol1 = find(f_col1);
        Idxcol2 = find(f_col2);

        if ReadPerTrial
            sum1 = zeros(numel(idxPlot),1); n1 = zeros(numel(idxPlot),1);
            for j = 1:numel(Idxcol1)
                x = squeeze(m1.normMUA(chIdx, Idxcol1(j), idxPlot));
                x = x(:); ok = ~isnan(x);
                sum1(ok) = sum1(ok) + x(ok); n1(ok) = n1(ok) + 1;
            end
            mu1 = sum1 ./ max(n1,1); mu1(n1==0)=NaN;

            sum2 = zeros(numel(idxPlot),1); n2 = zeros(numel(idxPlot),1);
            for j = 1:numel(Idxcol2)
                x = squeeze(m1.normMUA(chIdx, Idxcol2(j), idxPlot));
                x = x(:); ok = ~isnan(x);
                sum2(ok) = sum2(ok) + x(ok); n2(ok) = n2(ok) + 1;
            end
            mu2 = sum2 ./ max(n2,1); mu2(n2==0)=NaN;
        else
            mu1 = squeeze(nanmean(normMUA(chIdx, f_col1, idxPlot), 2));
            mu2 = squeeze(nanmean(normMUA(chIdx, f_col2, idxPlot), 2));
        end

        nexttile; hold on
        plot(tPlot, mu1, 'LineWidth', 1.2);
        plot(tPlot, mu2, 'LineWidth', 1.2);
        xline(0,'--'); yline(0,':');
        xlim(PlotWindow); box off
        title(sprintf('Ch %d | Arr %d | p=%.2g', chIdx, Arr, p_vals_col(chIdx)), ...
              'Interpreter','none','FontSize',9);
        if ii==1, legend({lab1,lab2},'Location','best','FontSize',8); end
        xlabel('Time (ms)'); ylabel('normMUA');
    end

    title(tl, sprintf('Color-tuned channels: page %d/%d (n=%d)', pg, nPages, nCh), ...
        'FontWeight','bold');
end

%% ===================== Time-resolved SVM: COLOR decoding + exGauss fit =====================
% Trains a fixed linear SVM on ResponseWindow, then tests with a sliding window.
% Fits an exGauss_mod() curve to the time-resolved accuracy and computes t33.

DoColorSVM = 1;
if DoColorSVM && ~isempty(col_tune_chans)

    ARR_COL   = 2;
    COLOR_COL = 8;
    colA = 2;   % color label A in ALLMAT(:,8)
    colB = 1;   % color label B in ALLMAT(:,8)

    SVM_Chans_use = col_tune_chans(:);
    SVM_Arrs_use  = ceil(SVM_Chans_use / ChansPerArray);

    TrainWindow = ResponseWindow;
    tidx_train  = (tb >= TrainWindow(1)) & (tb <= TrainWindow(2));
    tidx2_train = find(tidx_train);
    tRangeTrain = tidx2_train(1):tidx2_train(end);

    nRepeats        = 50;
    nPseudoPerClass = 200;

    winWidthMs = 50;
    stepMs     = 10;
    centers    = tb(1):stepMs:tb(end);
    nT         = numel(centers);

    XcolA = cell(numel(RelArray),1);
    XcolB = cell(numel(RelArray),1);

    fprintf('\n=== COLOR SVM: building TRAIN data [%g %g] ms ===\n', TrainWindow(1), TrainWindow(2));

    for a = 1:numel(RelArray)
        Arr = RelArray(a);
        ch = SVM_Chans_use(SVM_Arrs_use == Arr);
        ch = ch(:);
        if isempty(ch), continue; end

        fA = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colA);
        fB = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colB);
        idxA = find(fA); idxB = find(fB);

        if numel(idxA) < 5 || numel(idxB) < 5, continue; end

        if ReadPerTrial
            ch = sort(unique(ch));
            nCh = numel(ch);
            cuts = [1; find(diff(ch)~=1)+1; numel(ch)+1];

            A_A = NaN(nCh, numel(idxA));
            A_B = NaN(nCh, numel(idxB));

            for j = 1:numel(idxA)
                row0 = 0;
                for rruns = 1:numel(cuts)-1
                    chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                    rr    = (row0+1):(row0+numel(chRun));
                    row0  = row0 + numel(chRun);
                    X  = m1.normMUA(chRun(1):chRun(end), idxA(j), tRangeTrain);
                    X2 = reshape(X, numel(chRun), []);
                    A_A(rr,j) = mean(X2, 2, 'omitnan');
                end
            end

            for j = 1:numel(idxB)
                row0 = 0;
                for rruns = 1:numel(cuts)-1
                    chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                    rr    = (row0+1):(row0+numel(chRun));
                    row0  = row0 + numel(chRun);
                    X  = m1.normMUA(chRun(1):chRun(end), idxB(j), tRangeTrain);
                    X2 = reshape(X, numel(chRun), []);
                    A_B(rr,j) = mean(X2, 2, 'omitnan');
                end
            end
        else
            A_A = squeeze(nanmean(normMUA(ch, fA, tidx_train), 3));
            A_B = squeeze(nanmean(normMUA(ch, fB, tidx_train), 3));
        end

        XcolA{a} = A_A.'; % trials x chans
        XcolB{a} = A_B.'; % trials x chans
    end

    [Xtrain, ytrain] = make_pseudotrials_concat_binary(XcolA, XcolB, nPseudoPerClass);
    MdlFixed = fitcsvm(Xtrain, ytrain, 'KernelFunction','linear', 'Standardize', true);

    accTime = nan(nT,1);
    accSE_t = nan(nT,1);

    fprintf('=== COLOR SVM: time-resolved test (win=%dms, step=%dms) ===\n', winWidthMs, stepMs);

    for ti = 1:nT
        tc = centers(ti);
        tidx_t = (tb >= tc - winWidthMs/2) & (tb <= tc + winWidthMs/2);

        accRep = nan(nRepeats,1);
        for r = 1:nRepeats
            XcolA_t = cell(numel(RelArray),1);
            XcolB_t = cell(numel(RelArray),1);

            for a = 1:numel(RelArray)
                Arr = RelArray(a);
                ch = SVM_Chans_use(SVM_Arrs_use == Arr);
                ch = ch(:);
                if isempty(ch), continue; end

                fA = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colA);
                fB = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colB);

                if ReadPerTrial
                    idxA = find(fA); idxB = find(fB);
                    if numel(idxA)<2 || numel(idxB)<2, continue; end
                    A_A_t = NaN(numel(ch), numel(idxA));
                    A_B_t = NaN(numel(ch), numel(idxB));
                    for j = 1:numel(idxA)
                        x = squeeze(m1.normMUA(ch, idxA(j), tidx_t)); % [nCh x T]
                        A_A_t(:,j) = mean(x, 2, 'omitnan');
                    end
                    for j = 1:numel(idxB)
                        x = squeeze(m1.normMUA(ch, idxB(j), tidx_t));
                        A_B_t(:,j) = mean(x, 2, 'omitnan');
                    end
                else
                    A_A_t = squeeze(nanmean(normMUA(ch, fA, tidx_t), 3));
                    A_B_t = squeeze(nanmean(normMUA(ch, fB, tidx_t), 3));
                end

                XcolA_t{a} = A_A_t.';
                XcolB_t{a} = A_B_t.';
            end

            [Xtest, ytest] = make_pseudotrials_concat_binary(XcolA_t, XcolB_t, nPseudoPerClass);
            yhat = predict(MdlFixed, Xtest);
            accRep(r) = mean(yhat == ytest);
        end

        accTime(ti) = mean(accRep, 'omitnan');
        accSE_t(ti) = std(accRep, 'omitnan') / sqrt(sum(~isnan(accRep)));
    end

    SVM_Color = struct();
    SVM_Color.centers = centers(:);
    SVM_Color.accTime = accTime(:);
    SVM_Color.accSE_t = accSE_t(:);
    SVM_Color.MdlFixed = MdlFixed;
    SVM_Color.trainWindow = TrainWindow;
    SVM_Color.testWinWidthMs = winWidthMs;
    SVM_Color.testStepMs = stepMs;

    assignin('base','SVM_Color',SVM_Color);

    % ---------- Plot data ----------
    figure('Color','w'); hold on
    hData_col = plot(SVM_Color.centers, SVM_Color.accTime, '-', 'LineWidth', 2);
    col_fit = hData_col.Color;
    yline(0.5,'--k'); xline(0,'--k');
    grid on; box off
    xlabel('Time (ms) (window center)')
    ylabel('Color decoding accuracy')
    title(sprintf('Color SVM (fixed) trained on [%g %g] ms; tested with %d ms window', ...
        TrainWindow(1), TrainWindow(2), winWidthMs));

    % ---------- Fit (exGauss) + t33 ----------
    t   = SVM_Color.centers(:);
    acc = SVM_Color.accTime(:);
    y   = acc - 0.5;

    ok   = ~isnan(t) & ~isnan(y);
    tfit = t(ok);
    yfit = y(ok);

    y_smooth = smoothdata(yfit, 'movmean', SmoothW);

    if ~exist('MinSigma','var') || isempty(MinSigma), MinSigma = 5; end
    if ~exist('MaxSigma','var') || isempty(MaxSigma), MaxSigma = 400; end

    modfun = @(p,tt) exGauss_mod(p,tt);

    p0 = [200, 20, 0.01, 0.1, 0];
    lb = [min(tfit), MinSigma, 0, 0, 0];
    ub = [max(tfit), MaxSigma, 10, 1, 0.5];

    opts = optimoptions('lsqcurvefit','Display','off');
    pHat = lsqcurvefit(modfun, p0, tfit, y_smooth, lb, ub, opts);

    yFit_corr = modfun(pHat, tfit);
    accFit    = yFit_corr + 0.5;

    mx  = max(yFit_corr, [], 'omitnan');
    thr = 0.33 * mx;

    idx33 = find(yFit_corr >= thr, 1, 'first');
    if isempty(idx33)
        t33 = NaN;
    elseif idx33 == 1
        t33 = tfit(1);
    else
        t33 = interp1(yFit_corr(idx33-1:idx33), tfit(idx33-1:idx33), thr, 'linear');
    end

    plot(tfit, accFit, '--', 'LineWidth', 3, 'Color', col_fit);
    if isfinite(t33)
        xline(t33, ':k');
        text(t33, 0.52, sprintf('t_{33}=%.0f ms', t33), ...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end

    SVM_Color.fit = struct('pHat',pHat,'tfit',tfit,'accFit',accFit, ...
                           'yFit_corr',yFit_corr,'t33',t33,'thr',thr, ...
                           'SmoothW',SmoothW,'MinSigma',MinSigma,'MaxSigma',MaxSigma);
    assignin('base','SVM_Color',SVM_Color);

    fprintf('COLOR SVM exGauss fit: t33 = %.1f ms\\n', t33);

end

%% ===================== Helper: build concatenated pseudotrials (binary) =====================
function [X, y] = make_pseudotrials_concat_binary(XA_cell, XB_cell, nPseudoPerClass)
% XA_cell / XB_cell: cell arrays, one per array:
%   each cell is [nTrials x nCh] for that array (can be empty)
% Builds nPseudoPerClass pseudotrials for class A and B per array (bootstrap),
% concatenates arrays horizontally -> final X is [2*nPseudoPerClass x sum(nCh_arr)].

    nArr = numel(XA_cell);
    useArr = false(nArr,1);

    for a = 1:nArr
        XA = XA_cell{a}; XB = XB_cell{a};
        if isempty(XA) || isempty(XB), continue; end
        if size(XA,2) ~= size(XB,2), continue; end
        if size(XA,1) < 2 || size(XB,1) < 2, continue; end
        useArr(a) = true;
    end

    if ~any(useArr)
        X = []; y = [];
        return
    end

    X_A_all = [];
    X_B_all = [];

    for a = 1:nArr
        if ~useArr(a), continue; end
        XA = XA_cell{a};
        XB = XB_cell{a};

        nA = size(XA,1);
        nB = size(XB,1);

        PA = nan(nPseudoPerClass, size(XA,2));
        PB = nan(nPseudoPerClass, size(XB,2));

        for k = 1:nPseudoPerClass
            ia = randi(nA, [nA,1]);
            ib = randi(nB, [nB,1]);
            PA(k,:) = mean(XA(ia,:), 1, 'omitnan');
            PB(k,:) = mean(XB(ib,:), 1, 'omitnan');
        end

        X_A_all = [X_A_all, PA]; %#ok<AGROW>
        X_B_all = [X_B_all, PB]; %#ok<AGROW>
    end

    X = [X_A_all; X_B_all];
    y = [ones(size(X_A_all,1),1); -ones(size(X_B_all,1),1)];
end
