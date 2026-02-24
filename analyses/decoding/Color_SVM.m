%% ===================== Time-resolved SVM: COLOR decoding =====================
% Requires from above script:
%   tb, ALLMAT, RelArray, ChansPerArray, ResponseWindow, ReadPerTrial, m1 (if ReadPerTrial)
%   col_tune_chans (global channel indices)

if ~exist('col_tune_chans','var') || isempty(col_tune_chans)
    warning('col_tune_chans empty/missing. Skipping COLOR SVM.');
else
    ARR_COL   = 2;
    COLOR_COL = 8;

    % Define labels in ALLMAT(:,8)
    colA = 2;   % e.g. green
    colB = 1;   % e.g. purple

    % Which channels to use for decoding
    SVM_Chans_use = col_tune_chans(:);

    % Map channels to arrays (assumes contiguous blocks of 64 channels per array)
    SVM_Arrs_use  = ceil(SVM_Chans_use / ChansPerArray);

    % ---- Training window (collapsed in time) ----
    TrainWindow = ResponseWindow;  % e.g. [0 500]
    tidx_train  = (tb >= TrainWindow(1)) & (tb <= TrainWindow(2));
    tidx2_train = find(tidx_train);
    tRangeTrain = tidx2_train(1):tidx2_train(end); % contiguous for faster reading

    % ---- Pseudotrial settings ----
    nRepeats        = 50;    % repeats for time-resolved evaluation
    nPseudoPerClass = 200;   % pseudotrials per class (per repeat)
    nDrawPerPseudo = 25;
    
    % ---- Time-resolved test window ----
    winWidthMs = 50;
    stepMs     = 10;
    centers    = tb(1):stepMs:tb(end);
    nT         = numel(centers);

    % ---- Build per-array training feature matrices ----
    XcolA = cell(numel(RelArray),1);
    XcolB = cell(numel(RelArray),1);

    fprintf('\n=== COLOR SVM: building TRAIN data [%g %g] ms ===\n', TrainWindow(1), TrainWindow(2));

    for a = 1:numel(RelArray)
        Arr = RelArray(a);

        ch = SVM_Chans_use(SVM_Arrs_use == Arr);
        ch = ch(:);

        if isempty(ch)
            warning('No SVM channels for Arr=%d. Skipping array.', Arr);
            continue;
        end

        fA = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colA);
        fB = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colB);

        idxA = find(fA);
        idxB = find(fB);

        if numel(idxA) < 5 || numel(idxB) < 5
            warning('Too few trials for Arr=%d (nA=%d, nB=%d). Skipping array.', Arr, numel(idxA), numel(idxB));
            continue;
        end

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

                    X  = m1.normMUA(chRun(1):chRun(end), idxA(j), tRangeTrain); % [nRun×1×nT]
                    X2 = reshape(X, numel(chRun), []);                          % [nRun×nT]
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
            % normMUA: [channels x trials x time]
            A_A = squeeze(nanmean(normMUA(ch, fA, tidx_train), 3)); % [nCh x nTrialsA]
            A_B = squeeze(nanmean(normMUA(ch, fB, tidx_train), 3)); % [nCh x nTrialsB]
        end

        XcolA{a} = A_A.'; % trials × channels
        XcolB{a} = A_B.'; % trials × channels
    end

    % ---- Train ONE fixed classifier on the long window ----
    [Xtrain, ytrain] = make_pseudotrials_concat_binary(XcolA, XcolB, nPseudoPerClass, nDrawPerPseudo);

    if isempty(Xtrain)
        warning('COLOR SVM: training set empty. Aborting.');
    else
        MdlFixed = fitcsvm(Xtrain, ytrain, ...
            'KernelFunction','linear', 'Standardize', true);

        % ---- Optional speed cache for ReadPerTrial (prefix sums over time) ----
        cache = [];
        if ReadPerTrial
            fprintf('=== COLOR SVM: caching per-trial prefix sums for fast time-resolved test ===\n');
            nTime = numel(tb);
            cache = cell(numel(RelArray),1);

            for a = 1:numel(RelArray)
                Arr = RelArray(a);
                fprintf('  reading Array %d\n', Arr); drawnow;

                ch = SVM_Chans_use(SVM_Arrs_use == Arr);
                ch = ch(:);
                if isempty(ch), continue; end

                fA = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colA);
                fB = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colB);
                idxA = find(fA);
                idxB = find(fB);
                if numel(idxA) < 5 || numel(idxB) < 5, continue; end

                ch = sort(unique(ch));
                nCh = numel(ch);
                cuts = [1; find(diff(ch)~=1)+1; numel(ch)+1];

                sumA = zeros(nCh, numel(idxA), nTime, 'single');
                cntA = zeros(nCh, numel(idxA), nTime, 'single');
                sumB = zeros(nCh, numel(idxB), nTime, 'single');
                cntB = zeros(nCh, numel(idxB), nTime, 'single');

                for j = 1:numel(idxA)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), idxA(j), :);     % [nRun×1×nTime]
                        X2 = reshape(single(X), numel(chRun), 1, nTime);
                        M  = ~isnan(X2);
                        sumA(rr,j,:) = cumsum(X2 .* single(M), 3);
                        cntA(rr,j,:) = cumsum(single(M), 3);
                    end
                end

                for j = 1:numel(idxB)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), idxB(j), :);
                        X2 = reshape(single(X), numel(chRun), 1, nTime);
                        M  = ~isnan(X2);
                        sumB(rr,j,:) = cumsum(X2 .* single(M), 3);
                        cntB(rr,j,:) = cumsum(single(M), 3);
                    end
                end

                cache{a} = struct('Arr',Arr,'ch',ch,'idxA',idxA,'idxB',idxB, ...
                                  'sumA',sumA,'cntA',cntA,'sumB',sumB,'cntB',cntB);
            end
        end

        % ---- Time-resolved evaluation of FIXED classifier ----
        accTime = nan(nT,1);
        accSE_t = nan(nT,1);

        fprintf('=== COLOR SVM: time-resolved test (win=%dms, step=%dms) ===\n', winWidthMs, stepMs);

        for ti = 1:nT
            tc = centers(ti);
            tidx_t = (tb >= tc - winWidthMs/2) & (tb <= tc + winWidthMs/2);
            tIdx = find(tidx_t);
            if numel(tIdx) < 3, continue; end
            t1 = tIdx(1);
            t2 = tIdx(end);

            accRep = nan(nRepeats,1);

            for r = 1:nRepeats
                XcolA_t = cell(numel(RelArray),1);
                XcolB_t = cell(numel(RelArray),1);

                for a = 1:numel(RelArray)
                    Arr = RelArray(a);

                    if ReadPerTrial
                        C = cache{a};
                        if isempty(C), continue; end

                        if t1 == 1
                            sA = C.sumA(:,:,t2);  nA = C.cntA(:,:,t2);
                            sB = C.sumB(:,:,t2);  nB = C.cntB(:,:,t2);
                        else
                            sA = C.sumA(:,:,t2) - C.sumA(:,:,t1-1);
                            nA = C.cntA(:,:,t2) - C.cntA(:,:,t1-1);
                            sB = C.sumB(:,:,t2) - C.sumB(:,:,t1-1);
                            nB = C.cntB(:,:,t2) - C.cntB(:,:,t1-1);
                        end

                        A_A_t = double(sA) ./ double(nA);
                        A_B_t = double(sB) ./ double(nB);
                        A_A_t(nA==0) = NaN;
                        A_B_t(nB==0) = NaN;

                    else
                        ch = SVM_Chans_use(SVM_Arrs_use == Arr);
                        ch = ch(:);
                        if isempty(ch), continue; end

                        fA = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colA);
                        fB = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,COLOR_COL)==colB);

                        A_A_t = squeeze(nanmean(normMUA(ch, fA, tidx_t), 3));
                        A_B_t = squeeze(nanmean(normMUA(ch, fB, tidx_t), 3));
                    end

                    XcolA_t{a} = A_A_t.'; % trials × channels
                    XcolB_t{a} = A_B_t.'; % trials × channels
                end
                
                [Xtest, ytest] = make_pseudotrials_concat_binary(XcolA_t, XcolB_t, nPseudoPerClass, nDrawPerPseudo);
                if isempty(Xtest)
                    accRep(r) = NaN;
                else
                    yhat = predict(MdlFixed, Xtest);
                    accRep(r) = mean(yhat == ytest);
                end
            end

            accTime(ti) = mean(accRep, 'omitnan');
            accSE_t(ti) = std(accRep, 'omitnan') / sqrt(sum(~isnan(accRep)));
        end

        
        if ~exist('MinSigma','var') || isempty(MinSigma), MinSigma = 5; end
        if ~exist('MaxSigma','var') || isempty(MaxSigma), MaxSigma = 400; end

 
% ---------- fit (chance-corrected + smoothing, like previous script) ----------
        t   = SVM_Color.centers(:);
        acc = SVM_Color.accTime(:);

        y = acc - 0.5;  % chance-correct for fitting

        ok   = ~isnan(t) & ~isnan(y);
        tfit = t(ok);
        yfit = y(ok);

        modfun = @(p,tt) exGauss_mod(p,tt);

        p0 = [200, 20, 0.01, 0.1, 0];
        lb = [min(tfit), MinSigma, 0, 0, 0];
        ub = [max(tfit), MaxSigma, 10, 1, 0.5];

        opts = optimoptions('lsqcurvefit','Display','off');
        pHat = lsqcurvefit(modfun, p0, tfit, yfit, lb, ub, opts);

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

        % ---- Store + plot ----
        SVM_Color = struct();
        SVM_Color.centers = centers(:);
        SVM_Color.accTime = accTime(:);
        SVM_Color.accSE_t = accSE_t(:);
        SVM_Color.MdlFixed = MdlFixed;
        SVM_Color.trainWindow = TrainWindow;
        SVM_Color.testWinWidthMs = winWidthMs;
        SVM_Color.testStepMs = stepMs;
        assignin('base','SVM_Color',SVM_Color);
        nChTotal = numel(SVM_Chans_use);

        figure('Color','w'); hold on
        plot(SVM_Color.centers, SVM_Color.accTime, 'LineWidth', 2);
        plot(tfit, accFit, 'LineWidth', 2);
        yline(0.5,'--k'); xline(0,'--k');
        grid on; box off
        xlabel('Time (ms) (window center)')
        ylabel('Color decoding accuracy')
        title(sprintf('Color SVM [%g %g] ms; %d ms window, Nchans= %d, 33%% point = %.1f ms', ...
            TrainWindow(1), TrainWindow(2), winWidthMs, nChTotal, t33));
    end
end


function [X, y] = make_pseudotrials_concat_binary(XA_cell, XB_cell, nPseudoPerClass, nDrawPerPseudo)

    if nargin < 4 || isempty(nDrawPerPseudo)
        nDrawPerPseudo = 25;
    end

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

        dA = min(nDrawPerPseudo, nA);
        dB = min(nDrawPerPseudo, nB);

        PA = nan(nPseudoPerClass, size(XA,2));
        PB = nan(nPseudoPerClass, size(XB,2));

        for k = 1:nPseudoPerClass
            ia = randi(nA, [dA,1]);
            ib = randi(nB, [dB,1]);
            PA(k,:) = mean(XA(ia,:), 1, 'omitnan');
            PB(k,:) = mean(XB(ib,:), 1, 'omitnan');
        end

        X_A_all = [X_A_all, PA]; %#ok<AGROW>
        X_B_all = [X_B_all, PB]; %#ok<AGROW>
    end

    X = [X_A_all; X_B_all];
    y = [ones(size(X_A_all,1),1); -ones(size(X_B_all,1),1)];
end
