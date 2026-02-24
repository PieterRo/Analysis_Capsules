% runs this after running Latency_analysis_SVM
%
% FAST version: avoids repeated disk reads of normMUA when ReadPerTrial==true
% by caching per-trial prefix sums per distance condition.

%% ===================== SVM pseudotrials across arrays =====================
% Requires variables already in workspace:
%   tb [1 x NTime], ALLMAT [NTrials x ...], RelArray, RelCones, AttentionWindow
%   SVM_Chans (global channel indices), SVM_Arrs (array id per channel)
%   If p_val_local==true: rel_chans_per_cone{d}, rel_arr_per_cone{d}
%   If ReadPerTrial==true: m1 = matfile(...); with variable m1.normMUA on disk
%
% Assumed ALLMAT columns:
%   ALLMAT(:,2) = array id
%   ALLMAT(:,3) = distance condition (1..3)
%   ALLMAT(:,4) = attention (2 = attend RF, 1 = elsewhere)

% --- sanity for inputs ---
if ~exist('p_val_local','var'), p_val_local = false; end
if ~exist('ReadPerTrial','var'), ReadPerTrial = false; end

if (~p_val_local && (~exist('SVM_Chans','var') || isempty(SVM_Chans))) || ...
   ( p_val_local && (~exist('rel_chans_per_cone','var') || isempty(rel_chans_per_cone)) )
    warning('SVM channels not found / empty. Skipping SVM section.');
else
    ARR_COL  = 2;
    DIST_COL = 3;
    ATT_COL  = 4;
    attVal   = 2;
    unattVal = 1;

    distLevels = 1:3;

    % Collapse time from within the attention window (training)
    tidx_train = (tb >= AttentionWindow(1)) & (tb <= AttentionWindow(2));

    % SVM/pseudotrial settings
    nRepeats        = 50;   % resample pseudotrials this many times
    nPseudoPerClass = 200;  % pseudotrials per class per repeat (sample w/ replacement)

    % Results
    SVM_Results = struct();

    % Ensure column vectors (global lists)
    if exist('SVM_Chans','var'), SVM_Chans = SVM_Chans(:); end
    if exist('SVM_Arrs','var'),  SVM_Arrs  = SVM_Arrs(:);  end

    % matfile requirement for large data
    if ReadPerTrial
        if ~exist('m1','var')
            error('ReadPerTrial==true but matfile handle m1 does not exist. Create it with: m1 = matfile(''yourfile.mat'');');
        end
        if ~ismember('normMUA', who(m1))
            error('matfile m1 does not contain variable ''normMUA''.');
        end
        nTime = numel(tb);
    else
        if ~exist('normMUA','var')
            error('ReadPerTrial==false but variable normMUA is not in workspace.');
        end
    end

    fprintf('\n=== Running pseudotrial SVM (trained on AttentionWindow) ===\n');

    for d = distLevels
        fprintf('Distance condition %d\n', d);

        d_val = RelCones(d);

        % --- Select channel set for this cone/distance condition ---
        if p_val_local
            SVM_Chans_use = rel_chans_per_cone{d}(:);
            SVM_Arrs_use  = rel_arr_per_cone{d}(:);
        else
            SVM_Chans_use = SVM_Chans(:);
            SVM_Arrs_use  = SVM_Arrs(:);
        end

        % ---------- Build per-array trial feature matrices (TRAIN) ----------
        Xatt = cell(numel(RelArray),1);
        Xun  = cell(numel(RelArray),1);

        for a = 1:numel(RelArray)
            Arr = RelArray(a);

            ch = SVM_Chans_use(SVM_Arrs_use == Arr);
            ch = ch(:);

            if isempty(ch)
                warning('No SVM channels found for Arr=%d (dist=%d). Skipping this array.', Arr, d);
                continue;
            end

            f_att = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==attVal);
            f_un  = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==unattVal);

            attIdx = find(f_att);
            unIdx  = find(f_un);

            if numel(attIdx) < 5 || numel(unIdx) < 5
                warning('Too few trials for Arr=%d dist=%d (nAtt=%d, nUn=%d). Skipping this array.', Arr, d, numel(attIdx), numel(unIdx));
                continue;
            end

            if ReadPerTrial
                % Read each trial ONCE for the training window only
                % (channels may be non-contiguous -> split into contiguous runs)
                ch = sort(unique(ch));
                nCh = numel(ch);
                cuts = [1; find(diff(ch)~=1)+1; numel(ch)+1];

                A_att = NaN(nCh, numel(attIdx));
                A_un  = NaN(nCh, numel(unIdx));

                tidx2 = find(tidx_train);
                tRange = tidx2(1):tidx2(end); % contiguous

                for j = 1:numel(attIdx)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), attIdx(j), tRange);  % [nRun×1×nT]
                        X2 = reshape(X, numel(chRun), []);                        % [nRun×nT]
                        A_att(rr,j) = mean(X2, 2, 'omitnan');
                    end
                end

                for j = 1:numel(unIdx)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), unIdx(j), tRange);
                        X2 = reshape(X, numel(chRun), []);
                        A_un(rr,j) = mean(X2, 2, 'omitnan');
                    end
                end
            else
                A_att = squeeze(nanmean(normMUA(ch, f_att, tidx_train), 3));
                A_un  = squeeze(nanmean(normMUA(ch, f_un,  tidx_train), 3));
            end

            Xatt{a} = A_att.';  % trials × channels
            Xun{a}  = A_un.';   % trials × channels
        end

        % ---------- Train ONE fixed classifier on the long window ----------
        [Xtrain, ytrain] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass);

        if isempty(Xtrain)
            warning('No training data for distance %d. Skipping.', d);
            continue
        end

        MdlFixed = fitcsvm(Xtrain, ytrain, 'KernelFunction','linear', 'Standardize', true);

        % ---------- Cache per-trial prefix sums ONCE per distance (ReadPerTrial) ----------
        % This removes repeated disk reads inside the time-window loop.
        cache = [];
        if ReadPerTrial
            cache = cell(numel(RelArray),1);

            for a = 1:numel(RelArray)
                Arr = RelArray(a);
                fprintf(' reading data Array %d\n', Arr);
                drawnow;
                
                ch = SVM_Chans_use(SVM_Arrs_use == Arr);
                ch = ch(:);
                if isempty(ch), continue; end

                f_att = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==attVal);
                f_un  = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==unattVal);
                attIdx = find(f_att);
                unIdx  = find(f_un);
                if numel(attIdx) < 5 || numel(unIdx) < 5, continue; end

                ch = sort(unique(ch));
                nCh = numel(ch);
                cuts = [1; find(diff(ch)~=1)+1; numel(ch)+1];

                % prefix sums/counts over time for each trial
                sumAtt = zeros(nCh, numel(attIdx), nTime, 'single');
                cntAtt = zeros(nCh, numel(attIdx), nTime, 'single');
                sumUn  = zeros(nCh, numel(unIdx),  nTime, 'single');
                cntUn  = zeros(nCh, numel(unIdx),  nTime, 'single');

                for j = 1:numel(attIdx)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), attIdx(j), :);       % [nRun×1×nTime]
                        X2 = reshape(single(X), numel(chRun), 1, nTime);          % [nRun×1×nTime]
                        M  = ~isnan(X2);

                        sumAtt(rr,j,:) = cumsum(X2 .* single(M), 3);
                        cntAtt(rr,j,:) = cumsum(single(M), 3);
                    end
                end

                for j = 1:numel(unIdx)
                    row0 = 0;
                    for rruns = 1:numel(cuts)-1
                        chRun = ch(cuts(rruns):cuts(rruns+1)-1);
                        rr    = (row0+1):(row0+numel(chRun));
                        row0  = row0 + numel(chRun);

                        X  = m1.normMUA(chRun(1):chRun(end), unIdx(j), :);
                        X2 = reshape(single(X), numel(chRun), 1, nTime);
                        M  = ~isnan(X2);

                        sumUn(rr,j,:)  = cumsum(X2 .* single(M), 3);
                        cntUn(rr,j,:)  = cumsum(single(M), 3);
                    end
                end

                cache{a} = struct('Arr',Arr,'ch',ch,'attIdx',attIdx,'unIdx',unIdx, ...
                                  'sumAtt',sumAtt,'cntAtt',cntAtt,'sumUn',sumUn,'cntUn',cntUn);
            end
        end

        % ---------- Evaluate the FIXED classifier as a function of time ----------
        winWidthMs = 50;
        stepMs     = 10;
        centers    = tb(1):stepMs:tb(end);
        nT = numel(centers);

        accTime = nan(nT,1);
        accSE_t = nan(nT,1);

        for ti = 1:nT
            tc = centers(ti);
            tidx_t = (tb >= tc - winWidthMs/2) & (tb <= tc + winWidthMs/2);
            tIdx = find(tidx_t);
            if numel(tIdx) < 3
                continue
            end
            t1 = tIdx(1);
            t2 = tIdx(end);

            accRep_t = nan(nRepeats,1);

            for r = 1:nRepeats
                Xatt_t = cell(numel(RelArray),1);
                Xun_t  = cell(numel(RelArray),1);

                for a = 1:numel(RelArray)
                    Arr = RelArray(a);

                    if ReadPerTrial
                        C = cache{a};
                        if isempty(C), continue; end

                        if t1 == 1
                            sA = C.sumAtt(:,:,t2);  nA = C.cntAtt(:,:,t2);
                            sU = C.sumUn(:,:,t2);   nU = C.cntUn(:,:,t2);
                        else
                            sA = C.sumAtt(:,:,t2) - C.sumAtt(:,:,t1-1);
                            nA = C.cntAtt(:,:,t2) - C.cntAtt(:,:,t1-1);
                            sU = C.sumUn(:,:,t2)  - C.sumUn(:,:,t1-1);
                            nU = C.cntUn(:,:,t2)  - C.cntUn(:,:,t1-1);
                        end

                        A_att_t = double(sA) ./ double(nA);
                        A_un_t  = double(sU) ./ double(nU);
                        A_att_t(nA==0) = NaN;
                        A_un_t(nU==0)  = NaN;

                    else
                        ch = SVM_Chans_use(SVM_Arrs_use == Arr);
                        ch = ch(:);

                        f_att = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==attVal);
                        f_un  = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d_val) & (ALLMAT(:,ATT_COL)==unattVal);

                        A_att_t = squeeze(nanmean(normMUA(ch, f_att, tidx_t), 3));
                        A_un_t  = squeeze(nanmean(normMUA(ch, f_un,  tidx_t), 3));
                    end

                    Xatt_t{a} = A_att_t.';  % trials × channels
                    Xun_t{a}  = A_un_t.';   % trials × channels
                end

                [Xtest, ytest] = make_pseudotrials_concat(Xatt_t, Xun_t, nPseudoPerClass);
                if isempty(Xtest)
                    accRep_t(r) = NaN;
                    continue
                end

                yhat = predict(MdlFixed, Xtest);
                accRep_t(r) = mean(yhat == ytest);
            end

            accTime(ti) = mean(accRep_t, 'omitnan');
            accSE_t(ti) = std(accRep_t, 'omitnan') / sqrt(sum(~isnan(accRep_t)));
        end

        % Store timecourse per distance
        SVM_Results.dist(d).centers  = centers;
        SVM_Results.dist(d).accTime  = accTime;
        SVM_Results.dist(d).accSE_t  = accSE_t;
        SVM_Results.dist(d).MdlFixed = MdlFixed;
        SVM_Results.dist(d).trainWindow = AttentionWindow;
        SVM_Results.dist(d).testWinWidthMs = winWidthMs;
        SVM_Results.dist(d).testStepMs = stepMs;
    end

    assignin('base','SVM_Results',SVM_Results);
    fprintf('=== Done. Results in SVM_Results ===\n');
end

figure('Color','w'); hold on

hData = gobjects(numel(distLevels),1);
hFit  = gobjects(numel(distLevels),1);

t33_all = nan(numel(distLevels),1);
col_all = nan(numel(distLevels),3);

for ii = 1:numel(distLevels)
    d = distLevels(ii);
    if ~isfield(SVM_Results.dist(d),'centers'), continue; end

    % --- data ---
    t   = SVM_Results.dist(d).centers(:);
    acc = SVM_Results.dist(d).accTime(:);

    % chance-correct for fitting
    y = acc - 0.5;

    % ---------- plot SVM result (SOLID) ----------
    hData(ii) = plot(t, acc, '-', 'LineWidth', 2);
    col = hData(ii).Color;
    col_all(ii,:) = col;

    % ---------- fit ----------
    ok = ~isnan(t) & ~isnan(y);
    tfit = t(ok);
    yfit = y(ok);

    y_smooth = smoothdata(yfit, 'movmean', SmoothW);

    modfun = @(p,t) exGauss_mod(p,t);

    p0 = [200, 20, 0.01, 0.1, 0];
    lb = [min(tfit), MinSigma, 0, 0, 0];
    ub = [max(tfit), MaxSigma, 10, 1, 0.5];

    opts = optimoptions('lsqcurvefit','Display','off');
    pHat = lsqcurvefit(modfun, p0, tfit, y_smooth, lb, ub, opts);

    yFit_corr = modfun(pHat, tfit);
    accFit    = yFit_corr + 0.5;

    % ---------- compute 33% of max (on corrected fit) ----------
    mx = max(yFit_corr, [], 'omitnan');
    thr = 0.33 * mx;

idx33 = find(yFit_corr >= thr, 1, 'first');

    if isempty(idx33)
        t33 = NaN;
    elseif idx33 == 1
        t33 = tfit(1);
    else
        t33 = interp1( ...
            yFit_corr(idx33-1:idx33), ...
            tfit(idx33-1:idx33), ...
            thr, 'linear');
    end
    t33_all(ii) = t33;

    % ---------- plot FIT (DASHED, SAME COLOR) ----------
    hFit(ii) = plot(tfit, accFit, '--', 'LineWidth', 3, 'Color', col);

    % store (optional)
    SVM_Results.dist(d).fit.pHat = pHat;
    SVM_Results.dist(d).fit.t33  = t33;
end

% Cosmetics (draw AFTER lines; keep out of legend)
hChance = yline(0.5,'--k'); hChance.HandleVisibility = 'off';
hZero   = xline(0,'--k');   hZero.HandleVisibility   = 'off';

grid on; box off
xlabel('Time (ms) (window center)')
ylabel('Accuracy (fixed classifier)')

% Build legend from handles (this is the key fix)
legH = [hData; hFit];
legH = legH(isgraphics(legH));  % remove empties if any d skipped

labels = [ ...
    arrayfun(@(d)sprintf('Dist %d SVM', d), distLevels(:), 'UniformOutput', false); ...
    arrayfun(@(d)sprintf('Dist %d fit', d), distLevels(:), 'UniformOutput', false) ];
labels = labels(1:numel(legH)); % match if some skipped

legend(legH, labels, 'Location','southwest');

title(sprintf('Fixed SVM trained on [%g %g] ms; tested with %d ms window', ...
    AttentionWindow(1), AttentionWindow(2), SVM_Results.dist(distLevels(1)).testWinWidthMs))

% ---------- Add small arrows at the 33% time points ----------
ax = gca;
xlim_ = ax.XLim;
ylim_ = ax.YLim;

yArrow = 0.48;      % where arrow points to (in accuracy units)
arrowLen = 0.06;    % arrow length in normalized units (vertical)

for ii = 1:numel(distLevels)
    t33 = t33_all(ii);
    if isnan(t33), continue; end

    col = col_all(ii,:);

    % data -> normalized coordinates
    xNorm = (t33 - xlim_(1)) / diff(xlim_);
    yNorm = (0.5 - ylim_(1)) / diff(ylim_);

    % draw a short downward arrow
    annotation('arrow', ...
        [xNorm xNorm], ...
        [yNorm-0.05 yNorm], ...
        'Color', col, ...
        'LineWidth', 2);
end

% ---------- Add textbox with t33 values ----------
txtLines = cell(numel(distLevels),1);
for ii = 1:numel(distLevels)
    if isnan(t33_all(ii))
        txtLines{ii} = sprintf('Dist %d: n/a', distLevels(ii));
    else
        txtLines{ii} = sprintf('Dist %d: t_{33} = %.0f ms', ...
            distLevels(ii), t33_all(ii));
    end
end

txtStr = strjoin(txtLines, '\n');

annotation('textbox', ...
    [0.75 0.15 0.30 0.18], ...   % [x y w h] in normalized figure coords
    'String', txtStr, ...
    'FitBoxToText','on', ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.3 0.3 0.3], ...
    'FontSize', 10);



