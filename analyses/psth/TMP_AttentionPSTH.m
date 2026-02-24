
for ConeCond = 1:1           % 1, 2 or 3  % run only for one cone condition
    RelCone = RelCones(ConeCond); 

    % Loop through the arrays
    selArray   = [];   % array number (Arr)
    selChanG   = [];   % global channel index (1..NChansTotal)
    p_vals    = [];    % array with p-values
    t_vals   = [];     % array with t-values
    allRelChn = [];    % global indices of the included channels with enough activity

    for ArLoop=1 : NArrays
        Arr = RelArray(ArLoop);
% exclude day 3 here
        f_att = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 2 & ALLMAT(:,11) ~= 3; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
        f_unatt = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 1 & ALLMAT(:,11) ~= 3; % unattended trials: same, but attention==1 ("not on RF")


        % --- Select channels belonging to the chosen array -----------------------

        first_chan = 100;  % first channel index in that array
        last_chan  = 100;        % last channel index in that array
        
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
                    Xatt(j) = mean(x, 3, 'omitnan');  % array with all activities in the 55 selected trials per cone
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
