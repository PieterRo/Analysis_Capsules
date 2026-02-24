%% ================= Plot color-tuned channels split by cone condition (2/4/6) =================
% One row per significant channel; 3 columns per row (cone conditions 2,4,6).
% Each subplot: mean time course for the two colors (ALLMAT(:,8)==2 vs ==1).
% New figure every 5 channels.
% Y-limits are shared across the 3 cone panels within each row.

if ~exist('col_tune_chans','var') || isempty(col_tune_chans)
    warning('col_tune_chans is empty or missing. Nothing to plot.');
else
    % --- settings ---
    coneVals    = [2 4 6];   % ALLMAT(:,3)
    colA        = 2;         % ALLMAT(:,8) value for color A
    colB        = 1;         % ALLMAT(:,8) value for color B
    PlotWindow  = [-200 500];% ms relative to stim onset
    chansPerFig = 5;         % new figure every 5 channels
    yPadFrac    = 0.05;      % add 5% padding to row-wise y-lims

    % --- time indices (use contiguous range for matfile safety) ---
    tIdx = find(tb >= PlotWindow(1) & tb <= PlotWindow(2));
    if isempty(tIdx)
        error('No tb samples within PlotWindow [%g %g] ms.', PlotWindow(1), PlotWindow(2));
    end
    tRange = tIdx(1):tIdx(end);   % contiguous indices (matfile-safe)
    tPlot  = tb(tRange);

    % --- total channels/pages ---
    chList  = col_tune_chans(:);
    nCh     = numel(chList);
    nFigs   = ceil(nCh / chansPerFig);

    for f = 1:nFigs
        i1 = (f-1)*chansPerFig + 1;
        i2 = min(f*chansPerFig, nCh);
        chThis = chList(i1:i2);
        nRows  = numel(chThis);

        figure('Color','w');
        tl = tiledlayout(nRows, 3, 'Padding','compact', 'TileSpacing','compact');
        title(tl, sprintf('Color-tuned channels by cone condition (page %d/%d)', f, nFigs), ...
            'FontWeight','bold');

        % keep axes handles per row so we can apply shared ylim
        axRow = gobjects(nRows,3);

        for rr = 1:nRows
            chIdx = chThis(rr);
            Arr   = ceil(chIdx / ChansPerArray);

            % Precompute all 3 cone-condition traces for this row
            muA_all = cell(1,3);
            muB_all = cell(1,3);
            nA_all  = zeros(1,3);
            nB_all  = zeros(1,3);

            for cc = 1:3
                cone = coneVals(cc);

                fA = (ALLMAT(:,2)==Arr) & (ALLMAT(:,3)==cone) & (ALLMAT(:,8)==colA);
                fB = (ALLMAT(:,2)==Arr) & (ALLMAT(:,3)==cone) & (ALLMAT(:,8)==colB);

                idxA = find(fA);
                idxB = find(fB);
                nA_all(cc) = numel(idxA);
                nB_all(cc) = numel(idxB);

                if ReadPerTrial
                    muA = nan(numel(tRange),1);
                    muB = nan(numel(tRange),1);

                    if ~isempty(idxA)
                        s = zeros(numel(tRange),1); n = zeros(numel(tRange),1);
                        for j = 1:numel(idxA)
                            x = squeeze(m1.normMUA(chIdx, idxA(j), tRange)); x = x(:);
                            ok = ~isnan(x);
                            s(ok) = s(ok) + x(ok);
                            n(ok) = n(ok) + 1;
                        end
                        muA = s ./ max(n,1); muA(n==0) = NaN;
                    end

                    if ~isempty(idxB)
                        s = zeros(numel(tRange),1); n = zeros(numel(tRange),1);
                        for j = 1:numel(idxB)
                            x = squeeze(m1.normMUA(chIdx, idxB(j), tRange)); x = x(:);
                            ok = ~isnan(x);
                            s(ok) = s(ok) + x(ok);
                            n(ok) = n(ok) + 1;
                        end
                        muB = s ./ max(n,1); muB(n==0) = NaN;
                    end
                else
                    if any(fA)
                        muA = squeeze(nanmean(normMUA(chIdx, fA, tRange), 2));
                    else
                        muA = nan(numel(tRange),1);
                    end
                    if any(fB)
                        muB = squeeze(nanmean(normMUA(chIdx, fB, tRange), 2));
                    else
                        muB = nan(numel(tRange),1);
                    end
                end

                muA_all{cc} = muA;
                muB_all{cc} = muB;
            end

            % Row-wise shared y-limits across all 3 cones and both colors
            allY = [];
            for cc = 1:3
                allY = [allY; muA_all{cc}(:); muB_all{cc}(:)]; %#ok<AGROW>
            end
            yMin = min(allY, [], 'omitnan');
            yMax = max(allY, [], 'omitnan');

            if ~isfinite(yMin) || ~isfinite(yMax) || yMin == yMax
                % fallback
                yMin = -0.1; yMax = 0.1;
            else
                pad = yPadFrac * (yMax - yMin);
                yMin = yMin - pad;
                yMax = yMax + pad;
            end

            % Now plot the 3 panels
            for cc = 1:3
                cone = coneVals(cc);

                axRow(rr,cc) = nexttile;
                hold on;
                plot(tPlot, muA_all{cc}, 'LineWidth', 1.2);
                plot(tPlot, muB_all{cc}, 'LineWidth', 1.2);
                xline(0,'--');
                yline(0,':');
                xlim([tPlot(1) tPlot(end)]);
                ylim([yMin yMax]);
                box off;

                if rr == 1
                    title(sprintf('Cone %d', cone), 'FontWeight','bold');
                end
                if cc == 1
                    ylabel(sprintf('Ch %d (Arr %d)', chIdx, Arr));
                end
                if rr == nRows
                    xlabel('Time (ms)');
                end

                if rr == 1 && cc == 1
                    legend({sprintf('Color %d', colA), sprintf('Color %d', colB)}, ...
                        'Location','best', 'FontSize', 8);
                end

                txt = sprintf('nA=%d  nB=%d', nA_all(cc), nB_all(cc));
                text(0.02, 0.95, txt, 'Units','normalized', ...
                    'VerticalAlignment','top', 'FontSize', 8);
            end
        end
    end
end