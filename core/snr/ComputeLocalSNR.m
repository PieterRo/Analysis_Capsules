% look for SNR in the data if required for the analysis, note that this
% section has to run in case of switching to the other monkey
if SelectLocalSNR 
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

