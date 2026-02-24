
% Note that there also is a folder with RFs of arrays 6 and 7
% '/Users/pieter/Library/CloudStorage/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Extras/monkeyN/RFs/array9_good_rfs.mat'

% See also
% '/Users/pieter/Library/CloudStorage/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Extras/monkeyN/exp/comp1/array_2_class1_complexity1_cones2_stim_1.png'

Monkey = 1; % 1 for Nilson, 2 for Figaro
SNR_th = 0.6;                 % activity threshold for good enough SNR 
% 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.075;
% 
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
cfg = config();
Dir = cfg.rootDir;
StimDir = cfg.stimDir;
cd(Dir);

% if Monkey == 1
%     load('/Users/pieter/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_normMUA.mat');
%     load('/Users/pieter/Dropbox/Pieter/data/Mr Nilson/ObjAtt_lines_MUA_trials.mat');
% else
%     load('/Users/pieter/Dropbox/Pieter/data/Figaro/ObjAtt_lines_normMUA.mat');
%     load('/Users/pieter/Dropbox/Pieter/data/Figaro/ObjAtt_lines_MUA_trials.mat');
% end

if Monkey == 1
    m1 = matfile(fullfile(cfg.dataDir, 'Mr Nilson', 'ObjAtt_lines_normMUA.mat'));
    m2 = matfile(fullfile(cfg.dataDir, 'Mr Nilson', 'ObjAtt_lines_MUA_trials.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
else
    m1 = matfile(fullfile(cfg.dataDir, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
    m2 = matfile(fullfile(cfg.dataDir, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
    SNR = m1.SNR;
    [NChansGlob,NTrialsGlob,NTimesGlob] = size(m1, 'normMUA');
    ALLMAT = m2.ALLMAT;
    tb=m2.tb;
    ReadPerTrial = 1; % necessary because I can't load normMUA at once. 
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

% look for SNR of the visual response in the data if required for the analysis
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

% compute significances of attention effect in the three cone conditions and store in p_vals_global
if ~exist('p_vals_global','var')                % if true, looks at the p-value in the actual distance condition
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

pmin = min(p_vals_global, [], 2, 'omitnan');


% select channels from the first array with sufficient SNR and significant
% attentional modulation

Arr = RelArray(3);
RelCone = RelCones(1);
f_att = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 2; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
u_att = ALLMAT(:,3) == RelCone & ALLMAT(:,2) == Arr & ALLMAT(:,4) == 1; % attended trials: correct RelCones + RelArray and attention==2 ("on RF")
the_trialsA = find(f_att);
the_trialsU = find(u_att);

stimsA = unique(ALLMAT(the_trialsA,1));
stimsU = unique(ALLMAT(the_trialsU,1));

%% channels of a specific array
% first_chan = (Arr-1)*ChansPerArray + 1;  % first channel index in that array
% last_chan  = (Arr)*ChansPerArray;        % last channel index in that array
% 
% RelChans = ones(1, 1024);        % start with all
% RelChans(1:first_chan-1) = 0;     % exclude channels in arrays before RelArray
% RelChans(last_chan+1:end) = 0;    % exclude channels in arrays after RelArray
% 
% RelChans(SNR_local<SNR_th) = 0;   % exclude channels with weak visual response
% RelChans(pmin > p_th) = 0;   % exclude channels with weak visual response
% 
% goodGlobal = find(RelChans);      % global channel indices selected this iteration
% nGoodChans = numel(goodGlobal);

%% all V1 channels
first_chan = 1;  
last_chan  = 512;

RelChans = ones(1, 1024);        % start with all
RelChans(last_chan+1:end) = 0;     % exclude channels in arrays before RelArray
goodGlobal = find(RelChans);


RFs % run this first 

stimIdx = 1;                    % choose 1..328
cx = 512; cy = 384;             % bitmap pixel corresponding to RF (0,0)

% If RF y is positive upward (visual coordinates), set this true:
flipY = true;   % try true first; if it's mirrored vertically, set false

%% --- Load bitmap ---
bmpFile = fullfile(StimDir, sprintf('%03d.bmp', stimIdx));
I = imread(bmpFile);
[H,W,~] = size(I);

%% --- RF vectors (assumes you already have x,y for the current monkey) ---
x_rf = x(goodGlobal);   % RF centers in "pixel offsets" relative to (0,0)
y_rf = y(goodGlobal);

% Convert to bitmap pixels
x_pix = cx + x_rf;
if flipY
    y_pix = cy - y_rf;
else
    y_pix = cy + y_rf;
end

% Valid points: finite + inside image bounds
valid = isfinite(x_pix) & isfinite(y_pix) & x_pix>=1 & x_pix<=W & y_pix>=1 & y_pix<=H;

% ---- Create figure + axes ----
figure('Color','w');
ax = axes();                         % explicit axes handle
imshow(I, 'Parent', ax);
hold(ax, 'on');                      % IMPORTANT: after imshow, and on that axes
axis(ax, 'image');
title(ax, sprintf('Stim %03d with RF centers', stimIdx));


% Plot centers by area (use your existing colors cV1/cV4/cIT)
scatter(x_pix, y_pix, 18, cV1, 'filled', ...
    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha', 0.25);

% Mark the RF-origin / image center
plot(cx, cy, 'wx', 'MarkerSize', 12, 'LineWidth', 5); % white X
plot(cx, cy, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % black X on top

legend({'V1','(0,0) at (512,384)'}, 'Location','best');
