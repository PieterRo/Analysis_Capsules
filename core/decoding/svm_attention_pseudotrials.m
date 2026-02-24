function Results = svm_attention_pseudotrials(normMUA, tb, ALLMAT, RelArray, ChansPerArray, ...
                                              SigChansByTiming, timingNames)
% normMUA: [NChansGlob x NTrials x NTime]
% tb:      [1 x NTime] or [NTime x 1] in ms
% ALLMAT:  [NTrials x ...] metadata
% RelArray: e.g. [1 2 3] (array IDs)
% ChansPerArray: 64
% SigChansByTiming: cell array {nTiming x 1}, each is vector of GLOBAL channel indices to include
% timingNames: cell array of names, same length as SigChansByTiming

% ---- USER: set columns in ALLMAT ----------------------------------------
ARR_COL  = 2;          % array id
ATT_COL  = 4;          % attention: 2=attend, 1=unattend (matches your earlier code)
DIST_COL = 3;          % GC condition 
attVal   = 2;
unattVal = 1;

distLevels = 1:3;      % your three distance conditions

% ---- analysis time window -----------------------------------------------
tidx = tb >= 200 & tb <= 500;   % collapse 200–500 ms

% ---- pseudotrial / SVM params -------------------------------------------
nPseudoPerClass = 200;     % increase if you have enough trials
nFolds          = 10;      % cross-val folds
nRepeats        = 50;      % repeat pseudotrial sampling for stability
doPermutation   = true;
nPerm           = 200;     % permutations per (timing, distance)

Results = struct();
tb = tb(:)'; %#ok<NASGU>

for tIdx = 1:numel(SigChansByTiming)
    sigGlobal = SigChansByTiming{tIdx}(:);     % global chan indices
    timingLabel = timingNames{tIdx};

    fprintf('\n=== Timing condition: %s | #sig chans (global) = %d ===\n', ...
        timingLabel, numel(sigGlobal));

    for d = distLevels
        fprintf('Distance %d\n', d);

        % ------------------------------------------------------------
        % 1) Build per-array trial feature matrices (rows = trials)
        % ------------------------------------------------------------
        Xatt  = cell(numel(RelArray),1);
        Xun   = cell(numel(RelArray),1);

        for a = 1:numel(RelArray)
            Arr = RelArray(a);

            % trials for this array + distance + attention
            f_att = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d) & (ALLMAT(:,ATT_COL)==attVal);
            f_un  = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d) & (ALLMAT(:,ATT_COL)==unattVal);
            fprintf('Arr %d dist %d: nAtt=%d nUnatt=%d\n', Arr, d, nnz(f_att), nnz(f_un));

            % global channel range for this array
            first_chan = (Arr-1)*ChansPerArray + 1;
            last_chan  = Arr*ChansPerArray;

            % significant channels restricted to this array
            sigThisArray = sigGlobal(sigGlobal >= first_chan & sigGlobal <= last_chan);

            if isempty(sigThisArray)
                warning('No significant channels for Arr=%d (timing=%s, dist=%d). Skipping.', Arr, timingLabel, d);
                Xatt{a} = [];
                Xun{a}  = [];
                continue
            end

            % --- Extract features: collapse time 200-500, keep trial-by-trial ---
            % raw: [nCh x nTrials x nTime] -> mean over time -> [nCh x nTrials]
            A_att = squeeze(nanmean(normMUA(sigThisArray, f_att, tidx), 3));  % [nCh x nAttTrials]
            A_un  = squeeze(nanmean(normMUA(sigThisArray, f_un,  tidx), 3));  % [nCh x nUnTrials]

            % transpose so rows are trials: [nTrials x nCh]
            Xatt{a} = A_att.';   % [nAttTrials x nCh]
            Xun{a}  = A_un.';    % [nUnTrials  x nCh]
        end

        % require all arrays to contribute
        if any(cellfun(@isempty, Xatt)) || any(cellfun(@isempty, Xun))
            fprintf('  -> Skipping dist %d (missing sig chans/trials in at least one array)\n', d);
            continue
        end

        % ------------------------------------------------------------
        % 2) Build pseudotrials by pairing trials across arrays
        %    (concatenate features from each array)
        % ------------------------------------------------------------
        accRep = nan(nRepeats,1);

        for r = 1:nRepeats
            [Xp, yp] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass);

            % --------------------------------------------------------
            % 3) Linear SVM with cross-validation
            % --------------------------------------------------------
            Mdl = fitcsvm(Xp, yp, 'KernelFunction','linear', 'Standardize', true);

            CVMdl = crossval(Mdl, 'KFold', nFolds);
            yhat = kfoldPredict(CVMdl);

            accRep(r) = mean(yhat == yp);
        end

        accMean = mean(accRep);
        accSE   = std(accRep)/sqrt(nRepeats);

        % ------------------------------------------------------------
        % 4) Optional permutation test (shuffle labels)
        % ------------------------------------------------------------
        p_perm = NaN;
        if doPermutation
            permAcc = nan(nPerm,1);

            % build one pseudotrial dataset to permute labels on
            [Xp0, yp0] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass);

            for k = 1:nPerm
                yperm = yp0(randperm(numel(yp0)));

                MdlP = fitcsvm(Xp0, yperm, 'KernelFunction','linear', 'Standardize', true);
                CVP  = crossval(MdlP, 'KFold', nFolds);
                yhatP = kfoldPredict(CVP);
                permAcc(k) = mean(yhatP == yperm);
            end

            % p-value: fraction perm >= observed mean accuracy
            p_perm = mean(permAcc >= accMean);
        end

        % store results
        Results.(timingLabel).dist(d).acc_mean = accMean;
        Results.(timingLabel).dist(d).acc_se   = accSE;
        Results.(timingLabel).dist(d).p_perm   = p_perm;
        Results.(timingLabel).dist(d).acc_repeats = accRep;

        fprintf('  Acc=%.3f ± %.3f (SE), perm p=%.4f\n', accMean, accSE, p_perm);
    end
end
end


% ======================================================================
% Helper: build pseudotrials by concatenating arrays
% ======================================================================
function [X, y] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass)
% Xatt/Xun: cell array over arrays, each [nTrials x nFeatArray]
% returns:
%   X: [2*nPseudoPerClass x sum(nFeatArray)]  (concat across arrays)
%   y: labels (+1 attended, -1 unattended)

nArr = numel(Xatt);

% sample with replacement independently per array
Xa = cell(nArr,1);
Xu = cell(nArr,1);

for a = 1:nArr
    nA = size(Xatt{a},1);
    nU = size(Xun{a},1);

    ia = randi(nA, [nPseudoPerClass,1]);
    iu = randi(nU, [nPseudoPerClass,1]);

    Xa{a} = Xatt{a}(ia, :);
    Xu{a} = Xun{a}(iu, :);
end

X_att_pseudo = cat(2, Xa{:});   % concat features across arrays
X_un_pseudo  = cat(2, Xu{:});

X = [X_att_pseudo; X_un_pseudo];
y = [ones(nPseudoPerClass,1); -ones(nPseudoPerClass,1)];
end
