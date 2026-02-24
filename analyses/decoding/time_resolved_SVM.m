%% ===== Time-resolved SVM: classification accuracy vs time (per distance) =====
% Requires: normMUA, tb, ALLMAT, RelArray, SVM_Chans, SVM_Arrs
% ALLMAT(:,2)=array, ALLMAT(:,3)=distance(1..3), ALLMAT(:,4)=attention (2 attend RF, 1 elsewhere)

ARR_COL  = 2;
DIST_COL = 3;
ATT_COL  = 4;
attVal   = 2;
unattVal = 1;

distLevels = 1:3;

% --- Sliding window params (edit to taste) ---
winWidthMs = 50;   % window width in ms
stepMs     = 10;   % step in ms

% Convert ms to indices using tb
% (assumes tb is in ms, roughly 1 ms resolution)
tStart = 0;        % optionally start at 0 to avoid prestim baseline triggering
tEnd   = 500;

centers = tStart:stepMs:tEnd;
nT = numel(centers);

% --- SVM / pseudotrial params ---
nRepeats        = 30;     % resample pseudotrials for stability
nFolds          = 10;
nPseudoPerClass = 200;    % pseudotrials per class
doPermutation   = false;  % time-resolved perms can be slow

% Output arrays: [nT x 3] for mean accuracy and SE
accMean = nan(nT, numel(distLevels));
accSE   = nan(nT, numel(distLevels));

fprintf('\n=== Time-resolved SVM (win=%d ms, step=%d ms) ===\n', winWidthMs, stepMs);

for di = 1:numel(distLevels)
    d = distLevels(di);
    fprintf('Distance %d\n', d);

    % --------- Precompute trial masks per array for this distance ----------
    % and store per-array channel lists
    arrCh = cell(numel(RelArray),1);
    attMask = cell(numel(RelArray),1);
    unMask  = cell(numel(RelArray),1);

    ok = true;
    for a = 1:numel(RelArray)
        Arr = RelArray(a);

        ch = SVM_Chans(SVM_Arrs == Arr);
        ch = ch(:);
        if isempty(ch)
            warning('No SVM channels for Arr=%d (dist=%d).', Arr, d);
            ok = false; break;
        end

        f_att = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d) & (ALLMAT(:,ATT_COL)==attVal);
        f_un  = (ALLMAT(:,ARR_COL)==Arr) & (ALLMAT(:,DIST_COL)==d) & (ALLMAT(:,ATT_COL)==unattVal);

        if nnz(f_att) < 5 || nnz(f_un) < 5
            warning('Too few trials Arr=%d dist=%d (nAtt=%d, nUn=%d).', Arr, d, nnz(f_att), nnz(f_un));
            ok = false; break;
        end

        arrCh{a}   = ch;
        attMask{a} = f_att;
        unMask{a}  = f_un;
    end

    if ~ok
        fprintf('  -> Skipping distance %d\n', d);
        continue;
    end

    % --------- Loop over time windows ----------
    for ti = 1:nT
        tc = centers(ti);
        w1 = tc - winWidthMs/2;
        w2 = tc + winWidthMs/2;

        tidx = (tb >= w1) & (tb <= w2);
        if nnz(tidx) < 3
            continue
        end

        % Build per-array trial×feature matrices for this window
        Xatt = cell(numel(RelArray),1);
        Xun  = cell(numel(RelArray),1);

        for a = 1:numel(RelArray)
            ch    = arrCh{a};
            f_att = attMask{a};
            f_un  = unMask{a};

            % mean over time within window -> [nCh x nTrials]
            A_att = squeeze(nanmean(normMUA(ch, f_att, tidx), 3));
            A_un  = squeeze(nanmean(normMUA(ch, f_un,  tidx), 3));

            % rows=trials, cols=channels
            Xatt{a} = A_att.';
            Xun{a}  = A_un.';
        end

        % Repeat pseudotrial sampling + CV
        accRep = nan(nRepeats,1);
        for r = 1:nRepeats
            [Xp, yp] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass);

            Mdl   = fitcsvm(Xp, yp, 'KernelFunction','linear', 'Standardize', true);
            CVMdl = crossval(Mdl, 'KFold', nFolds);
            yhat  = kfoldPredict(CVMdl);

            accRep(r) = mean(yhat == yp);
        end

        accMean(ti, di) = mean(accRep);
        accSE(ti, di)   = std(accRep)/sqrt(nRepeats);
    end
end

% --------- Plot accuracy vs time ----------
figure('Color','w'); hold on
for di = 1:numel(distLevels)
    plot(centers, accMean(:,di), 'LineWidth', 2);
end
yline(0.5,'--k'); % chance for balanced classes
xline(0,'--k');
grid on; box off
xlabel('Time (ms) (window center)')
ylabel('SVM accuracy')
legend({'Dist 1','Dist 2','Dist 3'}, 'Location','best')
title(sprintf('Time-resolved SVM (win=%d ms, step=%d ms)', winWidthMs, stepMs))

% Save results
SVM_TimeCourse.centers = centers;
SVM_TimeCourse.accMean = accMean;
SVM_TimeCourse.accSE   = accSE;
assignin('base','SVM_TimeCourse',SVM_TimeCourse);

%% --- helper: pseudotrials by concatenating across arrays ---
function [X, y] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass)
nArr = numel(Xatt);
Xa = cell(nArr,1);
Xu = cell(nArr,1);

for a = 1:nArr
    nA = size(Xatt{a},1);
    nU = size(Xun{a},1);

    ia = randi(nA, [nPseudoPerClass,1]);
    iu = randi(nU, [nPseudoPerClass,1]);

    Xa{a} = Xatt{a}(ia,:);
    Xu{a} = Xun{a}(iu,:);
end

X = [cat(2, Xa{:}); cat(2, Xu{:})];
y = [ones(nPseudoPerClass,1); -ones(nPseudoPerClass,1)];
end
