function [X, y] = make_pseudotrials_concat(Xatt, Xun, nPseudoPerClass)
% Build pseudotrials by sampling trials independently per array and concatenating features.
% Xatt/Xun are cell arrays {nArr x 1}, each [nTrials x nFeatArr].
    nArr = numel(Xatt);

    Xa = cell(nArr,1);
    Xu = cell(nArr,1);

    for a = 1:nArr
        nA = size(Xatt{a},1);
        nU = size(Xun{a},1);

        % Skip this array entirely if one class is missing
        if nA == 0 || nU == 0
            Xa{a} = [];
            Xu{a} = [];
        else 
            ia = randi(nA, [nPseudoPerClass,1]);  % sample w/ replacement
            iu = randi(nU, [nPseudoPerClass,1]);
            Xa{a} = Xatt{a}(ia, :);
            Xu{a} = Xun{a}(iu, :);
        end
    end

    keep = ~cellfun(@isempty, Xa);   % only use arrays that have channels
    X_att = cat(2, Xa{keep});
    X_un = cat(2, Xu{keep});

    X = [X_att; X_un];
    y = [ones(nPseudoPerClass,1); -ones(nPseudoPerClass,1)]; % +1 attend, -1 unattend
end