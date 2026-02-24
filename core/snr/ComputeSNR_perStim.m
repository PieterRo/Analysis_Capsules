% variance of the spontaneous activity
w0 = 1;  % spontaneous window
n = double(R.nTrials(:));                 % [Nstim x 1]
validStim = (n > 0);
n = n(validStim);

% Extract spontaneous moments for ALL channels
% R.meanAct: [Nchan x Nstim x Nwin]
mu_sp_byStim = squeeze(R.meanAct(:, validStim, w0));     % [Nchan x Nstim]
m2_sp_byStim = squeeze(R.meanSqAct(:, validStim, w0));   % [Nchan x Nstim]

% Weighted sums across stimuli (i.e., across all trials pooled)
N = sum(n);                                               % scalar
sumX  = mu_sp_byStim  * n;                                % [Nchan x 1]  (since n is [Nstim x 1])
sumX2 = m2_sp_byStim  * n;                                % [Nchan x 1]

% Pooled spontaneous mean per channel
mu_spont = sumX / N;                                      % [Nchan x 1]

% Pooled unbiased variance per channel:
% Var = (sumX2 - sumX.^2 / N) / (N - 1)
var_spont = (sumX2 - (sumX.^2) / N) / (N - 1);            % [Nchan x 1]
sd_spont  = sqrt(var_spont);                              % [Nchan x 1]

% Z-score all stimuli and windows: Z(ch, stim, win) = (meanAct - mu_spont) / sd_spont
meanAct_all = R.meanAct;                                  % [Nchan x Nstim x Nwin]

% Broadcast subtract/divide (works in recent MATLAB). If your MATLAB is older, use bsxfun.
Z = (meanAct_all - mu_spont) ./ sd_spont;                 % [Nchan x Nstim x Nwin]

