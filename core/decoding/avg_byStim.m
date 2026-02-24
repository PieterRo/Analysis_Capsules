function [meanAct, meanSqAct, nTrialsPerStim, stimList] = avg_byStim(m1, m2, timeWindows, varargin)
%AVG_BYSTIM  Average activity by stimulus (stimulus ID is the 2nd dimension).
%
% Usage:
%   [meanAct, meanSqAct, nTrials, stimList] = avg_byStim(m1, m2, timeWindows, 'days', [1 2]);
%
% Required inputs:
%   m1: matfile with variable normMUA  [NChans x NTrials x NTimes]
%   m2: matfile with variables ALLMAT [NTrials x ...] and tb [1 x NTimes]
%   timeWindows: [Nwin x 2] in ms (default), or indices if windowsIn='idx'
%
% Guarantees:
%   - meanAct(:, stimID, win) corresponds directly to stimulus number stimID.
%   - If stimuli are 1..384 contiguous, output second dimension is 384.
%   - Per-channel non-NaN denominators (no NaN/denominator bias).
%   - MatFile-safe indexing (contiguous trial chunks).
%
% Name-value options:
%   'stimCol'     : ALLMAT column holding stimulus number (default 1)
%   'onlyCorrect' : whether to include only correct trials (default true)
%   'correctCol'  : ALLMAT column holding correctness flag (default 9)
%   'correctVal'  : value indicating correct (default 1)
%   'days'        : allowed day ids (default [])
%   'dayCol'      : ALLMAT column holding day id (default 11)
%   'chunkTrials' : chunk size for reading trials from matfile (default 200)
%   'windowsIn'   : 'ms' (default) or 'idx'
%   'verbose'     : true/false (default true)
%   'expectedMaxStim' : expected maximum stimulus ID (default 384). If empty, uses max found.
%
% Outputs:
%   meanAct       : [NChans x NStim x Nwin]
%   meanSqAct     : [NChans x NStim x Nwin]
%   nTrialsPerStim: [NStim x 1] included trials per stim (after filters), regardless of NaNs
%   stimList      : [NStim x 1] (1:NStim)'

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('m1', @(x) isa(x,'matlab.io.MatFile'));
p.addRequired('m2', @(x) isa(x,'matlab.io.MatFile'));
p.addRequired('timeWindows', @(x) isnumeric(x) && size(x,2)==2);

p.addParameter('stimCol',     1,    @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('onlyCorrect', true, @(x)islogical(x) && isscalar(x));
p.addParameter('correctCol',  9,    @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('correctVal',  1,    @(x)isnumeric(x) && isscalar(x));
p.addParameter('days',        [],   @(x)isnumeric(x) && (isempty(x) || isvector(x)));
p.addParameter('dayCol',      11,   @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('chunkTrials', 200,  @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('windowsIn',  'ms',  @(s)ischar(s) || isstring(s));
p.addParameter('verbose',     true, @(x)islogical(x) && isscalar(x));
p.addParameter('expectedMaxStim', 384, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));

p.parse(m1, m2, timeWindows, varargin{:});
opt = p.Results;
windowsIn = lower(string(opt.windowsIn));

% -------------------- sizes from m1 --------------------
[NChans, NTrials, NTimes] = size(m1, 'normMUA');

% -------------------- load ALLMAT and tb --------------------
ALLMAT = m2.ALLMAT;
tb     = m2.tb;

if size(ALLMAT,1) ~= NTrials
    error('Mismatch: m1.normMUA has NTrials=%d but m2.ALLMAT has %d rows.', ...
        NTrials, size(ALLMAT,1));
end

tb = tb(:)'; % force row
if numel(tb) ~= NTimes
    error('Mismatch: m1.normMUA has NTimes=%d but m2.tb has length %d.', ...
        NTimes, numel(tb));
end

% -------------------- determine included trials --------------------
if size(ALLMAT,2) < opt.stimCol
    error('ALLMAT has only %d columns; stimCol=%d is out of range.', size(ALLMAT,2), opt.stimCol);
end
stimPerTrial = ALLMAT(:, opt.stimCol);

isCorrect = true(NTrials,1);
if opt.onlyCorrect
    if size(ALLMAT,2) < opt.correctCol
        error('ALLMAT has only %d columns; correctCol=%d is out of range.', size(ALLMAT,2), opt.correctCol);
    end
    isCorrect = (ALLMAT(:, opt.correctCol) == opt.correctVal);
end

isDay = true(NTrials,1);
if ~isempty(opt.days)
    if size(ALLMAT,2) < opt.dayCol
        error('ALLMAT has only %d columns; dayCol=%d is out of range.', size(ALLMAT,2), opt.dayCol);
    end
    isDay = ismember(ALLMAT(:, opt.dayCol), opt.days(:));
end

isInclude = isCorrect & isDay;

% -------------------- check stimulus IDs --------------------
stimIDsIncl = unique(stimPerTrial(isInclude));
stimIDsIncl = stimIDsIncl(:);
stimIDsIncl = stimIDsIncl(isfinite(stimIDsIncl));

if isempty(stimIDsIncl)
    error('No trials left after filters (onlyCorrect/days).');
end

minStim = min(stimIDsIncl);
maxStimFound = max(stimIDsIncl);
nUnique = numel(unique(stimIDsIncl));

expectedMax = opt.expectedMaxStim;
if isempty(expectedMax)
    NStim = maxStimFound;
else
    NStim = expectedMax;
end

% Warn if not contiguous 1..NStim (based on included trials)
expectedSet = (1:NStim)';
if ~(minStim == 1 && maxStimFound <= NStim && nUnique == numel(intersect(stimIDsIncl, expectedSet)))
    warning(['avg_byStim:StimulusIDsNotContiguous: Included-trial stimulus IDs are not exactly ' ...
             'a contiguous subset of 1..%d. min=%g max=%g nUnique=%d. ' ...
             'The output will still use stimID indexing, but some stim columns may be all-NaN.'], ...
             NStim, minStim, maxStimFound, nUnique);
end

stimList = (1:NStim)';

% Map each trial to stimID position directly (0 for excluded or out-of-range)
stimPosPerTrial = zeros(NTrials,1);
stimOK = isfinite(stimPerTrial) & stimPerTrial>=1 & stimPerTrial<=NStim & (floor(stimPerTrial)==stimPerTrial);
stimPosPerTrial(stimOK) = stimPerTrial(stimOK);
stimPosPerTrial(~isInclude) = 0;

% Included trials per stim (regardless of NaNs later)
nTrialsPerStim = accumarray(stimPosPerTrial(stimPosPerTrial>0), 1, [NStim 1], @sum, 0);

% -------------------- resolve windows to indices --------------------
Nt = size(timeWindows,1);
winIdx = zeros(Nt,2);

if windowsIn == "ms"
    for w = 1:Nt
        t1 = timeWindows(w,1);
        t2 = timeWindows(w,2);
        i1 = find(tb >= t1, 1, 'first');
        i2 = find(tb <= t2, 1, 'last');
        if isempty(i1) || isempty(i2) || i2 < i1
            error('Bad ms window [%g %g]: does not overlap tb range [%g %g].', t1, t2, tb(1), tb(end));
        end
        winIdx(w,:) = [i1 i2];
    end
elseif windowsIn == "idx"
    winIdx = round(timeWindows);
    if any(winIdx(:) < 1) || any(winIdx(:) > NTimes)
        error('Index windows must be within [1..%d].', NTimes);
    end
else
    error("windowsIn must be 'ms' or 'idx'.");
end

% -------------------- accumulators --------------------
sumAct    = zeros(NChans, NStim, Nt, 'double');
sumSqAct  = zeros(NChans, NStim, Nt, 'double');
nNonNaN   = zeros(NChans, NStim, Nt, 'double'); % per-channel denominators

chunk = opt.chunkTrials;
tStarts = 1:chunk:NTrials;

% -------------------- main loop: contiguous read per chunk --------------------
for k = 1:numel(tStarts)
    t0 = tStarts(k);
    t1 = min(NTrials, t0 + chunk - 1);
    trials = t0:t1; % contiguous -> MatFile-safe

    if opt.verbose
        fprintf('avg_byStim: trials %d-%d of %d (%.1f%%)\n', t0, t1, NTrials, 100*t1/NTrials);
    end

    stimPosChunk = stimPosPerTrial(trials); % [nChunk x 1]
    good = stimPosChunk > 0;
    if ~any(good)
        continue;
    end

    % Read contiguous chunk once, then subselect in memory
    Xchunk = m1.normMUA(:, trials, :);  % [NChans x nChunk x NTimes]
    X = Xchunk(:, good, :);             % [NChans x nGood x NTimes]
    stimPosG = stimPosChunk(good);      % [nGood x 1]

    uStim = unique(stimPosG);

    for w = 1:Nt
        ii = winIdx(w,1):winIdx(w,2);

        % Per-trial window moments (mean over time, per trial)
        Xw  = mean(X(:,:,ii), 3, 'omitnan');      % [NChans x nGood]
        Xw2 = mean(X(:,:,ii).^2, 3, 'omitnan');   % [NChans x nGood]

        for us = 1:numel(uStim)
            s = uStim(us);  % s is stimID (1..NStim)
            idxS = (stimPosG == s);
            if ~any(idxS), continue; end

            Xs  = Xw(:, idxS);
            Xs2 = Xw2(:, idxS);

            sumAct(:,s,w)   = sumAct(:,s,w)   + sum(Xs,  2, 'omitnan');
            sumSqAct(:,s,w) = sumSqAct(:,s,w) + sum(Xs2, 2, 'omitnan');

            % Per-channel non-NaN trial counts
            nNonNaN(:,s,w)  = nNonNaN(:,s,w)  + sum(~isnan(Xs), 2);
        end
    end
end

% -------------------- finalize means (correct denominators) --------------------
meanAct   = sumAct   ./ nNonNaN;
meanSqAct = sumSqAct ./ nNonNaN;

meanAct(nNonNaN==0)   = NaN;
meanSqAct(nNonNaN==0) = NaN;

end