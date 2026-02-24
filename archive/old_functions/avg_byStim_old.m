function [meanAct, meanSqAct, nTrialsPerStim, stimList, winIdx] = ...
    avg_byStim(m1, m2, timeWindows, varargin)
%AVG_NORMMUA_BYSTIM_TIMEWINDOWS_TWOFILES
% Chunked averaging of m1.normMUA (chan x trial x time) by stimulus and time windows.
%
% Inputs
%   m1          : matfile object for file containing normMUA
%                 (accessed as m1.normMUA(ch, trial, t))
%   m2          : matfile object for file containing ALLMAT and tb
%                 ALLMAT = m2.ALLMAT; tb = m2.tb;
%   timeWindows : [Nwin x 2] matrix of windows. Interpreted as ms by default.
%
% Name-value options
%   'stimCol'     : column in ALLMAT with stimulus id (default 1)
%   'correctCol'  : column in ALLMAT coding correctness (default 9)
%   'dayCol'      : column in ALLMAT with recording day id (default 11)
%   'days'        : numeric vector of allowed day ids; empty => include all days (default [])
%   'correctVal'  : value in correctCol indicating correct (default 1)
%   'onlyCorrect' : if true, include only correct trials (default true)
%   'stimList'    : explicit list of stimulus ids to average (default [])
%   'chunkTrials' : number of trials per chunk read from disk (default 200)
%   'windowsIn'   : 'ms' or 'idx' (default 'ms')
%
% Outputs
%   meanAct        : [NChans x NStim x Nwin] mean activity per stimulus and window
%   meanSqAct      : [NChans x NStim x Nwin] mean squared activity per stimulus and window
%   nTrialsPerStim : [NStim x 1] number of included trials per stimulus (across all windows)
%   stimList       : [NStim x 1] list of stimuli included
%   winIdx         : [Nwin x 2] resolved indices into tb for each window
%
% Notes
% - This function reads m1.normMUA in trial-chunks to reduce memory usage.
% - Averaging uses simple sums / sums of squares; NaNs are omitted per window.

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('m1', @(x)isobject(x) || isa(x,'matlab.io.MatFile'));
p.addRequired('m2', @(x)isobject(x) || isa(x,'matlab.io.MatFile'));
p.addRequired('timeWindows', @(x)isnumeric(x) && size(x,2)==2);

p.addParameter('stimCol',     1,    @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('correctCol',  9,    @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('correctVal',  1,    @(x)isnumeric(x) && isscalar(x));
p.addParameter('onlyCorrect', true, @(x)islogical(x) && isscalar(x));
p.addParameter('dayCol',      11,   @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('days',        [],   @(x)isnumeric(x) && (isempty(x) || isvector(x)));
p.addParameter('stimList',    [],   @(x)isnumeric(x) && (isempty(x) || isvector(x)));
p.addParameter('chunkTrials', 200,  @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('windowsIn',  'ms',  @(s)ischar(s) || isstring(s));

p.parse(m1, m2, timeWindows, varargin{:});
opt = p.Results;

windowsIn = lower(string(opt.windowsIn));

% -------------------- sizes from m1 --------------------
[NChans, NTrials, NTimes] = size(m1, 'normMUA');

% -------------------- load ALLMAT and tb from m2 --------------------
ALLMAT = m2.ALLMAT;     % should be [NTrials x ...]
tb     = m2.tb;         % should be [1 x NTimes] or [NTimes x 1] in ms

if size(ALLMAT,1) ~= NTrials
    error('Mismatch: m1.normMUA has NTrials=%d but m2.ALLMAT has %d rows.', NTrials, size(ALLMAT,1));
end

tb = tb(:)'; % row
if numel(tb) ~= NTimes
    error('Mismatch: m1.normMUA has NTimes=%d but m2.tb has length %d.', NTimes, numel(tb));
end

% -------------------- trial -> stimulus mapping --------------------
stimPerTrial = ALLMAT(:, opt.stimCol);

% correctness filter
isCorrect = true(NTrials,1);
if opt.onlyCorrect
    isCorrect = (ALLMAT(:, opt.correctCol) == opt.correctVal);
end

% day filter (recording days)
isDay = true(NTrials,1);
if ~isempty(opt.days)
    if size(ALLMAT,2) < opt.dayCol
        error('ALLMAT has only %d columns; cannot use dayCol=%d.', size(ALLMAT,2), opt.dayCol);
    end
    isDay = ismember(ALLMAT(:, opt.dayCol), opt.days(:));
end

% combined inclusion mask
isInclude = isCorrect & isDay;

% choose stimList
if isempty(opt.stimList)
    stimList = unique(stimPerTrial(isInclude), 'stable');
else
    stimList = opt.stimList(:);
end
NStim = numel(stimList);

% map each trial to 1..NStim, set excluded trials to 0
[~, stimPosPerTrial] = ismember(stimPerTrial, stimList);
stimPosPerTrial(~isInclude) = 0;  % drop excluded trials

% -------------------- resolve time windows -> indices --------------------
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
meanAct   = nan(NChans, NStim, Nt);
meanSqAct = nan(NChans, NStim, Nt);

sumAct    = zeros(NChans, NStim, Nt, 'double');
sumSqAct  = zeros(NChans, NStim, Nt, 'double');
sumN      = zeros(1,      NStim, Nt, 'double');

% count trials per stim (independent of window)
nTrialsPerStim = accumarray(stimPosPerTrial(stimPosPerTrial>0), 1, [NStim 1], @sum, 0);

% -------------------- chunked read over trials --------------------
chunk = opt.chunkTrials;
tStarts = 1:chunk:NTrials;

for k = 1:numel(tStarts)
    t0 = tStarts(k);
    t1 = min(NTrials, t0 + chunk - 1);
    trials = t0:t1;
    
    fprintf('Processed trials %d-%d of %d (%.1f%%)\n', ...
        t0, t1, NTrials, 100*t1/NTrials);

    % read chunk: [NChans x NChunk x NTimes]
    X = m1.normMUA(:, trials, :);

    % loop windows, compute sums per stim
    for w = 1:Nt
        ii = winIdx(w,1):winIdx(w,2);

        % average within window per trial (omitnan over time)
        Xw  = mean(X(:,:,ii), 3, 'omitnan');        % [NChans x NChunk]
        Xw2 = mean(X(:,:,ii).^2, 3, 'omitnan');     % [NChans x NChunk]

        % accumulate per stimulus (within this chunk)
        stimPos = stimPosPerTrial(trials);          % [NChunk x 1]
        good = stimPos > 0;

        if any(good)
            stimPosG = stimPos(good);
            XwG  = Xw(:, good);
            Xw2G = Xw2(:, good);

            for s = 1:NStim
                idxS = (stimPosG == s);
                if any(idxS)
                    sumAct(:,s,w)   = sumAct(:,s,w)   + sum(XwG(:,idxS), 2, 'omitnan');
                    sumSqAct(:,s,w) = sumSqAct(:,s,w) + sum(Xw2G(:,idxS), 2, 'omitnan');
                    sumN(1,s,w)     = sumN(1,s,w)     + sum(~isnan(XwG(1,idxS))); %#ok<NASGU>
                end
            end
        end
    end
end

% -------------------- finalize means --------------------
for w = 1:Nt
    for s = 1:NStim
        n = nTrialsPerStim(s);
        if n > 0
            meanAct(:,s,w)   = sumAct(:,s,w)   ./ n;
            meanSqAct(:,s,w) = sumSqAct(:,s,w) ./ n;
        end
    end
end
end
