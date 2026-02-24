function OUT = attention_modulation_V1_3bin(DATA, Tall_V1, SNR, opts)
% ATTENTION_MODULATION_V1_3BIN
% Per-site object-based attention modulation in V1 (sites 1:512),
% using 3 time bins and focusing on timeIdx=3 (300–500 ms) by default.
%
% DATA must contain:
%   DATA.meanAct    [nSites x 384 x 3]  or [1024 x 384 x 3]
%   DATA.meanSqAct  [nSites x 384 x 3]
%   DATA.nTrials    [1 x 384] or [384 x 1] or [nSites x 384]
%   DATA.timeWindows [3 x 2] optional
%
% Tall_V1: 1x384 struct array. For each stimulus:
%   Tall_V1(stim).T is a table with >=512 rows and fields:
%     - assignment   : 'target'/'distractor'/'background'
%     - overlap      : 0/1  (optional; if missing assumed 0)
%     - center_color : e.g. 'yellowArm' or 'purple'
%
% SNR fields used for normalization:
%   muSpont, muYellowEarly, muYellowLate, muPurpleEarly, muPurpleLate
%
% Normalization per site:
%   Y = (X - muSpont) / scale, where scale = max(mu* fields) - muSpont
% For second moment:
%   E[Y^2] = (E[X^2] - 2*b*E[X] + b^2) / scale^2
%
% Color-balancing across conditions:
%   For each color c in {Y,P}, define effective weight w_c = min(N_Tc, N_Dc).
%   Then combine colors using weights wY,wP (so target and distractor are matched per color).
%
% Outputs (all [512 x 1] unless stated):
%   OUT.index, OUT.dprime
%   OUT.muT, OUT.muD, OUT.varT, OUT.varD
%   OUT.pValueTD, OUT.tStatTD, OUT.dofTD  (Target vs Distractor, matched sets)
%   OUT.pValueY,  OUT.pValueP             (color-specific matched sets)
%   OUT.N_TY, OUT.N_TP, OUT.N_DY, OUT.N_DP (raw trial totals)
%   OUT.wY, OUT.wP (effective matched weights per color)
%   OUT.validSite (logical)
%   OUT.timeWindowUsed [1x2]
%   OUT.stimMaskUsed [1x384]

%% Options
if nargin < 4, opts = struct(); end
if ~isfield(opts,'v1Sites'),        opts.v1Sites = 1:512; end
if ~isfield(opts,'timeIdx'),        opts.timeIdx = 3; end
if ~isfield(opts,'excludeOverlap'), opts.excludeOverlap = true; end
if ~isfield(opts,'epsDen'),         opts.epsDen = 1e-6; end
if ~isfield(opts,'verbose'),        opts.verbose = true; end
if ~isfield(opts,'debugSite'),      opts.debugSite = []; end
if ~isfield(opts,'stimMask'),       opts.stimMask = []; end
if ~isfield(opts,'stimIdx'),        opts.stimIdx = []; end

v1Sites = opts.v1Sites(:)';
nV1     = numel(v1Sites);
timeIdx = opts.timeIdx;
dbgRow = [];
if ~isempty(opts.debugSite)
    dbgRow = find(v1Sites == opts.debugSite, 1, 'first');
end

% Stimulus subset selection: either stimMask (logical 1x384) or stimIdx (indices).
stimMask = true(1,384);
if ~isempty(opts.stimMask)
    if ~(islogical(opts.stimMask) && numel(opts.stimMask) == 384)
        error('opts.stimMask must be a logical vector with 384 elements.');
    end
    stimMask = reshape(opts.stimMask, 1, 384);
end
if ~isempty(opts.stimIdx)
    if ~(isnumeric(opts.stimIdx) && isvector(opts.stimIdx))
        error('opts.stimIdx must be a numeric vector of stimulus indices (1..384).');
    end
    idx = opts.stimIdx(:)';
    if any(idx < 1 | idx > 384 | idx ~= floor(idx))
        error('opts.stimIdx contains invalid stimulus indices.');
    end
    tmp = false(1,384);
    tmp(idx) = true;
    stimMask = stimMask & tmp;
end

DBG = struct();
if ~isempty(dbgRow)
    DBG.site = opts.debugSite;
    DBG.rowInV1 = dbgRow;
    DBG.TY_stim = []; DBG.TY_w = [];
    DBG.TP_stim = []; DBG.TP_w = [];
    DBG.DY_stim = []; DBG.DY_w = [];
    DBG.DP_stim = []; DBG.DP_w = [];
end

%% Validate inputs
assert(isfield(DATA,'meanAct') && isfield(DATA,'meanSqAct') && isfield(DATA,'nTrials'), ...
    'DATA must contain meanAct, meanSqAct, nTrials.');
assert(numel(Tall_V1) == 384, 'Tall_V1 must be 1x384.');

meanAct   = DATA.meanAct;
meanSqAct = DATA.meanSqAct;

assert(size(meanAct,2) == 384 && size(meanSqAct,2) == 384, 'Expected 384 stimuli.');
assert(size(meanAct,3) >= timeIdx && size(meanSqAct,3) >= timeIdx, 'opts.timeIdx out of range.');

nTrials = DATA.nTrials;
if isvector(nTrials)
    nTrials = nTrials(:)'; % 1x384
    assert(numel(nTrials) == 384, 'DATA.nTrials vector must have 384 elements.');
    perSiteTrials = false;
elseif ismatrix(nTrials) && size(nTrials,2)==384
    perSiteTrials = true;
else
    error('DATA.nTrials must be vector(384) or matrix(nSites x 384).');
end

%% Build normalization constants (V1 sites)
b = SNR.muSpont(:);
topMat = [SNR.muYellowEarly(:), SNR.muYellowLate(:), SNR.muPurpleEarly(:), SNR.muPurpleLate(:)];
muTop = max(topMat, [], 2);

if numel(b) < max(v1Sites) || numel(muTop) < max(v1Sites)
    error('SNR fields are smaller than requested V1 site indices.');
end

b     = b(v1Sites);
muTop = muTop(v1Sites);

scale = muTop - b;
scale(~isfinite(scale) | scale <= 0) = NaN;

%% Accumulators for trial-weighted normalized moments per bucket
sumY_TY  = zeros(nV1,1); sumY2_TY = zeros(nV1,1); N_TY = zeros(nV1,1);
sumY_TP  = zeros(nV1,1); sumY2_TP = zeros(nV1,1); N_TP = zeros(nV1,1);
sumY_DY  = zeros(nV1,1); sumY2_DY = zeros(nV1,1); N_DY = zeros(nV1,1);
sumY_DP  = zeros(nV1,1); sumY2_DP = zeros(nV1,1); N_DP = zeros(nV1,1);

if opts.verbose
    fprintf('V1 attention modulation: timeIdx=%d over 384 stimuli...\n', timeIdx);
end

%% Main loop over stimuli
for stim = 1:384
    if ~stimMask(stim)
        continue;
    end

    % Trials for this stimulus (per-site or scalar)
    if ~perSiteTrials
        nTr_vec = repmat(nTrials(stim), nV1, 1);
    else
        if size(nTrials,1) < max(v1Sites)
            error('nTrials matrix has insufficient rows for V1 sites.');
        end
        nTr_vec = nTrials(v1Sites, stim);
    end

    if all(~isfinite(nTr_vec) | nTr_vec<=0)
        continue;
    end

    % Tall table
    if ~isfield(Tall_V1(stim),'T') || ~istable(Tall_V1(stim).T) || isempty(Tall_V1(stim).T)
        continue;
    end
    Ttab = Tall_V1(stim).T;

    % Required fields
    if ~ismember('assignment',   Ttab.Properties.VariableNames) || ...
       ~ismember('center_color', Ttab.Properties.VariableNames)
        error('Tall_V1(%d).T missing assignment and/or center_color.', stim);
    end
    hasOverlap = ismember('overlap', Ttab.Properties.VariableNames);

    % Ensure enough rows
    if height(Ttab) < max(v1Sites)
        continue;
    end

    % Extract metadata for V1 sites
    asg = Ttab.assignment(v1Sites);
    cc  = Ttab.center_color(v1Sites);
    if hasOverlap
        ov = Ttab.overlap(v1Sites);
    else
        ov = zeros(nV1,1);
    end

    asgStr = toCellStr(asg);
    colStr = toCellStr(cc);

    isBG = cellfun(@(x) strcmpi(x,'background'), asgStr);
    isT  = cellfun(@(x) strcmpi(x,'target'), asgStr);
    isD  = cellfun(@(x) strcmpi(x,'distractor'), asgStr);

    if opts.excludeOverlap
        isOV = (ov ~= 0);
    else
        isOV = false(nV1,1);
    end

    % Determine RF-center color
    isY = cellfun(@(x) contains(lower(x),'yellow'), colStr);
    isP = cellfun(@(x) contains(lower(x),'purple'), colStr);
    goodColor = isY | isP;

    % Pull moments
    EX  = squeeze(meanAct(v1Sites, stim, timeIdx));   % E[X]
    EX2 = squeeze(meanSqAct(v1Sites, stim, timeIdx)); % E[X^2]

    % Normalize moments
    bi = b(:); si = scale(:);
    EY  = (EX - bi) ./ si;
    EY2 = (EX2 - 2.*bi.*EX + (bi.^2)) ./ (si.^2);

    goodNorm = isfinite(EY) & isfinite(EY2) & isfinite(nTr_vec) & (nTr_vec > 0);

    useBase = goodNorm & ~isBG(:) & ~isOV(:) & goodColor(:);

    % --- Accumulate into 4 buckets (trial-weighted) ---
    m = (isT(:) & isY(:) & useBase);
    if any(m)
        w = nTr_vec(m);
        sumY_TY(m)  = sumY_TY(m)  + w .* EY(m);
        sumY2_TY(m) = sumY2_TY(m) + w .* EY2(m);
        N_TY(m)     = N_TY(m)     + w;
        if ~isempty(dbgRow) && m(dbgRow)
            DBG.TY_stim(end+1,1) = stim; %#ok<AGROW>
            DBG.TY_w(end+1,1) = nTr_vec(dbgRow); %#ok<AGROW>
        end
    end

    m = (isT(:) & isP(:) & useBase);
    if any(m)
        w = nTr_vec(m);
        sumY_TP(m)  = sumY_TP(m)  + w .* EY(m);
        sumY2_TP(m) = sumY2_TP(m) + w .* EY2(m);
        N_TP(m)     = N_TP(m)     + w;
        if ~isempty(dbgRow) && m(dbgRow)
            DBG.TP_stim(end+1,1) = stim; %#ok<AGROW>
            DBG.TP_w(end+1,1) = nTr_vec(dbgRow); %#ok<AGROW>
        end
    end

    m = (isD(:) & isY(:) & useBase);
    if any(m)
        w = nTr_vec(m);
        sumY_DY(m)  = sumY_DY(m)  + w .* EY(m);
        sumY2_DY(m) = sumY2_DY(m) + w .* EY2(m);
        N_DY(m)     = N_DY(m)     + w;
        if ~isempty(dbgRow) && m(dbgRow)
            DBG.DY_stim(end+1,1) = stim; %#ok<AGROW>
            DBG.DY_w(end+1,1) = nTr_vec(dbgRow); %#ok<AGROW>
        end
    end

    m = (isD(:) & isP(:) & useBase);
    if any(m)
        w = nTr_vec(m);
        sumY_DP(m)  = sumY_DP(m)  + w .* EY(m);
        sumY2_DP(m) = sumY2_DP(m) + w .* EY2(m);
        N_DP(m)     = N_DP(m)     + w;
        if ~isempty(dbgRow) && m(dbgRow)
            DBG.DP_stim(end+1,1) = stim; %#ok<AGROW>
            DBG.DP_w(end+1,1) = nTr_vec(dbgRow); %#ok<AGROW>
        end
    end

end

%% Bucket moments
[mu_TY, m2_TY] = momentsFromSums(sumY_TY, sumY2_TY, N_TY);
[mu_TP, m2_TP] = momentsFromSums(sumY_TP, sumY2_TP, N_TP);
[mu_DY, m2_DY] = momentsFromSums(sumY_DY, sumY2_DY, N_DY);
[mu_DP, m2_DP] = momentsFromSums(sumY_DP, sumY2_DP, N_DP);

% Bucket variances (may be NaN when N==0)
var_TY = m2_TY - mu_TY.^2;
var_TP = m2_TP - mu_TP.^2;
var_DY = m2_DY - mu_DY.^2;
var_DP = m2_DP - mu_DP.^2;

% Numerical guards
[var_TY, var_TP, var_DY, var_DP] = guardVars(var_TY, var_TP, var_DY, var_DP);

%% Color matching across conditions (pure attention modulation)
% Effective matched weights per color:
wY = min(N_TY, N_DY);
wP = min(N_TP, N_DP);

% At least one color must have matched trials to define T vs D
validSite = (wY > 0) | (wP > 0);

% Combine colors for Target and Distractor using matched weights
% Weighted mean across colors:
muT = weighted2(mu_TY, mu_TP, wY, wP);
muD = weighted2(mu_DY, mu_DP, wY, wP);

% Weighted second moment across colors:
m2T = weighted2(m2_TY, m2_TP, wY, wP);
m2D = weighted2(m2_DY, m2_DP, wY, wP);

varT = m2T - muT.^2;
varD = m2D - muD.^2;
[varT, varD] = guardVars(varT, varD);

%% Index and d-prime
den = abs((muT + muD)/2);
den(den < opts.epsDen) = opts.epsDen;
index = (muT - muD) ./ den;

sdP = sqrt(0.5*(varT + varD));
sdP(sdP < opts.epsDen) = opts.epsDen;
dprime = (muT - muD) ./ sdP;

%% p-values (Target vs Distractor) on matched sets
Nmatch = wY + wP;  % equal matched counts in T and D by construction
[pValueTD, tStatTD, dofTD] = welchPvalue(muT, varT, Nmatch, muD, varD, Nmatch, opts.epsDen);

% Color-specific matched-set p-values
[pValueY, ~, ~] = welchPvalue(mu_TY, var_TY, wY, mu_DY, var_DY, wY, opts.epsDen);
[pValueP, ~, ~] = welchPvalue(mu_TP, var_TP, wP, mu_DP, var_DP, wP, opts.epsDen);

% Apply validity
index(~validSite)  = NaN;
dprime(~validSite) = NaN;
muT(~validSite) = NaN; muD(~validSite) = NaN;
varT(~validSite) = NaN; varD(~validSite) = NaN;
pValueTD(~validSite) = NaN;
tStatTD(~validSite)  = NaN;
dofTD(~validSite)    = NaN;
pValueY(~validSite)  = NaN;
pValueP(~validSite)  = NaN;

%% Package output
OUT = struct();
OUT.index  = index;
OUT.dprime = dprime;
OUT.muT    = muT;
OUT.muD    = muD;
OUT.varT   = varT;
OUT.varD   = varD;

OUT.N_TY = N_TY; OUT.N_TP = N_TP;
OUT.N_DY = N_DY; OUT.N_DP = N_DP;

OUT.wY = wY; OUT.wP = wP;
OUT.validSite = validSite;
OUT.pValueTD = pValueTD;
OUT.tStatTD  = tStatTD;
OUT.dofTD    = dofTD;
OUT.pValueY  = pValueY;
OUT.pValueP  = pValueP;
OUT.stimMaskUsed = stimMask;
if ~isempty(dbgRow)
    OUT.debug = DBG;
end

if isfield(DATA,'timeWindows') && size(DATA.timeWindows,1) >= timeIdx
    OUT.timeWindowUsed = DATA.timeWindows(timeIdx,:);
else
    OUT.timeWindowUsed = [NaN NaN];
end

if opts.verbose
    fprintf('Done. Valid sites (>=1 matched color): %d / %d\n', nnz(validSite), nV1);
end

end

%% ===================== Local helpers =====================

function C = toCellStr(x)
% Convert table column to cell array of char vectors robustly
if iscell(x)
    C = x;
elseif isstring(x)
    C = cellstr(x);
elseif iscategorical(x)
    C = cellstr(x);
elseif ischar(x)
    C = cellstr(x);
else
    C = cellstr(string(x));
end
C = C(:);
for i = 1:numel(C)
    if iscell(C{i}) && numel(C{i})==1
        C{i} = C{i}{1};
    end
end
end

function [mu, m2] = momentsFromSums(sumY, sumY2, N)
mu = sumY ./ N;
m2 = sumY2 ./ N;
mu(N<=0) = NaN;
m2(N<=0) = NaN;
end

function out = weighted2(a, b, wa, wb)
% Weighted average of two quantities a and b with weights wa and wb.
% Handles wa+wb==0 -> NaN
w = wa + wb;
out = (wa .* a + wb .* b) ./ w;
out(w<=0) = NaN;
end

function [p, t, dof] = welchPvalue(mu1, var1, n1, mu2, var2, n2, epsDen)
% Welch t-test from summary moments/counts (site-wise vectorized).
% Returns NaN when inputs are insufficient (e.g., n<=1).

mu1 = mu1(:); var1 = var1(:); n1 = n1(:);
mu2 = mu2(:); var2 = var2(:); n2 = n2(:);

bad = ~isfinite(mu1) | ~isfinite(mu2) | ~isfinite(var1) | ~isfinite(var2) | ...
      ~isfinite(n1)  | ~isfinite(n2)  | (n1 <= 1)       | (n2 <= 1);

var1 = max(var1, 0);
var2 = max(var2, 0);

se2 = var1 ./ n1 + var2 ./ n2;
se2(se2 < epsDen^2) = epsDen^2;
se = sqrt(se2);

t = (mu1 - mu2) ./ se;

num = se2.^2;
den = (var1.^2) ./ (n1.^2 .* max(n1-1,1)) + (var2.^2) ./ (n2.^2 .* max(n2-1,1));
dof = num ./ den;

p = NaN(size(t));
use = ~bad & isfinite(t) & isfinite(dof) & (dof > 0);
if any(use)
    if exist('tcdf', 'file') == 2 || exist('tcdf', 'builtin') == 5
        p(use) = 2 * tcdf(-abs(t(use)), dof(use));
    else
        % Fallback: large-sample normal approximation
        p(use) = erfc(abs(t(use)) ./ sqrt(2));
    end
end

t(bad) = NaN;
dof(bad) = NaN;
end

function varargout = guardVars(varargin)
% Clamp tiny negative variances due to numerical precision to 0; else NaN.
varargout = varargin;
for k = 1:nargin
    v = varargin{k};
    v(v < 0 & v > -1e-10) = 0;
    v(v < 0) = NaN;
    varargout{k} = v;
end
end
