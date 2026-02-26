function OUT = compute_projected_delta_points_allbins( ...
    stimID_example, Tall_V1, ALLCOORDS, RTAB384, R_resp, SNRnorm, varargin)
% COMPUTE_PROJECTED_DELTA_POINTS_ALLBINS
% Build point-wise signed T-D values for all response bins and project them
% to post-affine coordinates on one example stimulus geometry.
%
% Output is intended for numeric analysis (not plotting).
%
% Required inputs:
%   stimID_example : scalar stimulus index used as canonical geometry
%   Tall_V1        : 1x384 struct array with Tall_V1(stim).T tables
%   ALLCOORDS      : stimulus geometry struct (fields stim_1..stim_384)
%   RTAB384        : stimulus metadata matrix/table (col 7 = width px)
%   R_resp         : response struct (many bins), used by compute_deltaQ_quartet
%   SNRnorm        : normalization struct for compute_deltaQ_quartet
%
% Name/value options:
%   'siteRange'       (default 1:512)
%   'excludeOverlap'  (default true)
%   'stimIdx'         (default 1:384)  % used only to choose which quartets to include
%   'sigSiteMask'     (default [])  logical mask over siteRange; true sites kept
%   'saveFile'        (default '')  full path to save OUT (-v7.3)
%   'verbose'         (default true)
%
% OUT fields:
%   OUT.meta
%   OUT.bins(tb).timeBin
%   OUT.bins(tb).timeWindow
%   OUT.bins(tb).x_px
%   OUT.bins(tb).y_px
%   OUT.bins(tb).delta
%   OUT.bins(tb).siteIdx
%   OUT.bins(tb).quartetIdx
%   OUT.bins(tb).stimIdxSource   % stimulus table used for geometry in this quartet
%   OUT.bins(tb).stimMembersUsed % [N x 4] quartet members for each point
%   OUT.bins(tb).stream          % 1=target-assigned point, 2=distractor-assigned point

p = inputParser;
p.addParameter('siteRange', 1:512, @(x) isnumeric(x) && isvector(x));
p.addParameter('excludeOverlap', true, @(x) islogical(x) && isscalar(x));
p.addParameter('stimIdx', 1:384, @(x) isnumeric(x) && isvector(x));
p.addParameter('sigSiteMask', [], @(x) isempty(x) || (islogical(x) && isvector(x)));
p.addParameter('saveFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

siteRange = double(opt.siteRange(:)');
stimList = double(opt.stimIdx(:)');
nFrames = size(R_resp.timeWindows, 1);

if isempty(stimList)
    error('compute_projected_delta_points_allbins:EmptyStimList', ...
        'stimIdx cannot be empty.');
end
if any(stimList < 1 | stimList > 384 | stimList ~= floor(stimList))
    error('compute_projected_delta_points_allbins:BadStimIdx', ...
        'stimIdx must contain integer values in [1..384].');
end

if isempty(opt.sigSiteMask)
    sigMask = true(numel(siteRange),1);
else
    sigMask = opt.sigSiteMask(:);
    if numel(sigMask) ~= numel(siteRange)
        error('compute_projected_delta_points_allbins:BadSigMaskSize', ...
            'sigSiteMask must have %d elements (matching siteRange).', numel(siteRange));
    end
end

% Canonical geometry from example stimulus.
W = 1024; H = 768;
toPx = @(q) [q(1) + W/2, H/2 - q(2)];
fName = sprintf('stim_%d', stimID_example);
if ~isfield(ALLCOORDS, fName)
    error('compute_projected_delta_points_allbins:MissingStimField', ...
        'ALLCOORDS missing field %s.', fName);
end
s = double(ALLCOORDS.(fName).s(:))';
tFig = double(ALLCOORDS.(fName).t_fig(:))';
tBack = double(ALLCOORDS.(fName).t_back(:))';
s_px = toPx(s);
tT = toPx(tFig);
tD = toPx(tBack);
uT = (tT - s_px) / norm(tT - s_px);
uD = (tD - s_px) / norm(tD - s_px);
nT = perp_toward_other(s_px, tT, tD);
nD = perp_toward_other(s_px, tD, tT);
if istable(RTAB384)
    widthEx = double(RTAB384{stimID_example, 7});
else
    widthEx = double(RTAB384(stimID_example, 7));
end

% Quartet map: [1 2 5 6], [3 4 7 8] per block of 8.
nBlocks = 384 / 8;
nQuartets = nBlocks * 2;
stimToQuartet = zeros(1, 384);
quartetMembers = zeros(nQuartets, 4);
q = 0;
for bIdx = 0:(nBlocks-1)
    base = 8*bIdx;
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    stimToQuartet(base + [1 2 5 6]) = q;
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
    stimToQuartet(base + [3 4 7 8]) = q;
end

selectedQuartets = unique(stimToQuartet(stimList));
selectedQuartets = selectedQuartets(selectedQuartets > 0);
if isempty(selectedQuartets)
    error('compute_projected_delta_points_allbins:NoQuartetsSelected', ...
        'stimIdx selection produced no valid quartets.');
end

% Geometry for projection uses canonical axes from example stimulus, but
% per-quartet site coordinates from the quartet stimulus tables.
requiredVars = ["assignment","overlap","along_GC","perp_signed_GC"];
qGeom = repmat(struct( ...
    'xT', nan(numel(siteRange),1), ...
    'yT', nan(numel(siteRange),1), ...
    'xD', nan(numel(siteRange),1), ...
    'yD', nan(numel(siteRange),1), ...
    'valid', false(numel(siteRange),1), ...
    'stimIdxSource', uint16(0), ...
    'stimMembers', uint16(zeros(1,4))), nQuartets, 1);

for qIdx = selectedQuartets(:)'
    stimsQ = quartetMembers(qIdx,:);
    Tq = [];
    stimGeom = 0;
    for stimNum = stimsQ
        if isfield(Tall_V1(stimNum), 'T') && istable(Tall_V1(stimNum).T)
            Ttmp = Tall_V1(stimNum).T;
            if height(Ttmp) >= max(siteRange) && ...
                    all(ismember(requiredVars, string(Ttmp.Properties.VariableNames)))
                Tq = Ttmp;
                stimGeom = stimNum;
                break;
            end
        end
    end
    if isempty(Tq)
        continue;
    end

    asg = string(Tq.assignment(siteRange));
    isBG = (asg == "background");
    isOV = Tq.overlap(siteRange) ~= 0;

    validGeom = ~isBG;
    if opt.excludeOverlap
        validGeom = validGeom & ~isOV;
    end
    validGeom = validGeom & sigMask(:);

    along = double(Tq.along_GC(siteRange)) * widthEx;
    perp = double(Tq.perp_signed_GC(siteRange)) * widthEx;

    xT = nan(numel(siteRange),1);
    yT = nan(numel(siteRange),1);
    xD = nan(numel(siteRange),1);
    yD = nan(numel(siteRange),1);
    if any(validGeom)
        xT(validGeom) = s_px(1) + along(validGeom).*uT(1) + perp(validGeom).*nT(1);
        yT(validGeom) = s_px(2) + along(validGeom).*uT(2) + perp(validGeom).*nT(2);
        xD(validGeom) = s_px(1) + along(validGeom).*uD(1) + perp(validGeom).*nD(1);
        yD(validGeom) = s_px(2) + along(validGeom).*uD(2) + perp(validGeom).*nD(2);
    end

    qGeom(qIdx).xT = xT;
    qGeom(qIdx).yT = yT;
    qGeom(qIdx).xD = xD;
    qGeom(qIdx).yD = yD;
    qGeom(qIdx).valid = validGeom;
    qGeom(qIdx).stimIdxSource = uint16(stimGeom);
    qGeom(qIdx).stimMembers = uint16(stimsQ);
end

% Prepare output struct.
bins = repmat(struct( ...
    'timeBin', [], ...
    'timeWindow', [], ...
    'x_px', single([]), ...
    'y_px', single([]), ...
    'delta', single([]), ...
    'siteIdx', uint16([]), ...
    'quartetIdx', uint16([]), ...
    'stimIdxSource', uint16([]), ...
    'stimMembersUsed', uint16([]), ...
    'stream', uint8([])), nFrames, 1);

if opt.verbose
    fprintf('Computing projected signed deltas for %d bins...\n', nFrames);
end

for tb = 1:nFrames
    deltaQ = compute_deltaQ_quartet(tb, R_resp, Tall_V1, SNRnorm, siteRange, quartetMembers, opt.excludeOverlap);

    X = [];
    Y = [];
    D = [];
    Site = [];
    Qid = [];
    StimSrc = [];
    StimMembers = zeros(0,4, 'uint16');
    Stream = zeros(0,1, 'uint8');

    for qIdx = selectedQuartets(:)'
        g = qGeom(qIdx);
        if g.stimIdxSource == 0
            continue;
        end

        dHere = deltaQ(:, qIdx);
        idxT = g.valid(:) & isfinite(dHere) & isfinite(g.xT) & isfinite(g.yT);
        if any(idxT)
            nAdd = nnz(idxT);
            X = [X; g.xT(idxT)]; %#ok<AGROW>
            Y = [Y; g.yT(idxT)]; %#ok<AGROW>
            D = [D; dHere(idxT)]; %#ok<AGROW>
            Site = [Site; siteRange(idxT(:)).']; %#ok<AGROW>
            Qid = [Qid; repmat(qIdx, nAdd, 1)]; %#ok<AGROW>
            StimSrc = [StimSrc; repmat(double(g.stimIdxSource), nAdd, 1)]; %#ok<AGROW>
            StimMembers = [StimMembers; repmat(g.stimMembers, nAdd, 1)]; %#ok<AGROW>
            Stream = [Stream; uint8(ones(nAdd,1))]; %#ok<AGROW>
        end

        idxD = g.valid(:) & isfinite(dHere) & isfinite(g.xD) & isfinite(g.yD);
        if any(idxD)
            nAdd = nnz(idxD);
            X = [X; g.xD(idxD)]; %#ok<AGROW>
            Y = [Y; g.yD(idxD)]; %#ok<AGROW>
            D = [D; dHere(idxD)]; %#ok<AGROW>
            Site = [Site; siteRange(idxD(:)).']; %#ok<AGROW>
            Qid = [Qid; repmat(qIdx, nAdd, 1)]; %#ok<AGROW>
            StimSrc = [StimSrc; repmat(double(g.stimIdxSource), nAdd, 1)]; %#ok<AGROW>
            StimMembers = [StimMembers; repmat(g.stimMembers, nAdd, 1)]; %#ok<AGROW>
            Stream = [Stream; uint8(2*ones(nAdd,1))]; %#ok<AGROW>
        end
    end

    bins(tb).timeBin = tb;
    bins(tb).timeWindow = double(R_resp.timeWindows(tb,:));
    bins(tb).x_px = single(X);
    bins(tb).y_px = single(Y);
    bins(tb).delta = single(D);
    bins(tb).siteIdx = uint16(Site);
    bins(tb).quartetIdx = uint16(Qid);
    bins(tb).stimIdxSource = uint16(StimSrc);
    bins(tb).stimMembersUsed = StimMembers;
    bins(tb).stream = Stream;

    if opt.verbose
        fprintf('  bin %2d/%2d [%4d %4d] ms -> %d points\n', ...
            tb, nFrames, R_resp.timeWindows(tb,1), R_resp.timeWindows(tb,2), numel(D));
    end
end

OUT = struct();
OUT.meta = struct();
OUT.meta.created = datestr(now, 30);
OUT.meta.stimID_example = stimID_example;
OUT.meta.siteRange = siteRange;
OUT.meta.excludeOverlap = opt.excludeOverlap;
OUT.meta.stimList = stimList;
OUT.meta.timeWindows = double(R_resp.timeWindows);
OUT.meta.nFrames = nFrames;
OUT.meta.nQuartets = nQuartets;
OUT.meta.selectedQuartets = selectedQuartets(:)';
OUT.meta.projectionReferenceStim = stimID_example;
OUT.meta.nSigSitesRequested = nnz(sigMask);
OUT.meta.qGeomStimIdxSource = arrayfun(@(qq) double(qq.stimIdxSource), qGeom(:))';
OUT.bins = bins;

saveFile = char(opt.saveFile);
if ~isempty(saveFile)
    save(saveFile, 'OUT', '-v7.3');
    if opt.verbose
        fprintf('Saved OUT to: %s\n', saveFile);
    end
end

end

function n = perp_toward_other(s_px, t_arm, t_other)
v = t_arm - s_px;
u = v / norm(v);
w = t_other - s_px;
w_perp = w - dot(w,u)*u;
n = w_perp / norm(w_perp);
end
