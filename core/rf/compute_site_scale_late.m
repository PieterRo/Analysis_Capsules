function siteScale = compute_site_scale_late(Rdata, Tall_V1, SNRn, siteRange, lateWindow, excludeOverlap, statMode)
% COMPUTE_SITE_SCALE_LATE
% Compute per-site scaling from late-phase windows using quartet-pooled deltas.
%
% statMode:
%   'mean_plus_sd' (default) -> mean(abs(deltaQ)) + std(abs(deltaQ))

if nargin < 6 || isempty(excludeOverlap), excludeOverlap = true; end
if nargin < 7 || isempty(statMode), statMode = 'mean_plus_sd'; end

tw = Rdata.timeWindows;
lateIdx = find(tw(:,1) >= lateWindow(1) & tw(:,2) <= lateWindow(2));
if isempty(lateIdx)
    error('No time windows fall within lateWindow [%g %g].', lateWindow(1), lateWindow(2));
end

% Quartets per block of 8 stimuli: [1 2 5 6] and [3 4 7 8]
nStim = 384;
nBlocks = nStim / 8;
nQuartets = nBlocks * 2;
quartetMembers = zeros(nQuartets, 4);
q = 0;
for bIdx = 0:(nBlocks-1)
    base = 8*bIdx;
    q = q + 1;
    quartetMembers(q,:) = base + [1 2 5 6];
    q = q + 1;
    quartetMembers(q,:) = base + [3 4 7 8];
end

nSites = numel(siteRange);
vals = cell(nSites,1);
for tb = lateIdx(:)'
    deltaQ = compute_deltaQ_quartet(tb, Rdata, Tall_V1, SNRn, siteRange, quartetMembers, excludeOverlap);
    d = abs(deltaQ);
    for s = 1:nSites
        v = d(s, :);
        v = v(isfinite(v));
        if ~isempty(v)
            vals{s} = [vals{s}; v(:)]; %#ok<AGROW>
        end
    end
end

siteScale = nan(max(siteRange),1);
for i = 1:nSites
    v = vals{i};
    if isempty(v)
        continue;
    end
    switch lower(statMode)
        case 'mean_plus_sd'
            siteScale(siteRange(i)) = mean(v) + std(v);
        otherwise
            error('Unknown statMode: %s', statMode);
    end
end
end
