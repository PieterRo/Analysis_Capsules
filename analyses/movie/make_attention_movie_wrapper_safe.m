function make_attention_movie_wrapper_safe(outFile, ...
    Tall_V1, ALLCOORDS, RTAB384, stimID_example, R_resp, SNRnorm, OUT3, varargin)
% MAKE_ATTENTION_MOVIE_WRAPPER_SAFE
% Render a time-resolved movie of attention effect (T-D) in 10 ms bins.
%
% Hard significance rule:
%   site is included iff OUT3.pValueTD < pThresh (default 0.05)
%
% Inputs:
%   outFile         : output .mp4 path
%   Tall_V1         : stimulus/site assignment table struct (1x384)
%   ALLCOORDS       : stimulus geometry struct
%   RTAB384         : stimulus metadata table
%   stimID_example  : stimulus ID for rendering frame (1..384)
%   R_resp          : response struct with many time windows (Resp_capsules_*)
%   SNRnorm         : normalization struct (muSpont/mu* fields)
%   OUT3            : baseline OUT struct providing pValueTD mask
%
% Optional name/value pairs:
%   'pThresh'         (0.05)
%   'frameRate'       (10)
%   'excludeOverlap'  (true)
%   'siteRange'       (1:512)
%   'alpha'           (0.50)
%   'markerSize'      (5)
%   'robustPct'       (95)
%   'bgColor'         ([0.5 0.5 0.5])
%   'cLow'            ([0.55 0.55 0.55])
%   'cHigh'           ([0.85 0.05 0.05])
%   'stimIdx'         ([])
%   'hideStimPreZero' (true)
%   'verbose'         (true)
%   'useWaitbar'      (true)

% ---- Defaults ----
optsMovie = struct();
optsMovie.pThresh = 0.05;
optsMovie.frameRate = 10;
optsMovie.excludeOverlap = true;
optsMovie.siteRange = 1:512;
optsMovie.alpha = 0.50;
optsMovie.markerSize = 5;
optsMovie.robustPct = 95;
optsMovie.bgColor = [0.5 0.5 0.5];
optsMovie.cLow = [0.55 0.55 0.55];
optsMovie.cHigh = [0.85 0.05 0.05];
optsMovie.stimIdx = [];
optsMovie.hideStimPreZero = true;
optsMovie.verbose = true;
optsMovie.useWaitbar = true;

% ---- Parse name/value ----
if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k+1};
    if ~ischar(name)
        error('Optional argument names must be char.');
    end
    if ~isfield(optsMovie, name)
        error('Unknown option: %s', name);
    end
    optsMovie.(name) = value;
end

% ---- Validate inputs ----
if ~isstruct(R_resp) || ~isfield(R_resp, 'timeWindows') || ~isfield(R_resp, 'meanAct')
    error(['R_resp must be a response struct with fields timeWindows and meanAct. ' ...
           'Use the struct loaded from Resp_capsules_* (many windows).']);
end
nFrames = size(R_resp.timeWindows,1);
if nFrames <= 3
    error(['R_resp.timeWindows has %d bins. ' ...
           'For this movie use the response struct from Resp_capsules_* (e.g., ~70 bins).'], nFrames);
end
if size(R_resp.meanAct,3) ~= nFrames
    error('Mismatch: size(R_resp.meanAct,3)=%d but size(R_resp.timeWindows,1)=%d.', ...
        size(R_resp.meanAct,3), nFrames);
end

if ~isstruct(OUT3) || ~isfield(OUT3,'pValueTD')
    error('OUT3 must contain pValueTD for fixed significance masking.');
end
if numel(OUT3.pValueTD) < max(optsMovie.siteRange)
    error('OUT3.pValueTD is smaller than max(siteRange).');
end

siteRange = optsMovie.siteRange(:)';
sigMask = isfinite(OUT3.pValueTD(siteRange)) & (OUT3.pValueTD(siteRange) < optsMovie.pThresh);
nSigSites = nnz(sigMask);

if optsMovie.verbose
    fprintf('\n============================================\n');
    fprintf('Starting attention-effect movie rendering (%d frames)\n', nFrames);
    fprintf('Output file: %s\n', outFile);
    fprintf('Example stimulus: %d\n', stimID_example);
    fprintf('Significant sites (fixed mask): %d / %d\n', nSigSites, numel(siteRange));
    fprintf('============================================\n');
end

wb = [];
if optsMovie.useWaitbar
    try
        wb = waitbar(0, 'Precomputing attention bins...');
    catch
        wb = [];
    end
end

% ---- Precompute OUT per frame and shared color scale ----
OUT_byFrame = cell(nFrames,1);
deltaAbs = [];
tp = tic;
for tb = 1:nFrames
    optsBin = struct('timeIdx',tb, ...
                     'excludeOverlap',optsMovie.excludeOverlap, ...
                     'verbose',false, ...
                     'v1Sites',siteRange);
    if ~isempty(optsMovie.stimIdx)
        optsBin.stimIdx = optsMovie.stimIdx;
    end

    OUTtb = attention_modulation_V1_3bin(R_resp, Tall_V1, SNRnorm, optsBin);
    OUTtb.pValueTD = OUT3.pValueTD;  % hard significance rule
    OUT_byFrame{tb} = OUTtb;

    d = OUTtb.muT - OUTtb.muD;
    d = d(:);
    d = d(sigMask);
    d = abs(d(isfinite(d)));
    if ~isempty(d)
        deltaAbs = [deltaAbs; d]; %#ok<AGROW>
    end

    if optsMovie.verbose
        tElapsed = toc(tp);
        pct = 100 * tb / nFrames;
        tPer = tElapsed / tb;
        eta = tPer * (nFrames - tb);
        fprintf('[Precompute] %2d/%2d | %5.1f%% | elapsed %.1fs | ETA %.1fs\n', ...
            tb, nFrames, pct, tElapsed, eta);
    end
    if ~isempty(wb)
        waitbar(0.45 * tb / nFrames, wb, sprintf('Precomputing bins... %d/%d', tb, nFrames));
    end
end

if isempty(deltaAbs)
    cShared = 1;
else
    cShared = prctile(deltaAbs, optsMovie.robustPct);
    if ~isfinite(cShared) || cShared <= 0
        cShared = max(deltaAbs);
    end
    if ~isfinite(cShared) || cShared <= 0
        cShared = 1;
    end
end

if optsMovie.verbose
    fprintf('Shared color scale cMax: %.6g\n', cShared);
end

% ---- Video setup ----
v = VideoWriter(outFile, 'MPEG-4');
v.FrameRate = optsMovie.frameRate;
open(v);

tStart = tic;

for tb = 1:nFrames
    if optsMovie.verbose
        tElapsed = toc(tStart);
        pct = 100 * tb / nFrames;
        tPer = tElapsed / max(tb-1,1);
        if tb == 1
            eta = NaN;
        else
            eta = tPer * (nFrames - tb + 1);
        end
        fprintf('Frame %2d / %2d  |  %5.1f%%  |  %4d-%4d ms  |  elapsed: %.1fs\n', ...
            tb, nFrames, pct, ...
            R_resp.timeWindows(tb,1), R_resp.timeWindows(tb,2), tElapsed);
        if isfinite(eta)
            fprintf('                 ETA %.1fs\n', eta);
        end
    end

    optsPlot = struct();
    optsPlot.RTAB384 = RTAB384;
    optsPlot.pThresh = optsMovie.pThresh;
    optsPlot.siteRange = siteRange;
    optsPlot.markerSize = optsMovie.markerSize;
    optsPlot.alpha = optsMovie.alpha;
    optsPlot.robustPct = optsMovie.robustPct;
    optsPlot.stimIdx = optsMovie.stimIdx;
    optsPlot.cMaxFixed = cShared;
    optsPlot.bgColor = optsMovie.bgColor;
    optsPlot.cLow = optsMovie.cLow;
    optsPlot.cHigh = optsMovie.cHigh;
    optsPlot.projectionMode = 'quartet_pooled';
    optsPlot.Rdata = R_resp;
    optsPlot.SNRnorm = SNRnorm;
    optsPlot.timeIdx = tb;

    h = plot_projected_attentiondiff_on_example_stim( ...
        stimID_example, OUT_byFrame{tb}, Tall_V1, ALLCOORDS, optsPlot);

    % Optional: hide stimulus outlines in pre-zero frames
    isPreZero = (R_resp.timeWindows(tb,2) <= 0);
    if optsMovie.hideStimPreZero && isPreZero
        ax = h.ax;
        hImg = findobj(ax, 'Type', 'image');
        if ~isempty(hImg)
            hImg = hImg(1);
            C = hImg.CData;
            [hh, ww, ~] = size(C);
            g = uint8(round(255 * optsMovie.bgColor(1)));
            hImg.CData = repmat(g, [hh, ww, 3]);
        end
        set(findobj(ax, 'Type', 'line'),  'Visible', 'off');
        set(findobj(ax, 'Type', 'patch'), 'Visible', 'off');
    end

    figure(h.fig);
    set(h.fig, 'Name', 'Attention TD Movie (10 ms bins)', 'NumberTitle', 'off');
    title(h.ax, sprintf('Attention effect (T-D) | frame %d/%d | %d-%d ms | stim %d', ...
        tb, nFrames, R_resp.timeWindows(tb,1), R_resp.timeWindows(tb,2), stimID_example), ...
        'Color','w','FontSize',14);

    text(h.ax, 0.02, 0.98, sprintf('t = %.0f ms', mean(R_resp.timeWindows(tb,:))), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize', 14, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'BackgroundColor','k', ...
        'Margin', 4);

    drawnow;
    frame = getframe(h.fig);
    writeVideo(v, frame);
    close(h.fig);
    if ~isempty(wb)
        waitbar(0.45 + 0.55 * tb / nFrames, wb, sprintf('Rendering frames... %d/%d', tb, nFrames));
    end
end

close(v);
if ~isempty(wb)
    try
        close(wb);
    catch
    end
end

if optsMovie.verbose
    fprintf('============================================\n');
    fprintf('Attention movie finished in %.1f seconds\n', toc(tStart));
    fprintf('Saved to: %s\n', outFile);
    fprintf('============================================\n\n');
end
end
