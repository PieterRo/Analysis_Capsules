function make_activity_movie_wrapper_safe(outFile, ...
    Tall, ALLCOORDS, RTAB384, exampleStimNum, R, SNR, varargin)

if ~isstruct(R) || ~isfield(R, 'timeWindows') || ~isfield(R, 'meanAct')
    error(['R must be a response struct with fields timeWindows and meanAct. ' ...
           'Use the struct loaded from Resp_capsules_* (70 windows).']);
end

nFrames = size(R.timeWindows,1);
if nFrames == 3
    error(['R.timeWindows has 3 bins (SNR-style struct). ' ...
           'For movies use the response struct from Resp_capsules_* (70 bins).']);
end
if size(R.meanAct,3) ~= nFrames
    error('Mismatch: size(R.meanAct,3)=%d but size(R.timeWindows,1)=%d.', ...
        size(R.meanAct,3), nFrames);
end

fprintf('\n============================================\n');
fprintf('Starting movie rendering (%d frames)\n', nFrames);
fprintf('Output file: %s\n', outFile);
fprintf('============================================\n');

v = VideoWriter(outFile, 'MPEG-4');
v.FrameRate = 10;
open(v);

tStart = tic;

for tb = 1:nFrames

    % ---- Progress reporting ----
    tElapsed = toc(tStart);
    pct = 100 * tb / nFrames;

    fprintf('Frame %2d / %2d  |  %5.1f%%  |  %4d–%4d ms  |  elapsed: %.1fs\n', ...
        tb, nFrames, pct, ...
        R.timeWindows(tb,1), R.timeWindows(tb,2), ...
        tElapsed);

    % ---- Call your existing plotting function ----
    h = plot_projected_activity_on_example_stim( ...
        Tall, ALLCOORDS, RTAB384, exampleStimNum, ...
        R, SNR, ...
        'TimeBin', tb, ...
        varargin{:});

    
    
 % --- If pre time-zero: show only grey background + activity ---
    isPreZero = (R.timeWindows(tb,2) <= 0);
    
    if isPreZero
        ax = gca;
    
        % 1) Replace underlying stimulus image by flat grey (keep size from current image)
        hImg = findobj(ax, 'Type', 'image');
        if ~isempty(hImg)
            % Use the first image object (imshow creates one)
            hImg = hImg(1);
            C = hImg.CData;              % existing image
            [hh, ww, ~] = size(C);
            greyVal = uint8(128);        % mid grey
            hImg.CData = repmat(greyVal, [hh, ww, 3]);
        end
    
        % 2) Hide any curves/outlines if they exist as line/patch objects
        set(findobj(ax, 'Type', 'line'),  'Visible', 'off');
        set(findobj(ax, 'Type', 'patch'), 'Visible', 'off');
    end   
    
    
    
    drawnow;

    % --- Time counter overlay (upper-left) ---
    ax = gca;  % assumes your plotting function leaves the right axes current
    tw = R.timeWindows(tb,:);

    % choose what to display:
    tMid = mean(tw); % midpoint in ms
    label = sprintf('t = %.0f ms', tMid);                 % or: sprintf('%d–%d ms', tw(1), tw(2));

    % place in normalized axes coordinates (0..1)
    text(ax, 0.02, 0.98, label, ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize', 16, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'BackgroundColor',[0 0 0 0.35], ...   % semi-transparent black (newer MATLAB)
        'Margin', 4);

    % ---- Capture frame safely ----
    if isstruct(h) && isfield(h,'fig')
        frame = getframe(h.fig);
        close(h.fig);   % prevents graphics accumulation
    else
        frame = getframe(gcf);
        close(gcf);
    end

    writeVideo(v, frame);

    drawnow;  % ensures console updates
end

close(v);

fprintf('============================================\n');
fprintf('Movie finished in %.1f seconds\n', toc(tStart));
fprintf('Saved to: %s\n', outFile);
fprintf('============================================\n\n');

end
