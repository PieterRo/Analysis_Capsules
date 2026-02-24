
% ---- Assemble figures 1-4 into one PDF page ----
outPdf = fullfile(pwd,SaveFilePdf);

% Layout figure: 4 rows × 3 columns = 12 slots (we use 10)
F  = figure('Color','w','Units','normalized','Position',[0.02 0.05 0.96 0.88]);
TL = tiledlayout(F, 4, 3, 'TileSpacing','compact', 'Padding','compact');
set(F,'PaperPositionMode','auto');

tile = 1;

if p_val_local
    ptxt = sprintf('p_{local} = true, p_{th} = %.3f', p_th);
else
    ptxt = sprintf('p_{local} = false, p_{th,arr} = [%g %g %g]', p_th_arr);
end

headerText = sprintf([ ...
    'SNR_{th} = %.2f | %s | MinSigma = %d | MaxSigma = %d'], ...
    SNR_th, ptxt, MinSigma, MaxSigma);


% Helper: copy *all* axes from a figure into successive tiles
for srcFigNum = 1:4
    figSrc = figure(srcFigNum);
    drawnow;

    % Find axes that are not legends/colorbars
    axList = findall(figSrc,'Type','axes');
    axList = axList(~strcmp(get(axList,'Tag'),'legend'));
    axList = axList(~strcmp(get(axList,'Tag'),'Colorbar'));

    % Keep only axes that actually contain plotted objects
    axList = axList(arrayfun(@(h) ~isempty(allchild(h)), axList));

    % If the source figure is a 3-panel figure, we want a stable order.
    % Sort by vertical position (top->bottom), then horizontal (left->right).
    pos = zeros(numel(axList),4);
    for k = 1:numel(axList)
        pos(k,:) = get(axList(k),'Position');
    end
    % sort rows by y descending, then x ascending
    [~,ord] = sortrows([-pos(:,2), pos(:,1)]);
    axList = axList(ord);

    % Copy each axes into a new tile
    for k = 1:numel(axList)
        if tile > 12
            warning('Not enough tiles. Increase tiledlayout size.');
            break
        end

        axDst = nexttile(TL, tile);
        axDst.Units = 'normalized';
        dstPos = axDst.Position;   % [x y w h] in destination figure

        cla(axDst);  % make sure it's empty

        axSrc = axList(k);

        % Copy plotted objects into destination axes
        copyobj(allchild(axSrc), axDst);
        
        % --- thin all line objects in this panel ---
        lines = findall(axDst, 'Type','line');
        for h = lines'
            if isprop(h,'LineWidth')
                h.LineWidth = 0.75;   % try 0.5–1.0 depending on taste
            end
        end

        % Copy basic formatting
        set(axDst, ...
            'XLim', get(axSrc,'XLim'), ...
            'YLim', get(axSrc,'YLim'), ...
            'YDir', get(axSrc,'YDir'), ...
            'XScale', get(axSrc,'XScale'), ...
            'YScale', get(axSrc,'YScale'), ...
            'CLim', get(axSrc,'CLim'));

        % Labels and title (string only)
        xlabel(axDst, get(get(axSrc,'XLabel'),'String'));
        ylabel(axDst, get(get(axSrc,'YLabel'),'String'));
        title(axDst,  get(get(axSrc,'Title' ),'String'));

         % --- shrink all text to avoid overlap ---
        set(axDst, 'FontSize', 8);   % tick labels

        t = get(axDst,'Title');
        set(t, 'FontSize', 8, 'FontWeight','normal');

        xl = get(axDst,'XLabel');
        yl = get(axDst,'YLabel');
        set(xl, 'FontSize', 8);
        set(yl, 'FontSize', 8);

        % --- copy legend (smaller; last panel gets upper-left) ---
        legSrc = legend(axSrc);
        if ~isempty(legSrc) && isvalid(legSrc) && strcmp(legSrc.Visible,'on')

            % Get legend strings
            legStr = legSrc.String;
            if isempty(legStr)
                h = findobj(axDst, '-property','DisplayName');
                legStr = get(h,'DisplayName');
                if ischar(legStr), legStr = {legStr}; end
                legStr = legStr(~cellfun(@isempty, legStr));
            end

            if ~isempty(legStr)
                % Default: upper-right. Exception: last plot (tile 10) -> upper-left
                if tile == 10   % <-- use your tile counter / destination tile index here
                    loc = 'northwest';
                else
                    loc = 'northeast';
                end

                legDst = legend(axDst, legStr, 'Location', loc, 'Box','off');

                % Smaller + more compact legend
                legDst.FontSize      = 6;
                legDst.ItemTokenSize = [6 6];
                legDst.Interpreter   = legSrc.Interpreter;
                legDst.AutoUpdate    = 'off';
            end
        end
        
        % --- copy annotation boxes (per subplot) ---
        axSrc.Units = 'normalized';
        srcPos = axSrc.Position;        % [x y w h] in SOURCE figure (normalized)

        axDst.Units = 'normalized';
        dstPos = axDst.Position;        % [x y w h] in DEST figure (normalized)

        % Find annotation shapes in the source figure.
        % Depending on how the annotations were created, they may show up via
        % 'Type' (textboxshape/rectangleshape) OR via class filtering.
        annTB1   = findall(figSrc,'Type','textboxshape');
        annRect1 = findall(figSrc,'Type','rectangleshape');
        annTB2   = findall(figSrc,'-isa','matlab.graphics.shape.TextBox');
        annRect2 = findall(figSrc,'-isa','matlab.graphics.shape.Rectangle');
        ann = unique([annTB1; annRect1; annTB2; annRect2]);

        % Some annotations sit just outside the axes box; allow a small margin.
        pad = 0.03;  % fraction of axes size
        x0 = srcPos(1) - pad*srcPos(3);
        y0 = srcPos(2) - pad*srcPos(4);
        x1 = srcPos(1) + srcPos(3) + pad*srcPos(3);
        y1 = srcPos(2) + srcPos(4) + pad*srcPos(4);

        for ia = 1:numel(ann)
            a = ann(ia);
            aNew = [];  

            % Annotation position in normalized SOURCE-FIGURE coordinates
            aPos = hgconvertunits(figSrc, a.Position, a.Units, 'normalized', figSrc);  % [x y w h]

            % Assign annotation to THIS subplot if its CENTER lies inside the source axes box
            cx = aPos(1) + 0.5*aPos(3);
            cy = aPos(2) + 0.5*aPos(4);
            if cx < x0 || cx > x1 || cy < y0 || cy > y1
                continue
            end

            % Map annotation position from source-axes-relative -> destination tile
            relX = (aPos(1) - srcPos(1)) / srcPos(3);
            relY = (aPos(2) - srcPos(2)) / srcPos(4);
            relW = aPos(3) / srcPos(3);
            relH = aPos(4) / srcPos(4);

            newPos = [ ...
                dstPos(1) + relX*dstPos(3), ...
                dstPos(2) + relY*dstPos(4), ...
                relW*dstPos(3), ...
                relH*dstPos(4) ];

            % Recreate as figure annotations (copyobj is not allowed for TextBox objects)
            if strcmpi(a.Type,'textboxshape')
                aNew = annotation(F, 'textbox', newPos, ...
                    'Units','normalized', ...
                    'String', a.String, ...
                    'FitBoxToText', 'off', ...
                    'LineWidth', 0.5);

                % Best-effort style transfer
                if isprop(a,'EdgeColor'), aNew.EdgeColor = a.EdgeColor; end
                if isprop(a,'FaceColor')
                    try, aNew.BackgroundColor = a.FaceColor; end %#ok<TRYNC>
                end
                if isprop(a,'Interpreter'), aNew.Interpreter = a.Interpreter; end
                if isprop(a,'FontWeight'),  aNew.FontWeight = a.FontWeight; end
                if isprop(a,'FontName'),    aNew.FontName   = a.FontName; end
                if isprop(a,'HorizontalAlignment'), aNew.HorizontalAlignment = a.HorizontalAlignment; end
                if isprop(a,'VerticalAlignment'),   aNew.VerticalAlignment   = a.VerticalAlignment; end

                % Compact styling for dense multi-panel PDF
                if isprop(aNew,'FontSize'), aNew.FontSize = 5; end
                if isprop(aNew,'Margin'),   aNew.Margin   = 1; end

            elseif strcmpi(a.Type,'rectangleshape')
                aNew = annotation(F, 'rectangle', newPos, 'Units','normalized');
                if isprop(a,'EdgeColor'), aNew.Color = a.EdgeColor; end
                aNew.LineWidth = 0.5;

                if isprop(a,'FaceColor')
                    if isprop(aNew,'FaceColor'), aNew.FaceColor = a.FaceColor; end
                    if isprop(aNew,'BackgroundColor'), aNew.BackgroundColor = a.FaceColor; end
                end
            end
        end

        % Grids/box
        if strcmp(get(axSrc,'XGrid'),'on'), grid(axDst,'on'); else, grid(axDst,'off'); end
        box(axDst, get(axSrc,'Box'));

        tile = tile + 1;
    end
end

% Turn off any remaining unused tiles
for t = tile:12
    ax = nexttile(TL, t);
    axis(ax,'off');
end
fig = gcf;

annotation(fig, 'textbox', ...
    [0 0.99 1 0.04], ...   % [x y width height] in normalized figure units
    'String', headerText, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 9, ...
    'Interpreter', 'tex', ...
    'FontWeight', 'normal');


t33txt = sprintf('t_{33} = [%d  %d  %d] ms', round(t33_all));
annotation(gcf, 'textbox', ...
    [0 0.01 1 0.04], ...    % bottom of page
    'String', t33txt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 9, ...
    'Interpreter', 'tex', ...
    'Color', [0.3 0.3 0.3]);
drawnow;

% Export (prefer print for R2020b reliability; omit -painters if images/transparency exist)
print(F, outPdf, '-dpdf', '-bestfit');
fprintf('Saved: %s\n', outPdf);


