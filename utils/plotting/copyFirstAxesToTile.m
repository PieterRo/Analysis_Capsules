function copyFirstAxesToTile(srcFigNum, TL, tileIdx)
    nexttile(TL, tileIdx);

    srcFig = figure(srcFigNum);

    % Pick the first non-legend/non-colorbar axes (adjust if needed)
    ax = findobj(srcFig, 'Type','axes', '-not','Tag','legend', '-not','Tag','Colorbar');
    ax = ax(end);  % often the "main" axes is last; if wrong, try ax(1)

    % Copy the axes into the current tile
    newAx = copyobj(ax, gcf);

    % Make it fill the tile nicely
    set(newAx, 'Units','normalized', 'Position', get(gca,'Position'));
    delete(gca);              % remove the empty tile axes
    set(newAx,'Box','off');
end