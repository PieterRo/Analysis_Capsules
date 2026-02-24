function copyAxesContentsToTile(srcFigNum, TL, tileIdx)
    % Destination axes (tile)
    axDst = nexttile(TL, tileIdx);

    % Source axes: choose the "main" one (largest)
    figSrc = figure(srcFigNum);
    axList = findall(figSrc, 'Type','axes');
    axList = axList(~strcmp(get(axList,'Tag'),'legend'));   % drop legends
    axList = axList(~strcmp(get(axList,'Tag'),'Colorbar')); % drop colorbars

    if isempty(axList)
        warning('No axes found in figure %d.', srcFigNum);
        axis(axDst,'off');
        return
    end

    % pick largest axes by area (robust across scalar/vector handles)
    areas = zeros(numel(axList),1);
    for k = 1:numel(axList)
        p = get(axList(k),'Position');   % always numeric 1x4
        areas(k) = p(3)*p(4);
    end
    [~,i] = max(areas);
    axSrc = axList(i);

    % Copy plotted objects (children) into the destination axes
    copyobj(allchild(axSrc), axDst);

    % Copy key axes properties
    set(axDst, ...
        'XLim', get(axSrc,'XLim'), ...
        'YLim', get(axSrc,'YLim'), ...
        'YDir', get(axSrc,'YDir'), ...
        'XScale', get(axSrc,'XScale'), ...
        'YScale', get(axSrc,'YScale'), ...
        'CLim', get(axSrc,'CLim'), ...
        'Box',  get(axSrc,'Box'));

    % Labels + title
    xlabel(axDst, get(get(axSrc,'XLabel'),'String'));
    ylabel(axDst, get(get(axSrc,'YLabel'),'String'));
    title(axDst,  get(get(axSrc,'Title' ),'String'));

    % Grid
    grid(axDst, get(axSrc,'XGrid'));  % copies 'on'/'off' (good enough)

    % Make sure limits apply after children copy
    axis(axDst,'tight');
    set(axDst, 'XLim', get(axSrc,'XLim'), 'YLim', get(axSrc,'YLim'));
end