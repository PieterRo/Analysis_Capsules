
function draw_rf_rects(idx, x, y, w, h, col, edgeAlpha, faceAlpha, lw)
% Draw axis-aligned RF rectangles centered at (x,y) with width w and height h
% Skips NaNs to preserve original indexing.

for ii = 1:numel(idx)
    i = idx(ii);

    if i < 1 || i > numel(x)
        continue
    end
    if isnan(x(i)) || isnan(y(i)) || isnan(w(i)) || isnan(h(i))
        continue
    end

    rectangle('Position', [x(i)-w(i)/2, y(i)-h(i)/2, w(i), h(i)], ...
        'EdgeColor', [col edgeAlpha], ...
        'FaceColor', [col faceAlpha], ...
        'LineWidth', lw, ...
        'HandleVisibility','off');   % prevent rectangles from flooding legend
end
end
