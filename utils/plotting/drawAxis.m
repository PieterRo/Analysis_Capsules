function drawAxis(c, oriDeg, L, col)
% Draw regionprops Orientation on an image (y axis points DOWN)
a = deg2rad(oriDeg);
dx = cos(a)*L;
dy = -sin(a)*L;    % <-- IMPORTANT for image coordinates
plot([c(1)-dx c(1)+dx], [c(2)-dy c(2)+dy], '-', ...
     'Color', col, 'LineWidth', 2);
end
