function a = overlapArea(p, q)
% p,q are [x y w h] in normalized coordinates
x1 = max(p(1), q(1));
y1 = max(p(2), q(2));
x2 = min(p(1)+p(3), q(1)+q(3));
y2 = min(p(2)+p(4), q(2)+q(4));
a = max(0, x2-x1) * max(0, y2-y1);
end
