function best = estimate_rigid_from_axes(Fmov, Ffix, imSize)
% Estimate mirror + rotation + translation using yellow?purple vector

H = imSize(1); W = imSize(2);
modes = ["none","fliplr","flipud"];

best.err = inf;

for m = modes
    [cy, oy] = mirrorFeature(Fmov.yellow.Centroid, Fmov.yellow.Orientation, m, H, W);
    [cp, op] = mirrorFeature(Fmov.purple.Centroid, Fmov.purple.Orientation, m, H, W);

    vM = cp - cy;
    vF = Ffix.purple.Centroid - Ffix.yellow.Centroid;

    rot = atan2d(vF(2),vF(1)) - atan2d(vM(2),vM(1));
    R = [cosd(rot) -sind(rot); sind(rot) cosd(rot)];

    cy2 = (R*cy')';
    t = Ffix.yellow.Centroid - cy2;

    cp2 = (R*cp')' + t;
    cy2 = cy2 + t;

    err = norm(cp2 - Ffix.purple.Centroid) + ...
          norm(cy2 - Ffix.yellow.Centroid);

    if err < best.err
        best.err = err;
        best.mirror = m;
        best.rotDeg = rot;
        best.R = R;
        best.t = t;
    end
end
end

function [c2,o2] = mirrorFeature(c,o,mode,H,W)
c2 = c; o2 = o;
switch mode
    case "fliplr"
        c2(1) = (W+1)-c(1); o2 = 180-o;
    case "flipud"
        c2(2) = (H+1)-c(2); o2 = -o;
end
o2 = wrapTo180(o2);
end
