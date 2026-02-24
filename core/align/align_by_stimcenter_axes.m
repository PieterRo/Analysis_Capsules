function best = align_by_stimcenter_axes(Fmov, Ffix, imSize)
% Choose mirror + (optional) scale + rotation + translation from extracted features.
% Uses stimCenter (axis intersection) as anchor and yellow axis as rotation reference.
% Picks mirror mode by maximizing IoU of union masks after warping.

H = imSize(1); W = imSize(2);

modes = ["none","fliplr","flipud"];

BWfix = Ffix.yellowMask | Ffix.purpleMask;
Rfix  = imref2d(size(BWfix));

best.iou = -inf;

for mm = modes
    % mirror masks & features
    [BWm, Fm] = mirrorAll(Fmov, mm, H, W);

    % --- optional uniform scale estimate (often ~1) ---
    % ratio based on mean major axis length of the two capsules
    s = mean([Ffix.yellow.MajorAxisLength, Ffix.purple.MajorAxisLength]) / ...
        mean([Fm.yellow.MajorAxisLength,  Fm.purple.MajorAxisLength]);

    % rotation: align yellow directed axis vectors
    vMov = directedAxis(Fm.yellow.Centroid, Fm.yellow.Orientation, Fm.purple.Centroid);
    vFix = directedAxis(Ffix.yellow.Centroid, Ffix.yellow.Orientation, Ffix.purple.Centroid);
    rotDeg = atan2d(vFix(2),vFix(1)) - atan2d(vMov(2),vMov(1));

    % build similarity transform about stimCenter
    Cmov = Fm.stimCenter_xy;
    Cfix = Ffix.stimCenter_xy;

    T = makeSimilarityAboutCenters(s, rotDeg, Cmov, Cfix);

    BWw = imwarp(BWm, T, 'OutputView', Rfix, 'Interp','nearest');

    iou = nnz(BWw & BWfix) / nnz(BWw | BWfix);

    if iou > best.iou
        best.iou = iou;
        best.mirror = mm;
        best.scale = s;
        best.rotDeg = rotDeg;
        best.tform = T;
    end
end
end

% ---------- helpers ----------
function [BWm, Fm] = mirrorAll(F, mode, H, W)
Fm = F;
BW = F.yellowMask | F.purpleMask;

switch mode
    case "none"
        BWm = BW;
    case "fliplr"
        BWm = fliplr(BW);
        Fm.yellow.Centroid(1) = (W+1) - Fm.yellow.Centroid(1);
        Fm.purple.Centroid(1) = (W+1) - Fm.purple.Centroid(1);
        Fm.stimCenter_xy(1)   = (W+1) - Fm.stimCenter_xy(1);
        Fm.yellow.Orientation = wrapTo180(180 - Fm.yellow.Orientation);
        Fm.purple.Orientation = wrapTo180(180 - Fm.purple.Orientation);
    case "flipud"
        BWm = flipud(BW);
        Fm.yellow.Centroid(2) = (H+1) - Fm.yellow.Centroid(2);
        Fm.purple.Centroid(2) = (H+1) - Fm.purple.Centroid(2);
        Fm.stimCenter_xy(2)   = (H+1) - Fm.stimCenter_xy(2);
        Fm.yellow.Orientation = wrapTo180(-Fm.yellow.Orientation);
        Fm.purple.Orientation = wrapTo180(-Fm.purple.Orientation);
end

Fm.yellowMask = BWm & (BWm); % not used further here
Fm.purpleMask = BWm & (BWm);
end

function v = directedAxis(centroid, orientationDeg, otherCentroid)
a = deg2rad(orientationDeg);
v1 = [cos(a), -sin(a)];  % image coords
v2 = -v1;
toOther = otherCentroid - centroid;
if dot(v1,toOther) >= dot(v2,toOther)
    v = v1 / norm(v1);
else
    v = v2 / norm(v2);
end
end

function T = makeSimilarityAboutCenters(s, rotDeg, Cmov, Cfix)
% Map moving -> fixed:
%   p' = (p - Cmov) * (sR) + Cfix
% using MATLAB affine2d row-vector convention: [x y 1]*T

R = [cosd(rotDeg) -sind(rotDeg);
     sind(rotDeg)  cosd(rotDeg)];
A = s * R;  % 2x2

% Row-vector affine2d expects:
% x' = x*T(1,1) + y*T(2,1) + T(3,1)
% y' = x*T(1,2) + y*T(2,2) + T(3,2)
Tlin = [ A(1,1)  A(1,2)  0;
         A(2,1)  A(2,2)  0;
         0       0       1];

Tto0  = [1 0 0;
         0 1 0;
        -Cmov(1) -Cmov(2) 1];

Tback = [1 0 0;
         0 1 0;
          Cfix(1)  Cfix(2) 1];

% IMPORTANT: for row-vectors, composition is left-to-right:
% p * Tto0 * Tlin * Tback
T = affine2d(Tto0 * Tlin * Tback);
end
