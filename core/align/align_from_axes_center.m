function out = align_from_axes_center(Fmov, Ffix, imSize, tolDeg)
% out = align_from_axes_center(Fmov, Ffix, imSize, tolDeg)
%
% Implements your logic:
%  1) compute rotation that aligns yellow axis
%  2) check if purple axis aligns under same rotation (within tolDeg)
%     -> if yes: pure rotation (about center) + translation
%  3) else: apply a vertical-line mirror through Cmov, mirror the orientations,
%     compute rotation from yellow again, check purple again
%     -> if yes: mirror + rotation + translation
%  4) else: flag failure
%
% Requires: F.stimCenter_xy, F.yellow.Orientation, F.purple.Orientation
% Optionally: masks for IoU diagnostics.

if nargin < 4, tolDeg = 3; end
H = imSize(1); W = imSize(2); %#ok<NASGU>

C1 = Fmov.stimCenter_xy(:)';   % [x y]
C6 = Ffix.stimCenter_xy(:)';

y1 = wrapTo180(Fmov.yellow.Orientation);
p1 = wrapTo180(Fmov.purple.Orientation);
y6 = wrapTo180(Ffix.yellow.Orientation);
p6 = wrapTo180(Ffix.purple.Orientation);

% ---------- attempt 1: no mirror ----------
rot1 = axisAngleDeltaDeg(y1, y6);           % rotation to align yellow axis
p1r  = wrapTo180(p1 + rot1);                % purple after applying same rotation

ok1 = axisAngleDistDeg(p1r, p6) <= tolDeg;  % purple matches too?

if ok1
    out.mirror = false;
    out.rotDeg = rot1;
    out.tform  = tform_centered(true, false, rot1, C1, C6); % rotate about C1 then translate to C6
    out.ok = true;
    out.reason = "rotation+translation";
    out.diag.purpleErrDeg = axisAngleDistDeg(p1r, p6);
    out.diag.yellowErrDeg = axisAngleDistDeg(wrapTo180(y1+rot1), y6);
    return;
end

% ---------- attempt 2: mirror (vertical through center) then rotate ----------
% Mirror vertical: orientation transforms as theta -> 180 - theta
y1m = wrapTo180(180 - y1);
p1m = wrapTo180(180 - p1);

rot2 = axisAngleDeltaDeg(y1m, y6);
p1mr = wrapTo180(p1m + rot2);

ok2 = axisAngleDistDeg(p1mr, p6) <= tolDeg;

out.mirror = true;
out.rotDeg = rot2;
out.tform  = tform_centered(true, true, rot2, C1, C6);  % mirror about vertical through C1, then rotate, then translate
out.ok = ok2;
out.reason = "mirror(vertical@center)+rotation+translation";
out.diag.purpleErrDeg = axisAngleDistDeg(p1mr, p6);
out.diag.yellowErrDeg = axisAngleDistDeg(wrapTo180(y1m+rot2), y6);

if ~ok2
    out.reason = "FAILED: purple axis does not align even after mirroring";
end
end

% ================= helper functions =================

function d = axisAngleDistDeg(a, b)
% distance between undirected axes (mod 180)
d0 = abs(wrapTo180(a - b));
d  = min(d0, abs(d0 - 180));   % because axis has 180° symmetry
end

function rot = axisAngleDeltaDeg(from, to)
% rotation to align axis 'from' to axis 'to' (choose minimal magnitude)
% i.e. find rot such that from+rot matches to modulo 180
cands = [wrapTo180(to - from), wrapTo180(to - from + 180), wrapTo180(to - from - 180)];
[~,k] = min(abs(cands));
rot = cands(k);
end

function T = tform_centered(doTranslate, doMirrorVert, rotDeg, Cmov, Cfix)
% Build affine2d in row-vector convention: [x y 1] * T
% Steps in centered coords:
%   translate by -Cmov
%   (optional) mirror about vertical axis: (dx,dy)->(-dx,dy)
%   rotate by rotDeg
%   (optional) translate by +Cfix (if doTranslate)

T_to0 = [1 0 0;
         0 1 0;
        -Cmov(1) -Cmov(2) 1];

if doMirrorVert
    M = [-1 0 0;
          0 1 0;
          0 0 1];
else
    M = eye(3);
end

R = [cosd(rotDeg) -sind(rotDeg) 0;
     sind(rotDeg)  cosd(rotDeg) 0;
     0            0            1];

if doTranslate
    T_back = [1 0 0;
              0 1 0;
              Cfix(1) Cfix(2) 1];
else
    T_back = eye(3);
end

T = affine2d(T_to0 * M * R * T_back);
end
