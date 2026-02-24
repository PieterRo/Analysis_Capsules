function out = align_centroids_about_center(F1, F6, varargin)
%ALIGN_CENTROIDS_ABOUT_CENTER
% Deterministic alignment using center + yellow/purple centroids.
%
% Steps:
%  1) Form vectors uY/uP from stimCenter -> centroids in moving & fixed.
%  2) Compute rotY and rotP. If they agree within rotTolDeg, take circular mean.
%  3) Else mirror moving vectors about vertical line through center (dx -> -dx),
%     recompute rotY/rotP, and accept if they agree.
%  4) Build affine2d that maps moving->fixed:
%       translate by -C1, optional mirror, rotate, translate by +C6
%
% Requires fields:
%   F.stimCenter_xy
%   F.yellow.Centroid
%   F.purple.Centroid
%
% Name/value options:
%   'rotTolDeg'  (default 8)   tolerance between rotY and rotP
%   'centTolPx'  (default 3)   tolerance for centroid check after transform

p = inputParser;
p.addParameter('rotTolDeg', 8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('centTolPx', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
rotTolDeg = p.Results.rotTolDeg;
centTolPx = p.Results.centTolPx;

C1 = F1.stimCenter_xy(:)';   C6 = F6.stimCenter_xy(:)';
Y1 = F1.yellow.Centroid(:)'; Y6 = F6.yellow.Centroid(:)';
P1 = F1.purple.Centroid(:)'; P6 = F6.purple.Centroid(:)';

% Centered vectors
uY1 = Y1 - C1;   uP1 = P1 - C1;
uY6 = Y6 - C6;   uP6 = P6 - C6;

% ---- attempt A: no mirror ----
[okA, rotDegA, diagA] = solveRot(uY1,uP1,uY6,uP6, rotTolDeg);
if okA
    doMirror = false;
    rotDeg   = rotDegA;
    diagRot  = diagA;
else
    % ---- attempt B: mirror about vertical line through center (dx -> -dx) ----
    uY1m = [-uY1(1), uY1(2)];
    uP1m = [-uP1(1), uP1(2)];
    [okB, rotDegB, diagB] = solveRot(uY1m,uP1m,uY6,uP6, rotTolDeg);
    if ~okB
        error("Alignment failed. rot disagreement: no-mirror=%.2f deg, mirror=%.2f deg.", ...
            diagA.deltaDeg, diagB.deltaDeg);
    end
    doMirror = true;
    rotDeg   = rotDegB;
    diagRot  = diagB;
end

% Build transform
T = makeTform(C1, C6, rotDeg, doMirror);

% Centroid check (after transform)
[yx, yy] = transformPointsForward(T, Y1(1), Y1(2));
[px, py] = transformPointsForward(T, P1(1), P1(2));
errY = norm([yx yy] - Y6);
errP = norm([px py] - P6);

out = struct();
out.ok = (errY <= centTolPx) && (errP <= centTolPx);
out.mirror = doMirror;
out.rotDeg = rotDeg;
out.tform  = T;
out.errY = errY;
out.errP = errP;

% diagnostics
out.diag.rotY = diagRot.rotY;
out.diag.rotP = diagRot.rotP;
out.diag.deltaDeg = diagRot.deltaDeg;

if ~out.ok
    out.reason = sprintf("Centroid check failed: errY=%.2f px, errP=%.2f px.", errY, errP);
else
    out.reason = "OK";
end
end

% ===================== helpers =====================

function [ok, rotDeg, diag] = solveRot(uY1,uP1,uY6,uP6, rotTolDeg)
% Compute rotY/rotP and test agreement
rotY = angleDeg(uY6) - angleDeg(uY1);
rotP = angleDeg(uP6) - angleDeg(uP1);
delta = circDiffDeg(rotY, rotP);

ok = (delta <= rotTolDeg);
rotDeg = circMeanDeg([rotY rotP]);

diag = struct('rotY',rotY,'rotP',rotP,'deltaDeg',delta);
end

function ang = angleDeg(v)
ang = atan2d(v(2), v(1));
end

function d = circDiffDeg(a, b)
% smallest absolute difference on circle (degrees)
d = abs(a - b);
d = mod(d, 360);
d = min(d, 360 - d);
end

function m = circMeanDeg(angles)
% circular mean in degrees
m = atan2d(mean(sind(angles)), mean(cosd(angles)));
end

function T = makeTform(C1, C6, rotDeg, doMirror)
% affine2d uses ROW-vector convention: [x y 1] * T

T_to0 = [1 0 0;
         0 1 0;
        -C1(1) -C1(2) 1];

if doMirror
    % mirror about vertical axis through origin in centered coords: (dx,dy)->(-dx,dy)
    M = [-1 0 0;
          0 1 0;
          0 0 1];
else
    M = eye(3);
end

R = [ cosd(rotDeg)  sind(rotDeg) 0;
     -sind(rotDeg)  cosd(rotDeg) 0;
      0             0            1];

T_back = [1 0 0;
          0 1 0;
          C6(1) C6(2) 1];

T = affine2d(T_to0 * M * R * T_back);
end
