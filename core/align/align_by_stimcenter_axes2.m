function best = align_by_stimcenter_axes2(Fmov, Ffix, varargin)
%ALIGN_BY_STIMCENTER_AXES  Rigid alignment using stimulus center + handedness.
%
% Assumes extract_capsule_axes has produced:
%   F.stimCenter_xy              [x y]
%   F.yellow.Centroid            [x y]
%   F.purple.Centroid            [x y]
%   (optional) F.yellowMask, F.purpleMask for IoU evaluation
%
% Logic:
%   1) Determine whether a mirror is needed by comparing handedness
%      of (yellow, purple) around the center (sign of 2D cross product).
%   2) If mirror needed, reflect MOVING coordinates about a horizontal line
%      through the moving center (any line reflection would work up to rotation).
%   3) Compute unique rotation angle that maps (center->purple) onto (center->purple).
%      (Optionally average with center->yellow for robustness.)
%   4) Build ONE affine2d mapping moving->fixed: translate to origin, mirror(if),
%      rotate, translate to fixed center.
%
% Optional name/value:
%   'UseBothColors' (default true): average rotation from yellow+purple vectors.
%   'ComputeIoU'    (default true if masks exist): compute IoU after warping union masks.

p = inputParser;
p.addParameter('UseBothColors', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ComputeIoU', [], @(x)islogical(x)&&isscalar(x)); % auto if []
p.parse(varargin{:});
UseBoth = p.Results.UseBothColors;

haveMasks = isfield(Fmov,'yellowMask') && isfield(Fmov,'purpleMask') && ...
            isfield(Ffix,'yellowMask') && isfield(Ffix,'purpleMask');

if isempty(p.Results.ComputeIoU)
    doIoU = haveMasks;
else
    doIoU = p.Results.ComputeIoU;
end

C1 = Fmov.stimCenter_xy(:)';   % [x y]
C6 = Ffix.stimCenter_xy(:)';

Y1 = Fmov.yellow.Centroid(:)';  P1 = Fmov.purple.Centroid(:)';
Y6 = Ffix.yellow.Centroid(:)';  P6 = Ffix.purple.Centroid(:)';

% --- 1) handedness test: sign(cross(C->Y, C->P)) ---
s1 = sign(cross2(Y1 - C1, P1 - C1));
s6 = sign(cross2(Y6 - C6, P6 - C6));

% If one is zero (rare), fall back to no-mirror assumption
if s1 == 0 || s6 == 0
    mirrorNeeded = false;
else
    mirrorNeeded = (s1 ~= s6);
end

% --- 2) reflect moving vectors around horizontal line through C1 if needed ---
% Reflection chosen: (dx,dy) -> (dx,-dy) in coordinates centered at C1
% This is sufficient because any line reflection differs only by a rotation.
uP1 = P1 - C1;
uY1 = Y1 - C1;

if mirrorNeeded
    uP1(2) = -uP1(2);
    uY1(2) = -uY1(2);
end

uP6 = P6 - C6;
uY6 = Y6 - C6;

% --- 3) unique rotation (use purple; optionally average with yellow) ---
rotP = angleDeg(uP6) - angleDeg(uP1);

if UseBoth
    rotY = angleDeg(uY6) - angleDeg(uY1);
    rotDeg = circMeanDeg([rotP rotY]);   % robust average on circle
else
    rotDeg = rotP;
end
rotDeg = wrapTo180(rotDeg);

% --- 4) build affine2d in MATLAB row-vector convention ---
T = makeCenteredMirrorRotateTranslate(rotDeg, mirrorNeeded, C1, C6);

best = struct();
best.mirrorNeeded = mirrorNeeded;
best.rotDeg = rotDeg;
best.tform = T;

% --- optional IoU for debugging ---
best.iou = NaN;
if doIoU
    BW1 = (Fmov.yellowMask | Fmov.purpleMask);
    BW6 = (Ffix.yellowMask | Ffix.purpleMask);
    BWw = imwarp(BW1, T, 'OutputView', imref2d(size(BW6)), 'Interp','nearest');
    best.iou = nnz(BWw & BW6) / nnz(BWw | BW6);
end

end

% ===================== helpers =====================

function z = cross2(a,b)
% 2D cross product (signed area): a_x b_y - a_y b_x
z = a(1)*b(2) - a(2)*b(1);
end

function ang = angleDeg(v)
ang = atan2d(v(2), v(1));
end

function m = circMeanDeg(anglesDeg)
% circular mean of angles in degrees
x = mean(cosd(anglesDeg));
y = mean(sind(anglesDeg));
m = atan2d(y,x);
end

function T = makeCenteredMirrorRotateTranslate(rotDeg, mirrorNeeded, Cmov, Cfix)
% Row-vector convention: [x y 1] * T
T_to0 = [1 0 0;
         0 1 0;
        -Cmov(1) -Cmov(2) 1];

R = [cosd(rotDeg) -sind(rotDeg) 0;
     sind(rotDeg)  cosd(rotDeg) 0;
     0            0            1];

if mirrorNeeded
    % reflect about x-axis in centered coords: (dx,dy)->(dx,-dy)
    M = [1  0 0;
         0 -1 0;
         0  0 1];
    core = M * R;
else
    core = R;
end

T_back = [1 0 0;
          0 1 0;
          Cfix(1) Cfix(2) 1];

T = affine2d(T_to0 * core * T_back);
end
