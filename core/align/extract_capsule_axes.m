function F = extract_capsule_axes(I, minAreaPx)
% Extract centroid + major axis for yellow and purple capsules

if nargin < 2, minAreaPx = 200; end
if size(I,3)==1, I = repmat(I,1,1,3); end
Iu = uint8(I);

% --- estimate background from borders ---
b = 10;
border = cat(1, ...
    reshape(Iu(1:b,:,:),[],3), ...
    reshape(Iu(end-b+1:end,:,:),[],3), ...
    reshape(Iu(:,1:b,:),[],3), ...
    reshape(Iu(:,end-b+1:end,:),[],3));
bg = uint8(round(median(double(border),1)));

% --- get object colors ---
cols = unique(reshape(Iu,[],3),'rows');
isBg = all(cols==bg,2);
if ~any(isBg)
    d = sum((double(cols)-double(bg)).^2,2);
    [~,k] = min(d); isBg(k)=true;
end
objCols = cols(~isBg,:);
assert(size(objCols,1)==2,'Expected exactly 2 object colors.');

% purple has higher blue channel
[~,ip] = max(objCols(:,3));
iy = 3-ip;
cPurple = objCols(ip,:);
cYellow = objCols(iy,:);

purple = all(Iu==reshape(cPurple,1,1,3),3);
yellow = all(Iu==reshape(cYellow,1,1,3),3);

purple = keepLargestCC(bwareaopen(purple,minAreaPx));
yellow = keepLargestCC(bwareaopen(yellow,minAreaPx));

rp = regionprops(purple,'Centroid','Orientation','MajorAxisLength','Area');
ry = regionprops(yellow,'Centroid','Orientation','MajorAxisLength','Area');

assert(~isempty(rp)&&~isempty(ry),'Could not detect both objects.');

F.yellow = ry(1);
F.purple = rp(1);
F.yellowMask = yellow;
F.purpleMask = purple;

% --- major-axis unit vectors in IMAGE coordinates (y down!) ---
vy = ori2vec_img(F.yellow.Orientation);
vp = ori2vec_img(F.purple.Orientation);

F.yellow.vAxis = vy;  % 1x2 unit vector
F.purple.vAxis = vp;

% --- intersection of the two axis LINES (stimulus center) ---
cy = F.yellow.Centroid;  % [x y]
cp = F.purple.Centroid;

[F.stimCenter_xy, F.axisIntersect_ok] = lineIntersectionLS(cy, vy, cp, vp);


end

function BW = keepLargestCC(BW)
cc = bwconncomp(BW);
if cc.NumObjects<=1, return; end
[~,k] = max(cellfun(@numel,cc.PixelIdxList));
BW = false(size(BW));
BW(cc.PixelIdxList{k}) = true;
end

function v = ori2vec_img(oriDeg)
% Convert regionprops Orientation (deg) to a unit direction vector in image coords
a = deg2rad(oriDeg);
v = [cos(a), -sin(a)];   % NOTE the minus for image y-down coords
v = v / norm(v);
end

function [p, ok] = lineIntersectionLS(c1, v1, c2, v2)
% Intersection of two infinite lines:
%   L1: c1 + t*v1
%   L2: c2 + s*v2
%
% Solve c1 + t*v1 = c2 + s*v2  =>  [v1  -v2]*[t;s] = (c2-c1)
A = [v1(:), -v2(:)];          % 2x2
b = (c2(:) - c1(:));          % 2x1

% If nearly parallel, least-squares still returns something but it can be unstable.
if abs(det(A)) < 1e-6
    ok = false;
    ts = A\b;                 % least squares
else
    ok = true;
    ts = A\b;
end

t = ts(1);
p = c1 + t*v1;                % [x y]
end

