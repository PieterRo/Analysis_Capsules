function make_activity_movie(Tall, ALLCOORDS, RTAB384, exampleStimNum, R, SNR, outFile)

% ---- Settings ----
useOnlyV1 = true;
clipRange = [-1 2];
alphaMax  = 0.85;
alphaThresh = 0.10;
markerSize = 18;

W = 1024; H = 768;
toPx = @(q) [q(1) + W/2, H/2 - q(2)];

fprintf('Precomputing geometry...\n');

% ----------------------------------------------------------
% 1) Compute projected RF coordinates ONCE (example frame)
% ----------------------------------------------------------

fieldName = sprintf('stim_%d', exampleStimNum);
s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

s_px = toPx(s);
tT   = toPx(tFig);
tD   = toPx(tBack);

uT = (tT - s_px) / norm(tT - s_px);
uD = (tD - s_px) / norm(tD - s_px);

nT = perpTowardOther(s_px, tT, tD);
nD = perpTowardOther(s_px, tD, tT);

widthEx = double(RTAB384(exampleStimNum,7));

if useOnlyV1
    sites = 1:512;
else
    sites = 1:size(R.meanAct,1);
end

X = [];
Y = [];
siteID = [];

seenPair = false(384,1);

for stimNum = 1:384
    comp = complementaryStim(stimNum);
    pairKey = min(stimNum, comp);
    if seenPair(pairKey), continue; end
    seenPair(pairKey) = true;

    T = Tall(stimNum).T;
    T = T(sites,:);
    assign = string(T.assignment);

    idx = (assign=="target") | (assign=="distractor");
    if ~any(idx), continue; end

    along = T.along_GC(idx) * widthEx;
    perp  = T.perp_signed_GC(idx) * widthEx;

    isTarget = assign(idx)=="target";
    baseDir  = zeros(sum(idx),2);
    basePerp = zeros(sum(idx),2);

    baseDir(isTarget,:)  = repmat(uT,sum(isTarget),1);
    baseDir(~isTarget,:) = repmat(uD,sum(~isTarget),1);

    basePerp(isTarget,:)  = repmat(nT,sum(isTarget),1);
    basePerp(~isTarget,:) = repmat(nD,sum(~isTarget),1);

    p = s_px + along.*baseDir + perp.*basePerp;

    X = [X; p(:,1)];
    Y = [Y; p(:,2)];
    siteID = [siteID; sites(idx)];
end

fprintf('Geometry done: %d plotted RF points\n', numel(X));

% ----------------------------------------------------------
% 2) Setup figure + scatter once
% ----------------------------------------------------------

figure('Color',[0.4 0.4 0.4]);
ax = axes('Position',[0 0 1 1]); hold(ax,'on');

img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, exampleStimNum);
imshow(img,'Parent',ax,'InitialMagnification','fit');
set(ax,'YDir','reverse');
axis(ax,'equal');

hSc = scatter(ax, X, Y, markerSize, zeros(numel(X),3), 'filled');
hSc.MarkerEdgeAlpha = 0;
hSc.MarkerFaceAlpha = 'flat';
hSc.AlphaData = zeros(numel(X),1);
hSc.AlphaDataMapping = 'none';

drawnow;

% ----------------------------------------------------------
% 3) Setup normalization
% ----------------------------------------------------------

muSpont = double(SNR.muSpont(:));
muTop   = max([SNR.muYellowEarly(:), SNR.muYellowLate(:), ...
               SNR.muPurpleEarly(:), SNR.muPurpleLate(:)],[],2);

scale = muTop - muSpont;
scale(scale<=1e-9)=1e-9;

% ----------------------------------------------------------
% 4) Setup video writer
% ----------------------------------------------------------

v = VideoWriter(outFile,'MPEG-4');
v.FrameRate = 10;
open(v);

fprintf('Rendering movie...\n');

% ----------------------------------------------------------
% 5) Loop time bins
% ----------------------------------------------------------

for tb = 1:70

    fprintf('Time bin %d / 70\n', tb);

    act = squeeze(double(R.meanAct(siteID, :, tb)));

    % Average complementary pairs
    actPair = zeros(size(act,1),1);
    count = 0;
    seenPair(:)=false;

    for stimNum = 1:384
        comp = complementaryStim(stimNum);
        pairKey = min(stimNum, comp);
        if seenPair(pairKey), continue; end
        seenPair(pairKey)=true;
        count = count + 1;
        actPair = actPair + 0.5*(act(:,stimNum)+act(:,comp));
    end

    actPair = actPair / count;

    z = (actPair - muSpont(siteID)) ./ scale(siteID);
    z = min(max(z,clipRange(1)),clipRange(2));

    [rgb, alpha] = valueToColorAlpha(z, alphaMax, alphaThresh, clipRange);

    hSc.CData = rgb;
    hSc.AlphaData = alpha;

    title(sprintf('%d–%d ms', R.timeWindows(tb,1), R.timeWindows(tb,2)), ...
          'Color','w','FontSize',16);

    drawnow;

    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

fprintf('Movie saved to %s\n', outFile);
end