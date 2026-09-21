% Diagnose_Comp4_ALLCOORDS
% Determine which saved ALLCOORDS generation matches the current bitmaps.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);

cfg = config();
stimDir = fullfile(cfg.extrasRoot, 'monkeyN', 'comp4');
logDir = fullfile(cfg.extrasRoot, 'monkeyN', '_logs');
geometryFile = fullfile(logDir, 'temp_figs.mat');
sources = {
    fullfile(logDir, 'ObjAtt_Shapes_comp4_monkeyN_20220502_B1.mat')
    fullfile(logDir, '_old', 'RANDTAB_shapes_comp4_monkeyN.mat')
    };
sourceNames = {'recording log 2022-05-02 B1', 'regenerated set 2022-05-03'};

G = load(geometryFile, 'temp_figures');
objects = {G.temp_figures.class1_complexity4, ...
           G.temp_figures.class2_complexity4};
purple = uint8([153 153 179]);
HX = 512;
HY = 384;

for sourceIdx = 1:numel(sources)
    S = load(sources{sourceIdx}, 'ALLMAT', 'ALLCOORDS');
    scores = nan(80, 2);
    angleDifference = nan(80, 1);

    for quartet = 1:80
        firstStim = 4 * (quartet - 1) + 1;
        row = S.ALLMAT(firstStim,:);
        primaryClass = row(5);
        secondaryClass = 3 - primaryClass;
        primaryCue = row(7);

        coordinates = S.ALLCOORDS.(sprintf('stim_%d', firstStim));
        primaryTargets = [coordinates.s; coordinates.t_back];
        secondaryTargets = [coordinates.s; coordinates.t_fig];

        primaryMask = frontPurpleMask(stimDir, firstStim + 2, purple);
        secondaryMask = frontPurpleMask(stimDir, firstStim, purple);

        primaryTransform = twoPointSimilarity( ...
            double(objects{primaryClass}.cues([primaryCue, 3-primaryCue],:)), ...
            double(primaryTargets), false);
        scores(quartet,1) = maskFit(objects{primaryClass}, ...
            primaryTransform, primaryMask, HX, HY);

        rr = row(3) / 2;
        if primaryCue == 2
            rr = 4 - rr;
        end
        rfCenter = applySimilarity(double(objects{primaryClass}.RFs(rr,:)), ...
            primaryTransform);
        rfRadius = primaryTransform.scale * objects{primaryClass}.mean_size;
        secondaryTransform = chooseSecondaryTransform( ...
            objects{secondaryClass}, secondaryTargets, rfCenter, rfRadius);
        scores(quartet,2) = maskFit(objects{secondaryClass}, ...
            secondaryTransform, secondaryMask, HX, HY);

        bitmapVector = primaryTargets(2,:) - primaryTargets(1,:);
        loggedAngle = atan2d(bitmapVector(2), bitmapVector(1));
        angleDifference(quartet) = wrapTo180Local(loggedAngle - row(6));
    end

    fprintf('\n%s\n', sourceNames{sourceIdx});
    fprintf('  median IoU, primary object:   %.4f\n', median(scores(:,1)));
    fprintf('  median IoU, secondary object: %.4f\n', median(scores(:,2)));
    fprintf('  minimum IoU, either object:   %.4f\n', min(scores(:)));
    fprintf('  quartets with both IoU > .98: %d / 80\n', ...
        nnz(all(scores > .98, 2)));
    fprintf('  first five primary IoUs:      %s\n', ...
        mat2str(scores(1:5,1).', 4));
end

% The two coordinate tables should differ if the stimulus set was regenerated.
A = load(sources{1}, 'ALLMAT', 'ALLCOORDS');
B = load(sources{2}, 'ALLMAT', 'ALLCOORDS');
angleDelta = abs(wrapTo180Local(A.ALLMAT(:,6) - B.ALLMAT(:,6)));
cueDelta = nan(320, 1);
for stim = 1:320
    fieldName = sprintf('stim_%d', stim);
    cueDelta(stim) = norm(A.ALLCOORDS.(fieldName).s - ...
        B.ALLCOORDS.(fieldName).s);
end
fprintf('\nDifference between the two saved generations\n');
fprintf('  stimuli with different angle: %d / 320\n', nnz(angleDelta > 1e-9));
fprintf('  median absolute angle change: %.2f deg\n', median(angleDelta));
fprintf('  median central-cue displacement: %.2f px\n', median(cueDelta));
fprintf('  maximum central-cue displacement: %.2f px\n', max(cueDelta));


function mask = frontPurpleMask(stimDir, stimulus, purple)
image = imread(fullfile(stimDir, sprintf('%03d.bmp', stimulus)));
mask = all(image == reshape(purple, 1, 1, 3), 3);
end


function score = maskFit(object, transform, observed, HX, HY)
polygon = double([object.X(:), object.Y(:)]);
xy = applySimilarity(polygon, transform);
predicted = poly2mask(xy(:,1) + HX, HY - xy(:,2), HY * 2, HX * 2);
score = nnz(predicted & observed) / nnz(predicted | observed);
end


function transform = chooseSecondaryTransform(object, targets, rfCenter, rfRadius)
theta = 0:pi/180:2*pi;
rfX = rfRadius * cos(theta) + rfCenter(1);
rfY = rfRadius * sin(theta) + rfCenter(2);
source = double(object.cues([1 2],:));
polygon = double([object.X(:), object.Y(:)]);
cost = nan(1, 2);
transforms = cell(1, 2);
for idx = 1:2
    transforms{idx} = twoPointSimilarity(source, double(targets), idx == 2);
    xy = applySimilarity(polygon, transforms{idx});
    cost(idx) = nnz(inpolygon(rfX, rfY, xy(:,1), xy(:,2)));
end
[~, best] = min(cost);
transform = transforms{best};
end


function transform = twoPointSimilarity(source, target, reflected)
sourceVector = source(2,:) - source(1,:);
targetVector = target(2,:) - target(1,:);
sourceUnit = sourceVector / norm(sourceVector);
targetUnit = targetVector / norm(targetVector);
sourceBasis = [sourceUnit; -sourceUnit(2), sourceUnit(1)];
targetNormal = [-targetUnit(2), targetUnit(1)];
if reflected
    targetNormal = -targetNormal;
end
targetBasis = [targetUnit; targetNormal];
transform.scale = norm(targetVector) / norm(sourceVector);
transform.rotation = sourceBasis.' * targetBasis;
transform.sourceOrigin = source(1,:);
transform.targetOrigin = target(1,:);
end


function xy = applySimilarity(xy, transform)
xy = transform.scale * (xy - transform.sourceOrigin) * ...
    transform.rotation + transform.targetOrigin;
end


function value = wrapTo180Local(value)
value = mod(value + 180, 360) - 180;
end
