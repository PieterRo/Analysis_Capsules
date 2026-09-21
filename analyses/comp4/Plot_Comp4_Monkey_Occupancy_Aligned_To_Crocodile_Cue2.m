% Plot_Comp4_Monkey_Occupancy_Aligned_To_Crocodile_Cue2
% Align every recorded configuration with crocodile cue 2 to the canonical
% crocodile and show the percentage of configurations occupied by monkey.

%% Paths and recorded-session data
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);

cfg = config();
comp4Dir = fullfile(cfg.resultsRoot, 'Comp4');
metadataDir = fullfile(comp4Dir, 'metadata');
figuresDir = fullfile(comp4Dir, 'figures');
geometryFile = fullfile(cfg.extrasRoot, 'monkeyN', '_logs', 'temp_figs.mat');
sessionFile = fullfile(cfg.extrasRoot, 'monkeyN', '_logs', ...
    'ObjAtt_Shapes_comp4_monkeyN_20220502_B1.mat');
trialFile = fullfile(cfg.dataRoot, 'Mr Nilson', 'comp4', ...
    'ObjAtt_Shapes_comp4_MUA_trials.mat');
cueFile = fullfile(metadataDir, 'comp4_cue_geometry.mat');

assert(exist(cueFile, 'file') == 2, ...
    'Run Plot_Comp4_Aligned_Cue_Locations.m first.');
assert(exist(geometryFile, 'file') == 2, 'Missing %s.', geometryFile);
assert(exist(sessionFile, 'file') == 2, 'Missing %s.', sessionFile);
assert(exist(trialFile, 'file') == 2, 'Missing %s.', trialFile);
if exist(figuresDir, 'dir') ~= 7
    mkdir(figuresDir);
end

G = load(geometryFile, 'temp_figures');
S = load(sessionFile, 'ALLMAT', 'ALLCOORDS');
T = load(trialFile, 'ALLMAT');
Q = load(cueFile, 'cueGeometry');
crocodile = G.temp_figures.class1_complexity4;
monkey = G.temp_figures.class2_complexity4;
objects = {crocodile, monkey};
cueGeometry = Q.cueGeometry;

% Guard against accidentally mixing the recording with another generation.
for stimulus = 1:320
    trialRows = double(T.ALLMAT(T.ALLMAT(:,1) == stimulus, 1:8));
    assert(~isempty(trialRows), 'Stimulus %d is absent from trial ALLMAT.', stimulus);
    assert(all(max(abs(trialRows - trialRows(1,:)), [], 1) < 1e-9), ...
        'Stimulus %d has inconsistent trial metadata.', stimulus);
    assert(max(abs(trialRows(1,:) - S.ALLMAT(stimulus,:))) < 1e-9, ...
        'Stimulus %d does not match the session ALLMAT.', stimulus);
end

%% Recover every unique configuration with crocodile cue 2
cue2Rows = find(cueGeometry.crocodileCueIndex == 2);
cue2Quartets = unique(cueGeometry.quartet(cue2Rows), 'stable');
selectedStimuli = 4 * (cue2Quartets - 1) + 1;
nConfigurations = numel(selectedStimuli);
nStimulusImages = numel(cue2Rows);
assert(nStimulusImages == 4 * nConfigurations, ...
    'Expected four stimulus images per unique configuration.');

alignedMonkeys = cell(nConfigurations, 1);
secondaryReflection = false(nConfigurations, 1);
for configuration = 1:nConfigurations
    stimulus = selectedStimuli(configuration);
    row = S.ALLMAT(stimulus,:);
    coordinates = S.ALLCOORDS.(sprintf('stim_%d', stimulus));

    [transforms, reflected] = recoverObjectTransforms( ...
        objects, row, coordinates);
    crocodileTransform = transforms{1};
    monkeyTransform = transforms{2};
    secondaryReflection(configuration) = reflected;

    monkeyDisplay = applySimilarity( ...
        double([monkey.X(:), monkey.Y(:)]), monkeyTransform);
    alignedMonkeys{configuration} = invertSimilarity( ...
        monkeyDisplay, crocodileTransform);
end

%% Rasterize monkey occupancy in canonical-crocodile coordinates
allX = double(crocodile.X(:));
allY = double(crocodile.Y(:));
for configuration = 1:nConfigurations
    allX = [allX; alignedMonkeys{configuration}(:,1)]; %#ok<AGROW>
    allY = [allY; alignedMonkeys{configuration}(:,2)]; %#ok<AGROW>
end
margin = 20;
pixelSize = 2;
xMin = floor(min(allX)) - margin;
xMax = ceil(max(allX)) + margin;
yMin = floor(min(allY)) - margin;
yMax = ceil(max(allY)) + margin;
xValues = xMin:pixelSize:xMax;
yValues = yMin:pixelSize:yMax;
nColumns = numel(xValues);
nRows = numel(yValues);

occupancyCount = zeros(nRows, nColumns);
for configuration = 1:nConfigurations
    xy = alignedMonkeys{configuration};
    columns = (xy(:,1) - xMin) / pixelSize + 1;
    rows = (xy(:,2) - yMin) / pixelSize + 1;
    occupancyCount = occupancyCount + ...
        poly2mask(columns, rows, nRows, nColumns);
end
crocodileColumns = (double(crocodile.X) - xMin) / pixelSize + 1;
crocodileRows = (double(crocodile.Y) - yMin) / pixelSize + 1;
crocodileMask = poly2mask( ...
    crocodileColumns, crocodileRows, nRows, nColumns);
occupancyPct = 100 * occupancyCount / nConfigurations;
occupancyPct(crocodileMask | occupancyCount == 0) = NaN;

%% Plot occupancy outside the canonical crocodile
h = figure('Color', 'w', ...
    'Name', 'Comp4 monkey occupancy aligned to crocodile cue 2', ...
    'NumberTitle', 'off', 'Position', [100 100 1100 800]);
ax = axes(h);
hold(ax, 'on');

mapHandle = imagesc(ax, xValues, yValues, occupancyPct);
set(ax, 'YDir', 'normal');
set(mapHandle, 'AlphaData', 0.92 * double(~isnan(occupancyPct)));
colormap(ax, occupancyColormap(256));
caxis(ax, [0 100]);
cb = colorbar(ax);
cb.Label.String = sprintf('Monkey occupancy (%% of %d configurations)', ...
    nConfigurations);

fill(ax, crocodile.X, crocodile.Y, [0.43 0.46 0.47], ...
    'EdgeColor', [0.16 0.18 0.18], 'LineWidth', 2.0);
scatter(ax, crocodile.cues(2,1), crocodile.cues(2,2), 180, ...
    [0.10 0.45 0.82], 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
text(ax, crocodile.cues(2,1), crocodile.cues(2,2), 'Cue 2  ', ...
    'Color', 'w', 'FontWeight', 'bold', 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'right');

axis(ax, 'equal');
axis(ax, 'off');
title(ax, sprintf(['Monkey occupancy after aligning all crocodile-cue-2 ' ...
    'configurations (n = %d)'], nConfigurations), 'FontWeight', 'bold');

fileStem = 'comp4_monkey_occupancy_aligned_to_crocodile_cue2';
savefig(h, fullfile(figuresDir, [fileStem '.fig']));
print(h, fullfile(figuresDir, [fileStem '.png']), '-dpng', '-r300');
sourceSessionFile = sessionFile;
save(fullfile(metadataDir, [fileStem '.mat']), ...
    'cue2Quartets', 'selectedStimuli', 'alignedMonkeys', ...
    'secondaryReflection', 'xValues', 'yValues', 'occupancyCount', ...
    'occupancyPct', 'crocodileMask', 'nConfigurations', 'nStimulusImages', ...
    'pixelSize', 'sourceSessionFile');

fprintf('Crocodile cue 2: %d stimulus images in %d unique configurations.\n', ...
    nStimulusImages, nConfigurations);
fprintf('Maximum monkey occupancy outside crocodile: %.1f%%.\n', ...
    max(occupancyPct(:), [], 'omitnan'));
fprintf('Saved occupancy plot to: %s\n', figuresDir);


function [transforms, secondaryReflected] = recoverObjectTransforms( ...
        objects, row, coordinates)
primaryClass = row(5);
secondaryClass = 3 - primaryClass;
primaryCue = row(7);
foreground = row(4);

if foreground == 1
    primaryOuterCue = coordinates.t_back;
    secondaryOuterCue = coordinates.t_fig;
else
    primaryOuterCue = coordinates.t_fig;
    secondaryOuterCue = coordinates.t_back;
end

primaryTargets = double([coordinates.s; primaryOuterCue]);
primarySource = double(objects{primaryClass}.cues( ...
    [primaryCue, 3-primaryCue],:));
primaryTransform = twoPointSimilarity(primarySource, primaryTargets, false);

rfIndex = row(3) / 2;
if primaryCue == 2
    rfIndex = 4 - rfIndex;
end
rfCenter = applySimilarity( ...
    double(objects{primaryClass}.RFs(rfIndex,:)), primaryTransform);
rfRadius = primaryTransform.scale * objects{primaryClass}.mean_size;

secondaryTargets = double([coordinates.s; secondaryOuterCue]);
[secondaryTransform, secondaryReflected] = chooseSecondaryTransform( ...
    objects{secondaryClass}, secondaryTargets, rfCenter, rfRadius);

transforms = cell(1, 2);
transforms{primaryClass} = primaryTransform;
transforms{secondaryClass} = secondaryTransform;
end


function [transform, reflected] = chooseSecondaryTransform( ...
        object, targets, rfCenter, rfRadius)
theta = 0:pi/180:2*pi;
rfX = rfRadius * cos(theta) + rfCenter(1);
rfY = rfRadius * sin(theta) + rfCenter(2);
source = double(object.cues([1 2],:));
polygon = double([object.X(:), object.Y(:)]);
cost = nan(1, 2);
transforms = cell(1, 2);
for idx = 1:2
    transforms{idx} = twoPointSimilarity(source, targets, idx == 2);
    xy = applySimilarity(polygon, transforms{idx});
    cost(idx) = nnz(inpolygon(rfX, rfY, xy(:,1), xy(:,2)));
end
[~, best] = min(cost);
transform = transforms{best};
reflected = best == 2;
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


function xy = invertSimilarity(xy, transform)
xy = ((xy - transform.targetOrigin) / transform.scale) / ...
    transform.rotation + transform.sourceOrigin;
end


function map = occupancyColormap(nColors)
anchors = [0.93 0.95 0.94; ...
           0.51 0.72 0.72; ...
           0.14 0.43 0.49; ...
           0.04 0.17 0.22];
anchorX = linspace(0, 1, size(anchors, 1));
mapX = linspace(0, 1, nColors);
map = interp1(anchorX, anchors, mapX, 'linear');
end
