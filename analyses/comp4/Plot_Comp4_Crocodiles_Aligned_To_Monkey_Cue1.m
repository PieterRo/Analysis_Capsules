% Plot_Comp4_Crocodiles_Aligned_To_Monkey_Cue1
% Reconstruct ten recorded configurations from the session ALLCOORDS,
% align their monkeys to the canonical monkey, and plot the crocodiles.

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

% The recorded trials must refer to this exact stimulus generation.
for stimulus = 1:320
    trialRows = double(T.ALLMAT(T.ALLMAT(:,1) == stimulus, 1:8));
    assert(~isempty(trialRows), 'Stimulus %d is absent from trial ALLMAT.', stimulus);
    assert(all(max(abs(trialRows - trialRows(1,:)), [], 1) < 1e-9), ...
        'Stimulus %d has inconsistent trial metadata.', stimulus);
    assert(max(abs(trialRows(1,:) - S.ALLMAT(stimulus,:))) < 1e-9, ...
        'Stimulus %d does not match the session ALLMAT.', stimulus);
end

%% Select ten spread-out configurations with front monkey and monkey cue 1
isMonkeyFront = strcmp(cueGeometry.attendedObject, 'monkey');
isMonkeyCue1 = cueGeometry.monkeyCueIndex == 1;
isPurple = strcmp(cueGeometry.attendedColor, 'purple');
candidateRows = find(isMonkeyFront & isMonkeyCue1 & isPurple);
assert(numel(candidateRows) >= 10, 'Fewer than ten matching configurations.');

pick = unique(round(linspace(1, numel(candidateRows), 10)));
assert(numel(pick) == 10, 'Could not select ten distinct examples.');
selectedRows = candidateRows(pick);
selectedStimuli = cueGeometry.stimulus(selectedRows);

%% Reconstruct each recorded configuration and align it to the monkey
alignedCrocodiles = cell(numel(selectedRows), 1);
secondaryReflection = false(numel(selectedRows), 1);
for example = 1:numel(selectedRows)
    stimulus = selectedStimuli(example);
    row = S.ALLMAT(stimulus,:);
    coordinates = S.ALLCOORDS.(sprintf('stim_%d', stimulus));

    [transforms, reflected] = recoverObjectTransforms( ...
        objects, row, coordinates);
    monkeyTransform = transforms{2};
    crocodileTransform = transforms{1};
    secondaryReflection(example) = reflected;

    crocodileDisplay = applySimilarity( ...
        double([crocodile.X(:), crocodile.Y(:)]), crocodileTransform);
    alignedCrocodiles{example} = invertSimilarity( ...
        crocodileDisplay, monkeyTransform);
end

%% Plot outlines behind the canonical monkey
h = figure('Color', 'w', ...
    'Name', 'Recorded Comp4 crocodiles aligned to monkey cue 1', ...
    'NumberTitle', 'off', 'Position', [100 100 1050 780]);
ax = axes(h);
hold(ax, 'on');

crocColor = [0.72 0.75 0.76];
for example = 1:numel(alignedCrocodiles)
    xy = alignedCrocodiles{example};
    plot(ax, xy(:,1), xy(:,2), 'Color', crocColor, 'LineWidth', 1.4);
end

fill(ax, monkey.X, monkey.Y, [0.43 0.46 0.47], ...
    'EdgeColor', [0.16 0.18 0.18], 'LineWidth', 2.0);
scatter(ax, monkey.cues(1,1), monkey.cues(1,2), 180, ...
    [0.10 0.45 0.82], 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
text(ax, monkey.cues(1,1), monkey.cues(1,2), '  Cue 1', ...
    'Color', 'w', 'FontWeight', 'bold', 'VerticalAlignment', 'middle');

axis(ax, 'equal');
axis(ax, 'off');
title(ax, ['Recorded configurations: front monkey aligned to ' ...
    'canonical cue 1'], 'FontWeight', 'bold');

fileStem = 'comp4_10_crocodiles_aligned_to_front_monkey_cue1';
savefig(h, fullfile(figuresDir, [fileStem '.fig']));
print(h, fullfile(figuresDir, [fileStem '.png']), '-dpng', '-r300');
sourceSessionFile = sessionFile;
save(fullfile(metadataDir, [fileStem '.mat']), ...
    'selectedStimuli', 'alignedCrocodiles', 'secondaryReflection', ...
    'sourceSessionFile');

fprintf('Selected recorded stimuli: %s\n', num2str(selectedStimuli.'));
fprintf('Verified trial metadata against session ALLMAT: 320 / 320 stimuli\n');
fprintf('Saved aligned configuration plot to: %s\n', figuresDir);


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

% Recover the RF circle used by the generator from the primary transform.
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
