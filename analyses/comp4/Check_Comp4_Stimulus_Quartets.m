% Check_Comp4_Stimulus_Quartets
% Verify the quartet structure of all 320 Nilson comp4 stimuli.

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(repoRoot);

cfg = config();
bitmapDir = fullfile(cfg.extrasRoot, 'monkeyN', 'comp4');
namedStimDir = fullfile(cfg.extrasRoot, 'monkeyN', 'exp', 'comp4');
outputDir = fullfile(cfg.resultsRoot, 'Comp4');
metadataDir = fullfile(outputDir, 'metadata');

assert(exist(bitmapDir, 'dir') == 7, 'Missing bitmap directory: %s', bitmapDir);
assert(exist(namedStimDir, 'dir') == 7, ...
    'Missing named-stimulus directory: %s', namedStimDir);
if exist(metadataDir, 'dir') ~= 7
    mkdir(metadataDir);
end

%% Parse stimulus metadata from the informative PNG filenames
files = dir(fullfile(namedStimDir, '*.png'));
assert(numel(files) == 320, 'Expected 320 PNG stimuli, found %d.', numel(files));

stimulus = nan(320, 1);
array = nan(320, 1);
objectClass = nan(320, 1);
complexity = nan(320, 1);
cones = nan(320, 1);

pattern = ['^array_(\d+)_class(\d+)_complexity(\d+)_' ...
           'cones(\d+)_stim_(\d+)\.png$'];
for fileIdx = 1:numel(files)
    tokens = regexp(files(fileIdx).name, pattern, 'tokens', 'once');
    assert(~isempty(tokens), 'Unexpected comp4 filename: %s', files(fileIdx).name);
    values = cellfun(@str2double, tokens);
    row = values(5);
    assert(row >= 1 && row <= 320 && isnan(stimulus(row)), ...
        'Invalid or duplicate stimulus number in %s.', files(fileIdx).name);
    array(row) = values(1);
    objectClass(row) = values(2);
    complexity(row) = values(3);
    cones(row) = values(4);
    stimulus(row) = row;
end

assert(isequal(stimulus, (1:320).'), ...
    'Stimulus numbers must form the complete sequence 1:320.');

%% Assign the four conditions within each quartet
quartet = ceil(stimulus / 4);
positionInQuartet = mod(stimulus - 1, 4) + 1;
foregroundCode = ceil(positionInQuartet / 2);

arrayObject = strings(320, 1);
arrayObject(objectClass == 1) = "crocodile";
arrayObject(objectClass == 2) = "monkey";
assert(all(arrayObject ~= ""), 'Unexpected object class outside 1:2.');

otherObject = strings(320, 1);
otherObject(arrayObject == "crocodile") = "monkey";
otherObject(arrayObject == "monkey") = "crocodile";

frontObject = otherObject;
frontObject(foregroundCode == 2) = arrayObject(foregroundCode == 2);
attendedObject = frontObject;

attendedColor = repmat("purple", 320, 1);
attendedColor(mod(positionInQuartet, 2) == 0) = "yellow";

stimulusConditions = table(stimulus, quartet, positionInQuartet, array, ...
    objectClass, arrayObject, foregroundCode, frontObject, attendedObject, ...
    attendedColor, complexity, cones);

%% Verify metadata and image consistency within all 80 quartets
background = uint8([128 128 128]);
purple = uint8([153 153 179]);
yellow = uint8([153 158 120]);

metadataConsistent = false(80, 1);
unionMaskConsistent = false(80, 1);
firstPairColorSwap = false(80, 1);
secondPairColorSwap = false(80, 1);

for quartetIdx = 1:80
    rows = (quartetIdx - 1) * 4 + (1:4);
    metadata = [array(rows), objectClass(rows), complexity(rows), cones(rows)];
    metadataConsistent(quartetIdx) = size(unique(metadata, 'rows'), 1) == 1;

    images = cell(4, 1);
    masks = cell(4, 1);
    for position = 1:4
        imagePath = fullfile(bitmapDir, sprintf('%03d.bmp', rows(position)));
        assert(exist(imagePath, 'file') == 2, 'Missing bitmap: %s', imagePath);
        images{position} = imread(imagePath);
        masks{position} = any(images{position} ~= ...
            reshape(background, 1, 1, 3), 3);
    end

    unionMaskConsistent(quartetIdx) = all(cellfun( ...
        @(mask) isequal(mask, masks{1}), masks(2:4)));
    firstPairColorSwap(quartetIdx) = isequal( ...
        swapColors(images{1}, purple, yellow), images{2});
    secondPairColorSwap(quartetIdx) = isequal( ...
        swapColors(images{3}, purple, yellow), images{4});
end

assert(all(metadataConsistent), 'Metadata differs within one or more quartets.');
assert(all(unionMaskConsistent), ...
    'Object-union geometry differs within one or more quartets.');
assert(all(firstPairColorSwap), ...
    'Stimuli 1 and 2 are not exact color swaps in one or more quartets.');
assert(all(secondPairColorSwap), ...
    'Stimuli 3 and 4 are not exact color swaps in one or more quartets.');

%% Report and save the compact audit
fprintf('Verified %d stimuli forming %d complete quartets.\n', ...
    height(stimulusConditions), max(stimulusConditions.quartet));
fprintf('All quartets have consistent metadata and object-union geometry.\n');
fprintf('Positions 1/2 show the non-array object in front: purple/yellow.\n');
fprintf('Positions 3/4 show the array object in front: purple/yellow.\n');
fprintf('Class 1 (crocodile at array): monkey first, crocodile second.\n');
fprintf('Class 2 (monkey at array): crocodile first, monkey second.\n');

summary = struct();
summary.nStimuli = height(stimulusConditions);
summary.nQuartets = max(stimulusConditions.quartet);
summary.arrays = unique(stimulusConditions.array).';
summary.nQuartetsPerArray = arrayfun( ...
    @(value) nnz(stimulusConditions.array == value) / 4, summary.arrays);
summary.nClass1Quartets = nnz(stimulusConditions.objectClass == 1) / 4;
summary.nClass2Quartets = nnz(stimulusConditions.objectClass == 2) / 4;
summary.metadataConsistent = all(metadataConsistent);
summary.unionMaskConsistent = all(unionMaskConsistent);
summary.colorSwapsConsistent = ...
    all(firstPairColorSwap) && all(secondPairColorSwap);

writetable(stimulusConditions, ...
    fullfile(metadataDir, 'comp4_stimulus_conditions.csv'));
save(fullfile(metadataDir, 'comp4_stimulus_conditions.mat'), ...
    'stimulusConditions', 'summary');


function swapped = swapColors(image, colorA, colorB)
swapped = image;
maskA = all(image == reshape(colorA, 1, 1, 3), 3);
maskB = all(image == reshape(colorB, 1, 1, 3), 3);
for channel = 1:3
    plane = swapped(:,:,channel);
    plane(maskA) = colorB(channel);
    plane(maskB) = colorA(channel);
    swapped(:,:,channel) = plane;
end
end
