function test_attention_modulation_stimulus_count()
% Check 216-stimulus inputs against an equivalent zero-trial-padded input.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'core', 'decoding'));
for perSiteCounts = [false true]
    [data, geometry, snr] = fixture(216, perSiteCounts);
    opts = struct('v1Sites', 1:4, 'verbose', false);
    small = attention_modulation_V1_3bin(data, geometry, snr, opts);
    paddedData = data;
    paddedData.meanAct(:,217:384,:) = 0;
    paddedData.meanSqAct(:,217:384,:) = 0;
    paddedData.nTrials(:,217:384) = 0;
    paddedGeometry = geometry;
    for stimulus = 217:384
        paddedGeometry(stimulus).T = table();
    end
    padded = attention_modulation_V1_3bin(paddedData, paddedGeometry, snr, opts);
    names = fieldnames(small);
    for i = 1:numel(names)
        if strcmp(names{i}, 'stimMaskUsed'); continue; end
        assert(isequaln(small.(names{i}), padded.(names{i})), ...
            'Padding changed %s.', names{i});
    end
    assert(numel(small.stimMaskUsed)==216 && ...
        numel(padded.stimMaskUsed)==384);
    [data, geometry, snr] = fixture(384, perSiteCounts);
    full = attention_modulation_V1_3bin(data, geometry, snr, opts);
    assert(all(full.validSite) && all(isfinite(full.dprime)));
end
fprintf('Stimulus-count checks passed for 216/384 and vector/site trial counts.\n');
end

function [data, geometry, snr] = fixture(nStimuli, perSiteCounts)
data = struct();
data.meanAct = zeros(4,nStimuli,3);
geometry = struct('T', cell(1,nStimuli));
for stimulus = 1:nStimuli
    target = mod(stimulus,2)==1;
    assignment = repmat("distractor",4,1);
    if target; assignment(:) = "target"; end
    center_color = repmat("purple",4,1);
    if mod(floor((stimulus-1)/2),2)==0; center_color(:) = "yellow"; end
    overlap = false(4,1);
    if mod(stimulus,11)==0; overlap(4) = true; end
    geometry(stimulus).T = table(assignment,center_color,overlap);
    response = 2+(1:4)'/10+sin(stimulus/3)/8+double(target)/3;
    data.meanAct(:,stimulus,:) = repmat(response,1,1,3);
end
data.meanSqAct = data.meanAct.^2+0.3;
data.nTrials = 8+mod(1:nStimuli,3);
if perSiteCounts; data.nTrials = repmat(data.nTrials,4,1)+(0:3)'; end
snr = struct('muSpont',ones(4,1),'muYellowEarly',4*ones(4,1), ...
    'muYellowLate',4.5*ones(4,1),'muPurpleEarly',3*ones(4,1), ...
    'muPurpleLate',3.5*ones(4,1));
end
