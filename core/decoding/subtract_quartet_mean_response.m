function [residualResponse, quartetMeanResponse, quartetTrialCount] = ...
    subtract_quartet_mean_response(response, trialQuartet, includeTrial, nQuartets)
%SUBTRACT_QUARTET_MEAN_RESPONSE Center each site's response within quartet.

if nargin < 4
    nQuartets = 96;
end

trialQuartet = double(trialQuartet(:));
includeTrial = logical(includeTrial(:));
assert(size(response, 2) == numel(trialQuartet), ...
    'Response columns must match trialQuartet.');
assert(numel(includeTrial) == numel(trialQuartet), ...
    'includeTrial must match trialQuartet.');

nSites = size(response, 1);
residualResponse = nan(size(response));
quartetMeanResponse = nan(nSites, nQuartets);
quartetTrialCount = zeros(nQuartets, 1);

for quartet = 1:nQuartets
    useTrial = includeTrial & trialQuartet == quartet;
    quartetTrialCount(quartet) = nnz(useTrial);
    if ~any(useTrial)
        continue;
    end

    quartetMeanResponse(:, quartet) = ...
        mean(response(:, useTrial), 2, 'omitnan');
    residualResponse(:, useTrial) = ...
        response(:, useTrial) - quartetMeanResponse(:, quartet);
end

end
