function S = IT_functional_inclusion_N(varargin)
%IT_FUNCTIONAL_INCLUSION_N Functional IT-site inclusion for Mr Nilson.
% Retain a site when it is visually responsive, attention-modulated, or
% color-modulated. Only sites lacking all three signals are excluded.

p = inputParser;
p.addParameter('Monkey','N',@(x) any(strcmpi(string(x),["N" "F"])));
p.addParameter('Config',struct(),@isstruct);
p.addParameter('VisualSNRThreshold', 0.7, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.parse(varargin{:});
opt = p.Results;

cfg = config();
if ~isempty(fieldnames(opt.Config)); cfg=opt.Config; end
suffix=upper(char(opt.Monkey));
attentionFile = fullfile(cfg.matDir, ...
    ['Attention_TargetSide_DistanceControl_IT_' suffix '.mat']);
if suffix == 'F' && ~isfile(attentionFile)
    attentionFile = fullfile(cfg.resultsDir, ...
        'Attention_TargetSide_DistanceControl_IT_F.mat');
end
colorDirectionFile = fullfile(cfg.resultsDir, ...
    ['Color_direction_tuning_IT_respSNR0_' suffix '.mat']);
colorDistanceFile = fullfile(cfg.resultsDir, ...
    ['Color_tuning_IT_nearest_curve_respSNR0_' suffix '.mat']);
assert(isfile(attentionFile), 'Missing %s.', attentionFile);
assert(isfile(colorDirectionFile), 'Missing %s.', colorDirectionFile);
assert(isfile(colorDistanceFile), 'Missing %s.', colorDistanceFile);

A = load(attentionFile, 'OUT');
Cdir = load(colorDirectionFile, 'OUT');
Cdist = load(colorDistanceFile, 'OUT');
A = A.OUT;
Cdir = Cdir.OUT;
Cdist = Cdist.OUT;

siteIndexGlobal = double(A.RFrange(:));
assert(~isempty(siteIndexGlobal) && ...
    isequal(double(Cdir.siteIndexGlobal(:)), siteIndexGlobal), ...
    'The all-site IT attention and color analyses are misaligned.');
assert(numel(Cdir.visualResponseScoreAllSites) == numel(siteIndexGlobal) && ...
    numel(Cdist.early.significantFDR) == numel(siteIndexGlobal), ...
    'The IT functional masks have inconsistent lengths.');

visual = double(Cdir.visualResponseScoreAllSites(:)) > ...
    opt.VisualSNRThreshold;
attentionDirection = logical(A.ctrlSig(:));
attentionDistance = logical(A.distGivenDirSig(:));
attention = attentionDirection | attentionDistance;
colorDirection = logical(Cdir.early.significantFDR(:)) | ...
    logical(Cdir.late.significantFDR(:));
colorDistance = logical(Cdist.early.significantFDR(:)) | ...
    logical(Cdist.late.significantFDR(:));
color = colorDirection | colorDistance;
included = visual | attention | color;

S = struct();
S.siteIndexGlobal = siteIndexGlobal;
S.visual = visual;
S.attention = attention;
S.attentionDirection = attentionDirection;
S.attentionDistance = attentionDistance;
S.color = color;
S.colorDirection = colorDirection;
S.colorDistance = colorDistance;
S.included = included;
S.excludedTripleNull = ~included;
S.visualSNRThreshold = opt.VisualSNRThreshold;
S.attentionDefinition = ['p<0.05 for direction after distance control or ' ...
    'distance after direction control in the late IT attention model.'];
S.colorDefinition = ['FDR q<0.05 for color direction or nearest-curve ' ...
    'color distance in either the 40-240 or 300-500 ms window.'];
S.inclusionDefinition = ['Visual-response SNR > threshold OR attention ' ...
    'effect OR color effect; exclude only triple-null sites.'];
S.sourceFiles = {attentionFile, colorDirectionFile, colorDistanceFile};

fprintf(['IT functional inclusion: %d/%d sites retained; visual %d, ' ...
    'attention %d, color %d; triple-null %d.\n'], ...
    nnz(included), numel(included), nnz(visual), nnz(attention), ...
    nnz(color), nnz(~included));
end
