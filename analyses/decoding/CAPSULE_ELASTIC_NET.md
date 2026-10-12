# Capsule attention decoders

This directory contains the elastic-net attention decoders and the distance-dependent latency analyses for Nilson and Figaro. The code uses the existing capsule geometry, normalization and site-selection rules. Raw recordings, response summaries, geometry caches and fitted results are external inputs and are not included in the repository.

The current Fig. 1 comparison uses a common **full-data fitting rule** across the six monkey/area combinations. The earlier cross-validated decoders are also provided. These are distinct analyses: full-data accuracy is in-sample and must not be described as held-out performance.

## Entry points

| File | Purpose |
| --- | --- |
| `Attention_decoding_elastic_net_CV_V1_V4.m` | Existing Nilson V1/V4 decoder with repeated, quartet-balanced cross-validation. |
| `Attention_decoding_elastic_net_CV_IT.m` | Nilson/Figaro IT decoder combining target-side direction and nearest-capsule RF proximity, with one coefficient per site. |
| `IT_functional_inclusion_N.m` | Shared IT inclusion mask for either monkey: visual OR attention OR colour eligibility. |
| `Attention_decoding_Figaro_capsules.m` | Builds the historical Figaro RF-aligned candidate pool and d-prime reference readout used by the lower-area comparison. |
| `Compare_Figaro_Capsule_ElasticNet.m` | Full-data Figaro elastic-net comparison; V1/V4 are reused in the current figure. Its historical IT selection is superseded by the matched IT pipeline below. |
| `Compare_Figaro_IT_Nilson_Method.m` | Prepares Figaro IT geometry and functional masks and applies the same IT model used for Nilson. |
| `Compare_Capsule_IT_FullData.m` | Refits saved IT feature matrices with the common full-data rule. |
| `../../tmp/it_full_data/refit_lower.m` | Existing Nilson V1/V4 full-data refit script, with portable repository-path discovery. |
| `Compare_Figaro_IT_Matched_Distance.m` | Shared fixed-weight site bootstrap and latency fit; accepts either monkey and V1/V4/IT via `Monkey` and `Area`. |
| `Decode_Nilson_Capsule_Attention_By_Distance.m` | Earlier distance analysis using saved CV-selected Nilson coefficients. |
| `Decode_Figaro_Capsule_Attention_By_Distance.m` | Figaro distance analysis, including its original full-data V1/V4 fits; choose `ElasticNet=true` and `BootstrapUnit='site'`. |

The colour-direction and colour-proximity analyses under `analyses/color_tuning/`, and the target-side analyses under `analyses/line_stimuli/`, provide the IT functional masks. The shared attention-summary function supports Nilson's 384 stimuli and Figaro's 216 stimuli.

## Setup and external inputs

Use MATLAB R2020b or a compatible later release with the Statistics and Machine Learning Toolbox (`lassoglm`). MATLAB's dependency audit also identifies the Image Processing Toolbox for the stimulus-rendering helpers.

Configure machine-specific paths in the ignored `config_local.m`, following `config.m` and the repository setup instructions. Run from the repository root:

```matlab
addpath(pwd);
addpath(genpath(fullfile(pwd, 'core')));
addpath(genpath(fullfile(pwd, 'analyses')));
cfg = config();
```

The decoder inputs include:

- Trial recordings `ObjAtt_lines_normMUA.mat` and trial metadata/time base `ObjAtt_lines_MUA_trials.mat` under `cfg.dataRoot/Mr Nilson` and `cfg.dataRoot/Figaro`.
- Geometry caches `Tall_V1_lines_*.mat`, `Tall_V4_lines_*.mat`, `Tall_IT_lines_*.mat`, containing the existing `Tall_*`, `RFrange` and stimulus-coordinate fields.
- Response summaries `SNR_capsules_N_d12.mat` and `SNR_capsules_F_d12.mat`.
- IT attention and colour-analysis caches consumed by `IT_functional_inclusion_N`.

Figaro geometry can be stored in the existing sibling `Figaro files` folder. `Compare_Figaro_IT_Nilson_Method` copies its native stimulus table to the shared `RTAB384` field without changing the values. It writes isolated working geometry and results under its output directory.

The loaders retain the existing site/session exclusions. Cached candidate IDs, trial IDs, response normalization, quartet centering and geometry alignment are retained when refitting on saved features. Any user-maintained exclusion configuration must accompany the external data setup.

## Reproducing the full-data comparison

First prepare the existing feature sources. If they are already available, reuse them rather than rebuilding the recordings:

```matlab
% Nilson lower-area feature source, and optional earlier CV comparison.
Attention_decoding_sparse_linear_V1_V4('MakeFigure', false);
Attention_decoding_elastic_net_CV_V1_V4('MakeFigure', false);

% Nilson IT functional masks and CV feature/alignment source.
Attention_TargetSide_DistanceControl_IT(struct('makeSummaryFigure', false));
Color_direction_tuning_IT('ResponseSNRThreshold', 0, 'MakeFigure', false);
Color_tuning_IT_nearest_curve('ResponseSNRThreshold', 0, 'MakeFigure', false);
nITsource = Attention_decoding_elastic_net_CV_IT( ...
    'Monkey', 'N', 'ResponseSNRThreshold', 0.7, 'MakeFigure', false);

% Figaro lower areas and matched IT source.
Attention_decoding_Figaro_capsules('MakeFigure', false);
Compare_Figaro_Capsule_ElasticNet;
fITsource = Compare_Figaro_IT_Nilson_Method('RunFixedPool', false);
```

Nilson's visual responsiveness threshold remains 0.7 and Figaro's remains 0.5. IT sites can also qualify through attention or colour responsiveness. IT alignment includes the full-data direction/proximity model; sites without an RF centre can retain direction-only alignment. Functional selection is estimated on the full data, including for the CV comparison, and is not nested within cross-validation.

Apply the common fitting rule:

```matlab
paperOut = fullfile(pwd, 'tmp', 'it_full_data');
run(fullfile(pwd, 'tmp', 'it_full_data', 'refit_lower.m'));
nLower = OUT;
nIT = Compare_Capsule_IT_FullData(nITsource.resultFile, 'OutputDir', paperOut);
fIT = Compare_Capsule_IT_FullData(fITsource.matched.resultFile, 'OutputDir', paperOut);
```

The rule is binomial `lassoglm`, alpha 0.25, 12 lambda values and `LambdaRatio=0.05`, using the first returned (smallest) lambda. Predictor standardization is MATLAB's default. There are no class-balancing weights or held-out folds in this refit. It changes penalty selection relative to the CV model; it is not merely a refit at the previously CV-selected penalty. Figaro V1/V4 already use this rule.

The existing `tmp/it_full_data/setup.m`, `run_N.m`, `run_F.m` and `run_lower_N.m` scripts retain the filenames used for the saved figure comparison. If a source has a day-specific filename, use the returned `resultFile` as above instead of assuming the cache name in those scripts. They preserve their existing entrypoint names; generated MAT files remain ignored by Git.

## Distance groups, latency and uncertainty

```matlab
for area = {'V1', 'V4'}
    Compare_Figaro_IT_Matched_Distance(nLower.resultFile, ...
        'Monkey', 'N', 'Area', area{1}, ...
        'Config', cfg, 'OutputDir', paperOut, 'BootstrapIterations', 10000);
end
Compare_Figaro_IT_Matched_Distance(nIT.resultFile, ...
    'Monkey', 'N', 'Config', cfg, ...
    'OutputDir', paperOut, 'BootstrapIterations', 10000);

% Use the matched Figaro geometry for its updated IT fit.
fCfg = cfg;
fCfg.matDir = fullfile(pwd, 'tmp', 'figaro_nilson_matched', 'data_mat');
Compare_Figaro_IT_Matched_Distance(fIT.resultFile, ...
    'Monkey', 'F', 'Config', fCfg, ...
    'OutputDir', paperOut, 'BootstrapIterations', 10000);

% Figaro V1/V4 use this historical full-data distance output.
fLowerDistance = Decode_Figaro_Capsule_Attention_By_Distance( ...
    'ElasticNet', true, 'BootstrapUnit', 'site', ...
    'OutputDir', paperOut, 'BootstrapIterations', 10000);
```

Use `fLowerDistance.result{1}` and `{2}` for Figaro V1/V4. The historical IT result in `{3}` is not the updated matched IT result; use the separate `Compare_Figaro_IT_Matched_Distance` output for that panel.

Six rank-based distance groups are formed from non-overlapping site-quartet combinations using `along_GC`, in capsule-width units. Native 1-ms responses are averaged into 10-ms decoder bins. The shared-slope sigmoid fit uses 80-450 ms and subtracts the -150 to 0 ms baseline. The regression relates fitted t50 to median cue distance per group.

A site bootstrap carries a selected site's quartets together while keeping learned coefficients, alignments, distance groups and median distances fixed. It does not refit the decoder or repeat site selection. The two-sided slope p-value is twice the smaller bootstrap tail proportion, capped at one. A zero observed tail should be displayed as below the attainable resolution, not as p=0. Count finite successful bootstrap fits when reporting that resolution.

The number of decoder candidates, nonzero coefficients, spatially eligible sites and site-quartet combinations are different quantities. The MAT outputs store the IDs and group maps needed to distinguish them. IT sites without valid non-overlap RF distance can contribute to decoder training but not to the distance regression.

## Validation

```matlab
addpath(fullfile(pwd, 'tests'));
test_attention_modulation_stimulus_count;
checkcode('analyses/decoding/Attention_decoding_elastic_net_CV_IT.m');
checkcode('analyses/decoding/Compare_Capsule_IT_FullData.m');
checkcode('analyses/decoding/Compare_Figaro_IT_Matched_Distance.m');
```

The stimulus-count check compares 216-stimulus summaries with equivalent 384-stimulus inputs padded with zero-trial stimuli, for both vector and per-site trial counts. Publishing validation also compares the updated 384-stimulus summary with the previously tracked implementation and checks full-data weights/scores against saved decoder artifacts. No raw datasets or generated results are committed.
