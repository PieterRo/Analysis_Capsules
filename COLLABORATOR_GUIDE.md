# Collaborator Guide: Capsule Line-Stimulus Analysis

This guide is for a new collaborator who wants to understand the MATLAB
analysis step by step and eventually run the same kind of analysis for the
second monkey, Figaro.

The safest approach is to first reproduce the existing Nilson pipeline, then
move one analysis block at a time to Figaro. Do not change scientific
definitions such as time windows, SNR thresholds, RF/stimulus alignment,
trial inclusion, or normalization unless this is an explicit analysis decision.

## 1. Repository Map

- `RUN_PIPELINE.m` is the root startup script. It adds paths, loads `config.m`,
  and runs selected analysis scripts.
- `config.m` defines shared paths and requires a local, untracked
  `config_local.m`.
- `analyses/line_stimuli/` contains the main line/capsule task analyses.
- `analyses/color_tuning/`, `analyses/psth/`, `analyses/decoding/`, and
  `analyses/movie/` contain downstream analyses and figures.
- `core/` contains reusable functions for RF/stimulus geometry, response
  averaging, SNR, color tuning, decoding, rendering, and plotting.
- `data_mat/` and `results/` are generated/cached output locations and are
  ignored by git.

## 2. Local Setup

Start MATLAB in the repository root:

```matlab
cd('/path/to/Analysis_c')
RUN_PIPELINE
```

On Pieter's Mac, MATLAB R2020b is available through:

```sh
arch -x86_64 /Applications/MATLAB_R2020b.app/bin/matlab -batch "RUN_PIPELINE"
```

The machine-specific `config_local.m` must provide:

```matlab
cfg.dataRoot = '/path/to/data';
cfg.extrasRoot = '/path/to/Object-based attention/Extras';
cfg.resultsRoot = '/path/to/Object-based attention/Analysis';
```

`config.m` then derives:

- `cfg.logsDir`: stimulus logs, including `ObjAtt_lines_monkeyN_20220201_B1.mat`
  and `RTAB384.mat`
- `cfg.stimDir`: stimulus bitmap folder
- `cfg.matDir`: cached `.mat` analysis files
- `cfg.resultsDir`: figures, movies, and larger result files

## 3. Important Inputs

Raw response data:

- `cfg.dataRoot/Mr Nilson/ObjAtt_lines_normMUA.mat`
- `cfg.dataRoot/Mr Nilson/ObjAtt_lines_MUA_trials.mat`
- `cfg.dataRoot/Figaro/ObjAtt_lines_normMUA.mat`
- `cfg.dataRoot/Figaro/ObjAtt_lines_MUA_trials.mat`

RF fits:

- `THINGS_RF1s_N.mat` for Nilson
- `THINGS_RF1s_F.mat` for Figaro

Stimulus geometry:

- `ObjAtt_lines_monkeyN_20220201_B1.mat` provides `ALLCOORDS`
- `RTAB384.mat` provides the stimulus table for the 384 line-task stimuli

Generated response summaries expected by most scripts:

- `SNR_capsules_N_d12.mat`, `Resp_capsules_N_d12.mat`
- `SNR_capsules_F_d12.mat`, `Resp_capsules_F_d12.mat`

At the time this guide was written, the Nilson summaries were present in the
external `cfg.matDir`; the Figaro raw data and RF file were present, but the
Figaro `SNR_capsules_F_d12.mat` and `Resp_capsules_F_d12.mat` summaries were
not found and likely need to be built.

## 4. Suggested Reading Order

1. Read `config.m` and local `config_local.m`.
   Confirm where data, stimulus logs, cached `.mat` files, and results live.

2. Read `RUN_PIPELINE.m`.
   It currently runs the V1 Nilson line-stimulus scripts by default:
   `Line_Stimuli.m` and `Analyse_Line_Stimuli.m`.

3. Read the geometry builders:
   `analyses/line_stimuli/Line_Stimuli.m`,
   `Line_Stimuli_V4.m`, and `Line_Stimuli_IT.m`.
   These load `ALLCOORDS`, `RTAB384`, and RF centers, then build `Tall_*`
   geometry tables.

4. Read the geometry helpers:
   `core/rf/rf_table_target_distractor.m`,
   `core/metrics/add_arm_projection_metrics.m`,
   `core/metrics/add_polar_about_s.m`,
   `core/metrics/add_arc_about_s_edges.m`, and
   `core/metrics/add_GC_normalization.m`.

5. Read response averaging and exclusion-aware loading:
   `core/decoding/avg_byStim.m`,
   `core/decoding/load_capsules_day_split.m`,
   `core/decoding/load_capsules_struct_exclusion_aware.m`, and
   `core/decoding/site_session_exclusions.m`.

6. Read the first-pass analyses:
   `Analyse_Line_Stimuli.m` for the historical V1 workflow,
   `Analyse_Line_Stimuli_V4.m`, and `Analyse_Line_Stimuli_IT.m`.

7. Read downstream analyses only after the geometry and response summaries
   are clear: color tuning, PSTH, decoding, attention modulation, and movies.

## 5. Nilson Baseline Run

Before changing anything for Figaro, verify that the existing Nilson pipeline
runs:

```matlab
cd('/path/to/Analysis_c')
set(0, 'DefaultFigureVisible', 'off')  % optional for batch/smoke tests
RUN_PIPELINE
```

Useful sanity checks from the run:

- `ALLCOORDS` should contain 384 stimuli.
- `RTAB384` should be `384 x 8`.
- `Tall_V1` should contain 384 entries.
- `SNR_capsules_N_d12.mat` should load into `R.meanAct` with dimensions
  `[1024 x 384 x 3]`.
- Check one RF/stimulus overlay visually before trusting downstream metrics.

## 6. Building Figaro Response Summaries

Most Figaro-ready scripts expect `SNR_capsules_F_d12.mat` and
`Resp_capsules_F_d12.mat`. If they are missing, build them from Figaro's raw
trial files using the same helper and windows used in the Nilson workflow:

```matlab
cfg = config();
m1 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_normMUA.mat'));
m2 = matfile(fullfile(cfg.dataRoot, 'Figaro', 'ObjAtt_lines_MUA_trials.mat'));
tb = m2.tb;

timeWindowsSNR = [-200 0; 40 240; 300 500];
[meanAct, meanSqAct, nTrials, stimList] = avg_byStim( ...
    m1, m2, timeWindowsSNR, 'days', [1 2]);

R = struct();
R.meanAct = meanAct;
R.meanSqAct = meanSqAct;
R.nTrials = nTrials;
R.stimList = stimList;
R.timeWindows = timeWindowsSNR;
R.tb = tb;
R.file_m1 = m1;
R.file_m2 = m2;
save(fullfile(cfg.matDir, 'SNR_capsules_F_d12.mat'), 'R', '-v7.3');

edges = -200:10:500;
timeWindowsResp = [edges(1:end-1).' edges(2:end).'];
[meanAct, meanSqAct, nTrials, stimList] = avg_byStim( ...
    m1, m2, timeWindowsResp, 'days', [1 2]);

R = struct();
R.meanAct = meanAct;
R.meanSqAct = meanSqAct;
R.nTrials = nTrials;
R.stimList = stimList;
R.timeWindows = timeWindowsResp;
R.tb = tb;
R.file_m1 = m1;
R.file_m2 = m2;
save(fullfile(cfg.matDir, 'Resp_capsules_F_d12.mat'), 'R', '-v7.3');
```

After building these files, check:

```matlab
S = load(fullfile(cfg.matDir, 'SNR_capsules_F_d12.mat'), 'R');
size(S.R.meanAct)
size(S.R.nTrials)
```

Expected response dimensions are 1024 channels by 384 stimuli by 3 windows
for the SNR summary.

## 7. Building Figaro Geometry

The V4 and IT geometry scripts already have a monkey switch:

```matlab
Monkey = 2;        % Figaro
ForceRebuild = true;
```

Run:

```matlab
run(fullfile(cfg.repoRoot, 'analyses', 'line_stimuli', 'Line_Stimuli_V4.m'))
run(fullfile(cfg.repoRoot, 'analyses', 'line_stimuli', 'Line_Stimuli_IT.m'))
```

These should produce:

- `Tall_V4_lines_F.mat`
- `Tall_IT_lines_F.mat`

For V1, `Line_Stimuli.m` can load Figaro when `Monkey = 2`, but its full
`Tall_V1` save block is currently commented and historically Nilson-oriented.
If building Figaro V1 geometry, do not overwrite `Tall_V1_lines_N.mat`; save
the result explicitly as `Tall_V1_lines_F.mat`.

## 8. Running Figaro Analyses

After the Figaro response summaries and geometry files exist, start with the
scripts that already expose `Monkey = 1/2` near the top:

- `analyses/line_stimuli/Analyse_Line_Stimuli_V4.m`
- `analyses/line_stimuli/Analyse_Line_Stimuli_IT.m`
- `analyses/color_tuning/ColorTuning_Capsules_V4.m`
- `analyses/color_tuning/ColorTuning_Capsules_IT.m`
- `analyses/decoding/Color_decoding_V4.m`
- `analyses/decoding/Color_decoding_IT_RFm.m`
- `analyses/movie/Make_V4_activity_movie.m`
- `analyses/movie/Make_IT_activity_movie.m`

Set `Monkey = 2` and run one script at a time. Check each output before
continuing. For V1, treat `Analyse_Line_Stimuli.m`,
`ColorTuning_Capsules.m`, `PSTH_colorPref_V1_grayDist.m`, and
`Color_decoding_V1.m` as Nilson-first templates until they are made
suffix-aware.

## 9. Validation Checklist

For each new Figaro step, record these checks:

- Number of included trials per stimulus.
- `size(R.meanAct)`, `size(R.meanSqAct)`, and `size(R.nTrials)`.
- Number of RF sites per area: V1, V4, IT.
- Number of sites assigned to target, distractor, background, and overlap.
- Main summary counts such as SNR-selected sites and attention-rescue sites.
- One representative RF-on-stimulus overlay.
- One representative activity/PSTH/decoding plot, with axes and time windows
  matching the Nilson analysis.

If any count differs unexpectedly, stop there and inspect inputs before moving
downstream.

## 10. Practical Rules

- Keep Nilson and Figaro outputs suffix-separated: `_N` and `_F`.
- Do not commit `.mat`, movie, or generated result files.
- Do not change `ALLCOORDS` alignment or affine/homogeneous transform logic.
- Do not change time windows, thresholds, or inclusion criteria silently.
- Prefer small commits: one analysis step or one documentation update at a
  time.
- When in doubt, compare against the Nilson run first.
