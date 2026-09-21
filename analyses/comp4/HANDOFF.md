# Comp4 analysis handoff

This document summarizes the current Mr Nilson Comp4 analysis and the steps
needed to continue it on another computer.

## Current status

- The dataset contains 320 stimulus images forming 80 unique quartets.
- Each quartet contains two foreground-object conditions and two color
  conditions. The four images share the same spatial configuration.
- Only Mr Nilson is included. Figaro is deliberately out of scope.
- Analysis code is on GitHub `main` under `analyses/comp4/`.
- Large source data and generated results are synchronized through Dropbox.

The cue-condition counts are:

| Canonical object | Cue | Stimulus images | Unique configurations |
| --- | ---: | ---: | ---: |
| Monkey | 1 | 224 | 56 |
| Monkey | 2 | 96 | 24 |
| Crocodile | 1 | 224 | 56 |
| Crocodile | 2 | 96 | 24 |

The occupancy analyses use unique configurations. Counting all four quartet
members would give identical percentages while unnecessarily weighting each
geometry four times.

## Critical geometry provenance

The bitmaps currently stored under `cfg.extrasRoot/monkeyN/comp4` are not the
stimulus generation used for the neural recordings. They were regenerated on
3 May 2022 and overwrote the original bitmap filenames.

For analyses tied to the neural recordings, use this 2 May session file as the
geometry ground truth:

```text
cfg.extrasRoot/monkeyN/_logs/ObjAtt_Shapes_comp4_monkeyN_20220502_B1.mat
```

Key validation results:

- Neural trial `ALLMAT` matches the 2 May session generation for 320/320
  stimuli.
- Neural trial `ALLMAT` matches the 3 May regeneration for 0/320 stimuli.
- The `ALLCOORDS` structures in recording blocks B1, B2 and B3 are identical.
- Current bitmaps match the 3 May archived generation exactly (mask IoU 1.0
  for both objects in all 80 quartets).
- The two generations differ in all 320 rotation angles. Their median absolute
  angle difference is 34.38 degrees, and their median central-cue displacement
  is 73.15 pixels.

`Diagnose_Comp4_ALLCOORDS.m` reproduces this audit. Do not use the current
bitmaps for spatial analyses of the recorded neural responses.

## Required Dropbox inputs

The two large neural-data files are already in Dropbox:

```text
cfg.dataRoot/Mr Nilson/comp4/ObjAtt_Shapes_comp4_MUA_trials.mat
cfg.dataRoot/Mr Nilson/comp4/ObjAtt_Shapes_comp4_normMUA.mat
```

Their contents are:

| File | Variables |
| --- | --- |
| `ObjAtt_Shapes_comp4_MUA_trials.mat` | `ALLMAT [2079 x 11]`, `ALLMUA [1024 x 2079 x 700]`, `tb [1 x 700]` |
| `ObjAtt_Shapes_comp4_normMUA.mat` | `normMUA [1024 x 2079 x 700]`, `SNR [1024 x 1]`, `lats [1024 x 1]` |

Required geometry files are also already in Dropbox:

```text
cfg.extrasRoot/monkeyN/_logs/temp_figs.mat
cfg.extrasRoot/monkeyN/_logs/ObjAtt_Shapes_comp4_monkeyN_20220502_B1.mat
```

The Dropbox copies of the two large MAT files have the same byte sizes and
timestamps as their network originals and have been read successfully by
MATLAB. No additional network copy is required for the current analyses.

The RF files are already available under:

```text
cfg.extrasRoot/monkeyN/RFs/
```

These may become useful for later site- or RF-specific analyses.

## Results location

All generated results are written below:

```matlab
comp4ResultsDir = fullfile(cfg.resultsRoot, 'Comp4');
```

This contains:

```text
Comp4/figures/
Comp4/metadata/
```

Important compact metadata files include:

```text
Comp4/metadata/comp4_stimulus_conditions.csv
Comp4/metadata/comp4_stimulus_conditions.mat
Comp4/metadata/comp4_cue_geometry.csv
Comp4/metadata/comp4_cue_geometry.mat
```

## Analysis scripts

`Check_Comp4_Stimulus_Quartets.m`

Verifies the 320-image/80-quartet organization and writes the compact stimulus
condition table. Its bitmap checks describe the currently stored 3 May set;
the categorical quartet labels remain valid for the recorded set.

`Diagnose_Comp4_ALLCOORDS.m`

Checks the geometry provenance and demonstrates that the current bitmaps match
the 3 May regeneration rather than the 2 May recording session.

`Plot_Comp4_Aligned_Cue_Locations.m`

Creates canonical monkey and crocodile cue-location plots and writes
`comp4_cue_geometry.*`.

`Plot_Comp4_Crocodiles_Aligned_To_Monkey_Cue1.m`

Plots ten reconstructed crocodile outlines after aligning their front monkeys
to canonical monkey cue 1.

`Plot_Comp4_Crocodile_Occupancy_Aligned_To_Monkey_Cue1.m`

Aligns all 56 unique monkey-cue-1 configurations and plots crocodile occupancy
outside the canonical monkey. Maximum occupancy is 96.4%.

`Plot_Comp4_Crocodile_Occupancy_Aligned_To_Monkey_Cue2.m`

Aligns all 24 unique monkey-cue-2 configurations. Maximum crocodile occupancy
outside the monkey is 83.3%.

`Plot_Comp4_Monkey_Occupancy_Aligned_To_Crocodile_Cue1.m`

Aligns all 56 unique crocodile-cue-1 configurations. Maximum monkey occupancy
outside the crocodile is 100%.

`Plot_Comp4_Monkey_Occupancy_Aligned_To_Crocodile_Cue2.m`

Aligns all 24 unique crocodile-cue-2 configurations. Maximum monkey occupancy
outside the crocodile is 75.0%.

All four occupancy figures use the same fixed 0-100% color scale. Occupancy
inside the canonical aligned object is masked.

## Resume at home

1. Let Dropbox finish synchronizing before shutting down the current computer.
   The two raw Comp4 data files occupy approximately 21 GiB together.
2. On the home computer, update the repository:

   ```bash
   git switch main
   git pull origin main
   ```

3. Confirm that the home computer has its own `config_local.m`. This file is
   intentionally not tracked by Git. It should point to the home computer's
   Dropbox locations:

   ```matlab
   cfg.dataRoot = '/path/to/Dropbox/Pieter/data';
   cfg.extrasRoot = ['/path/to/Dropbox/Pieter/Text/Papers/Paolo/' ...
       'Object-based attention/Extras'];
   cfg.resultsRoot = ['/path/to/Dropbox/Pieter/Text/Papers/Paolo/' ...
       'Object-based attention/Analysis'];
   ```

4. Make the two large MAT files available offline in Dropbox before running
   neural-response analyses. Online-only placeholders may cause long delays.
5. Start MATLAB in the repository root and run a Comp4 script, for example:

   ```matlab
   run(fullfile('analyses', 'comp4', ...
       'Plot_Comp4_Crocodile_Occupancy_Aligned_To_Monkey_Cue1.m'));
   ```

Each occupancy script first verifies all 320 stimulus definitions against the
recorded trial `ALLMAT`. It will stop rather than silently use the wrong
stimulus generation.

## Logical next analysis

The current work establishes stimulus identities, cue identities and canonical
geometry. A logical next stage is to connect neural responses to these aligned
object coordinates. Before implementing that analysis, explicitly decide:

- Site inclusion criterion, including whether to use the stored `SNR` and what
  threshold to apply.
- Response and baseline time windows.
- Whether V1 and V4 are analyzed separately from the start.
- Whether responses are summarized per unique quartet before pooling.
- Whether RF position is represented relative to the attended object's cue,
  object boundary, or occupancy map.

Do not change these scientific definitions implicitly in plotting code.
