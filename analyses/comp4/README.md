# comp4 analysis

This folder contains the MATLAB code for analyses of the `comp4` stimulus set.
Large input and output files do not belong in this Git repository.

## Raw inputs

The current analysis includes Mr Nilson only. Its Dropbox-synchronized source
files are:

- `cfg.dataRoot/Mr Nilson/comp4/ObjAtt_Shapes_comp4_MUA_trials.mat`
- `cfg.dataRoot/Mr Nilson/comp4/ObjAtt_Shapes_comp4_normMUA.mat`

The original copies remain on VS03 under
`/Volumes/VS03-VandC-2/Object_attention/Shapes/monkeyN/`. Figaro is deliberately
out of scope for the current analysis.

## Generated files

Every comp4 script must initialize paths with `cfg = config()` and write only
below:

```matlab
comp4ResultsDir = fullfile(cfg.resultsRoot, 'Comp4');
```

On the current computer this resolves to the Dropbox-synchronized tree:

```text
/Users/pieter/Dropbox/Pieter/Text/Papers/Paolo/Object-based attention/Analysis/Comp4/
```

The current stimulus audit writes the canonical per-stimulus condition table
to:

```text
Comp4/metadata/comp4_stimulus_conditions.mat
Comp4/metadata/comp4_stimulus_conditions.csv
```

The table records the attended object and color for every stimulus. Later
scripts should load this table rather than infer quartet positions again.

## Stimulus geometry provenance

The bitmaps currently present in `cfg.extrasRoot/monkeyN/comp4` were regenerated
on 3 May 2022 and are not the images indexed by the neural recordings from
2 May 2022. Their geometry matches
`_logs/_old/RANDTAB_shapes_comp4_monkeyN.mat` exactly, but not the recording
session logs. The first eight columns of the neural `ALLMAT` match all 320 rows
of `ObjAtt_Shapes_comp4_monkeyN_20220502_B1.mat`; the `ALLCOORDS` structures in
blocks B1-B3 are identical. Analyses tied to the recordings must therefore use
the session `ALLCOORDS`, not the current bitmap geometry.

`Diagnose_Comp4_ALLCOORDS.m` reproduces this provenance check. The current
bitmaps have an exact mask overlap (IoU 1.0) with the 3 May regeneration for
all 80 quartets and do not match the 2 May session geometry.

`Plot_Comp4_Aligned_Cue_Locations.m` uses the original canonical object
contours and cue anchors from `temp_figs.mat`. It writes per-stimulus aligned
cue coordinates to `Comp4/metadata/comp4_cue_geometry.*` and a two-panel
monkey/crocodile cue plot to `Comp4/figures/`.

`Plot_Comp4_Crocodiles_Aligned_To_Monkey_Cue1.m` selects ten configurations
with the monkey in front at monkey cue 1, aligns the monkeys exactly, and plots
the associated crocodile contours underneath the canonical monkey. It
reconstructs the recorded shapes from the 2 May session `ALLCOORDS` and first
verifies all 320 stimulus definitions against the neural trial metadata.

`Plot_Comp4_Crocodile_Occupancy_Aligned_To_Monkey_Cue1.m` uses all 56 unique
recorded configurations (224 stimulus images) with monkey cue 1. It aligns the
monkeys, rasterizes the complete crocodile shapes, and plots crocodile pixel
occupancy outside the canonical monkey as a percentage of configurations.

`Plot_Comp4_Crocodile_Occupancy_Aligned_To_Monkey_Cue2.m` provides the same
analysis for the 24 unique configurations (96 stimulus images) with monkey
cue 2.

`Plot_Comp4_Monkey_Occupancy_Aligned_To_Crocodile_Cue1.m` reverses the
alignment: it uses all 56 unique configurations (224 stimulus images) with
crocodile cue 1 and plots monkey pixel occupancy outside the canonical
crocodile.

`Plot_Comp4_Monkey_Occupancy_Aligned_To_Crocodile_Cue2.m` provides the same
analysis for the 24 unique configurations (96 stimulus images) with crocodile
cue 2.

Do not copy the multi-gigabyte source files into the result tree. Read the
Dropbox data copies with `matfile` and save only compact quantities required to
reproduce figures and statistics.
