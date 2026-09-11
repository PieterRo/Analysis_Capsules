# Visual-response spread analysis

This folder starts an analysis of the spatial spread of visual responses in
Nilson V1 and V4. The current figures average the short response bins from
300 to 500 ms and project each RF/stimulus combination into the canonical
geometry of stimulus 38.

## Current entry scripts

- `Visual_Response_Spread_V1.m`: V1 channels 1:512.
- `Visual_Response_Spread_V4.m`: V4 channels 513:768, mapped explicitly to
  local rows 1:256 in `Tall_V4_lines_N.mat`.

Both scripts are standalone: they initialize paths with `config()`, load the
required data, reproduce the site-selection rules, average the response bins,
and make the projected activity figure.

## Current analysis choices

- Monkey: Nilson.
- Canonical stimulus: 38.
- Response window: 300-500 ms.
- Response estimate: unweighted mean of twenty equal-width 10-ms bins from
  `Resp_capsules_N_d12.mat`.
- RF/stimulus combinations on target, distractor, and background are shown
  (`onlyOnObjects = false`).
- Overlap exclusion and T-D statistics follow the existing attention code.

The site-inclusion mask is the union of:

- maximum color/window SNR greater than 0.7;
- T-D attention p-value below 0.05 with at least 20 matched trials;
- early color p-value below 0.05.

Activity is baseline-subtracted and divided by the established per-site
response scale:

```text
z = (activity - muSpont) / (max(muYellowEarly, muYellowLate,
                                muPurpleEarly, muPurpleLate) - muSpont)
```

Positive values are red/orange/yellow and negative values are blue. The
current plotting helper gives every non-zero value at least 0.08 opacity, so
small negative responses remain faintly visible. This is a visualization
choice to review before interpreting the spatial extent quantitatively.

## Validation on 2026-09-11

V1:

- 190 sites passed SNR > 0.7.
- 19 additional sites entered through T-D rescue.
- 0 additional sites entered through color rescue.
- 209/512 sites were included.
- 40,128 projected activity points were plotted.
- The median absolute difference between the short-bin mean and the broad
  300-500 ms bin was 0.00221151; the maximum was 0.444261.

V4:

- 137 sites passed SNR > 0.7.
- 8 additional sites entered through T-D rescue.
- 0 additional sites entered through color rescue.
- 145/256 sites were included.
- 27,840 projected activity points were plotted.
- Recomputed V4 SNR matched the stored V4 SNR to machine precision
  (maximum absolute difference 2.22e-15).
- The median absolute difference between the short-bin mean and the broad
  300-500 ms bin was 0.00162347; the maximum was 0.209941.

## Running the analysis

From the repository root in MATLAB:

```matlab
run(fullfile('analyses', 'spread', 'Visual_Response_Spread_V1.m'));
run(fullfile('analyses', 'spread', 'Visual_Response_Spread_V4.m'));
```

Set `P.saveFigure = true` in an entry script to save its `.fig` and `.png`
files under:

```text
fullfile(cfg.resultsDir, 'spread')
```

This resolves to the Dropbox results tree configured by `config_local.m`.
Generated figures and large result files should not be committed to Git.

## Suggested next steps

1. Define the quantitative spatial-spread metric before adding inferential
   statistics (for example area, radial profile, or along/perpendicular GC
   profiles).
2. Decide whether spread should be measured in pixels, visual degrees, GC
   units, or relative to each area's RF size.
3. Add the attention-effect spread while preserving the existing matched
   target/distractor and color-balancing definitions.
4. Decide whether values below `AlphaThresh` should become fully transparent
   for display-only figures.
5. Investigate the few large differences between the mean of short bins and
   the independently computed broad 300-500 ms bin.

## Continuing on another computer

The work is stored on branch `codex/spread-analysis`. On a computer where the
branch does not yet exist:

```bash
git fetch origin
git switch --track origin/codex/spread-analysis
```

If the local branch already exists:

```bash
git switch codex/spread-analysis
git pull
```

Machine-specific paths remain in the untracked `config_local.m`; ensure that
`cfg.dataRoot`, `cfg.extrasRoot`, `cfg.matRoot`, and `cfg.resultsRoot` point to
the corresponding Dropbox folders on that computer.
