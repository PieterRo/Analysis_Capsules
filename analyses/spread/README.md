# Visual-response spread analysis

This folder starts an analysis of the spatial spread of visual responses in
Nilson V1 and V4. The current figures average the short response bins from
300 to 500 ms and project each RF/stimulus combination into the canonical
geometry of stimulus 38.

## Current entry scripts

- `Visual_Response_Spread_V1_V4.m`: recommended entrypoint; runs both regions,
  combines them into one two-panel figure, and saves `.fig` and `.png` output.
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
- Projected samples are displayed as a density-normalized Gaussian average
  (`smoothSigmaPx = 15`) on a white background. The target capsule has a
  solid outline and the distractor capsule a dashed outline; neither object
  is filled. This smoothing changes only the visualization, not site
  inclusion, normalization, or the underlying response values.
- Neutral gray shows the projected coverage of all RFs in the selected
  cortical region. The colored response overlay still uses only sites that
  pass the established inclusion rule. Its display alpha is amplified by a
  factor of 3 to make positive responses more clearly red; response values
  themselves are unchanged. The smoothed display keeps positive values red
  over most of the range and shifts toward yellow only near the upper limit.
- The display includes 100 pixels of white padding around the original
  stimulus frame so peripheral RF coverage remains visible. The solid target
  and dashed distractor contours use a 4-pixel line width.
- The combined figure includes a shared color scale for the normalized
  response (`z`, clipped to -1 through 2). It is not an SNR scale.

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
run(fullfile('analyses', 'spread', 'Visual_Response_Spread_V1_V4.m'));
```

The individual panels can still be run separately:

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

## Attentional-modulation spread

`Attention_Modulation_Spread_V1_V4.m` is the recommended entrypoint for the
combined V1/V4 attentional-modulation figure. The individual scripts are
`Attention_Modulation_Spread_V1.m` and
`Attention_Modulation_Spread_V4.m`.

`Attention_Modulation_Spread_Activity_Inclusion_V1_V4.m` is a separate
comparison analysis using the broader visual-activity inclusion mask instead
of the attention-significance mask. This selects the same 209 V1 and 145 V4
sites as the visually driven activity figure and saves to a distinct
`*_activity_inclusion_*` output name. It does not overwrite the
attention-significant version.

The map uses the established 300-500 ms attention analysis and its fixed
`pValueTD < 0.05` site mask, as in the attention movie. Unlike that movie,
the spatial calculation includes target, distractor, and background RF
locations. It uses the 192 target-swap pairs `[1 6]`, `[2 5]`, `[3 8]`, and
`[4 7]` repeated in every block of eight. Only matched assignment pairs
(target/distractor, distractor/target, or background/background) with the
same RF-center color are included.

For each directed pair, target and distractor first and second response
moments are normalized per site using the established SNR response scale.
The contribution is weighted by the smaller trial count of the two stimuli.
After projection into the canonical stimulus-38 frame, the weighted moments
are Gaussian-averaged (`sigma = 15 px`) and converted locally to:

```text
d' = (muTarget - muDistractor) /
     sqrt(0.5 * (varTarget + varDistractor))
```

Red indicates target greater than distractor, blue indicates distractor
greater than target, and neutral gray shows all-RF sampling coverage. The
shared display range is `[-0.5 0.5]`; values near zero fade to the gray
coverage layer. Background values are direct comparisons from matched
background/background stimulus pairs, not copies of object-assigned d-prime.

Validation on 2026-09-12 included 129/512 V1 sites and 55/256 V4 sites. The
projected matched samples comprised 24,728 V1 and 13,512 V4 background
locations, in addition to equal target and distractor sample counts. The
supported field ranges were `[-0.2261 0.2266]` in V1 and
`[-0.3366 0.3391]` in V4, so no supported pixels were clipped by the shared
display range.

This is a descriptive smoothed map. Demonstrating statistically significant
modulation outside the object contours will require a separate spatial test
or confidence interval that accounts for smoothing and repeated site/stimulus
contributions.

Run the combined attention map from the repository root:

```matlab
run(fullfile('analyses', 'spread', ...
    'Attention_Modulation_Spread_V1_V4.m'));
```

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
