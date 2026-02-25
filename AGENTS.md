# AGENTS.md — Project operating rules for Codex (MATLAB)

This repository contains MATLAB analysis code for visual neuroscience experiments (stimulus alignment, receptive-field overlays, PSTHs/SNR metrics, attention modulation analyses, and figure generation). The code is used on multiple machines and must remain reproducible and backwards-compatible.

## 1. Primary goals
- Make small, reviewable changes that preserve scientific meaning.
- Prefer robustness and clarity over cleverness.
- Keep outputs stable unless the request explicitly changes analysis logic.

## 2. Golden rules (do not violate)
1) Do not change scientific definitions (normalization, SNR, windowing, inclusion criteria, stats) unless explicitly requested.
2) Do not change coordinate systems or stimulus geometry alignment unless explicitly requested.
   - In particular: treat ALLCOORDS and the per-stimulus affine/homogeneous transforms as ground truth for alignment to a canonical frame.
3) Do not rename or move entry scripts/functions without explicit instruction.
4) Do not introduce new external dependencies/toolboxes unless explicitly requested.
5) Never delete code or data-loading logic unless you can prove it is unused AND the user asked for cleanup.
6) Avoid broad rewrites. Keep diffs minimal and localized.

## 3. Repository structure (assumptions)
- analyses/ : main analysis scripts and pipelines
- plotting/ (or equivalent): plotting utilities
- stimuli/  (or equivalent): stimulus geometry helpers, templates
- data_mat/ and results/ may exist but large data should not be added to git

If paths differ, ask or infer cautiously from the repository tree.

## 4. MATLAB coding conventions
- Prefer functions over scripts when editing, but do not convert scripts to functions unless requested.
- Keep function signatures stable.
- Use clear variable names for derived quantities (e.g., snrNorm, pValueTD).
- Add comments only where they clarify non-obvious intent.
- Use `fullfile` for paths; avoid hard-coded absolute paths.
- Keep compatibility with the user’s current MATLAB version (do not require newest features).

## 5. Data and reproducibility policy
- Assume large data files are stored outside git (e.g., Dropbox) and referenced by paths.
- Do not add `.mat` datasets or large binaries to git.
- When adding ignores, ignore MATLAB build artifacts (slprj/, codegen/) and editor/OS temp files.
- Do not change default RNG behavior unless explicitly requested.

## 6. Safety checks before proposing changes
When modifying analysis code:
- Identify the entrypoint(s) and the minimal set of files to change.
- Explain what will change and why in 2–5 bullet points.
- Provide an explicit validation plan (what to run, what outputs should match).

## 7. Validation guidance (preferred)
If you need to propose commands:
- Prefer “dry” checks first (syntax, existence of fields, dimensions).
- If running MATLAB is requested/allowed, use a non-interactive run such as:
  - `matlab -batch "run('path/to/script.m')"`
- Always suggest checking key outputs:
  - number of included sites/trials
  - main summary metrics (e.g., meanAct, SNRnorm)
  - one representative figure/plot for sanity

## 8. How to work with git in this repo
- Always work on a feature branch named `codex/<topic>`.
- Keep commits focused: one logical change per commit.
- Do not commit generated results or large outputs.
- After changes, show the `git diff`-level summary and list modified files.

## 9. Common pitfalls (avoid)
- Accidental changes to selection criteria (e.g., “minDist” logic, overlap exclusion).
- Changing time windows or bins without explicit request.
- Implicitly changing statistics (e.g., p-value thresholds, multiple comparisons handling).
- Changing plotting code in a way that alters scientific interpretation (scales, normalization).

## 10. What I need from the user if unclear
Ask only if necessary; otherwise make a conservative choice and state it.
Examples:
- Which script is the canonical entrypoint for this analysis?
- Which MATLAB version/toolboxes are available?
- Which outputs should be treated as “golden” for regression testing?
