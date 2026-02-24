Here is the cleaned, ready-to-paste README text:

---

# Object-Based Attention Analysis (V1)

This repository contains analysis pipelines for object-based attention experiments in V1. It includes receptive-field (RF) alignment to stimuli, color tuning quantification, time-resolved decoding (SVM-based), attention modulation analyses, and visualization utilities such as RF overlays and activity movies.

The codebase is organized to separate **top-level scientific analyses** from **reusable computational modules**.

---

## Repository Structure

* **analyses/** – Top-level scripts organized by scientific question. These are the entry points for running analyses.
* **core/** – Reusable functions (alignment, rendering, RF mapping, SNR computation, decoding utilities, metrics, stimulus pairing).
* **utils/** – Shared plotting and helper utilities.
* **archive/** – Legacy or superseded scripts kept for provenance.
* **data_mat/** – Generated/intermediate `.mat` artifacts (not raw data).
* **config.m** – Central configuration file for paths and global parameters.
* **RUN_PIPELINE.m** – Root entry script to initialize paths and execute selected analyses.

---

## Analyses Overview

### Latency (analyses/latency)

**Goal:** Determine onset timing of attention and color effects in V1.

Main scripts:

* Latency_analysis_SVM.m
* Latency_analysis_color.m
* Latency_analysis_color_integrated.m

Outputs:

* Latency figures
* Time-resolved performance summaries
* Optional saved result structures

---

### Decoding (analyses/decoding)

**Goal:** Quantify how accurately stimulus or condition information can be decoded from V1 activity over time.

Main scripts:

* Color_SVM.m
* Color_decoding_V1.m
* time_resolved_SVM.m

Outputs:

* Decoding accuracy curves
* Condition-wise performance comparisons
* Optional saved decoding result files

---

### PSTH / Color Preference (analyses/psth)

**Goal:** Compare temporal dynamics of preferred vs non-preferred color responses, including RF–stimulus distance constraints.

Main scripts:

* PSTH_colorPref_V1.m
* PSTH_colorPref_V1_grayDist.m

Outputs:

* PSTH figures
* Summary traces and statistics

---

### Line Stimuli & Attention (analyses/line_stimuli)

**Goal:** Characterize representation of line-task stimuli relative to RF geometry and attention condition.

Main scripts:

* Line_Stimuli.m
* Analyse_Line_Stimuli.m
* Attention_Line_Stimuli.m

Outputs:

* Stimulus–RF summary tables
* Attention modulation plots
* Intermediate `.mat` datasets

---

### RF Alignment & Overlay (analyses/rf_overlay)

**Goal:** Align RF locations to rendered stimuli and validate geometric transformations.

Main scripts:

* RF_on_stim_V1.m
* RF_on_stim.m
* AlignBitmapsRFs.m
* AllRFs_AlignBitmaps.m

Outputs:

* RF-on-stimulus overlay figures
* Alignment diagnostics

---

### Color Tuning (analyses/color_tuning)

**Goal:** Identify and quantify color-tuned V1 sites.

Main scripts:

* ColorTuning_Capsules.m
* color_tuning_cone.m

Outputs:

* Color tuning metrics
* Tuning classification figures
* Optional ColorTune result structures

---

### Stimulus Generation (analyses/stimulus_generation)

**Goal:** Generate and validate RANDTAB/ALLCOORDS stimulus sets.

Main scripts:

* Eval_RANDTAB.m
* makle_randatab_shapes_compl.m

Outputs:

* RANDTAB-related `.mat` files
* Validation diagnostics

---

### Exploratory (analyses/exploratory)

Ad hoc diagnostics and exploratory checks (pairing logic, SNR comparisons, inclusion criteria, RF inspection).

---

### Activity Movies (analyses/movie)

**Goal:** Visualize time-resolved activity projected onto stimulus geometry.

Main scripts:

* make_activity_movie.m
* make_activity_movie_wrapper_safe.m

Outputs:

* Activity movies (e.g., MP4 files)
* Associated projection plots

---

## How to Run

1. Open MATLAB in the repository root.

2. Run:

   RUN_PIPELINE

3. In RUN_PIPELINE.m, uncomment the analysis blocks you want to execute.

4. Re-run RUN_PIPELINE after modifying settings in config.m if needed.

Recommendation: Validate one analysis at a time when first setting up paths or updating data.

---

## Reproducibility Notes

* Avoid relying on variables from the base workspace.
* Save result structures with timestamps and configuration parameters.
* Prefer modifying config.m over editing analysis scripts directly.
* Archived scripts are kept for provenance but are not guaranteed to reproduce current results.

---
