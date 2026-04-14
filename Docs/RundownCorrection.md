# Run-down Correction Procedure

## Problem

During repeated activations of skinned cardiac muscle, progressive force loss (run-down) occurs due to cumulative damage and ATP depletion.  A reference recording `F_ref` and a later run-down recording `F_rd` show a systematic force deficit that varies with sarcomere length (SL) and time within the activation.  The goal is to reconstruct what `F_rd` would have looked like without run-down, and to correct `F_ref` for its own within-recording decay.

## Data

| Recording | Dataset | Description |
|-----------|---------|-------------|
| F_ref     | data(2) | 8 mM ATP, active — first (reference) activation |
| F_rd      | data(4) | 8 mM ATP, active — repeat activation (impaired) |

Both recordings follow the same staircase protocol (SL = 2.0 um baseline, steps to 2.2 um and back) over a ~10 s window (wall-clock t = 68 -- 78 s).  Force is measured in kPa.

Two stable isometric zones at SL = 2.0 um are used for calibration:

- **Zone 1**: t = 71.5 -- 72.4 s  (reference time point)
- **Zone 2**: t = 77.4 -- 78.4 s

## Measured quantities

| Quantity | Value | Source |
|----------|-------|--------|
| F_ref at zone 1 | 68.4 kPa | Linear fit to zone 1+2 points of data(2) |
| F_rd at zone 1 | 56.4 kPa | Linear fit to zone 1+2 points of data(4) |
| slope_ref | -0.497 kPa/s | Within-recording force decay of F_ref |
| slope_rd | -0.467 kPa/s | Within-recording force decay of F_rd |
| Force gap at zone 1 | 12.0 kPa | F_ref - F_rd at zone 1 midpoint |

## T0: hypothetical decay distance

The key derived quantity is **T0** -- the hypothetical activated time it would take for force to decay from the reference level to the run-down level at a constant rate:

    T0 = (F_ref(zone1) - F_rd(zone1)) / |slope_ref| = 12.0 / 0.497 = 24.1 s

T0 is *not* wall-clock time between recordings.  It represents **equivalent activated time**: the cumulative hot-activated exposure (muscle under calcium) that separates the two recordings.  Inactive recovery periods between activations do not contribute to T0.

T0 establishes a common time axis for the decay model: zone 1 of F_ref is at T = 0, zone 1 of F_rd is at T = T0.

![T0 concept and zoomed correction detail](figures/fig47_zoom_and_concept.png)

*Panel A: Zoomed timecourse (t = 71–79.5 s) comparing F_ref correction methods. Zone 1 (t = 71.5–72.4 s) and Zone 2 (t = 77.4–78.4 s) are shaded. Slope correction (additive, SL-independent) and M1 model correction (multiplicative, SL-dependent) diverge at SL = 2.2 during the stretch step.  Panel B: Decay timeline concept. The dashed line shows force decay from zone 1 of F_ref (T = 0) at rate slope_ref. Zone 1 of F_rd falls at T = T0 = 24.1 s. Correction arrows illustrate the gap that the model must close.*

## Step 1: Spatial correction model (M1)

A linear correction surface is fitted from paired (F_rd, F_ref) data over the valid window (F > 3 kPa, t = 68--78 s).  Both signals are first slope-detrended to their zone-1 values before fitting, to remove within-recording trends.

**Model M1** (SL-linear, 2 parameters):

    F_corrected = F_rd * (a + b * (SL - 1.0))

    a = 1.214,  b = -0.746

    R^2 = 0.977,  RMSE = 1.61 kPa

The correction ratio `a + b*(SL - 1.0)` depends only on sarcomere length:

| SL (um) | L = SL/2 | Correction ratio |
|---------|----------|-----------------|
| 1.6     | 0.8      | 1.363 |
| 1.8     | 0.9      | 1.289 |
| 2.0     | 1.0      | 1.214 |
| 2.2     | 1.1      | 1.140 |

Higher correction at shorter SL, lower at longer SL -- consistent with the observation that run-down impairs force more at shorter sarcomere lengths.

A more complex model (M3: +F and F x SL interaction, 4 params, R^2 = 0.980) is available in the code but M1 is sufficient for practical data processing.

## Step 2: Model-based F_ref correction

The reference recording also decays within its own activation window.  Rather than correcting with a simple constant slope (-0.497 kPa/s), we reuse the spatial model M1 scaled by the fraction of T0:

    frac(t) = (t - t_zone_ref) / T0

    r_model(t) = a + b * (SL(t) - 1.0)      [same M1 correction ratio]

    F_ref_corrected(t) = F_ref(t) * (1 + (r_model(t) - 1) * frac(t))

At zone 1 (frac = 0): no correction (identity).  
At zone 2 (frac = 0.247): correction of ~5.3% at SL = 2.0.  
At T0 distance (frac = 1): full model correction (~21.4%).

### Comparison: slope vs model at zone 2 (SL = 2.0)

| Method | Correction | Effect on 68 kPa |
|--------|-----------|-------------------|
| None (original) | 0 | 68.4 kPa -> ~65.4 kPa at zone 2 |
| Slope only | +2.96 kPa (additive) | 65.4 + 3.0 = 68.4 kPa |
| M1 model | x 1.053 (multiplicative) | 65.4 x 1.053 = 68.9 kPa |

The model-based correction is **multiplicative and SL-dependent**: at SL = 2.2 the correction is smaller (~4.0%), at SL = 1.8 it is larger (~6.5%).  The slope correction is additive and constant across all SL levels, which is a less accurate approximation of the underlying decay process.

![Correction comparison overview](figures/fig46_correction_comparison.png)

*Top: F_ref with no correction, slope-only correction, and M1/M3 model-based corrections.  Bottom: F_rd with no correction, M1 static correction, M1 with calculated slope rate, and M1 with freely fitted rate.  Model-based F_ref correction is multiplicative and SL-dependent; slope correction is additive and constant.*

## Step 3: F_rd correction (run-down compensation)

To reconstruct what F_rd would have been without run-down:

    F_rd_corrected(t) = F_rd(t) * (a + b * (SL(t) - 1.0))

This is the standard M1 static correction applied to the (slope-detrended) run-down signal.  The result should approximate F_ref at zone 1 level.

Time-extended variants are also computed (using within-recording decay rate from the slope difference), but for practical data processing the static M1 correction is sufficient -- the within-window time effect is small (~1% over 10 s) compared to the spatial correction (~21%).

## Output files

All saved as tab-delimited `[time, SL, force]` in the data directory:

| File | Description |
|------|-------------|
| `Corr_Fref_original.txt` | F_ref as recorded (no correction) |
| `Corr_Fref_slopeOnly.txt` | F_ref with linear slope correction to zone 1 |
| `Corr_Fref_M1_model.txt` | F_ref with M1 model-based correction to zone 1 |
| `Corr_Fref_M3_model.txt` | F_ref with M3 model-based correction to zone 1 |
| `Corr_Frd_original.txt` | F_rd as recorded (no correction) |
| `Corr_Frd_M1_static.txt` | F_rd corrected by M1 static model |
| `Corr_Frd_M1t_slopeRate.txt` | F_rd corrected by M1 + time (rate from slope) |
| `Corr_Frd_M1t_fittedRate.txt` | F_rd corrected by M1 + time (rate freely fitted) |
| `Corr_Frd_M3_static.txt` | F_rd corrected by M3 static model |

## How to apply to new data

For a new run-down recording with a known reference:

1. Fit the within-recording slope at two SL = 2.0 zones
2. Compute T0 = force_gap / |slope|
3. Apply M1 correction: `F_corrected = F_rd * (1.214 - 0.746 * (SL - 1.0))`

For correcting F_ref itself:

1. Compute `frac = (t - t_zone1_midpoint) / T0`
2. Apply `F_ref_corrected = F_ref * (1 + (1.214 - 0.746*(SL-1.0) - 1) * frac)`

For multiple sequential recordings at known activation intervals: each recording is at a different T on the decay timeline.  The same model and T0 concept extends naturally -- the correction at each recording scales with its T distance from the reference.

## Script

All fitting, visualization, and export is in `Workbench/FitRundownCorrection.m`.  The script is self-contained (loads data from `AllDataMerged.mat`).

---

## Feature analysis: effect of WRR correction on ktr and slack

The corrected F_ref variants (original, slope-only, M1 model, M3 model) were processed through `CompareProtocols.m` to extract ktr and slack features.  Results below are from the 03/27/2026 8 mM dataset.

### ktr experiment

| Dataset | ktr (s⁻¹) | ΔF (kPa) | RMSE (kPa) |
|---------|-----------|----------|------------|
| Original (no correction) | 45.7 | 68.7 | 1.58 |
| Slope corrected | 45.2 | 68.3 | 1.58 |
| M1 model corrected | 45.1 | 68.2 | 1.57 |
| M3 model corrected | 45.1 | 68.3 | 1.58 |

**ktr is correction-invariant.**  All three WRR correction methods agree to <1%.  This is expected: ktr is an exponential rate constant extracted from the restretch recovery transient.  The correction applies a slowly-varying amplitude scaling (~5% over 10 s) that does not meaningfully alter the fitted rate.

### Slack ktr curves

The figure below overlays ktr vs SL_slack curves (lines) from the slack protocol with single-point ktr values from the restretch experiment (filled circles).

![ktr vs SLslack: slack curves vs ktr experiment](figures/fig500_ktr_vs_SLslack_zoom.png)

Observations:

- **All four correction variants stack on top of each other** in the slack ktr curves.  The spread between corrected and uncorrected is ≤ 0.5 s⁻¹ (< 1%), well within noise.
- **A systematic gap of ~3–6 s⁻¹ separates the ktr experiment from the slack-derived ktr**, regardless of correction.  At SL ≈ 2.0 µm: slack ktr ≈ 49 s⁻¹ vs restretch ktr ≈ 45 s⁻¹.  This gap is a property of the two protocols (full detachment in the restretch experiment vs partial re-attachment during slack), not an artifact of rundown.

### Slack force features

The figure below shows all slack features across the five slack segments (SL_slack = 1.88–2.04 µm) for all four correction variants.

![Slack feature comparison](figures/fig403_slack_features.png)

Key feature changes after WRR correction:

| Feature | Original | Slope | M1 | M3 | Notes |
|---------|----------|-------|----|----|-------|
| A (force amplitude, kPa, at SL 2.04) | 69.6 | 70.9 | 71.2 | 70.9 | +1.3–1.6 kPa after correction |
| steady (kPa, SL 2.04) | 80.8 | 82.3 | 82.3 | 81.7 | +1.0–1.5 kPa |
| peak1_y (kPa, SL 2.04) | 93.6 | 95.0 | 95.7 | 94.4 | +0.8–2.1 kPa |
| ktr (s⁻¹, SL 2.04) | 50.5 | 50.3 | 50.2 | 50.8 | < 0.3 s⁻¹ change |
| t0 (onset time, s) | 0.0055 | 0.0052 | 0.0055 | 0.0054 | negligible |

The correction acts exclusively on force **amplitude**: steady-state isometric force, peak after restretch, valley, and fitted amplitude A all shift up by ~1.5–2.5 kPa.  The kinetic features (ktr rate, t0 onset, oscillation ω) are unchanged to within noise.

Slope and M1 corrections give near-identical steady-state force (+~1.5 kPa), while M3 gives a slightly smaller lift (+~0.9 kPa).  The SL-dependence of the multiplicative M1/M3 models introduces a fractionally larger correction at shorter SL (as designed) but the effect on the final features is minor — all three methods are interchangeable for practical data processing.

### Verdict

**Nailed for force; irrelevant for kinetics.**

The within-recording rundown correction successfully recovers ~1.5–2.5 kPa of force amplitude that would otherwise be lost to drift over the 10 s activation window.  The kinetic parameters that determine cross-bridge cycling rates (ktr, t0, ω) are indifferent to the correction at the level of precision achievable from these recordings.  The M1 model correction is recommended as it is SL-aware and physically principled, but slope correction is sufficient for most analyses.
