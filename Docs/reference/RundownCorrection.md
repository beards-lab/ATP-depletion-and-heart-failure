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
| F_ref at zone 1 | 68.3 kPa | Linear fit to zone 1+2 points of data(2) |
| F_rd at zone 1 | 56.4 kPa | Linear fit to zone 1+2 points of data(4) |
| slope_ref | -0.495 kPa/s | Within-recording force decay of F_ref |
| slope_rd | -0.467 kPa/s | Within-recording force decay of F_rd |
| Force gap at zone 1 | 11.9 kPa | F_ref - F_rd at zone 1 midpoint |

## T0: hypothetical decay distance

The key derived quantity is **T0** — the hypothetical activated time it would take for force to decay from the reference level to the run-down level at a constant rate:

    T0 = (F_ref(zone1) - F_rd(zone1)) / |slope_ref| = 11.9 / 0.495 = 24.0 s

T0 is *not* wall-clock time between recordings.  It represents **equivalent activated time**: the cumulative hot-activated exposure (muscle under calcium) that separates the two recordings.  Inactive recovery periods between activations do not contribute to T0.

T0 establishes a common time axis for the decay model: zone 1 of F_ref is at T = 0, zone 1 of F_rd is at T = T0.

![T0 concept and within-recording decay](figures/fig47_zoom_and_concept.png)

*Top panel: Zoomed timecourse (t = 71–79.5 s) showing F_ref original (black) and slope-corrected (blue dashed). Zone 1 (t = 71.5–72.4 s) and Zone 2 (t = 77.4–78.4 s) are shaded. The slope correction (+2.94 kPa additive, SL-independent) is shown between the two zones.  Bottom panel: Decay timeline concept. The dashed line shows force decay from zone 1 of F_ref (T = 0) at rate slope_ref. Zone 1 of F_rd falls at T = T0 = 24.0 s. Blue dashed arrow: slope correction at zone 2. Red arrow: M1 model correction at zone 2.*

---

## Step 1: Spatial correction model (M1) and F_rd correction

### M1 model

A linear correction surface is fitted from paired (F_rd, F_ref) data over the valid window (F > 3 kPa, t = 68--78 s).  Both signals are first slope-detrended to their zone-1 values before fitting, to remove within-recording trends.  F_ref is processed on the hi-fi smoothed grid (`t_fineShifted`); F_rd is processed on the original measured points (`t_shifted_orig`), preserving the raw signal without smoothing.

**Model M1** (SL-linear, 2 parameters):

    F_corrected = F_rd * (a + b * (L - L0))    where L = SL/2, L0 = 1.0

    a = 1.214,  b = -0.759

    R^2 = 0.962,  RMSE = 2.1 kPa  (999 fitting points)

![M1 correction curves at 4 SL levels](figures/fig50_M1_slopes.png)

*M1 correction curves at four sarcomere lengths.  Left axis: corrected force vs run-down force.  Right axis: correction ratio r = F_corr/F_rd.  Dashed horizontal lines indicate the mean corrected force at each SL level.  Higher correction at shorter SL, lower at longer SL — consistent with run-down impairing force more at shorter lengths.*

The correction ratio `a + b*(L - L0)` depends only on sarcomere length:

| SL (um) | L = SL/2 | Correction ratio |
|---------|----------|-----------------|
| 1.6     | 0.8      | 1.366 |
| 1.8     | 0.9      | 1.290 |
| 2.0     | 1.0      | 1.214 |
| 2.2     | 1.1      | 1.138 |

A more complex model (M3: +F×L interaction, 4 params, R²=0.963) is available in the code but M1 is sufficient for practical data processing.

### F_rd correction timecourse

Applying the M1 correction brings the run-down recording back up to the slope-detrended reference level (dashed goal line).

![F_rd correction timecourse](figures/fig46_bottom_frd.png)

*F_rd correction timecourse.  Black: F_ref original.  Grey: F_rd raw.  Dotted blue: M1 static correction.  Solid blue: M1 + time correction (rate from slope difference).  Purple: M1 + time (freely fitted rate).  Dashed black: F_ref slope-detrended (target).  The goal line is plotted last and sits on top.*

The zoom below confirms that all M1 variants converge tightly on the target:

![F_rd correction zoom 72–76 s](figures/fig48_zoom_frd.png)

*Zoom of the F_rd correction at t = 72–76 s, 40–100 kPa.  Static M1 (dotted) and time-extended variants (solid) land essentially on the slope-detrended F_ref target (dashed) throughout the slack and isometric segments.*

---

## Step 2: Model-based F_ref correction

The reference recording also decays within its own activation window.  Rather than correcting with a simple constant slope (-0.495 kPa/s), we reuse the spatial model M1 scaled by the fraction of T0:

    frac(t) = (t - t_zone_ref) / T0

    r_model(t) = a + b * (L(t) - L0)      [same M1 correction ratio]

    F_ref_corrected(t) = F_ref(t) * (1 + (r_model(t) - 1) * frac(t))

At zone 1 (frac = 0): no correction (identity).  
At zone 2 (frac ≈ 0.25): correction of ~5.3% at SL = 2.0 (~+3.5 kPa on 68 kPa).  
At T0 distance (frac = 1): full model correction (~21.4%).

### Comparison: slope vs model at zone 2 (SL = 2.0)

| Method | Correction | Effect on 68 kPa |
|--------|-----------|-------------------|
| None (original) | 0 | 68.3 kPa → ~65.3 kPa at zone 2 |
| Slope only | +2.94 kPa (additive) | 65.3 + 2.9 = 68.3 kPa |
| M1 model | ×1.053 (multiplicative) | 65.3 × 1.053 = 68.7 kPa |

The model-based correction is **multiplicative and SL-dependent**: at SL = 2.2 the correction is smaller (~4.0%), at SL = 1.8 it is larger (~6.5%).  The slope correction is additive and constant across all SL levels — a less accurate approximation of the underlying decay process.

![F_ref correction timecourse](figures/fig46_top_fref.png)

*F_ref correction comparison over the full activation window.  Black: original.  Blue dashed: slope-corrected.  Cyan: M1 model-based.  Green: M3 model-based.  All three correction methods diverge slightly at the SL = 2.2 stretch step where the M1/M3 corrections are smaller (as expected from the SL-dependence of run-down).*

The last slack segment (76.7–77.4 s) provides the best view of the correction magnitude:

![F_ref last slack zoom](figures/fig49_zoom_fref_slack.png)

*Zoom of the last slack segment (SL ≈ 1.88 µm, t = 76.7–77.4 s).  Slope correction (blue dashed) and M1/M3 model corrections (cyan/green) agree to within ~0.5 kPa.  The correction lifts the steady-state isometric level by ~1.5 kPa relative to the uncorrected signal.*

---

## Output files

All saved as tab-delimited `[time, SL, force]` in the data directory:

| File | Grid | Description |
|------|------|-------------|
| `Corr_Fref_original.txt` | hi-fi | F_ref as recorded (no correction) |
| `Corr_Fref_slopeOnly.txt` | hi-fi | F_ref with linear slope correction to zone 1 |
| `Corr_Fref_M1_model.txt` | hi-fi | F_ref with M1 model-based correction to zone 1 |
| `Corr_Fref_M3_model.txt` | hi-fi | F_ref with M3 model-based correction to zone 1 |
| `Fref_lofi.txt` | lo-fi | F_ref original, unsmoothed (t_shifted_orig) |
| `Corr_Frd_original.txt` | lo-fi | F_rd as recorded (no correction) |
| `Corr_Frd_M1_static.txt` | lo-fi | F_rd corrected by M1 static model |
| `Corr_Frd_M1t_slopeRate.txt` | lo-fi | F_rd corrected by M1 + time (rate from slope) |
| `Corr_Frd_M1t_fittedRate.txt` | lo-fi | F_rd corrected by M1 + time (rate freely fitted) |
| `Corr_Frd_M3_static.txt` | lo-fi | F_rd corrected by M3 static model |

F_ref files are saved on the smoothed hi-fi grid (`t_fineShifted`); F_rd files are saved on the original measured points (`t_shifted_orig`).

## How to apply to new data

For a new run-down recording with a known reference:

1. Fit the within-recording slope at two SL = 2.0 zones
2. Compute T0 = force_gap / |slope|
3. Apply M1 correction: `F_corrected = F_rd * (1.214 - 0.759 * (SL/2 - 1.0))`

For correcting F_ref itself:

1. Compute `frac = (t - t_zone1_midpoint) / T0`
2. Apply `F_ref_corrected = F_ref * (1 + (1.214 - 0.759*(SL/2-1.0) - 1) * frac)`

For multiple sequential recordings at known activation intervals: each recording is at a different T on the decay timeline.  The same model and T0 concept extends naturally — the correction at each recording scales with its T distance from the reference.

## Script

All fitting, visualization, and export is in `Workbench/FitRundownCorrection.m`.  The script is self-contained (loads data from `AllDataMerged.mat`).

---

## Feature analysis: effect of WRR correction on ktr and slack

The corrected F_ref variants (original, slope-only, M1 model, M3 model) were processed through `CompareProtocols.m` to extract ktr and slack features (Preset A).  Results below are from the 03/27/2026 8 mM dataset.

### ktr experiment

| Dataset | ktr (s⁻¹) | ΔF (kPa) | RMSE (kPa) |
|---------|-----------|----------|------------|
| Original (no correction) | 45.7 | 68.7 | 1.66 |
| Slope corrected | 45.2 | 68.3 | 1.65 |
| M1 model corrected | 45.1 | 68.2 | 1.63 |
| M3 model corrected | 45.1 | 68.3 | 1.64 |

**ktr is correction-invariant.**  All three WRR correction methods agree to <1%.  This is expected: ktr is an exponential rate constant extracted from the restretch recovery transient.  The correction applies a slowly-varying amplitude scaling (~5% over 10 s) that does not meaningfully alter the fitted rate.

### Slack ktr curves and feature comparison

![Slack feature comparison (WRR)](figures/fig403_WRR.png)

*All four WRR correction variants (original, slope, M1, M3) overlaid.  Top-left panel: ktr vs SLslack — all four traces stack on top of each other (spread ≤ 0.5 s⁻¹, < 1%).  Force amplitude panels (A, ktr2_A, peak1_y, steady): correction lifts all force features by ~1.5–2.5 kPa.  Kinetic panels (ktr2_k, ktr2_omega, t0): indistinguishable across correction variants.*

![ktr vs SLslack: slack curves vs ktr experiment (WRR)](figures/fig420_WRR.png)

*ktr vs SLslack from the slack protocol (lines) overlaid with single-point ktr values from the restretch experiment (filled markers).  All four correction variants stack on top of each other.  A systematic gap of ~3–5 s⁻¹ separates the ktr experiment from the slack-derived ktr — this is a protocol difference (full detachment in the restretch experiment vs partial re-attachment during slack), not an artifact of rundown.*

### Slack force features summary

Key feature changes after WRR correction (at SL_slack = 2.04 µm):

| Feature | Original | Slope | M1 | M3 | Notes |
|---------|----------|-------|----|----|-------|
| A (force amplitude, kPa) | 69.6 | 70.9 | 71.2 | 71.1 | +1.3–1.6 kPa after correction |
| steady (kPa) | 80.8 | 82.3 | 82.2 | 81.9 | +1.0–1.5 kPa |
| peak1_y (kPa) | 93.6 | 95.0 | 95.7 | 95.2 | +0.8–2.1 kPa |
| ktr (s⁻¹) | 50.5 | 50.3 | 50.2 | 50.4 | < 0.3 s⁻¹ change |
| t0 (onset time, s) | 0.0055 | 0.0052 | 0.0055 | 0.0055 | negligible |

The correction acts exclusively on force **amplitude**: steady-state isometric force, peak after restretch, valley, and fitted amplitude A all shift up by ~1.5–2.5 kPa.  The kinetic features (ktr rate, t0 onset, oscillation ω) are unchanged to within noise.

### Verdict

**Nailed for force; irrelevant for kinetics.**

The within-recording rundown correction successfully recovers ~1.5–2.5 kPa of force amplitude that would otherwise be lost to drift over the 10 s activation window.  The kinetic parameters that determine cross-bridge cycling rates (ktr, t0, ω) are indifferent to the correction at the level of precision achievable from these recordings.  The M1 model correction is recommended as it is SL-aware and physically principled, but slope correction is sufficient for most analyses.

---

## Feature analysis: compensating for rundown

This section asks whether M1/M3 spatial compensation can bring the run-down recording (`F_rd`) back to the reference level.  The comparison uses `CompareProtocols.m` Preset B:

| Dataset | Description |
|---------|-------------|
| 03/27 8mM Lo-Fi | Unsmoothed F_ref (`t_shifted_orig`) — the ground truth |
| 03/27 8mM repeat | F_rd as recorded (run-down, no correction) |
| 03/27 8mM M1 compensated | F_rd × M1 static correction |
| 03/27 8mM M3 compensated | F_rd × M3 static correction |

### ktr experiment

| Dataset | ktr (s⁻¹) | ΔF (kPa) |
|---------|-----------|----------|
| Lo-Fi reference | 46.9 | 68.7 |
| Repeat (run-down) | 37.2 | 56.4 |
| M1 compensated | 37.2 | 68.5 |
| M3 compensated | 37.5 | 68.5 |

The M1/M3 correction **fully restores force amplitude** (ΔF: 56.4 → ~68.5 kPa, matching the Lo-Fi reference at 68.7 kPa).  However, **ktr remains depressed**: the compensated recordings show ktr ≈ 37 s⁻¹ vs the Lo-Fi reference at 47 s⁻¹ — a persistent ~10 s⁻¹ deficit.

### Slack feature comparison

![Slack feature comparison (F_rd compensation)](figures/fig403_FrdComp.png)

*Preset B comparison: Lo-Fi reference (blue), repeat/run-down (red), M1 compensated (cyan), M3 compensated (green).  Top-left: ktr vs SLslack — the repeat and both compensated variants sit ~8–10 s⁻¹ below the Lo-Fi reference across all SL levels.  Force amplitude panels (A, ktr2_A, peak1_y, steady): M1 and M3 compensated lines land close to the Lo-Fi reference, substantially above the uncorrected repeat.  Kinetic panels: compensation does not recover ktr.*

![ktr vs SLslack + ktr experiment (F_rd compensation)](figures/fig420_FrdComp.png)

*ktr vs SLslack (lines) with ktr experiment single points (filled markers).  The Lo-Fi reference (blue) sits ~8–10 s⁻¹ above all run-down variants regardless of compensation.  The ktr experiment points (filled markers) confirm the pattern: Lo-Fi ktr ≈ 47 s⁻¹ vs repeat/compensated ktr ≈ 37 s⁻¹.*

### Summary

| Feature | Lo-Fi ref | Repeat (raw) | M1 comp | M3 comp | Recovery? |
|---------|-----------|--------------|---------|---------|-----------|
| A (kPa, SL 2.04) | ~70 | ~42 | ~67 | ~67 | ✓ amplitude recovered |
| steady (kPa, SL 2.04) | ~80 | ~71 | ~80 | ~80 | ✓ recovered |
| peak1_y (kPa, SL 2.04) | ~80 | ~65 | ~80 | ~80 | ✓ recovered |
| ktr (s⁻¹, SL 2.04) | ~47 | ~38 | ~38 | ~38 | ✗ not recovered |
| ktr2_k (s⁻¹) | ~47 | ~38 | ~38 | ~38 | ✗ not recovered |

### Verdict

**Force amplitude: fully recovered.  Kinetics: not recovered.**

The M1 spatial correction restores the amplitude of force generation from the run-down recording to near-reference levels across all SL values.  However, the cross-bridge cycling rate (ktr) remains depressed by ~10 s⁻¹ in the compensated signal.  This means the run-down recording captures a state of impaired kinetics that is **not** explained by amplitude scaling alone — the muscle has changed its kinetic state, not just its force capacity.

The spatial model corrects for the force-amplitude deficit caused by run-down but cannot restore the intrinsic rate of cross-bridge cycling.  When using M1/M3 compensated F_rd data, **force amplitudes are reliable** but **ktr values should be treated with caution** — they reflect the kinetic state of the impaired muscle, not the reference state.
