# Restretch Force Discrepancy Analysis

## Observation

During the re-stretch phase of the slack protocol (t ≈ 2.9–2.93 s, `RunSlackSegments = 'Last'`),
the experimental data shows a characteristic double-peak shape:

1. **First peak** — sharp force overshoot at restretch onset
2. **Valley** — brief dip as overshoot collapses
3. **Second peak** — force recovery as restretch velocity drops to zero

The model produces none of this: force rises weakly and gradually, with errors up to ~23 units
at the first peak and ~17 at the second.

## Diagnostic (figure 201)

Three panels were examined over t = 2.87–3.0 s:

### Force decomposition
Model force is far too small and too slow during restretch. Passive force (FXBPassive) is
negligible throughout, so the discrepancy is entirely in the active cross-bridge force.

### Attached fractions (p1, p2)
- Total attached fraction drops only modestly during the slack: ~0.35 → ~0.22.
  The bridges are *not* all detaching — there are plenty left to generate a restretch response.
- p2 (strongly attached) **decreases further during the restretch itself** (0.22 → 0.18)
  before recovering. The model is losing bridges while being actively stretched.
- p1 stays near zero throughout (fast p1→p2 transition keeps p1 depopulated — expected).

### Lengths (SL vs. LXB)
**Key finding:** LXB barely moves during restretch. SL jumps 1.89 → 2.0 µm in ~0.02 s,
but LXB lags ~0.1 µm behind and only catches up after SL has stopped moving.
The serial elastic element (kSE) is absorbing essentially all the rapid length change.
Cross-bridges never experience the strain increase — they cannot generate an elastic force spike.

## Root Causes

### 1. kSE too low — primary, parametric

The 0.1 µm SL–LXB gap during restretch means cross-bridges see almost no strain increase.
The first force peak in the data is the elastic response of attached bridges yanked to high
positive strain. With kSE too compliant, that strain never reaches the bridges.

**Current value:** `kSE = 2000`
**Suggested action:** Increase kSE by 3–5× and re-evaluate — this is the single most
impactful change and will reveal how much remaining discrepancy is structural.

### 2. p2 detaching during restretch — secondary, parametric

Even without a strain increase at LXB, p2 decreases during the restretch. Bridges that should
survive (and generate force once kSE is stiffened) are being lost. Likely cause: R1 or R1D at
near-zero or slightly negative strain is too fast. With stiffer kSE those same bridges would
instead be pushed to moderate positive strain — verify that R1 at s ≈ 0.01–0.02 is not too
aggressive.

### 3. UseA2AttachmentShift is the reattachment mechanism — but gated by kSE

`UseA2AttachmentShift` (ON in baseline) allows p2 heads at s > `s_threshold_R` (≈0.0046 µm)
to hop by one actin monomer step (`-d_actin` = −0.0055 µm) toward zero strain. This produces
a burst of force redistribution that generates the valley/second-peak pattern seen in the data.

However, with kSE = 2000 the serial elastic element absorbs the restretch — bridges never
reach the strain threshold. The mechanism exists but cannot fire. This is **not a structural
gap**: the physics is implemented. The fix is to increase kSE so that the restretch strain
actually reaches the cross-bridges, at which point A2Shift fires automatically.

Secondary tuning knobs: `slope` (hopping rate per unit excess strain) and `s_threshold_R`
(threshold for hopping onset). Lowering `s_threshold_R` and/or raising `slope` will sharpen
the valley and second peak.

## Summary Table

| Issue | Type | Suggested fix |
|---|---|---|
| LXB does not follow SL during restretch | Parametric | Increase `kSE` |
| First elastic peak absent | Consequence of low kSE | Same |
| p2 detaches during restretch | Parametric | Reduce R1 at moderate positive strain |
| Second peak too slow | Parametric + structural | Increase `ka`; add velocity-dependent ka |
| Double-peak shape not reproducible | Parametric (kSE gates A2Shift) | Increase `kSE`; tune `slope` and `s_threshold_R` |

## Recommended Next Step

Increase `kSE` (try 3–5×). `UseA2AttachmentShift` is already implemented and will fire once
bridges reach threshold strain. After kSE is tuned, adjust `slope` and `s_threshold_R` to
match the sharpness of the valley and second peak. See `TuneRestretch.m` for a sweep.
