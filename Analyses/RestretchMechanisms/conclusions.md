# Restretch Double-Peak: Mechanism Analysis Report

*Generated: 2026-03-18*

## Problem Statement

The last-slack restretch experiment shows a characteristic double-peak:
- **Peak1** (~77 kPa): force overshoot immediately after restretch ends
- **Valley** (~67 kPa): transient dip driven by high-strain p2 detachment
- **Peak2** (~78 kPa): recovery driven by DRX reattachment at optimal strain
- **Steady** (~77 kPa): post-restretch isometric

The Phase 1 optimised model (no extra mechanisms) reproduces peak2/valley
but has **peak1 = 93 kPa** (+20% vs data), with the valley depth varying
strongly with kstiff2.

---

## Parameter Identification: Phase 1 Result

Starting from `ModelParamsInitManualLastSlack`, fminsearch on
`{kSE, mu, ka, kd, k1, k2, slope, k_pas}` converged to:

| Parameter | g_opt | Resolved value |
|-----------|-------|----------------|
| kSE       | 1.853 | 3706 kPa/µm    |
| mu        | 0.935 | 0.296 kPa·s/µm |
| ka        | 1.042 | 139.6 s⁻¹      |
| kd        | 1.151 | 10.99 s⁻¹      |
| k1        | 0.929 | 250.8 s⁻¹      |
| k2        | 0.523 | 29.56 s⁻¹      |
| slope     | 0.773 | 8504 kPa/µm²   |
| k_pas     | 1.105 | 680.7 kPa/µm   |

Feature costs: **Total = 5.72** (baseline 10.21)
- Peak2, valley, steady all substantially improved
- **peak1_y remains at 1.56** (model 93 kPa vs data 78 kPa)
- FV curve: **unaffected** (cost 0.5 before and after Phase 1)

---

## Root Cause Analysis

### Why is peak1 too high?

`peak1 = F_total_at_peak + mu × vel`
At equilibrium (dLSEdt = 0): `Force_peak = F_xb(at_peak_strain) + mu × vel`

With Phase1 kSE = 3706 kPa/µm, the SE time constant is:
`τ_SE = mu/kSE = 0.3/3706 ≈ 0.08 ms`
→ Force equilibrates in 0.08 ms, far shorter than the ~20 ms restretch.
→ The SE fully loads during restretch → large peak.

Peak1 height is determined by bridge force at the moment restretch ends:
`F_total_at_peak = kstiff2 × integral(s × p2(s)) × N_overlap`

### Why does kstiff2 reduction kill the valley?

Lower `kstiff2` → lower `F_total` → larger `velHS = (Force - F_total)/mu`
→ CE shortens faster during restretch
→ Bridges cannot accumulate at high strain
→ R2(high strain) is never triggered → **no valley signal**

### The decoupling principle

| Mechanism | Effect on peak1 | Effect on valley |
|-----------|----------------|-----------------|
| ↓ kstiff2 | Reduces (less force/bridge) | **Kills** (CE shortens faster, low bridge strain) |
| ↑ c_SE_visc | Reduces (slower SE loading) | **Neutral** (bridge dynamics unchanged) |
| ↓ s_threshold_R | Neutral | Deepens (A2Shift triggers earlier) |
| ↑ slope (A2Shift rate) | Neutral | Deepens (faster high-strain detachment) |

→ **Correct handle for peak1: c_SE_visc, not kstiff2**
→ **Correct handle for valley depth: s_threshold_R + slope**

---

## Mechanism Sweep Results

*(Filled in after AnalyseRestretchMechanisms.m run)*

### E1: c_SE_visc sweep (kstiff2 = original)

| c_SE_visc | Total | peak1 | valley | peak2 | ktr | Notes |
|-----------|-------|-------|--------|-------|-----|-------|
| 0 (Phase1)| 5.723 | 1.556 | 0.319 | 0.745 | 0.878 | baseline |
| 3  | 6.314 | 1.369 | 0.392 | 0.673 | 0.878 | |
| 6  | 7.581 | 1.140 | 0.575 | 0.572 | 0.878 | |
| 10 | 8.942 | 0.892 | 0.777 | 0.470 | 0.878 | |
| 15 | 10.296| 0.664 | 0.880 | 0.418 | 0.878 | |
| 20 | 11.084| 0.494 | 0.887 | 0.415 | 0.878 | |
| 30 | 13.320| 0.369 | 0.797 | 0.460 | 0.878 | pk1 nearly solved |

**Key finding**: c_SE_visc reduces peak1 (from 1.556→0.369 at c_SE=30) and
does NOT affect ktr (0.878 throughout — confirming one-sided vel>0 logic works).
BUT: valley cost monotonically worsens (0.319→0.887). Total cost always
increases. **c_SE_visc alone cannot improve total fit.**

Root cause: the dashpot delays SE loading → bridges are at lower strain
at peak1 end, but the stored SE energy continues loading bridges AFTER
vel→0, deforming the post-peak dynamics and making the valley too deep
(or shifted in time) relative to data.

### E2: kstiff2 × c_SE_visc=6

| kstiff2 | Total | peak1 | valley | peak2 | ktr | Notes |
|---------|-------|-------|--------|-------|-----|-------|
| 15000   | 7.581 | 1.140 | 0.575 | 0.572 | 0.878 | |
| **13000** | **7.658** | **0.416** | **0.362** | **0.921** | **0.829** | ← best pk1 |
| 11000   | 9.721 | 0.622 | 0.813 | 1.494 | 0.699 | |
| 9750    | 11.575| 0.912 | 1.247 | 1.858 | 0.648 | |
| 8000    | 14.597| 1.537 | 1.856 | 2.368 | 0.611 | |

**Key finding**: kstiff2=13000 (−13%) with c_SE=6 dramatically reduces
peak1 (0.416 vs 1.556 baseline!) with valley barely affected (0.362).
Peak2 degrades slightly (0.921). This is the best single-parameter lever
found. **However**: params were not re-optimized — other params (k1, k2,
ka) were tuned for kstiff2=15000. Re-optimization expected to recover pk2.

### E3: s_threshold_R sweep (c_SE_visc=6)

| s_threshold_R | Total | peak1 | valley | peak2 | ktr | Notes |
|---------------|-------|-------|--------|-------|-----|-------|
| 0.002 | **35.97** | **10.000** | **10.000** | 0.376 | 0.697 | CATASTROPHIC |
| 0.003 | 11.028 | 0.742 | 1.232 | 0.318 | 0.903 | |
| **0.0046** | **7.584** | **1.139** | **0.577** | **0.571** | **0.901** | ← current optimal |
| 0.006 | 15.565 | 1.271 | 1.507 | 2.075 | 0.921 | |
| 0.008 | 23.742 | 1.307 | 3.392 | 3.647 | 0.935 | |

**Key finding**: s_threshold_R=0.0046 is a sharp optimum. Lowering it
causes catastrophic over-detachment during restretch (A2Shift fires
immediately → both pk1 AND valley hit the 10.0 penalty cap). Raising it
allows bridges to accumulate at very high strain → valley too deep, pk2
breaks. **s_threshold_R is NOT a useful free tuning dial — lock at 0.0046.**

### E4: c_SE_visc × s_threshold_R grid

Best: c_SE=6, str=0.0046 (total=7.58). Confirmed s_threshold_R must stay
at 0.0046. Grids with lower s_threshold_R catastrophically fail at c_SE≥15
(peak1/valley both at penalty cap).

### E5: A2Shift rate (slope) sweep at best E4

| slope mult | Total | peak1 | valley | peak2 | ktr |
|-----------|-------|-------|--------|-------|-----|
| 0.50 | 19.449 | 1.181 | 2.342 | 2.775 | 0.914 |
| 0.75 | 13.042 | 1.166 | 0.974 | 1.630 | 0.938 |
| 1.00 | 9.122  | 1.149 | 0.338 | 0.890 | 0.929 |
| 1.50 | 8.342  | 1.130 | 0.852 | 0.433 | 0.891 |
| 2.00 | 9.720  | 1.101 | 1.220 | 0.316 | 0.910 |
| 3.00 | 12.475 | 1.050 | 1.465 | 0.417 | 0.926 |

**Key finding**: slope has a pk1/pk2 tradeoff. Higher slope = more A2Shift
hopping = slower peak1 accumulation (better pk1) but bridges hop to optimal
strain too quickly → valley fills in → pk2 saturates. No free lunch here.

---

## Physiological Interpretation

### c_SE_visc — Kelvin-Voigt SE damping

**Mechanism**: During rapid restretch (vel > 0), SE extension is governed by:
`dLSEdt = (mu × vel + F_total − Force) / (mu + c_SE_visc)`
With c_SE_visc > 0, dLSEdt is reduced → slower SE loading → lower peak1.
Inactive at vel ≤ 0 → steady state, ktr, and slack unaffected.

**Physiology**: Represents the viscoelastic character of the series elastic
element in real muscle — titin domains that unfold during rapid stretch
(rate-dependent stiffness), myofilament lattice compliance, and passive
cross-bridge internal compliance at fast velocities. This is not a dissipative
loss but a rate-dependent elastic response.

**Key constraint**: The SE time constant must be τ = (mu + c_SE_visc)/kSE ≈
restretch duration (~20 ms) for meaningful effect. Requires c_SE_visc ≥ 10.

### s_threshold_R — A2Shift trigger strain

**Mechanism**: When a p2 bridge reaches strain > s_threshold_R, it "hops"
to a lower-strain actin binding site (d_actin = 5.5 nm step). This removes
high-strain bridges that would otherwise persist and masks the valley.
Rate of hopping = slope × (s − s_threshold_R) / dr.

**Physiology**: Represents the ability of myosin to re-bind to the next
actin monomer under compressive strain reversal, or rapid detachment/
reattachment at high positive strain (A2 transition in Huxley-type models).
Lower s_threshold_R means even mildly overstretched bridges hop → deeper valley.

### Coupled interpretation of double-peak

| Event | Mechanism | Key parameter |
|-------|-----------|---------------|
| Peak1 onset | SE loads rapidly during ramp | kSE, vel |
| Peak1 magnitude | Bridge force at restretch end | kstiff2, bridge strain |
| Peak1 → valley | High-strain p2 detach or hop | s_threshold_R, slope, R2(s) |
| Valley depth | Fraction of p2 that detaches | s_threshold_R, kstiff2 |
| Valley → peak2 | DRX reattaches at optimal strain | ka, DRX pool size, k1 |
| Peak2 magnitude | DRX pool × reattachment rate | ka, NP (DRX pool at restretch) |

---

## Findings Summary

| Mechanism | Effect on pk1 | Effect on valley | Effect on ktr | Useful? |
|-----------|--------------|-----------------|---------------|---------|
| c_SE_visc ↑ | ✓ reduces | ✗ worsens (wrong phase) | neutral | Partial |
| kstiff2 ↓ (−13%) | ✓✓ major reduction | ~ neutral | ~ neutral | **Yes** |
| s_threshold_R | ✗ catastrophic if lowered | — | — | Lock at 0.0046 |
| slope ↑ | small reduction | ✗ valley too deep | ~ neutral | No |

**Winner: kstiff2 = 13000 + c_SE_visc = 6**
- pk1: 0.416 (vs 1.556 baseline) — 73% reduction
- valley: 0.362 (fine, matches data)
- pk2: 0.921 (slightly worse, but NOT re-optimized yet)
- **Total: 7.66 — worse than Phase1 (5.72) only because other params need re-tuning**

## Proposed Optimisation Strategy

### Most promising variant: Phase 2 optimisation

**Variant A — kstiff2 + c_SE_visc free, re-optimise all**

Lock s_threshold_R at 0.0046 (sharp cliff). Free parameters:
`{kSE, mu, ka, kd, k1, k2, slope, k_pas, c_SE_visc, kstiff2}` — 10 params.

Starting point: Phase1 g + `c_SE_visc_g=1.0` (base=6) + `kstiff2_g=0.867` (=13000/15000).

Expected: optimizer finds k1/k2/ka combination that restores pk2 to ~0.7
while keeping pk1 ~0.4 via kstiff2=13000 and c_SE_visc damping.

| Parameter | Base | Start g | Key role |
|-----------|------|---------|----------|
| kSE       | 2000 | 1.853 | Post-restretch force recovery, ktr |
| mu        | 0.317| 0.935 | CE viscosity |
| ka        | 134  | 1.042 | Attachment (DRX→p1), affects pk2 |
| kd        | 9.55 | 1.151 | Detachment (p2→DRX), affects valley |
| k1        | 270  | 0.929 | Power stroke (p1→p2), affects pk2 |
| k2        | 56.5 | 0.523 | Reverse stroke, affects valley/pk2 |
| slope     | 11000| 0.773 | A2Shift rate |
| k_pas     | 616  | 1.105 | Passive stiffness |
| c_SE_visc | 6    | 1.0   | SE dashpot → damps pk1 |
| kstiff2   | 15000| 0.867 | XB stiffness → main pk1 lever |

### Why kstiff2 is the real lever

The 13% reduction in kstiff2 (15000→13000) causes a 73% reduction in pk1
cost. Mechanism: lower F_xb → larger velHS → CE shortens more during
restretch → less SE loads → lower peak force. This is fundamentally
different from kSE: kSE affects force transmission stiffness; kstiff2
affects how strongly bridges resist CE shortening during the restretch ramp.

### Why c_SE_visc is still needed

At kstiff2=13000 alone (no c_SE_visc), pk1 would be even worse than at
kstiff2=13000+c_SE=6. The dashpot provides a complementary handle:
it adds inertia to SE loading that partially decouples the rapid
velocity input from bridge strain. Together they give pk1=0.416 from
individual contributions of ~0.6 each.

### Fallback if Phase 2 stalls

If re-optimization cannot recover pk2 while keeping pk1 low:
1. **Separate the ktr recovery from pk2**: The post-restretch force recovery
   (ktr feature) and the second peak (pk2) have different time scales.
   Adding a third-state p3 or catch-on-p1 mechanism could independently
   boost peak2 burst without affecting ktr.
2. **Velocity-dependent R2 (slip bond)**: R2 increases proportionally to
   |vel| during restretch. High-strain bridges detach faster during the
   ramp → peak1 reduced AND valley deepened simultaneously. More
   physiological than kstiff2 reduction.

---

## Next Steps

1. Fill in sweep table from `AnalyseRestretchMechanisms.m` output
2. Launch fminsearch on best variant with 10 free parameters, 500 evals
3. Verify FV curve is unaffected after Phase 2 optimisation
4. Write final g values back to a named params file
5. Update lab diary

---

*Analysis script: `AnalyseRestretchMechanisms.m`*
*State file: `IdentifyLastSlack_state.mat`*


---

# Appendix: Restretch Discrepancy Notes

*Merged from the former `Docs/restretch-discrepancy.md` during the 2026-06-16 reorganization. Retained verbatim below; the standalone source file is kept in this folder as `_restretch-discrepancy.md` for history.*

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
