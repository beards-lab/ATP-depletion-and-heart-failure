# Lab Diary

## 2025-03-10 — Vernier mismatch / target zone saturation

### Motivation

The force-velocity (FV) data shows that the measured isometric force (F0) is lower than what you get by extrapolating the Hill curve fit to the shortening data back to v=0. We explored whether myosin crown spacing and the vernier mismatch between thick and thin filament repeats could explain this.

### Analysis

- **Myosin crowns** repeat every 14.3 nm (42.9 nm super-repeat). **Actin target zones** repeat every ~37 nm. The mismatch means not all crowns are aligned with binding sites at any given registration.
- In the current **mean-field model**, a simple vernier effect reduces to a constant scaling factor on ka — it cannot produce velocity-dependent modulation. The ergodic time-average availability is W/L regardless of velocity.
- However, there is an **asymmetry**: at isometric, a head that detaches into a "blocked" position (between target zones) is stuck forever. During shortening, filament sliding unblocks it. This makes the effective cycling pool velocity-dependent.
- More practically: the finite number of actin sites creates a **saturation effect** on the per-bin attached density. At isometric (high duty ratio, population piles up at s≈0), the local density saturates. During shortening (convection spreads population), saturation is relieved.

### Changes made

**Lattice spacing (one-sided Gaussian)** — corrected from symmetric to one-sided:
- `dPUdT_CombinedTransitions.m`: `f_lattice` now only penalizes when lattice is wider than optimal (shorter SL). At longer SL (narrower lattice), `f_lattice = 1` — consistent with Frank-Starling.

**Target zone saturation** — new mechanism (`UseTargetZoneSaturation`):
- `getParams.m`: Added `UseTargetZoneSaturation` (default false), `max_attached_per_bin` (default 0.01)
- `dPUdT_CombinedTransitions.m`: Per-bin saturation factor `f_sat(j) = max(0, 1 - (p1(j)+p2(j)+p3(j))*dS / max_attached_per_bin)` applied at each attachment bin in all three attachment branches (cached kernel, point, triangular)
- `evaluateModel.m`: Exposed `out.f_saturation` (value at attachment point s≈0) in output struct
- Output vectors extended by 1 element (f_saturation appended)

### Files modified

| File | Change |
|------|--------|
| `Model/getParams.m` | +2 params: `UseTargetZoneSaturation`, `max_attached_per_bin` |
| `Model/dPUdT_CombinedTransitions.m` | One-sided lattice factor; per-bin saturation on all attachment branches; f_saturation in outputs |
| `Model/evaluateModel.m` | `out.f_saturation` unpacking for 2-state and 3-state |

### Notes

- `max_attached_per_bin` default is 0.01 — set to be influential with current low attached fractions. Physical estimate is ~dS*163/300 (≈0.001 for dS=2nm). Will need tuning.
- The current parametrization (ModelParamsSlackKtrOpt) has very low attached fractions (~5-15%) due to xrate=0.435 and rate-limiting hydrolysis. Physiological is ~30-40%.
- After enabling saturation, ka should be re-optimized upward to restore isometric force. The increased ka then gives relatively more force during shortening (where saturation is weaker), producing the F0 deficit.

---

## 2026-03-17 — Restretch double-peak: catch bonds

### Observation

The slack/restretch experiment shows a characteristic double-peak during restretch (~t=2.91–2.93 s): a sharp first peak, a valley, then a second peak. The model (baseline) cannot reproduce it — force rises slowly and monotonically.

### Diagnostic (from figure 201)

- **LXB barely moves during restretch**: SL jumps 1.89→2.0 µm in 0.02 s but LXB follows sluggishly due to serial element compliance (kSE=2000). Cross-bridges never see the rapid strain increase → no elastic force spike.
- **p2 decreases during restretch** (0.22→0.18): bridges detach *during* the stretch instead of being held by it.
- **p1 near zero throughout**: p1→p2 transition is fast, consistent with normal kinetics.

### Hypotheses tested

| ID | Mechanism | Physical basis |
|---|---|---|
| H1 | One-sided ViscoSE (vel>0 only) | SE resists rapid extension, transmits more strain to XBs |
| H2 | Dynamic passive / titin viscoelasticity | Titin stiffer at high stretch rates; re-engages after slack |
| H3 | **Catch bonds** — reduced R1D+R2 during vel>0 | Actomyosin bonds strengthen under rapid stretch |
| H4 | Stretch activation — ka boost during vel>0 | Rapid lengthening opens binding sites |

### Results

| Variant | E_slack |
|---|---|
| Baseline | 486.8 |
| kSE x3 + A2Shift (prev best) | 421.6 |
| H1 ViscoSE c=0.15 | 483.8 |
| **H3 Catch k=0.05** | **217.7** |
| H3 k=0.06 | 223.7 |
| H3 + ViscoSE 0.15 | 230.8 |
| H4 SA k=0.5 | 462.9 |

**Catch bonds (H3, k=0.05) win decisively** — 55% reduction vs baseline, 48% better than previous best approach.

### Interpretation

During rapid restretch, `R1D` (p1 detachment) and `R2` (p2 detachment) are suppressed by factor `max(0, 1 - k_catch*vel)`. At vel=6 µm/s and k=0.05: suppression factor = 1 - 0.05*6 = 0.7 (30% rate reduction). Bridges survive longer → generate force spike (first peak). As velocity drops to zero, catch effect vanishes and bridges detach at normal rates → valley. Fresh reattachment from DRX pool → second peak. This is exactly the observed double-peak structure.

**Physiological basis**: Actomyosin is a catch bond. Single-molecule experiments show the bond lifetime *increases* under applied load/rapid stretch. This is also consistent with "residual force enhancement" after stretch (post-stretch muscles are stronger than pre-stretch at same length).

### Implementation

`dPUdT_CombinedTransitions.m` — after `UseStrictDetachmentAt` block:
```matlab
if params.UseCatchBond && vel > 0
    catch_factor = max(0, 1 - params.k_catch_bond * vel);
    R1D = R1D * catch_factor;
    R2  = R2  * catch_factor;
end
```
`getParams.m`: `'UseCatchBond', false` and `'k_catch_bond', 0` as defaults.

### Feature cost analysis (evalFeatureCost)

Full feature-by-feature breakdown revealed hidden tradeoffs:

| Feature | Baseline | Catch only | kSE=3000+Catch |
|---|---|---|---|
| ktr | 0.595 | 0.603 (+1%) | 1.196 (+101%) |
| peak1_y | 0.669 | 0.491 (−27%) | 0.491 |
| peak1_dSL | 1.507 | 3.697 (+145%) | 1.766 (+17%) |
| peak2 | 1.687 | 0.577 (−66%) | 0.359 |
| vall_y | 1.787 | 0.246 (−86%) | 0.194 |
| vall2_dy | 3.347 | 1.324 (−60%) | 1.205 |
| **TOTAL** | **16.99** | **14.14** | **12.52** |

**Key finding**: kSE=3000 *alone* doubles ktr error (0.595→1.197). kSE controls both restretch peak timing AND ktr dynamics — they cannot be independently tuned. The catch bond alone barely touches ktr (+1.3%).

**FV unaffected** by catch bond: FV shortening uses vel<0, catch bond activates only for vel>0.

**VisSE + Catch**: does not fix peak1_dSL (makes it worse, 3.697→4.619). Not useful.

### Recommendation

**Catch bond k=0.05 alone** is the safe improvement: 17% reduction in total feature cost, ktr preserved. The peak1_dSL tradeoff (+145%) means the first peak arrives slightly late in SL-space — acceptable given the large gains in valley and second-peak shape.

If pursuing kSE=3000+Catch for the extra 12% gain, must re-optimize k1/kd to restore ktr.

### p1-burst hypothesis — tested and rejected

Hypothesis: second peak driven by p1 doing synchronised power strokes after vel→0.

**Disproved**: p1=0.035 at vel=0, p2=0.189. Even full p1→p2 conversion yields ~3 force units. Second peak in data is ~20 units above valley. p1 is too small to matter.

**Real mechanism (confirmed)**: R2-catch on p2 drives both peaks.
1. R2 catch keeps p2 alive at positive strain during restretch → **first peak**
2. vel→0: catch gone → high-strain p2 detaches via R2 → **valley**
3. Large DRX pool (PD grew as bridges detached) → rapid reattachment → **second peak**

### Strain-limited catch bond — tested and failed

Hypothesis: apply R2 catch only at s < smax to fix peak1_dSL timing (catch at low strain, release at high strain).

**Failed**: All cutoffs (0.004–0.010 µm) give identical results. Bridges during restretch are ALL at low positive strain — LXB barely moves (kSE absorbs all velocity), so no bridge ever reaches the slip threshold. The strain distribution during restretch stays near s≈0 regardless of smax.

### peak1_dSL root cause: the kSE coupling

The first peak occurs at wrong dSL because R2 catch keeps p2 alive across the full restretch before the valley forms. Fixing this requires bridges to reach high strain quickly → stiffer SE → kSE=3000. But kSE=3000 also speeds up ktr force redevelopment → ktr cost doubles.

**Resolution path**: kSE=3000 + catch k=0.05 (total feature cost −26%), then re-tune k1/kd to restore ktr. The strain-limited catch becomes useful once kSE is stiffened (bridges will then actually reach the slip threshold during restretch).

### Files modified

| File | Change |
|---|---|
| `Model/dPUdT_CombinedTransitions.m` | H1 one-sided ViscoSE; H2 DynPassive; H3 CatchBond (with R1DOnly and CatchBondStrainMax variants); H4 StretchActivation |
| `Model/getParams.m` | Defaults for above + `UseCatchBondR1DOnly`, `CatchBondStrainMax` |

---

## 2026-03-17 — Kelvin-Voigt SE: implementation fix + full comparison

### Bug fix

The initial `UseViscoelasticSE` implementation was incorrect. The old code:
```matlab
mu_eff = params.mu + params.c_SE_visc;
dLSEdt = vel - (Force - F_total) / mu_eff;
```
expands to `dLSEdt = ((mu + c_SE)*vel + F_total - Force) / (mu + c_SE)` — the velocity term is scaled by `(mu + c_SE)`. With `c_SE=1000` and `mu≈0.3`, this gives `dLSEdt ≈ vel`, meaning all velocity dumps into SE extension (opposite of intended).

The correct Kelvin-Voigt implicit solve is `(mu + c_SE)*dLSEdt = mu*vel + F_total - Force`, giving:
```matlab
dLSEdt = (params.mu * vel + F_total - Force) / (params.mu + params.c_SE_visc);
```
With `c_SE=1000`: `dLSEdt ≈ 0.3*vel / 1000 ≈ 0` — SE barely extends, full velocity reaches cross-bridges. Fixed in `dPUdT_CombinedTransitions.m`.

### Results after fix

| Variant | E_slack |
|---|---|
| Baseline | 486.8 |
| kSE x3 + A2Shift | 421.6 |
| ViscoSE c=1000 | 454.6 |
| ViscoSE c=5000 | 459.7 |
| **ViscoSE c=5000 + kSE x2** | **391.3** |
| **Catch k=0.05** | **217.6** |
| Catch k=0.05 + kSE=3000 | 217.9 |
| Catch k=0.05 + ViscoSE c=1000 | 303.5 |
| Catch k=0.05 + ViscoSE c=5000 | 311.2 |

**ViscoSE alone** provides modest improvement (391 best), but is far inferior to catch bonds (217.6). ViscoSE **hurts** when combined with catch — the two mechanisms antagonize each other. Catch bonds alone remain the clear winner.

### Feature cost analysis (6-feature subset)

| Feature | Baseline | Catch k=0.05 | Catch+kSE=3000 | ViscoSE+kSEx2 |
|---|---|---|---|---|
| ktr | 1.175 | 1.175 | **0.602** | 0.670 |
| peak1_y | 0.602 | 0.430 | 0.399 | 10.000 |
| peak1_dSL | 0.342 | 1.751 | **0.401** | 10.000 |
| peak2 | 2.066 | 1.244 | 1.163 | 1.810 |
| vall_y | 1.498 | 0.599 | 0.541 | 10.000 |
| vall2_dy | 0.999 | 3.335 | 3.399 | 0.931 |
| **TOTAL** | **6.681** | **8.534** | **6.503** | **33.4** |

**Key finding**: `Catch k=0.05 + kSE=3000` beats baseline on total feature cost (6.503 < 6.681) AND beats catch-alone (8.534):
- ktr **improves** (1.175 → 0.602): stiffer SE speeds force redevelopment in ktr protocol as well
- peak1_dSL **fixed** (1.751 → 0.401): stiffer SE transmits restretch strain to bridges faster
- **Remaining issue**: `vall2_dy` = 3.399 (second peak rise rate) — worsens in all catch bond variants

### Updated recommendation

**Use Catch k=0.05 + kSE=3000** as the base configuration going forward. No k1/kd re-tuning needed (ktr already improved vs baseline). The remaining issue is the `vall2_dy` — investigate mechanisms that control the second peak rise rate (DRX pool size, ka rate).

### Files modified

| File | Change |
|---|---|
| `Model/dPUdT_CombinedTransitions.m` | Fix KV formula (direct dLSEdt); remove intermediate mu_eff variables |
| `TuneRestretch.m` | Added catch bond variants to sweep |

---

## Running TODOs

1. **Implement velocity-dependent vernier effect as a test** — add `f_vernier(|v_hs|) = 1 + alpha * |v|/(|v| + v_ref)` on attachment rate as an alternative/complementary test of the vernier hypothesis
2. **Set the occupation of attached states (p1+p2) to ~30% at steady state** — adjust rate constants or xrate to reach physiological attached fraction
3. **Fit the time-course of a single slack** — start from best parametrization (ModelParamsSlackKtrOpt), focus on matching the full transient shape
4. **Debug why steady-state force depends so much on the right side of the rate transition curve** — investigate piecewise strain-dependent detachment (`PieceWiseStrainDepParams/X`), understand why right-side breakpoints disproportionately affect isometric force
5. **Fix vall2_dy**: `Catch k=0.05 + kSE=3000` is the best config (total feature cost −3% vs baseline, E_slack −55%). Remaining issue is the second peak rise rate (`vall2_dy` = 3.4 vs baseline 1.0). Investigate DRX pool dynamics and ka to control reattachment speed after the valley.
