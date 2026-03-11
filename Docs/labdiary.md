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

## Running TODOs

1. **Implement velocity-dependent vernier effect as a test** — add `f_vernier(|v_hs|) = 1 + alpha * |v|/(|v| + v_ref)` on attachment rate as an alternative/complementary test of the vernier hypothesis
2. **Set the occupation of attached states (p1+p2) to ~30% at steady state** — adjust rate constants or xrate to reach physiological attached fraction
3. **Fit the time-course of a single slack** — start from best parametrization (ModelParamsSlackKtrOpt), focus on matching the full transient shape
4. **Debug why steady-state force depends so much on the right side of the rate transition curve** — investigate piecewise strain-dependent detachment (`PieceWiseStrainDepParams/X`), understand why right-side breakpoints disproportionately affect isometric force
