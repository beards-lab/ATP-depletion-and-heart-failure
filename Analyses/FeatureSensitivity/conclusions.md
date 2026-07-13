# Feature-resolved sensitivity & unused-mechanism screen (2-state)

**Question (user).** Without the 3rd state, go through each feature discrepancy and find
improvement: reweighting, unused mechanisms, or parameters. Base = `params_2state_a2hop.m`
(2-state, cost ~2.95). Tools: `featureSensitivity.m` (elasticity of every feature to 47
params), `unusedMechScreen.m` (OFF mechanisms), `probeP1dSL.m` (kstiff1/k_1 test).

## 1. The two master levers: `dr`, `dr2`

The power-stroke strain offsets `dr` (p2 force center) and `dr2` (R2 detachment center)
dominate **almost every feature** with elasticities 2–4.5 — far above any other parameter.
They slide the whole strain distribution relative to force and detachment, so they move
everything at once. High leverage, but the most coupled knobs (can't move one feature alone).
They are the key optimization DOF and are in `RunOptimFull`'s pool.

## 2. Per-feature levers (sign = toward data)

| feature | model / data | best levers (direction toward data) | free or walled |
|---|---|---|---|
| ktr | 52.4 / 49.2 | `k2`↓, `kmsrd`↓, `dr2`↓ | **walled** (↔ FV via k2/kmsrd); near-matched |
| peak1_y | 95.5 / 89.3 (too high) | `kstiff1`↓, `kstiff2`↓, `ka`↓, `dr`↓ | partly free via `kstiff1` (see §4) |
| peak1_dSL | 0.0345 / 0.0256 (too high) | `kstiff1`↓ (2.69), `k_1`↓ (2.32), `kSE`↑ | see §4 (kSE↑ was walled ↔ ktr) |
| FV2/FV3/FV4 | .55/.22/.09 vs .66/.32/.11 | `dr`↑ (1.2–2.7), `k2`↑, `kstiff2`↑, `kmsrd`↑ | **walled** (↑ raises peaks &/or ktr) |
| vall_y | 70.0 / 71.1 | (A2 hop already fixed) | matched |
| peak2 | 81.0 / 78.4 | `dr2`↑, `kstiff1`↓, `ka`↓ | small, coupled to peak1 |
| ovrsht_dy | 1.26 / 1.33 | `eta_M`↑, `kSE_M`↑, `k2`↑; **viscoSE** | free (small) |
| restretchSlope | 1506 / 1588 | `kSE`↑, `kstiff2`↑ | small, coupled |
| A | 60.7 / 62.2 | `kstiff2`↑, `ka`↑ | near-matched |
| steady, t0, FV1 | — | matched | — |

**The dominant residual (FV mid-tail) is walled**: its master lever `dr`↑ raises the
already-too-high peaks, and its secondary levers `k2`/`kmsrd`↑ raise ktr. So FV-tail is welded
to the peaks (via `dr`) and to ktr (via `k2`) — reweighting FV up only trades peaks/ktr for
tail. This is the structural residual that the 3rd state targets.

## 3. Unused mechanisms (screen)

| mechanism | effect | verdict |
|---|---|---|
| **Kelvin–Voigt series viscosity** (`c_SE_visc`=0.02) | cost 2.952→**2.890** via `ovrsht`+ktr; no peak1_dSL change | keep (small free win, vel-gated) |
| velocity-dependent titin (`c_titin_visc`) | inflates all restretch force → cost 15–248 | reject |
| stretch activation (`k_SA`) | inflates force → cost 173–442 | reject |
| built-in A2 attachment-shift | inflates force → cost 62 (ungated/untuned) | reject as-is |
| W-detachment, vernier | inert / unstable | reject |

**Key point:** every mechanism that *adds restretch force* overshoots, because the restretch
amplitude is already correct (A2 hop fixed valley/peak2). The residual is **shape/timing**
(`peak1_dSL`) and **FV-tail**, not force magnitude. Only the *dissipative* Kelvin–Voigt element
(which shapes without adding net force) helps, and only marginally.

## 4. peak1_dSL via kstiff1 / k_1

_[probeP1dSL pending — does kstiff1↓ / k_1↓ move peak1_dSL toward 0.0256 without the ktr
penalty kSE carried?]_

## 5. Reweighting verdict

_[to finalize with §4]_ Reweighting helps only the **free** features (ovrsht; peaks IF kstiff1
is clean). The **walled** dominant residual (FV mid-tail) cannot be improved by reweighting —
it just relocates error onto the peaks/ktr. For literal timecourse fidelity, add a direct
time-series residual term rather than reweighting landmarks.
