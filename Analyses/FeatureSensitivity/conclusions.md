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

## 4. peak1_dSL is doubly-walled (probeP1dSL, probeP1dSLcomp)

`kstiff1`↓ DOES fix `peak1_dSL` (0.0345→0.028) and `peak1_y`, and — unlike `kSE`↑ — it HOLDS ktr
(52.2). So it escapes the `peak1_dSL`↔ktr wall. But it hits a *second* wall: it collapses `vall_y`
(70→55) and `peak2` (81→66), because lowering p1 stiffness removes the force that holds the
restretch trough. Compensation fails: even `kA2hop`×2.2 + `kstiff2`×1.12 lifts `vall_y` only to
~59 (data 71) while `peak1_dSL` sits at 0.027 — best combined cost 6.8 vs base 2.95. (`k_1`↓ is
worse: it inflates `peak1_y`.)

**Verdict: `peak1_dSL` is doubly-welded** — to ktr (via `kSE`) AND to the valley/peak2 (via
`kstiff1`/p1-force). No 2-state lever or compensation reduces the excursion without breaking a
coupled feature. This is the cleanest structural case for the 3rd state.

## 5. Answers to the user's questions

**Reweighting the fn:** helps only *free* features. The two dominant residuals — FV mid-tail and
`peak1_dSL` — are **walled**, so up-weighting them only relocates error (buy FV-tail with
peaks/ktr; buy `peak1_dSL` with the valley/ktr). It changes which features look good in a figure
(a legitimate presentation choice) but cannot lower the achievable cost. Only `ovrsht` (via
Kelvin–Voigt `c_SE_visc`) is genuinely free, and it's already small.

**Unused mechanisms:** only the dissipative Kelvin–Voigt series element helps (2.95→2.89); every
force-adding mechanism overshoots because restretch amplitude is already correct.

**Parameters:** the master levers are `dr`/`dr2` (strain offsets) — highest leverage but move
everything. Fresh single-feature levers exist (`kstiff1` for the peaks, `eta_M`/`kSE_M` for
`ovrsht`) but all the big ones are coupled.

**Bottom line:** the 2-state fit is near its floor (~2.9); the two dominant residuals are
structural 2-state walls. The one caveat is that these are *point-probe* couplings — the true
floor test is a full high-dimensional re-optimization over `dr`/`dr2` + the PCHIP knot *positions*
(not just their rate values), which can align the strain distribution in ways pairwise probes
can't see. For literal timecourse fidelity (vs landmark features), add a direct time-series
residual term. The durable fix for the walled residuals remains the 3rd state.
