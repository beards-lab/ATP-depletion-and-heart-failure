# What the model is actually made of right now — and how well it fits

State at `params/realign2_opt.m` (the completed reoptim, snapshot #7 of 7), 2-state,
`@dPUdT_CombinedTransitions`, scored on **FV + slack + passive** (`RunKtr = 0`,
`RunStairs = 0`) with the optimiser's objective (25 feature terms, `ovrsht_dy`
dropped).

**Headline: L2 6.103 → 5.060, −17.1 %.** Reproduce with
[`PlotReportFigures.m`](PlotReportFigures.m).

The run accepted an improvement in **all 7 rounds and never stalled** before hitting
its 4 h budget, so this is a direction of travel, not a converged optimum.

---

## 1. The mechanisms in play

### Core cycle — 2 states, strain-resolved

`p1` (weakly attached) and `p2` (strongly attached, post-stroke) are probability
densities over cross-bridge strain `s`; the detached side is the scalar pools
`PT` (ATP-cocked, "DRX"), `PD` (primed, hydrolysed, attachment-ready) and the two
super-relaxed pools. Every rate that matters is a function of `s`.

| mechanism | flag | what it does | why it is there |
|---|---|---|---|
| **Piecewise strain dependence** | `UsePieceWiseStrainDep` | `R12`, `R21`, `R1D`, `R2` are pchip splines in `s` rather than exponentials | exponentials could not fit attachment and detachment simultaneously; the knots are fitted. `R2(s)` reaches 700–4600 s⁻¹ at re-stretch strains — pinned at both ends by the slack re-stretch peak and the ktr overstretch peak ([ApoPoolDetachment](../ApoPoolDetachment/conclusions.md)) |
| **Space extension** | `UseSpaceExtension` | grows the strain grid when the distribution hits a boundary | large re-stretches push `p2` past the static grid |
| **Filament overlap** | `UseOverlap`, `UseOverlapFactor`, `UseCalculatedN` | attachment ∝ thick/thin overlap from `L_thick`, `L_hbare`, `L_thin` | the slack protocol sweeps 0.94–1.10 ML, so overlap is not constant |
| **Lattice spacing** | `UseLatticeSpacing` | `d10` interfilament spacing modulates attachment | length change alters lattice spacing at constant volume |

### Regulation of attachment — three independent gates on `ka`

This is where most of the recent work has gone, and it is worth being blunt: **all
three are phenomenological.**

| mechanism | flag / params | what it does | status |
|---|---|---|---|
| **Global occupancy saturation** | `UseGlobalOccupancySaturation`, `P_bound_max` | `RD1 *= (1 − P_bound/P_bound_max)` — a mean-field cap on how many heads can be bound | the *only* faithful occupancy form available here: the strain axis is not a site axis, so per-bin caps are a category error ([BindingSiteOccupancy](../BindingSiteOccupancy/conclusions.md)). `P_bound_max` = 1.02 is above any physical bound, i.e. the term is nearly inert at this fit |
| **Registration availability** | `UseRegistrationAvailability`, `A0`, `v_ref_reg`, `tau_reg`, `RegAvailShorteningOnly` | a state `A_reg ∈ [A0,1]` gating `ka`, charged by *shortening* speed, relaxing with `tau_reg` | the static part (`A0` < 1, target-zone mismatch of actin's 5.5 nm monomer against myosin's 14.3 nm crown) is well founded. The **velocity-dependent part is admitted double-counting** — advection already carries heads past sites — and `RegAvailShorteningOnly` is *"a description of a fit, not of a mechanism"* ([RestretchVsKtrRecovery](../RestretchVsKtrRecovery/conclusions.md) §C) |
| **Compliant realignment (NEW)** | `UseLoadRealign`, `k_cr`, `tau_cr`, `F_cr`, `Adef_max` | a state `A_def` gating `ka` by `(1−A_def)`, driven by **bound mass stripped while the filament is loaded**: `dA_def/dt = k_cr·max(0,−d(bound)/dt)·F/(F+F_cr) − A_def/tau_cr` | §8. Zero at any steady state by construction, so it cannot move the operating point. The force weighting is what separates post-re-stretch (strips `p2` at 0.84 F_iso) from post-slack/ktr (strips `p2` at zero force) — proven by the load-blind control, which gets the same `rsK` gain but breaks `ktr` 1.07 → 0.83 |

### Super-relaxed / OFF states

`PT ⇄ P_SR ⇄ P_SRD → PD`. With `kmsr = 0` the SRX is a **one-way force-gated trap**:
entry `ksr0·e^(−F/σ₂)`, exit `kmsrd·e^(+F/σ_srd1)` — *both gates point the same way*,
so the in/out ratio swings **60×** across the force range (16.6 at F = 0 → 0.28 at
0.87 F_iso). Consequence, established in §4: **the SRX brake auto-disengages under
load**, which structurally obliges the model to make the post-re-stretch recovery its
fastest window. `UseSuperRelaxed`, `UseSuperRelaxedADP`.

`SRXFromR2HighStrain = 0.2` routes a fifth of the high-strain (rupture) detachment
flux into SRX instead of `PT`. **It is carried, not earned** — it was seeded at 0.2 so
the optimiser could explore it and was never drawn in the rounds that ran (§11). One-at-a-time
it was mildly harmful (§2). Treat as untested at this point, not as supported.

### Mechanics

| mechanism | flag / params | what it does |
|---|---|---|
| **Series elastic element** | `UseSerialStiffness`, `kSE`, `ekSE`, `mu`, `mu_neg` | a nonlinear spring + dashpot between the cross-bridges and the motor. `LSE` is a state; `dSLc = vel − dLSEdt` is the *contractile* velocity |
| **Titin passive** | `UsePassive`, `UseDynamicPassive`, `k_pas`, `Lsc0`, `gamma` | nonlinear passive force from sarcomere length |
| **Maxwell dashpot** | `UseMaxwellDashpot`, `kSE_M`, `eta_M` | a viscoelastic stress `X` that **charges on lengthening only** and decays with τ = `eta_M`. Carries 11–17 % of total force at the end of each re-stretch ramp. Held responsible for ~half the `rsK` excess ([RestretchRecoveryFit](../RestretchRecoveryFit/conclusions.md) §6b) |
| **Negative-strain stiffness** | `UseNegativeKstiff`, `kstiff1_n`, `kstiff2_n` | separate stiffness for negatively strained bridges |
| **A2 re-attachment** | `UseA2Reattaching`, `kA2re`, `kA2hop`, `sA2hop` | during lengthening only: `kA2hop` moves high-strain `p2` heads back to low strain ("jump to the next site"), `kA2re` binds `PT` heads directly into `p2` at low strain. `kA2re` was 0 at baseline and is now ~30 — it is the compensator for the realignment mechanism's valley damage (§8) |

**Not active:** `UseD0State`, `UseR2Ceiling` (both falsified, [ApoPoolDetachment](../ApoPoolDetachment/conclusions.md)), `UseCatchBond`, `UseA2AttachmentShift`, `UseA2MechanicalRecocking`, `UseAtpDetachR2`, `UseAtpKmsrd` (low-ATP only), `SRDOutflowK`/`D0OutflowK` (fixed-dwell, falsified §9), **`UseMaxwellTensionOnly`** (§10b — should be on; worth ≈ −0.99 by itself and is not in this fit).

---

## 2. The fit, critically

### Slack protocol

![full slack](results/fig_slack_full.png)

Plateaus, the five descending steady levels, and the redevelopment from slack are all
essentially exact in both fits — that part of the model has been right for a while and
is not what is being tested.

![first cycle](results/fig_slack_cycle1.png)

The zoom is where the honest reading is. **Good:** from ~340 ms the reoptim (red)
tracks the data where the baseline (blue) sits visibly above it — that is the whole
point of the exercise, and `coolDownLS` confirms it (0.449 → 0.181). `rsK` goes 102.6 → 75.9
against a data 43.7 (×2.35 → ×1.74).

**Bad, and it is the same picture the features give:** the re-stretch spike is *worse*.
The model peaks at ~100 kPa where the data peaks at ~86, and the valley that follows is
deeper than the data. `peak1_y` 0.166 → 0.277, `vall2_dy` 0.081 → 0.371. **The model is
buying recovery kinetics partly by exaggerating the transient that precedes it** —
a bigger excursion recovering at the same absolute rate reads as a slower normalised
rate. That is a mechanism-shaped concern, not a cosmetic one, and it is the first thing
to check if `rsK` improves further.

### Force–velocity — the term that must not break

![force-velocity](results/fig_fv.png)

**This is the weakest part of the fit, the reoptim made it worse, and it is getting
worse with every round.** The model sits *below* the data at every non-zero velocity,
in both fits:

| v (µm/s) | data | baseline | reoptim final |
|---|---|---|---|
| 0.5 | 0.918 | 0.951 | 0.914 |
| **1** | **0.663** | **0.566** | **0.513** |
| **2** | **0.316** | **0.226** | **0.212** |
| 4 | 0.114 | 0.107 | 0.090 |

`FV_fnorm` **0.518 → 0.922**, the single largest degradation in the refit
(+0.404 of a −1.042 net) — and it **grew again over the last round** (round 6: 0.881
→ round 7: 0.922) while `rsK` fell 0.687 → 0.633. The optimiser is now buying `rsK`
with FV at a poor exchange rate.

**`FV_fnorm` has overtaken `rsK` as the largest single term in the objective**
(0.922 vs 0.633). At 1 µm/s the model gives 85 % of the measured force at baseline and
77 % after the refit. The shortfall at 1–2 µm/s is a **pre-existing defect** — the
baseline is already 15–28 % low — but the trend is the thing to act on: this is what a
longer run will keep doing unless FV is protected.

### Passive

![passive](results/fig_passive.png)

Peaks and decay shape are tracked; the data is very noisy (±3 kPa on a 0–20 kPa
signal). Two visible changes, both mild and both the wrong way: the reoptim's
re-stretch peaks are *lower* than baseline (16.6 vs 18.8 kPa on cycle 1, against a
noisy data peak ~21) and its inter-cycle holds sit ~1 kPa high where baseline sits at
~0.2. `PS_rampupRMSE` 0.148 → 0.285 and `PS_holdDecayRMSE` 0.187 → 0.247. The
series-viscoelastic parameters the refit used (`eta_M`, `kSE_M`, `mu_neg`) are shared
with the passive experiment, so this is the expected bill and it was correctly paid
rather than hidden.

### Feature ledger

![features](results/fig_features.png)

| improved | Δ | | degraded | Δ |
|---|---|---|---|---|
| **rsK** | **−1.685** | | **FV_fnorm** | **+0.404** |
| coolDownLS | −0.274 | | vall2_dy | +0.280 |
| peak1_dSL | −0.123 | | PS_rampupRMSE | +0.152 |
| ktr | −0.104 | | A | +0.123 |
| vall_y | −0.047 | | peak1_y | +0.110 |
| steady | −0.045 | | PS_holdDecayRMSE | +0.063 |
| peak2, vall_t | −0.013 | | doublePeak, t0_crossing | +0.106 |
| **sum** | **−2.291** | | **sum** | **+1.249** |
| | | | **TOTAL** | **6.103 → 5.060** |

`rsK` alone (−1.685) is larger than the whole net gain (−1.042). Strip `rsK` out and
the refit is a **net loss of 0.64** on everything else — i.e. the mechanism is not
improving the fit broadly, it is trading a large gain on one feature against a spread
of smaller losses, with FV the biggest.

---

## 3. What this does and does not establish

**Does:**
- A load-engaged attachment gate can slow the post-re-stretch recovery by 26 % while
  *improving* `ktr` and isometric force — something no previous candidate managed.
- It survives a joint refit that could have switched it off (`k_cr` stayed on, at 1.56).
- The net objective improves 15.5 % with FV and passive both scored.

**Does not:**
- `rsK` is still **×1.70** of data (74.6 vs 43.7). The defect is reduced, not solved.
- **The gain is one feature.** `rsK` −1.685 against a −1.042 net: everything else is a
  combined **+0.64**. This is not a broadly better fit, it is a concentrated trade.
- The refit is **7 rounds of a 29-parameter pool**, accepted in all 7 and stopped by
  the clock, not by convergence. A direction, not an optimum.
- `k_cr` has **no independent measurement or bound**. Until it is bounded in
  `parameterBounds.m` this is a working hypothesis with a free gain, not physiology.
- **FV is deteriorating monotonically** (0.518 → 0.881 → 0.922 over rounds 6→7) and is
  now the largest term in the objective. A longer run will keep spending it.

### Immediate actions, in order

1. **Protect FV before the next run.** Either raise its weight or add the missing
   1–2 µm/s force as an explicit constraint. As it stands the objective permits the
   trade that is being made, and the trade is getting worse per round.
2. **Set `UseMaxwellTensionOnly = 1`** (§10b) — worth ≈ −0.99 on its own, about the
   same as this entire 4 h refit bought, and it acts on `rsK`/`coolDownLS`, i.e. it
   should *reduce* the pressure to pay for `rsK` out of FV.
3. **Resolve `SRXFromR2HighStrain`** — currently carried at 0.2 with no evidence
   (never drawn in any round). Test it or zero it.
4. **Bound `k_cr`** in `parameterBounds.m`.
5. Only then run long (8–18 h).
