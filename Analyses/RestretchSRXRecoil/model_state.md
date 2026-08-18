# What the model is actually made of right now — and how well it fits

> **CURRENT STATE: `params/realign4_opt.m`, L2 = 3.8513** (baseline 7.0815 on the
> same objective, **−45.6 %**). Figures: [`PlotFinalReport.m`](PlotFinalReport.m) →
> `results/fin_*.png`. Sections 1–3 below describe the mechanism set (unchanged);
> §5 has the realign4 parameter values.

State at `params/realign3_opt.m`, 2-state, `@dPUdT_CombinedTransitions`, scored on
**FV + slack + passive** (`RunKtr = 0`, `RunStairs = 0`), objective = 25 feature
terms with `ovrsht_dy` dropped and **`FV_fnorm` at weight 30** (FV protection).

**Headline: L2 7.1392 → 5.0819, −28.8 %, with FV held to +0.9 %.**
Reproduce with [`PlotReportFigures.m`](PlotReportFigures.m).

| | total | FV term | rsK cost |
|---|---|---|---|
| baseline `rskR2_w025_opt` | 7.1392 | 1.5550 | 2.3175 |
| seed (realign2 + `UseMaxwellTensionOnly`) | 6.6219 | 2.7675 | — |
| **realign3 final** | **5.0819** | **1.5686** | **0.3191** |

Two things made this run different from realign2, which had spent FV monotonically
(0.518 → 0.881 → 0.922) to buy `rsK`:

1. **`UseMaxwellTensionOnly = 1`** (§10b) — worth −0.28 at this point on its own.
2. **The FV levers were put in the pool.** realign2's pool contained *no FV
   mechanism at all*: `A0`, `v_ref_reg`, `tau_reg` and the `R2(s)` shape knots were
   all absent, so the optimiser could only ever pay FV, never repair it. That was a
   setup error, not a property of the mechanism. With them available FV came back to
   +0.9 % of baseline.

And the gain is now **broad, not concentrated**. In realign2, removing `rsK` left a
net loss of 0.64 on everything else. Here the non-`rsK` features net to ≈ −0.06:

| improved | Δ | | degraded | Δ |
|---|---|---|---|---|
| **rsK** | **−1.998** | | PS_rampupRMSE | +0.190 |
| coolDownLS | −0.306 | | A | +0.138 |
| **ktr** | **−0.085** | | peak1_y | +0.062 |
| vall_y | −0.082 | | PS_holdDecayRMSE | +0.059 |
| peak1_dSL | −0.058 | | doublePeak | +0.030 |
| PS_restretchPeak | −0.028 | | vall2_dy | +0.030 |
| steady | −0.023 | | t0_crossing | +0.019 |
| | | | **FV_fnorm** | **+0.014** |

`ktr` is now *better* than baseline as well as `rsK`. What it still costs is the
passive ramp and `peak1_y` — the re-stretch peak still overshoots.

**Not converged.** The run stopped on its 3 h clock at round 16 having just made its
**largest** gain (5.4870 → 5.0693, −7.6 %) after three stalled rounds including a
failed ±10 % kick. Block-coordinate search on a 39-parameter pool clearly still had
unexplored directions. (An earlier round-13 crash — a parpool worker died and took
the pool down — cost nothing: the state file resumed exactly.)

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

---

## 4. The FV power metric (new)

**Why.** `FV_fnorm` is a normalised *force* curve whose first point is the isometric
force — fixed at 1 by construction and carrying no kinetic information. Power
(`force × velocity`) is the observable that actually distinguishes 8 from 2 mM ATP:
both the peak and the velocity at which it occurs move. It is also where the model's
residual defect lives.

**The two-dataset problem.** The FV reference comes from two different sources with
incompatible absolute forces (isometric 56.4 vs 67.4 kPa) but comparable *shapes*.
Each is therefore normalised by its **own** isometric force before power is formed,
the two are averaged, and their per-point variance becomes the weight.

New fields, emitted by `runFVExperiment` (data) and `extractForceVelocityAttributes`
(model): `FV_power1/2`, `FV_normpower1/2`, `FV_normpowerAvg`, `FV_normpowerVar`,
plus `FV_powerv` (the velocities they are on). On the model side there is one curve,
so the `*1`/`*2` fields are aliases of it.

**`v = 0` is excluded.** Power there is identically zero in every dataset and every
model, by definition. Included, it contributed exactly zero error while absorbing
43 % of the weight mass, making the weights on the informative points much weaker
than they looked.

| v (µm/s) | Baker2022 | Isovel2021 | avg | variance | **model** | error |
|---|---|---|---|---|---|---|
| 0.5 | 0.4593 | 0.4525 | 0.4559 | 2.3e-5 | 0.4923 | **+8.0 %** |
| 1 | 0.6639 | 0.6337 | 0.6488 | 4.6e-4 | 0.5746 | **−11.4 %** |
| 2 | 0.6313 | 0.5595 | 0.5954 | 2.6e-3 | 0.4699 | **−21.1 %** |
| 4 | 0.4443 | 0.3824 | 0.4133 | 1.9e-3 | 0.4047 | −2.1 % |

Inverse-variance weights (mean 1): **[2.889, 0.784, 0.139, 0.188]**. The two sources
agree closely at 0.5 µm/s and diverge most at 2 µm/s — which is also where the model
is worst. Weighting therefore *lowers* the term (0.0811 → 0.0371 at weight 1), which
is the correct behaviour: **do not chase a discrepancy the data itself cannot pin
down.**

### Vector weights in the objective

`evalFeatureCost` previously supported only a scalar weight per feature. Added a
`w:FIELD` token:

```
'FV_normpowerAvg|w:FV_normpowerVar|20'
```

`FIELD` is read from **feats_data** (a weight describes how well a measurement is
determined, so it must not depend on the model), treated as a variance, converted to
inverse-variance, floored at 10 % of its own mean (otherwise a point where the sources
agree exactly gets infinite weight), and renormalised to mean 1 so the scalar weight
keeps its usual meaning. Verified: identical to the unweighted path when the variance
is uniform, degenerate cases fall back to uniform, errors on a missing or
wrong-length field.

### Cost-neutral introduction

The term was added at weight 20 and the passive weights scaled by **0.4466** so the
total is unchanged at this parameter point — **5.0819 before, 5.0819 after** (Δ = 7e-7),
so optimisation can continue and any further drop is real progress rather than a
change of yardstick.

| term | weight before | weight after | cost before | cost after |
|---|---|---|---|---|
| `FV_normpowerAvg` | — | **20** (var-weighted) | — | **0.7421** |
| PS_steady20 | 0.5 | 0.2233 | 0.4205 | 0.1878 |
| PS_rampupRMSE | 0.005 | 0.002233 | 0.3385 | 0.1512 |
| PS_steady22 | 0.5 | 0.2233 | 0.2628 | 0.1174 |
| PS_holdDecayRMSE | 0.01 | 0.004466 | 0.2457 | 0.1097 |
| PS_restretchPeak | 1 | 0.4466 | 0.0735 | 0.0328 |
| **passive total** | | | **1.3410** | **0.5989** |

Passive drops from 26.4 % to 11.8 % of the objective. It is deliberately **not**
removed: `eta_M`, `kSE_M` and `mu_neg` are shared between the re-stretch transient and
the passive experiment, so with no passive term the dashpot would be unconstrained.

⚠️ The neutrality is exact **only at this parameter point**. As the optimiser moves,
the two groups scale differently — that is intended; what matters is that the starting
value is unchanged.

---

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

---

## 5. realign4 — the current fit, mechanisms and parameter values

`params/realign4_opt.m`, scored on the **power objective** (FV + slack + passive;
`ovrsht_dy` and `FV_fnorm` dropped, `FV_normpowerAvg|w:FV_normpowerVar|62.2779`,
passive weights ×0.446621).

| snapshot | L2 | FV power term | rsK cost | ktr cost |
|---|---|---|---|---|
| baseline `rskR2_w025_opt` | 7.0815 | 2.1142 | 2.3175 | 0.3206 |
| realign3 | 5.0819 | 2.3107 | 0.3191 | 0.2357 |
| **realign4 (current)** | **3.8513** | **1.3420** | 0.3810 | **0.1430** |

*(realign3 scores 5.0819 on **both** objectives — that is the cost-neutral swap
working exactly as designed.)*

**Note realign3 was slightly WORSE than baseline on power** (2.3107 vs 2.1142) even
though it was much better overall: it was optimising `FV_fnorm`, and the force-normalised
curve and the power curve are not the same target. Changing the metric exposed that.

### Power–velocity, the fitted quantity

| v (µm/s) | data (mean of 2 sources) | weight | baseline | realign3 | **realign4** |
|---|---|---|---|---|---|
| 0.5 | 0.4559 | 2.89 | 0.4759 | 0.4923 | **0.4870** |
| 1 | 0.6488 | 0.78 | 0.5660 | 0.5746 | **0.5986** (−7.7 %) |
| 2 | 0.5954 | 0.14 | 0.4539 | 0.4699 | **0.5011** (−15.8 %) |
| 4 | 0.4133 | 0.19 | 0.3849 | 0.4047 | **0.4150** (+0.4 %) |

The peak is still too low and still slightly too fast: the model's power maximum sits
at v = 1 with 0.599 against a measured 0.649, and the model falls away faster than the
data through 2 µm/s. **This is the remaining structural defect** — and note it is now
the largest single term in the objective (1.342 of 3.851, 35 %).

### Other observables

| observable | data | baseline | realign3 | **realign4** |
|---|---|---|---|---|
| `rsK` (s⁻¹) | 43.74 | 102.6 | 65.3 | **67.1** |
| `ktr` (s⁻¹) | 49.21 | 52.49 | 51.09 | **50.49** |
| steady (kPa) | 77.35 | 76.19 | 76.31 | **76.99** |
| `peak1_y` (kPa) | 89.28 | 93.57 | 95.09 | **93.26** |
| `vall_y` (kPa) | 71.11 | 67.56 | 69.23 | **69.50** |
| `peak2` (kPa) | 78.42 | 77.32 | 79.48 | 80.16 |

Improved over realign3: `ktr` (0.236→0.143), `peak1_y` (0.227→0.114), `steady`
(0.056→0.012), `A` (0.224→0.103), `vall_y`, `doublePeak`.
Worse: `peak1_dSL` (0.241→0.328), `restretchSlopeStart` (0.072→0.115), `peak2`,
and `rsK` drifted back slightly (0.319→0.381).

### Active mechanisms

Unchanged from §1 — no mechanism was added or removed, only retuned:

`UseSuperRelaxed` + `UseSuperRelaxedADP` (SRX/SRXD, one-way force-gated trap) ·
`UsePieceWiseStrainDep` (spline `R12`/`R21`/`R1D`/`R2`) · `UseSpaceExtension` ·
`UseOverlap` + `UseOverlapFactor` + `UseCalculatedN` · `UseLatticeSpacing` ·
`UseGlobalOccupancySaturation` · `UseRegistrationAvailability` +
`RegAvailShorteningOnly` · **`UseLoadRealign`** (compliant realignment, §8) ·
`UseSerialStiffness` · `UsePassive` + `UseDynamicPassive` · `UseMaxwellDashpot` +
**`UseMaxwellTensionOnly`** · `UseNegativeKstiff` · `UseA2Reattaching`.

### Parameters that moved from baseline (|ln ratio| > 0.02)

| parameter | baseline | realign4 | ratio | what it is |
|---|---|---|---|---|
| `k_cr` | 3 | **1.519** | 0.51 | realignment gain (seeded 6, settled ~1.5 across three refits) |
| `eta_M` | 0.02700 | **0.01561** | 0.58 | Maxwell decay time (s) |
| `sigma2` | 41.37 | **28.64** | 0.69 | SRX entry force sensitivity |
| `sigma_srd1` | 26.70 | **34.12** | 1.28 | SRX exit force sensitivity |
| `v_ref_reg` | 0.2659 | **0.3381** | 1.27 | FV shoulder half-speed |
| `kSE_M` | 77.01 | 95.37 | 1.24 | Maxwell charging stiffness |
| `kmsrd` | 19.78 | 23.50 | 1.19 | SRXD → PD rate |
| `ka` | 216.2 | 253.3 | 1.17 | attachment rate |
| `kah` | 101.6 | 118.7 | 1.17 | hydrolysis / priming rate |
| `mu` | 0.01564 | 0.01786 | 1.14 | series dashpot |
| `k2` | 100.5 | 111.0 | 1.10 | master detachment rate |
| `kSE` | 4912 | 5392 | 1.10 | series elastic stiffness |
| `k_1` | 33.16 | 30.17 | 0.91 | p2 → p1 reversal |
| `sSRXrip` | 0.0150 | 0.01366 | 0.91 | strain threshold of the SRX routing |
| `sA2hop` | 0.005126 | 0.004750 | 0.93 | A2 hop strain threshold |
| `mu_neg` | 0.01544 | 0.01452 | 0.94 | negative-velocity dashpot |
| `F_cr` | 25 | 23.60 | 0.94 | realignment force half-point (kPa) |
| `tau_cr` | 0.0200 | 0.02125 | 1.06 | realignment relaxation (s) |
| `kA2hop`, `kstiff2_n`, `ksr0`, `x_M_slack`, `kstiff1_n` | | | 1.03–1.05 | minor |
| `ekSE`, `kstiff2`, `R_thick` | | | 0.96–0.97 | minor |

Notable: the two SRX force sensitivities moved in **opposite** directions —
`sigma2` down 31 % (entry more force-sensitive) and `sigma_srd1` up 28 % (exit less
force-sensitive). Both changes *weaken* the auto-disengagement described in §4, which
is the structural reason post-restretch was the fastest window. The optimiser found
that on its own.

`SRXFromR2HighStrain` = **0.2094**, essentially untouched from the 0.2 it was seeded
at across three campaigns. It has still never been shown to earn its place.

### ⚠️ How this run ended

The 12 h budget **did not apply**. `optimizeFeatures` checks elapsed time only
*between* rounds, and round 13 entered a region where every evaluation hit
`MaxRunTime` and returned the 1e6 failure penalty. It stayed there for **four days**
before I killed it. The best snapshot (#10, L2 3.8449 internal / 3.8513 re-scored) was
written on day one and is unaffected.

**Fix needed:** move the budget check inside the simplex callback, or add a
consecutive-penalty circuit breaker that abandons a round once N evaluations in a row
return ≥1e6.
