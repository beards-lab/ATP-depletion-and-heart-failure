# Restretch feature fit — conclusions

**Question.** Improve the feature-cost match (evalFeatureCost) on the slack-restretch
**restretch features** (`Am`, `peak1_y`, `peak1_dSL`, `vall_y`, `peak2`, `vall2_dy`)
without degrading force-velocity or the physiology output bounds. Base =
`params/params_reseeded_regavail_opt2.m` (8 mM, protocol_03_27_2026).

## Headline result

| Config | change | feature cost |
|---|---|---|
| Baseline | symmetric registration availability (current best-fit) | **15.86** |
| E01 | **directional availability** (shortening-only) | 11.38 (−28%) |
| C3 | R1D⁺ clamp↓, k2↓, kSE↑, kstiff2↓ | 9.93 |
| D1 | + kstiff1↓ | 8.23 |
| D2 | + v_ref_reg↓ (FV shoulder width) | 7.34 (−54%) |
| J1/K1 | + kSE↑, kstiff↑ (dSL, force level) | 6.23 |
| Q2 | + **A0↓ with ka-compensation** (FV shoulder depth) | 4.91 |
| **S3** | + A0→0.42, kstiff fine-trim | **4.56 (−71%)** |

Best config: `results/params_best.m` (= `params/params_restretch_best.m`, S3).
At S3: `steady` [80 80 81 81 64] = data, `vall2_dy` [−13..−11] vs [−13..−10] = data, `peak2`
seg1 81 vs 80, `Am`/`restretchSlopeStart` in range, and **FV cost 2.09→0.44** (fnorm now
[1 .85 .60 .34 .13] vs data [1 .92 .66 .32 .11]). Residual: `vall_y` (2-state wall, 0.97),
physiology bounds (0.74), `peak1_dSL` (0.69), t0 (0.64), `peak1_y` (~105 vs 96, residual
force/peak coupling, 0.54). Prior milestones: `results/params_D2_best2state.m` (7.34).

## Three mechanistic findings

**1. Registration availability must be shortening-directional.** The mechanism added for
the FV shoulder drove attachment off `abs(vel)`, so a re-stretch (lengthening) transiently
*boosted* attachment — inflating `peak2` (+18 kPa) and erasing the post-restretch undershoot.
Changing `A_inf` to track `max(0,-vel)` (new flag `RegAvailShorteningOnly`,
`Model/dPUdT_CombinedTransitions.m` + `Model/getParams.m`) restored the undershoot and halved
the peak2 excess **while leaving the FV shoulder untouched** (shortening branch unchanged).
This single change is the largest single step and is the recommended default for this config.

**1b. FV shoulder is fixable via A0-depth + attachment compensation (not stiffness).** FV
steepness (initially the dominant residual, cost ~2.0) is NOT structural. Deepening the
isometric availability floor `A0` (0.6→0.42) widens the velocity-driven availability swing and
flattens FV, but depresses isometric force. The compensation lever is decisive: compensating
with `kstiff` restores force but inflates the P1 `peak1_y` and overshoots — whereas
compensating with **`ka` (attachment rate)** restores force by attaching *more* bridges, not
stiffer ones, so peaks stay put. `A0↓ + ka↑` cut FV cost 2.09→0.44 with `peak1_y`/`peak2` held.

**2. `vall_y` vs `vall2_dy` is a 2-state wall.** The fast shallow post-peak valley (`vall_y`,
data only ~4 kPa below peak2) and the deep slow undershoot (`vall2_dy`, ~13 kPa) are both set
by the *same* positive-strain R1D detachment (confirmed via both R1D value and knot-position
sweeps). Any change that lifts `vall_y` to data (≈77) simultaneously flattens `vall2` to ≈−4
and inflates the peaks — cost rises to ~24. A single strain-dependent rate cannot separate a
fast-shallow from a slow-deep response; they are timescale-separated. **The durable fix is a
3rd, low-force slow-detaching cross-bridge state** (`NumberOfStates=3`), which also relieves
the ktr iron law (ktr was dropped from the cost as unfittable in 2-state).

## Lever reference (validated sensitivities)

- **Post-restretch undershoot depth** (`vall2_dy`): `k2`↓ deepens it; positive-strain R1D
  (`PieceWiseStrainDepR1DParams__4` clamp, `__3`) shapes it.
- **Fast valley** (`vall_y`): positive-strain R1D — but coupled oppositely to `vall2` (the wall).
- **Peaks** (`peak1_y`, `peak2`): `kstiff1`/`kstiff2`↓; occupancy ceiling `P_bound_max`↓.
- **Peak length excursion** (`peak1_dSL`): `kSE`↑ (series stiffness); still ~30% high at D2.
- **FV shoulder** (`FV_fnorm`): `v_ref_reg`↓ sets shoulder WIDTH cheaply; `A0`↓ sets shoulder
  DEPTH but must be paired with **`ka`↑ force compensation** (see finding 1b) — do NOT compensate
  with `kstiff` (inflates peaks). R2 negative knots move FV counter to naive expectation and are
  coupled to the slack phase; leave R2 to the optimizer.

## Recommendations / next steps

1. `RegAvailShorteningOnly=true` is now the default in `params/params_restretch_best.m` and is the
   seed for `Workbench/RunSurrogateSimplex_Optim.m` (its pool includes A0, v_ref_reg, ka, kstiff1/2,
   k2, kSE, R1D knots). Run it on the grid to polish the S3 basin.
2. Remaining genuine misses to target: `peak1_dSL` (~0.032 vs 0.027), t0_crossing, `peak1_y` (residual
   force/peak coupling). `vall_y` is the 2-state wall — do not over-chase.
3. Physiology bounds (XTOR_vmax ~18 vs 15; SRX_ss 0.07 vs 0.10): candidates for JUSTIFIED relaxation
   (cardiac Fenn factor can exceed 3×; biochemical SRX at max activation can be <0.10, Ma 2022) —
   or address SRX via SR kinetics (ksr0/kmsr), not SRXT_0 (initial condition only, no steady-state effect).
4. Structural: `NumberOfStates=3` to break the `vall_y`/`vall2` wall and the ktr iron law together.

Reproduce: `cd(root); addpath(genpath('.'))`, load `results/params_best.m` (S3), set the fn
(see labdiary/`showFit.m`), and run `Analyses/RestretchFeatureFit/showFit.m`.
