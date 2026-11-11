# How to refit — the objective audit that preceded the ktrfence campaign (2026-09-18)

Reproduce: `RunCrossScore.m` → `RunRefitStrategy.m` → `RunDeadbandAudit.m` →
`FixDeadbandArithmetic.m` → `BuildKtrFenceSeed.m` → `RunKtrFenceOptim.m`.
Outcome: `params/ktrfence_seed.m` (0.6997) and the `ktrfence` campaign.

---

## 0. The three questions, and their answers

| question | answer |
|---|---|
| Which lineage gets the budget? | **`lean_opt` + both ktr fixes.** Not a trade — the fixes are *free* on the lean objective and take the held-out ktr from 13.9× to 1.27× the data noise. |
| How far can the deadband narrow? | **To 0.25 before F-V detonates; 1/3 is the principled value and sits inside that.** The floor is not one number — it is s²/n with n the number of sources averaged into each term's target. |
| Score the ktr protocol? | **As an inert fence, `ktrP_k63\|fl:0.0174\|0.610`.** It pays nothing at the seed and only bites if `tau_reg` drifts back. The ktr *trace* stays fully held out. |

---

## 1. There were two lineages, and only one of them was floored

The premise going in was "we have a floored strategy at 1 SD". Half true. `params/estiff125_opt.m`
and its fast-gate descendant `params/fastgate_best.m` carry **no `fl:` tokens at all**. The floored
objective lives in a different line: `rigaware_a_opt` → `floored_seed` → `lean_seed` →
**`params/lean_opt.m`** (the MechanismAblation output).

| | fastgate_best | lean_opt |
|---|---|---|
| objective | 24 terms, no floors, `coolDownLS`, `rsK`, FV w = 125, single-prep 03/27 target | 13 floored terms @ 1 prep-width, `rsK80`, no `coolDownLS`, FV share 0.25, **pooled 3-prep target** |
| E_t under A | **0.9390** | 34.89 |
| E_t under B | 1.1481 | **0.5320** |
| held-out ktr | 1.66× noise | **13.94× noise** |
| mechanisms | A2Reattaching, MaxwellTensionOnly, DynamicPassive all on | all three deleted |

Each is far better under its own objective, so the columns are not comparable — but the last row is,
and it is damning for the lean line.

## 2. The ablation deleted the one flag the ktr protocol needs — for the same reason as last time

`lean_opt` has `UseMaxwellTensionOnly = 0`. The ablation measured deleting it at **+0.037, "free",
and `rsK80` 0.06 better**. That flag exists for exactly one purpose: to stop the titin dashpot
leaving a force floor across the ktr redevelopment window
([dPUdT_CombinedTransitions.m:308](../../Model/dPUdT_CombinedTransitions.m)). It measured as free
because **the ktr protocol was not in the objective that measured it** — the identical blind spot
documented in [`../KtrProtocolVsSlackKtr`](../KtrProtocolVsSlackKtr/conclusions.md), one campaign later.

Its real price, on data nobody was scoring: ktr `k63` **231.9 s⁻¹ against a data 43.2**, a trace
error **13.9× the measurement's own noise floor**.

### Restoring it is free, but only together with the fast gate

| variant | E_t (lean) | Δ | ktr k63 | ktr trace ÷ noise |
|---|---|---|---|---|
| L0 `lean_opt` as is | 0.5311 | — | 231.9 | 13.94× |
| L1 + `UseMaxwellTensionOnly` | 0.5324 | +0.0013 | 96.4 | 6.15× |
| L2 + fast gate (`tau_reg` 0.2 ms) | 0.6309 | +0.0999 | 105.2 | 8.43× |
| **L3 + both** | **0.5169** | **−0.0141** | **53.6** | **1.27×** |
| *data 03/27* | | | *49.4* | *1.00×* |

Either fix alone is mediocre — L2 costs 0.10. **Together they are better than the incumbent on its
own objective** and land the held-out protocol at 1.27× noise. The two defects had been propping
each other up: a dashpot that could not discharge left a standing force floor, and a 27 ms
availability memory left extra attachment drive, and the fit had balanced one against the other.

> **Trap, recorded.** `BuildLeanSeed` sets `x_M_slack = 0` when it deletes the mechanism, and
> `dPUdT_CombinedTransitions.m:322` *divides* by `x_M_slack`. Flipping the flag back on its own
> yields `Inf` with no error — it cost one wasted batch run here. Guard added at
> [getParams.m:641](../../Model/getParams.m).

---

## 3. The deadband: one number was never the right shape

Every data term carries `fl = FLOOR_CHI2 · nEl · s²` at `FLOOR_CHI2 = 1` — "being within one
between-preparation SD is free". At the L3 seed that leaves **10 of 16 data terms paying nothing**:

| feature | pooled | χ² | floor | prep-widths | pays |
|---|---|---|---|---|---|
| peak1_dSL | ✓ | 5.5100 | 1.3710 | 2.00 | 0.1604 |
| PS_steady20 | | 0.8183 | 0.2253 | 1.91 | 0.0280 |
| PS_steady22 | | 0.6385 | 0.2253 | 1.68 | 0.0195 |
| vall_t | ✓ | 4.2615 | 1.6324 | 1.62 | 0.0856 |
| peak1_y | ✓ | 0.2121 | 0.1354 | 1.25 | 0.0301 |
| vall2_dy | ✓ | 0.2394 | 0.2261 | 1.03 | 0.0031 |
| restretchSlopeStart | ✓ | 0.1558 | 0.1994 | 0.88 | — |
| t0_crossing | ✓ | 0.2753 | 0.4209 | 0.81 | — |
| vall_y | ✓ | 0.0452 | 0.1343 | 0.58 | — |
| rsK80[1:4] | ✓ | 0.0265 | 0.0997 | 0.52 | — |
| ktr | ✓ | 0.0254 | 0.1398 | 0.43 | — |
| PS_restretchPeak | | 0.1287 | 1.0564 | 0.35 | — |
| peak2 | ✓ | 0.0148 | 0.1239 | 0.35 | — |
| **FV_normpowerAvg** | | **0.0062** | **0.0515** | **0.35** | — |
| steady | ✓ | 0.0130 | 0.1149 | 0.34 | — |
| A | ✓ | 0.0069 | 0.1359 | 0.23 | — |

### Why it cannot simply be narrowed: F-V is a 1749× lever

For a term weighted `w = λ/s²` with `floor = nEl·s²`, the product `w × floor = λ·nEl` — **one
preparation-width costs the same everywhere**. That is the whole point of a chi-square. It holds for
every term here except one:

| term | w × floor | × median |
|---|---|---|
| **FV_normpowerAvg** | **92.97** | **1749×** |
| rsK80[1:4] | 0.0850 | 1.6× (the deliberate `RSK_BOOST = 2`) |
| the twelve 5-element slack terms | 0.0531 | 1.0× |
| PS_steady22 / PS_steady20 | 0.0106 | 0.2× |

`FV_SHARE` pins F-V's weight ~2187× above its precision-based value. **The deadband is the only
thing standing between the objective and that lever** — with no floor at all, F-V is **92.5% of the
objective**, which is precisely the rigaware failure the floor was invented to stop
([`../RecoveryOscillation`](../RecoveryOscillation/conclusions.md) §7). Uniform narrowing:

| FLOOR_CHI2 | E_t | n live | FV share | newly live |
|---|---|---|---|---|
| 1.000 | 0.5169 | 6 | 0 % | — |
| 0.500 | 0.6570 | 8 | 0 % | t0_crossing, restretchSlopeStart |
| **0.333** | **0.7138** | **9** | **0 %** | vall_y |
| 0.250 | 0.7480 | 10 | 0 % | rsK80 |
| 0.100 | 2.6505 | 15 | **68.8 %** | ktr, peak2, steady, PS_restretchPeak, **FV** |
| 0.000 | 12.0197 | 16 | **92.5 %** | — |

**The safe zone ends at 0.25.** F-V is inert down to ≈ 0.12 and then takes the objective over.

### The principled value, which happens to sit inside the safe zone

The floor was set to one *preparation* SD when the target was one preparation. It is not any more:
`UsePooledTarget = 1` replaces the slack features with the **mean of three** preparations
([pooledTargetFeatures.m](../../Model/pooledTargetFeatures.m)), and the sampling uncertainty on a
mean of three is `s/√3`. A model inside that is indistinguishable from a correct one, which is
exactly what a deadband should encode. So the floor is not one number — it is **s²/n, with n the
number of sources averaged into that term's target**:

| terms | target | n | FLOOR_CHI2 |
|---|---|---|---|
| the twelve pooled slack features + rsK80 | mean of 3 preparations | 3 | **1/3** |
| FV_normpowerAvg | mean of the 2 lab curves it averages | 2 | **1/2** |
| PS_restretchPeak, PS_steady22, PS_steady20 | single source (`pooledTargetFeatures` leaves them alone) | 1 | **1** |

This is a derivation, not a preference, and it defuses the F-V cliff without a special case: F-V
narrows only to 1/2, still comfortably inside its band at 0.35 prep-widths.

---

## 4. The ktr fence

`tau_reg` is in the optimiser pool, so a free refit can undo the fix. Rather than freeze it, the
protocol is fenced:

```
ktrP_k63|fl:0.0174|0.610       s = 0.132 (3-day spread 49.4 / 42.1 / 38.1)
                               floor = s² = one preparation-width
                               w = λ/s² = 0.610,  so w·fl = λ exactly
```

Commensurate with every other term by construction. The model sits **0.12 prep-widths** off, so the
term is **0.000000 at the seed** — a fence, not a fitting target. Cost: **+6 % per evaluation**
(20.6 s vs 19.5 s) for having `RunKtr = 1`.

**What this gives up, stated plainly.** The ktr *rate* is no longer a clean held-out quantity; it is
a guarded one. The ktr *trace* — `rmse` against the measurement's own noise floor, on all three days
— is untouched and remains the honest out-of-sample test. `validateKtrProtocol('<tag>')` reports it.

---

## 5. The seed

`params/ktrfence_seed.m`, **0.6997**, round-trips exactly (worst relative error 0). Composition:

| term | cost | share |
|---|---|---|
| peak1_dSL | 0.1959 | 28.0 % |
| vall_t | 0.1210 | 17.3 % |
| PS_holdDecayRMSE | 0.0992 | 14.2 % |
| PS_rampupRMSE | 0.0908 | 13.0 % |
| peak1_y | 0.0655 | 9.4 % |
| vall2_dy | 0.0385 | 5.5 % |
| PS_steady20 | 0.0280 | 4.0 % |
| restretchSlopeStart | 0.0238 | 3.4 % |
| PS_steady22 | 0.0195 | 2.8 % |
| t0_crossing | 0.0170 | 2.4 % |
| **ktrP_k63 (fence)** | **0.0000** | — |

The objective is now the **re-stretch transient shape** (peak1_dSL + vall_t + peak1_y + vall2_dy =
60 %) plus the **passive ramp** (28 %). That is the right target: it is what
[`../RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) identifies as the open structural
defect.

Held-out ktr at the seed: `k63` 53.6 / 54.2 / 53.9 against 49.4 / 42.1 / 38.1; ratio to post-slack
**1.04** (data 0.99); trace 1.27× / 1.60× / 2.62× noise.

**Pool: 77 parameters**, taken from `lean_seed.m` — `lean_opt.m`'s `leanPool` field was overwritten
with the *feature list* by the snapshot writer and is unusable. `x_M_slack` added back (live again)
and made **compulsory**: it is the only parameter whose current value is transplanted from a
different parametrisation and has never been fitted on this objective.

---

## 5b. The campaign result — objective down 23.7 %, held-out ktr up 40 %

`ktrfence`, 34 rounds, 14 improving snapshots, 8 h budget reached.
`RunFinalVerify.m`, `CheckTimecourses.m` → `figures/timecourse_check.png`.

**Quote the reloaded cost, not the header.** Campaign printed 0.5271; reloading
`params/ktrfence_opt.m` gives **0.5342** — the `%g` snapshot round-trip.

| | seed | opt |
|---|---|---|
| **objective** | 0.6997 | **0.5342 (−23.7 %)** |
| held-out ktr `k63` mean (data 43.2) | 53.9 | **58.4** |
| held-out ktr trace ÷ noise, 03/27 / 04/03 / 04/10 | 1.27 / 1.60 / 2.62 | **1.81 / 2.32 / 3.37** |
| ktr overstretch peak (data 0.724) | 0.925 = 1.28× | **0.970 = 1.34×** |
| fence term `ktrP_k63` | 0.000000 | **0.009653** |

### It bought exactly what it was pointed at

| bought | Δ | | spent | Δ |
|---|---|---|---|---|
| peak1_dSL | −0.0870 | | rsK80 | +0.0112 |
| peak1_y | −0.0589 | | restretchSlopeStart | +0.0099 |
| vall_t | −0.0324 | | **ktrP_k63** | +0.0097 |
| t0_crossing | −0.0153 | | PS_steady22 | +0.0035 |
| PS_holdDecayRMSE | −0.0047 | | ktr | +0.0026 |

The four terms bought are the **re-stretch transient shape** — the objective's whole point, and
[`../RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md)'s open structural defect.

### The fence detected the drift but did not stop it

`ktrP_k63` left its deadband (0 → 0.0097) and `ktr` left its own (0 → 0.0026). **Both redevelopment
rates rose together** — post-slack #2 51.6 → 56.0, ktr protocol 53.9 → 58.4 — against data of 50.0
and 43.2. The ratio held at **1.04** throughout, so the `tau_reg` fix survived exactly as designed:
what drifted is not the protocol split but the *global rate*, and the pooled `ktr` target falls
across the slack ladder (51.6 → 44.0) while the model rises.

At `w = 0.610` the fence is commensurate with the other terms by construction — which is precisely
why it could not outbid a 0.087 gain on `peak1_dSL`. **A fence sized as one term among twenty cannot
hold a quantity the rest of the objective is pushing on.** If the ktr rate is to be held, it needs
either a weight above the chi-square rule (i.e. an honest pin, argued as such) or `tau_reg` frozen.

### No bug — the state vector is clean

`CheckTimecourses.m`: probability conserved to **2.2e-16**, no negative-force excursions, all pools
in bounds (attached ≤ 0.332, A_reg ∈ [0.210, 0.995] snapping to A₀ between events, A_def 0.280
against a 0.90 ceiling), strain-grid edge mass 7e-4. The improvement is real, not a degenerate trace.

### But three features have the wrong SIGN of segment-dependence

| feature | model (opt) | pooled target |
|---|---|---|
| peak1_y | 82.1 → 52.0, monotone fall | 66.7, 72.6, **75.0**, 70.1, 51.9 — rises then falls |
| ktr | 50.1 → 53.0, rises | 51.6 → 44.0, falls |
| peak1_dSL | 0.031–0.034, flat | 0.007, 0.017, **0.030**, 0.025, 0.009 — an arc |

The *level* is fittable; the *shape across the slack ladder* is not. `peak1_dSL` is still the largest
term after the refit (0.109) and remains ~4× the target on segment 1. This is the structural
residual, and no further tuning on this pool will remove it.

### Parameters: 12 of 77 moved by more than 2 %

`kah` ×1.15, `k_cr` ×1.12, `k2` ×1.11, `ka` ×1.09, `v_ref_reg` ×1.09, `eta_M` ×0.94,
`L_thick` ×0.92, `x_M_slack` ×1.06, `P_bound_max` ×1.04, `kmsrd` ×1.03, `d_optimal` ×0.98,
`kamh` ×1.02. Two notes: **`x_M_slack` barely moved**, so the value transplanted from
`fastgate_best` was already near-optimal — the compulsory slot was cheap insurance, not a discovery.
And **`L_thick` moved 8.5 % to 1.17 µm**, already well below anatomical; it is acting as a fitting
device inside `f_lattice`, which the degeneracy analysis flagged and the lean pool did not remove.

### Verdict

Take it or not depending on what the next question is. **On the objective it is a clear win and the
gain is in the right place.** **On the one thing held out, it is a clear loss.** That is what
held-out data is for, and the honest reading is that this objective — with the re-stretch shape
weighted as it is and the rate terms deadbanded — prefers a globally faster model than the data
supports.

## 6. Left on the table

- **`PS_rampupRMSE` and `PS_holdDecayRMSE` have no deadband** and are 27 % of the seed cost. They
  are unbounded raw RMSE terms, not chi-squares — the same pathology the ablation removed
  `coolDownLS` for ("a guard term masquerading as a data term"). They need a noise model and a
  floor. Not changed here: one objective surgery at a time.
- **The `FV_SHARE` pin itself.** Left in place deliberately — with a floor it is unexploitable,
  which is the argument `BuildFlooredSeed` makes. But the objective is only a chi-square *because*
  the floor hides the one term that is not, and that should not be permanent.
- **`peak1_dSL` at 2.0 preparation-widths** is the single largest residual and the campaign's main
  target. `vall_t` at 1.6 is second. Both are re-stretch transient shape.

## 7. Process notes

- Two MATLAB sessions were lost to this analysis: the MCP session wedged on a script carrying the
  `x_M_slack` divide-by-zero, and burned **16,330 s of CPU** before being killed. It was starving
  the campaign's workers badly enough to trip `MaxRunTime = 150 s` on a chunk that takes ~15 s.
  **Check the machine is quiet before launching a long campaign** — the preflight
  `checkFeatureCoverage` call aborts the whole run if a chunk times out.
- `evaluateBakersExp`'s third output `featCost` is **already weighted**. Dividing by `w` again put
  F-V at 14.7 preparation-widths in the first pass of this audit when it is at 0.35.
