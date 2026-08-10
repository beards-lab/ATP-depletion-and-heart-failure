# Fitting the force redevelopment after re-stretch (2026-08-04)

**Question.** The worst part of the fit is the force redevelopment *after the
re-stretch is done*, compared against `ktr` and the post-slack redevelopment.
Would a **first-order fit** be a better metric than the least-squares term
currently in the objective, and what does it hint about the model?

**Short answer.** Yes — decisively. The rate metric turns a defect that was worth
1 % of the objective into one worth 78 %, and once it is measurable the defect
resolves into two halves with two different levers: **`k2` sets the level, the
Maxwell dashpot `eta_M` removes the amplitude dependence**, and together they take
the post-restretch rate from ×2.70 to **×1.07** of data. The bill is large
(feature total 28 → 176, mostly through `k2`'s force and `ktr` couplings) and
`kstiff2`↓ pays about half of it; a joint refit is running to see how much of the
rest it can find. The one thing still without a mechanism is the data's mild
*negative* amplitude dependence.

Base: `params/optfull7_opt_mov.m`, protocol 03/27 8 mM.
Reproduce, in order: [`Baseline.m`](Baseline.m) → [`RunAmplitudeTest.m`](RunAmplitudeTest.m)
→ [`RunAmplitudeVsData.m`](RunAmplitudeVsData.m) → [`RunR2CeilingSlope.m`](RunR2CeilingSlope.m)
→ **[`RunSlopeUnderLevers.m`](RunSlopeUnderLevers.m)** and **[`RunComboCost.m`](RunComboCost.m)**
(the two that carry the main result, §6). Full file list in §8; step-by-step log
in [labdiary.md](labdiary.md).

⚠️ Run these **one at a time**. Every evaluation is guarded by a wall-clock
`MaxRunTime`, so two heavy MATLAB jobs on this box make healthy parameter points
score the ~1e6 timeout penalty — which is how the first overnight launch spent an
hour optimising a failure. See [`TimeBattery.m`](TimeBattery.m).

---

## 1. The LSQE really was blind, and structurally so

The post-restretch window is ~280 ms. The rate information lives in its first
~30 ms; the remaining ~250 ms is plateau the model already gets right. A
trace-wise least-squares term therefore averages the defect away — this is not
a weight mistake but a property of the metric:

| | at the seed |
|---|---|
| `coolDownLS` (the LSQE over that window) | **0.32 of 28.18 = 1.1 %** of the feature cost |
| `rsK` (the new first-order rate) | **14.91 of 19.18 = 78 %** of what the optimiser minimises |

## 2. The defect is a *rate* defect, and nothing else

New per-cycle features in [`Model/extractSlackAttributes.m`](../../Model/extractSlackAttributes.m):
`rsK` (rate), `rsA` (amplitude), `rsT0`, `rsR2`, `rsK63`. Baseline vs data:

| observable | model / data |
|---|---|
| **post-restretch rate `rsK`** | **×2.70** |
| post-restretch amplitude `rsA` | ×1.11 |
| post-slack rate `ktr` | ×1.05 |
| isometric force | ×0.96 |

Everything else in that window is right. That is a far sharper statement than
"the restretch shape is off", and it is visible only because the metric is a rate.

The estimator is deliberately amplitude-fixed (only `k` and `t0` free), which
makes it span-insensitive — refitting over 30 ms vs the full window moves the
data value from 48.5 to 48.4 s⁻¹ — and keeps a force-scale error out of the rate.

![Recovery windows, data vs model](results/baseline_windows.png)

*Top row: post-restretch. The model (blue) is done by ~30 ms where the data
(black) is still climbing at 80 ms. Bottom row: post-slack, where the same model
tracks the data to within 10–20 %.*

## 3. The gap is half rate error, half nonlinearity

Shrinking the perturbation (synthetic release/re-stretch of matched depth) makes
the model's post-restretch rate fall toward its small-signal limit:

| depth (ML) | 0.080 | 0.060 | 0.040 | 0.030 | 0.020 | 0.015 | 0.010 |
|---|---|---|---|---|---|---|---|
| model `rsK` | 137.2 | 109.2 | 102.7 | 92.6 | 93.2 | 87.4 | 79.8 |

| | | | fixed by |
|---|---|---|---|
| model small-signal limit vs data | ~74–80 vs 44 s⁻¹ | **×1.8 — plain rate error** | `k2`↓ (§6) |
| large-stretch speed-up | 80 → 137 s⁻¹ | **×1.7 — nonlinearity** | `eta_M`↓ (§6) |

Roughly half the defect is that the model's *intrinsic* relaxation at its
operating point is too fast; the other half only appears when you stretch hard.
This split turned out to be the useful one — §6 shows each half has its own lever.

## 4. The decisive constraint: the data has no amplitude dependence

The protocol already sweeps re-stretch amplitude (cycles 1–5 re-stretch by
0.080, 0.101, 0.121, 0.141, 0.060 ML), so this is free on every protocol day:

| series | n | slope d`rsK`/d(amplitude), s⁻¹/ML | r |
|---|---|---|---|
| data 03/27 8 mM | 5 | −94 | −0.64 |
| data 04/03 8 mM | 5 | −95 | −0.72 |
| data 04/10 8 mM | 5 | −188 | −0.45 |
| **data pooled 8 mM** | 15 | **−126** | −0.31 |
| **model, synthetic sweep** | 7 | **+714** | **+0.97** |

**The model gets faster the harder you stretch it; the real muscle, if anything,
gets slower.** The sign is wrong, not merely the magnitude. This is a *shape*
constraint: it survives any per-prep force scaling and cannot be re-weighted away.

![Model vs data amplitude dependence](results/amplitude_vs_data.png)

*Left: the model (red) climbs with re-stretch size; every data point (circles,
three preps) sits on a mildly falling trend. Right: mean `rsK` per prep — note
the 2 mM bars at roughly a third of 8 mM.*

*Caveat.* In the synthetic sweep, release depth and re-stretch amplitude move
together; in the real protocol they are decoupled (cycle 5 releases 0.161 ML but
re-stretches 0.060). The synthetic sweep is the one that isolates a single knob,
and the model's response to it is near-deterministic (r = 0.97); the data slope
is negative on all three days independently.

## 5. R2-capping is falsified again, on a new observable

The natural suspect is the unbounded strain acceleration of the strong-bridge
detachment `R2` (2100–4250 s⁻¹ at re-stretch strains). `UseR2Ceiling` saturates it.

| `k2max` | slope | `rsK`@0.08 | undershoot A@0.08 (kPa) | F_iso |
|---|---|---|---|---|
| off | 714 | 137.2 | 12.95 | 76.3 |
| 1500 | 483 | 99.1 | 7.77 | 77.6 |
| 800 | 423 | 92.5 | 5.57 | 78.9 |
| 400 | 977 | 135.4 | 3.09 | 81.5 |
| 200 | 1607 | 209.1 | 0.93 | 85.3 |
| *data* | *−126* | *~52* | *~13* | *~80* |

![R2 ceiling vs the amplitude slope](results/r2ceiling_slope.png)

**Rejected.** The slope improves marginally then gets much worse, because the
ceiling keeps bridges attached through the stretch and *destroys the undershoot
itself* (12.95 → 0.93 kPa). It removes the phenomenon rather than slowing it.
`rsK` at small amplitude also stays 73–80 s⁻¹ at every setting, so it does not
touch the ×1.8 linear half either. This is a **second, independent falsification**
of R2-capping — [`ApoPoolDetachment`](../ApoPoolDetachment/conclusions.md) killed
it on the ktr overstretch peak; this kills it on the amplitude slope.

## 6. Two levers, and they fix both halves — the main result

The one-at-a-time scan (~100 evaluations) found nothing that reaches the needed
`d ln rsK` = **−0.995** cheaply, and I initially concluded the amplitude half was
structural. **That was wrong, and the combination test says so.** Two levers act
on the two halves *separately*:

| variant (synthetic amplitude sweep) | slope | `rsK`@0.080 | `rsK`@0.010 |
|---|---|---|---|
| baseline | **+714** | 137.2 | 79.8 |
| `eta_M`×0.3 | +83 | 80.6 | 70.3 |
| `k2`×0.6 | +243 | 65.4 | **44.7** |
| **`k2`×0.6 + `eta_M`×0.3** | **+21** | **48.7** | **43.7** |
| `k2`×0.7 + `eta_M`×0.3 | **−32** | 59.7 | 62.0 |
| *data* | *−126* | *~52* | — |

- **`k2` (master cycling rate) sets the LEVEL** — it alone brings the
  small-signal rate from 79.8 to 44.7 s⁻¹, on top of the data.
- **`eta_M` (the Maxwell dashpot viscosity) removes the AMPLITUDE DEPENDENCE** —
  slope 714 → 83 on its own, and it does so while leaving `ktr` and isometric
  force *completely untouched*.

The second point generalises to the whole **series-viscoelastic class**: `eta_M`,
`kSE_M`, `mu_neg`, `mu`, `kstiff1_n` and `kstiff2_n` all move `rsK` while leaving
`ktr` at 51.5–52.5 and `steady` at exactly 76.2 (the scan's `ktr`/`steady`
columns). `eta_M` is simply the strongest of them (−0.337 in `d ln rsK` vs
−0.247 for `kSE_M`, −0.110 for `mu_neg`, −0.103 for `kstiff1_n`). **This class is
the only route found in this repo that moves post-restretch kinetics without
paying the cross-bridge iron law** — and it is worth remembering that it is
*mechanical*, not kinetic: the natural instinct on a rate defect is to reach for
a rate constant, and every rate constant here is welded to force or `ktr`.

Together they give **`rsK` ≈ 49 s⁻¹ at every amplitude with slope +21** — the
level *and* the flatness. On the real protocol with the full battery,
`k2`×0.6 + `eta_M`×0.3 gives **`rsK` = 46.8, ×1.07 of data** (from ×2.70).

![Slope under levers](results/slope_under_levers.png)

### The bill, honestly

That combination costs **feature total 28.2 → 175.9** on the full battery, because
`k2`↓ carries its usual couplings: isometric force 76 → 92 kPa and `ktr` 52 → 39
(×0.79 of data). The damage is concentrated in `A` (+53), `steady` (+45),
`peak1_y` (+15), `ktr` (+7.7), `FV_fnorm` (+6.1).

Compensation direction matters and is *not* the obvious one:

| variant | `rsK` ratio | feature total |
|---|---|---|
| baseline | ×2.70 | 28.2 |
| `k2`.6 + `etaM`.3 | ×1.07 | 175.9 |
| + `ka`×1.4 | ×1.05 | **367.6** ← wrong way |
| + **`kstiff2`×0.85** | ×1.09 | **89.9** ← halves the bill |
| `k2`.7 + `etaM`.3 | ×1.26 | 124.7 |

`ka`↑ makes it much worse (it raises force further); `kstiff2`↓ pays for k2's
force excess and halves the cost. A 76-parameter joint refit has many more knobs
than these three, which is what the overnight run is for — and it now has a real
target to aim at rather than a hopeless one.

**`eta_M` alone is the cheap win:** ×2.70 → ×1.58 for a feature cost of +9.8,
localised in `ovrsht_dy` (+3.6), `vall2_dy` (+2.2), `vall_y` (+2.0) and
`PS_restretchPeak` (+1.2) — i.e. it does send part of its bill to the passive
experiment, as a series-element change should.

### What remains structural

Neither lever makes the slope *negative*; the best is +21 (`k2`.6+`etaM`.3) or
−32 at a level that is then too high (`k2`.7+`etaM`.3, `rsK` 60). So the residual
question is unchanged in kind, just much smaller: the data's mild *negative*
amplitude dependence still has no mechanism in this model. The candidate the
literature supports is **compliant realignment** (Daniel 1998; Tanner, Daniel &
Regnier 2007) — bound bridges set the register of neighbouring sites through
filament compliance, so stripping a large bound population *while the filament
stays loaded* leaves the register disturbed and re-equilibration takes time.
That is driven by **bound mass disrupted under load**: direction- and
strain-agnostic, and naturally larger after a bigger re-stretch, hence negative.
[`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) §C arrived
at the same candidate independently. This analysis supplies its target:
**slope ≈ −126 s⁻¹/ML at `rsK` ≈ 44–52 s⁻¹.**

## 6b. REVISION (2026-08-10) — the amplitude dependence was the Maxwell dashpot,
## not missing physiology

*Assessed at `params/rskR2_w025_opt.m` (the optimiser's point after the overnight
runs). Scripts: [`AssessSnapshot.m`](AssessSnapshot.m),
[`RunRecoveryAnatomy.m`](RunRecoveryAnatomy.m),
[`RunMechanismProbe.m`](RunMechanismProbe.m),
[`RunAmpNoMaxwell.m`](RunAmpNoMaxwell.m).*

**§6 above said the model's positive amplitude dependence was structural and needed
a new stretch-engaged mechanism. That was wrong.** Every amplitude sweep behind §6
ran with the Maxwell dashpot on. Switching it off:

| synthetic depth sweep | slope (s⁻¹/ML) | r | level |
|---|---|---|---|
| Maxwell **ON** | **+276** | 0.58 | 76–115 |
| Maxwell **OFF** | **−52** | −0.61 (flat) | 70–77 |
| *data* | *−94* | — | *~44* |

With the dashpot off the cross-bridge machinery already delivers a **nearly
amplitude-independent** rate — qualitatively the behaviour the data shows. The
`+714` slope quoted in §6 was a property of the series-viscoelastic element.

### Anatomy: what actually carries the recovery

Fraction of each carrier's own excursion completed when **force** is 63 % recovered
(1.00 = finished early, not rate-limiting; ~0.63 = tracks force):

| window | Force | F_active | F_passive | LSE | bound (p1+p2) | t63 |
|---|---|---|---|---|---|---|
| postRestretch c1 | 0.66 | 0.57 | 0.41 | 0.65 | 0.75 | 13.8 ms |
| postSlack c1 | 0.64 | 0.65 | **1.00** | 0.94 | **0.28** | 21.7 ms |

And the composition of the recovered force:

| window | ΔForce | ΔF_active | ΔF_passive |
|---|---|---|---|
| postRestretch | **+14.3** | +22.2 | **−7.9** |
| postSlack | +67.6 | +70.2 | −2.6 |

**The post-restretch net recovery is the small difference of two large opposing
terms** — a +22 kPa active rise against a −8 kPa passive decay. The post-slack
window has no such component (−2.6 kPa against +70). The two windows are
rate-limited by different things *in the model*: post-slack by strain accumulation
on a slowly-growing population (`bound` lags at 0.28), post-restretch by the bound
population itself (0.75) blended with a decaying dashpot.

### Revised decomposition of the gap

| | `rsK` | vs data |
|---|---|---|
| model as-is | 102.6 | ×2.35 |
| **Maxwell dashpot OFF** | 66.9 | **×1.53** |
| data | 43.7 | ×1.00 |

Roughly **half the excess is the dashpot** (ln 102.6/66.9 = 0.43 of ln 102.6/43.7
= 0.85) and half is genuinely fast cross-bridge kinetics. Turning it off is not the
fix — feature total 6.5 → 18.3, because the element is load-bearing for
`ovrsht_dy`, `vall2_dy` and the passive features. It has to be made *right*.

### The residual protocol dependence, correctly sized

| | ktr | rsK | rsK/ktr |
|---|---|---|---|
| data | ~49 | 43.7 | **0.89** |
| model, Maxwell ON | 52.5 | 102.6 | 1.95 |
| model, Maxwell OFF | 52.2 | 66.9 | **1.28** |

So a real but much smaller structural residual: the model is ~1.4× too
protocol-dependent after the dashpot is accounted for, not the ~2× §6 implied.

### Negative results (do not re-run these)

| probe | `rsK` | reading |
|---|---|---|
| `kah`×0.4 (slower priming) | 109.7 (**up**) | the existing hydrolysis step does **not** gate this window — fast alternative paths exist |
| `ka`×0.5 | 107.9, force → 52.8 kPa | attachment rate is not the lever, and it wrecks force |
| `UseMaxwellTensionOnly=1` | 103.5, featL2 6.5→8.0 | no help *here*; its benefit in RestretchVsKtrRecovery Part 2 was specific to the **ktr protocol**, which this objective does not run |

The `kah` result matters for mechanism design: **slowing an existing rate does not
clamp the recovery, because the model always has a faster parallel route.** Only an
*obligatory in-series* step can.

## 7. Unlooked-for: `rsK` is a strong ATP signal

Extracting `rsK` on every protocol file gave the 2 mM / 8 mM comparison for free:

| prep | 8 mM | 2 mM | ratio |
|---|---|---|---|
| 03/27 | 43.7 | 15.9 | **0.36** |
| 04/03 | 51.5 | 20.4 | **0.40** |
| 04/10 | 65.1 | 18.0 | **0.28** |

Low ATP slows the post-restretch recovery **~3×**, consistent across three
protocol days — far more than it slows `ktr` (×0.55). The 2 mM fits are also
*cleaner* first-order (`rsR2` 0.84–0.94 vs 0.51–0.80). This is a sharp new
constraint for the low-ATP work, recorded but not pursued here.

---

## ⚠️ Two issues found in passing, not fixed

1. **The optimiser minimises L2, the reports show L1.** `evaluateBakersExp`
   calls `evalFeatureCost` without `costExp`, so it defaults to 2, while
   `RunBakersExp` passes 1. Features whose normalised error exceeds 1 are
   amplified in the objective relative to how they are reported. Worth making
   explicit either way.
2. **The 03/27 data features were extracted by an older extractor than the one
   the model is scored through.** Re-extracting today reproduces the 04/10 files
   exactly, but not 03/27 (`peak1_dSL` 97 %, `peak1_t` 92 %, `vall_t` 96 %,
   `peak1_y` 26 %) or 04/03 (`ovrsht_dy` 91 %, `vall2_t` 53 %). `peak1_y` and
   `peak1_dSL` carry weights 10 and 1 in the current objective. Left untouched
   because correcting it moves the fit target and invalidates the current best
   parameters — but it means part of the current cost compares model and data
   through different conventions. **Needs a decision.**

## 8. Files

| script | what it does |
|---|---|
| [`Baseline.m`](Baseline.m) | reference run + feature-cost breakdown + rate tables |
| [`PlotWindows.m`](PlotWindows.m) | window overlays, span-insensitivity check |
| [`AddRestretchFeatsToData.m`](AddRestretchFeatsToData.m) | writes `rsK…` into the protocol `.mat` files, refusing to move existing fields |
| [`evalRs.m`](evalRs.m) | one evaluation → feature costs + rates (shared by the sweeps) |
| [`RunAmplitudeTest.m`](RunAmplitudeTest.m) | §3 — rate vs perturbation size, small-signal limit |
| [`RunAmplitudeVsData.m`](RunAmplitudeVsData.m) | §4 — the slope comparison against all 6 recordings |
| [`RunR2CeilingSlope.m`](RunR2CeilingSlope.m) | §5 — the R2-ceiling falsification |
| [`RunLeverScan.m`](RunLeverScan.m) / [`RunLeverScan2.m`](RunLeverScan2.m) / [`AnalyzeLevers.m`](AnalyzeLevers.m) | ~100-evaluation one-at-a-time sweep + ranking |
| [`RunSlopeUnderLevers.m`](RunSlopeUnderLevers.m) | §6 — which levers move the slope vs the level, and the combination |
| [`RunComboCost.m`](RunComboCost.m) | §6 — full-battery bill and the compensation direction |
| [`RunOptimRsK.m`](RunOptimRsK.m) | the joint refit driver |
| [`TimeBattery.m`](TimeBattery.m) | where an evaluation spends its time (and what `MaxRunTime` needs to be) |

## Status of the refit

An overnight block-coordinate refit with `rsK` in the objective is the remaining
step; see [labdiary.md](labdiary.md) for its configuration and outcome.

The needed slowing is `d ln rsK` = ln(43.7/118.2) = **−0.995**. The one-at-a-time
scan (~100 evaluations, 2 factors × ~50 params) found **no lever that reaches it at
tolerable cost**. Two exceed it, both catastrophically — `ekSE`×0.6 gives −2.18 but
sends `ktr` to 329 s⁻¹ and the feature cost to +546; `kstiff2`×1.7 gives −0.48 but
puts `steady` at 119 kPa. The cheapest useful movers are:

| lever | `d ln rsK` | % of needed | ΔfeatCost | side effect |
|---|---|---|---|---|
| `kSE`×0.6 | −0.247 | 25 % | **+13.1** | `ktr` 52 → 41 |
| `k1`×1.7 | −0.215 | 22 % | +31.2 | — |
| `k_1`×0.6 | −0.148 | 15 % | +15.9 | — |
| `kah`×0.6 | −0.141 | 14 % | +22.1 | — |
| `k2`×0.6 | −0.630 | 63 % | +144.2 | `steady` 76 → 93, `ktr` 52 → 39 |

So a *combination* could plausibly reach the linear half if the optimiser can find
compensations that cancel the side effects — which is exactly what the refit tests,
and is not a foregone conclusion. What it cannot do is change the **sign** of the
amplitude slope; that is the residual the next mechanism has to carry.
