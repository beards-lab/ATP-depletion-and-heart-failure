# Lab diary — fitting the force redevelopment after re-stretch

Running log. Conclusions are distilled in [conclusions.md](conclusions.md).

---

## 2026-08-04 — session start

**Question (from FJ).** The worst part of the fit is the force redevelopment
*after the re-stretch is done*. Compared against `ktr` and against the post-slack
redevelopment, the model misbehaves. A **first-order fit** may be a better metric
than the LSQE currently in the objective; that should hint at what the model needs.
Base: `params/optfull7_opt_mov.m`.

### Step 0 — what was already known

- [`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) established
  that the **data has ONE recovery rate** (~41–47 s⁻¹) shared by the ktr protocol,
  the post-slack redevelopment and the post-restretch recovery, on all three
  protocol days. The model splits it into three.
- Its Part 2 fixed two of three defects (`UseMaxwellTensionOnly`, `tau_reg`).
  **Defect 3 — post-restretch ~2× too fast — was left open**, and Part 4 concluded
  no parameter fix was in sight: the SRX force gates make it *worse*, and the
  pools that carry the effect are the *state*, not parameters.
- That analysis measured the rate with a standalone script
  (`recoveryWindows.m`). **The rate was never a term in the fitting objective.**

### Step 1 — why the current objective cannot see the defect

The post-restretch window `[ramp end, next release]` is ~280 ms. The only
terms scoring it are `vall2_dy` (undershoot depth), `ovrsht_dy` (overshoot) and
`coolDownLS` (a least-squares error over the first 150 ms).

At the seed, `coolDownLS` contributes **0.32 of a 28.18 total feature cost — 1.1 %**.

The reason is structural, not a weight mistake: the rate information lives in the
first ~30 ms, while the remaining ~250 ms is plateau that the model already gets
right. A trace-wise least-squares average is therefore dominated by the part that
is correct. **FJ's instinct was right** — a rate metric exposes what an LSQE
averages away.

### Step 2 — the estimator

Added `fitRestretchRecovery` to [`Model/extractSlackAttributes.m`](../../Model/extractSlackAttributes.m),
emitting five new per-cycle features:

| feature | meaning |
|---|---|
| `rsK` | first-order rate (s⁻¹) — **the fit target** |
| `rsA` | recovery amplitude P−B (kPa) |
| `rsT0` | dead time (s) |
| `rsR2` | fit r² (guardrail) |
| `rsK63` | model-free crossing rate (diagnostic / cross-check) |

Conventions deliberately match `RestretchVsKtrRecovery/recoveryWindows.m` so the
numbers are comparable to that analysis' published tables:

- **baseline `B`** = the `vall2` minimum within 80 ms of the ramp end. The recovery
  starts there, *not* at the ramp end — the first ~7 ms are still falling, and a
  free-offset fit absorbs that descent and reports a spuriously slow rate (this is
  the error that produced the earlier, wrong "post-restretch is 1.18× slower" claim).
- **amplitude `A` = P−B held FIXED** (P = median of the last 15 % of the window).
  Only `k` and `t0` are free. Fixing `A` is what makes the estimate
  span-insensitive and keeps a force-scale error out of the rate — the amplitude
  is already scored separately by `vall2_dy`/`steady`.
- the anchor minimum is located on the lightly-smoothed trace (data path) so a
  single noise spike cannot set the baseline; the fit itself is on the raw trace.

**Span-insensitivity check** (refitting over progressively shorter spans):

| fit span | 30 ms | 50 ms | 80 ms | 120 ms | 200 ms | full |
|---|---|---|---|---|---|---|
| data, cycle 1 | 48.5 | 47.2 | 47.4 | 48.2 | 48.4 | 48.4 |
| model, cycle 1 | 113.0 | 121.8 | 122.1 | 122.1 | 122.1 | 122.1 |

Stable on both sides ⇒ the 2.7× gap is not an artefact of where the window ends.

### Step 3 — data features re-extracted, safely

`runSlackExperiment` reads `features_data` from the protocol `.mat` files, so the
new fields had to be written there.
[`AddRestretchFeatsToData.m`](AddRestretchFeatsToData.m) re-extracts, then
**compares every pre-existing field and refuses to overwrite any of them** — the
stored features are live fit targets and must not move silently.

⚠️ **Pre-existing inconsistency found and NOT changed.** Three of the seven files
no longer reproduce their own stored features under the current extractor:

| file | fields that differ (max rel. change) |
|---|---|
| `protocol_03_27_2026_8mM_slack` | `peak1_y` (26 %), `peak1_t` (92 %), `peak1_dSL` (**97 %**), `vall_y` (9 %), `vall_t` (96 %) |
| `protocol_03_27_2026_2mM_slack` | `peak1_y` (17 %), `peak1_t` (64 %), `peak1_dSL` (67 %) |
| `protocol_04_03_2026_*_slack` | `vall2_dy` (23 %), `vall2_t` (53 %), `ovrsht_dy` (**91 %**), `ovrsht_t` (17 %) |

The 04/10 files reproduce exactly. So the 03/27 targets — **the active fit set** —
were extracted by an *older* `extractSlackAttributes` than the one the model is
scored through. `peak1_dSL` and `peak1_y` carry weights 1 and 10 in the objective.
This is a real apples-to-oranges term in the current cost, independent of this
analysis. Left untouched here because correcting it moves the fit target and would
invalidate the current best parameters; **flagged for a separate decision.**

### Step 4 — the baseline gap

`params/optfull7_opt_mov.m`, protocol 03/27 8 mM:

| | c1 | c2 | c3 | c4 | c5 | mean |
|---|---|---|---|---|---|---|
| `rsK` data | 46.3 | 37.8 | 41.6 | 42.9 | 50.0 | **43.7** |
| `rsK` model | 118.4 | 127.7 | 134.0 | 101.9 | 106.9 | **117.8** |
| ratio | 2.56 | 3.38 | 3.22 | 2.38 | 2.14 | **2.70** |

For contrast, on the same run:

| observable | model / data |
|---|---|
| post-restretch rate `rsK` | **×2.70** |
| post-restretch amplitude `rsA` | ×1.11 |
| post-slack rate `ktr` | ×1.05 |

**The defect is clean and specific: a rate error in one window.** Amplitude,
plateau and the other two redevelopment rates are all essentially right. That is
a much sharper statement than "the restretch shape is off", and it is only
visible because the metric is now a rate.

Figure: [`results/baseline_windows.png`](results/baseline_windows.png) — the model
is done by ~30 ms where the data is still climbing at 80 ms.

### Step 5 — a strong, unlooked-for ATP result

Extracting `rsK` on every protocol file gave the 2 mM / 8 mM comparison for free:

| prep | 8 mM `rsK` | 2 mM `rsK` | ratio 2/8 |
|---|---|---|---|
| 03/27 | 43.7 | 15.9 | **0.36** |
| 04/03 | 51.5 | 20.4 | **0.40** |
| 04/10 | 65.1 | 18.0 | **0.28** |

Low ATP slows the post-restretch recovery by **~3×**, consistently across three
protocol days — far more than it slows `ktr` (×0.55, per
[`ATPEffectReconciliation`](../ATPEffectReconciliation/conclusions.md)). The
low-ATP fits are also *cleaner* first-order (`rsR2` 0.84–0.94 vs 0.51–0.80 at
8 mM). This is a sharp new constraint for the low-ATP phase — recorded here, not
pursued yet, per FJ's ordering.

### Step 6 — objective and runs

Added to the objective: `rsK|1` plus a guardrail `rsR2|0.4-1.0|0.5` (so the
optimiser cannot win on `rsK` by making the window non-exponential). At the seed:

```
rsK|1                8.4652   <- largest single term
FV_fnorm|FV_v|10     4.0182
A|50                 4.2300
steady|50            3.8132
...
TOTAL               36.8056   (was 28.18 without rsK)
```

Two jobs launched:
- `RunLeverScan.m` — one-at-a-time sweep over ~40 candidate levers, slack-only.
- `RunOptimRsK.m` (tag `rsk_w1`) — 14 h block-coordinate refit, full battery.

*(A machine restart killed the first launch of both; relaunched, no results lost.)*

**Two objective-plumbing facts found while wiring this up, both worth knowing:**

1. `evaluateBakersExp` calls `evalFeatureCost` **without** the `costExp` argument,
   so the optimiser minimises the **L2** (squared normalised) feature cost, while
   `RunBakersExp` reports **L1**. The two rank features very differently: a
   feature whose normalised error exceeds 1 is *amplified* by L2. `rsK`'s errors
   are 1.3–2.1, so `rsK|1` is 8.47 of 36.81 (23 %) in the reported L1 breakdown
   but **14.91 of 19.18 (78 %)** in what the optimiser actually minimises.
2. The first overnight launch silently optimised a **failed** evaluation: its
   baseline came back 1016441, the `1e6 + 1e3·unsimulated` ODE-timeout penalty,
   because `MaxRunTime = 60 s` left no margin once the lever scan was competing
   for cores (the slack init run alone is ~22 s — [`TimeBattery.m`](TimeBattery.m)).
   `RunOptimRsK.m` now refuses to start if the baseline exceeds 1e4.

### A dead end worth recording — don't linearise this model naively

My first instinct for "which relaxation modes does each window excite?" was to
build the Jacobian at the isometric steady state by finite differences and look
at the eigenvalue spectrum and modal force weights. **It does not work here**, for
three reasons that are properties of the model, not of the coding:

- most strain bins hold ~0 probability, so a relative finite-difference step is
  ~1e-12 there and the columns come back as exact zeros — the spectrum fills up
  with spurious zero eigenvalues and `inv(V)` is singular (RCOND ~1e-292);
- the serial element has τ = mu/kSE ≈ 3 µs, so the true spectrum spans ~9 orders
  of magnitude and is hopelessly stiff for a dense eigendecomposition;
- `f(PU_ss)` is not small (~1e15) when `out.params` is replayed with
  `Velocity = 0`, i.e. the "steady state" taken from a trajectory is not a fixed
  point of the RHS as re-assembled.

The empirical route below (shrink the perturbation, watch the rate) answers the
same question robustly using only forward integration, and is what all the
results here rest on. Script deleted rather than left broken.

### Step 7 — is the fast recovery nonlinear, or structural?

[`RunAmplitudeTest.m`](RunAmplitudeTest.m). Each window is a relaxation back to
the *same* isometric steady state from a different displaced state. For an
approximately linear relaxation the rate is a property of the system; only the
mode amplitudes depend on the starting point. The data behaves that way — ktr,
post-slack and post-restretch all give ~41–47 s⁻¹ at amplitudes of 1.00, 0.74
and 0.22 F_iso. The model does not. So: shrink the perturbation and watch.

Synthetic protocols, release and re-stretch by the same depth so the muscle
always returns to the same length:

| depth (ML) | 0.080 | 0.060 | 0.040 | 0.030 | 0.020 | 0.015 | 0.010 |
|---|---|---|---|---|---|---|---|
| model `rsK` | 137.2 | 109.2 | 102.7 | 92.6 | 93.2 | 87.4 | 79.8 |
| amplitude (kPa) | 12.9 | 15.2 | 16.0 | 14.4 | 13.2 | 12.7 | 6.8 |

(Below 0.010 ML the excursion falls under 3 kPa and the fitted "rate" is
estimator noise — those points are excluded, not read as convergence.)

**The gap decomposes into two halves of almost equal size:**

| | | |
|---|---|---|
| model small-signal limit | ~74–80 s⁻¹ vs data 44 | **×1.8 — a plain rate error** |
| large-stretch speed-up | 80 → 137 s⁻¹ | **×1.7 — a nonlinearity** |
| total at protocol depth | 137 vs 44 | ×3.1 |

In logs: ln 3.1 = 1.13 ≈ 0.59 + 0.54. So roughly half the defect is that the
model's *intrinsic* relaxation is too fast at its operating point — which
parameters could in principle fix — and half is a speed-up that only appears
when you stretch hard, which they cannot.

### Step 8 — the falsification: the data has no such amplitude dependence

[`RunAmplitudeVsData.m`](RunAmplitudeVsData.m). The protocol already sweeps
re-stretch amplitude — cycles 1–5 re-stretch by 0.080, 0.101, 0.121, 0.141 and
0.060 ML — so the answer is in the data for free, on every protocol day.

| series | n | slope d`rsK`/d(amplitude), s⁻¹ per ML | r |
|---|---|---|---|
| data 03/27 8 mM | 5 | **−94** | −0.64 |
| data 04/03 8 mM | 5 | **−95** | −0.72 |
| data 04/10 8 mM | 5 | **−188** | −0.45 |
| **data pooled 8 mM** | 15 | **−126** | −0.31 |
| model, synthetic sweep | 7 | **+714** | **+0.97** |

**The model gets faster the harder you stretch it (r = 0.97). The real muscle,
if anything, gets slower.** The sign is wrong, not just the magnitude. This is
the sharpest constraint the analysis produced, and it is a *shape* constraint —
it survives any per-prep force scaling, and it cannot be absorbed by re-weighting.

*Caveat, stated plainly:* in the synthetic sweep the release depth and the
re-stretch amplitude move together (by construction), whereas in the real
protocol they are decoupled (cycle 5 releases 0.161 ML but re-stretches only
0.060). So the two slopes are not measuring an identical quantity. The synthetic
sweep is the one that isolates a single knob, and the model's response to it is
near-deterministic; the data's slope over its own amplitude range is negative on
all three days independently. Fitting `rsK` per cycle scores this automatically.

### Step 9 — the leading candidate mechanism, tested and rejected

[`RunR2CeilingSlope.m`](RunR2CeilingSlope.m). The obvious suspect for an
amplitude-dependent speed-up is the unbounded strain acceleration of the
strong-bridge detachment `R2`, which reaches 2100–4250 s⁻¹ at re-stretch strains
(RestretchVsKtrRecovery Part 3): stretch harder → higher strain → faster
detachment → faster turnover. `UseR2Ceiling` (wired by
[`ApoPoolDetachment`](../ApoPoolDetachment/conclusions.md)) saturates it via
`1/R2_eff = 1/R2_strain + 1/R2max`. That analysis falsified the ceiling on the
ktr overstretch peak; it never tested it against the amplitude *slope*.

| `k2max` | slope | r | `rsK`@0.08 | `rsK`@0.01 | undershoot A@0.08 | F_iso |
|---|---|---|---|---|---|---|
| off | 714 | 0.97 | 137.2 | 79.8 | 12.95 | 76.3 |
| 1500 | 483 | 0.78 | 99.1 | 73.7 | 7.77 | 77.6 |
| 800 | 423 | 0.77 | 92.5 | 74.1 | 5.57 | 78.9 |
| 400 | 977 | 0.93 | 135.4 | 74.6 | 3.09 | 81.5 |
| 200 | 1607 | 0.89 | 209.1 | 79.7 | 0.93 | 85.3 |
| *data* | *−126* | *−0.31* | *~52* | | *~13* | *~80* |

**Rejected.** The slope improves marginally at `k2max` 800–1500 and then gets
*much worse*. The reason is visible in the amplitude column: the ceiling keeps
bridges attached through the stretch, so the undershoot itself collapses
(12.95 → 0.93 kPa). It removes the phenomenon rather than slowing it, and what
is left is a noise-dominated rate on a ~1 kPa excursion. Note also that
`rsK` at small amplitude stays 73–80 s⁻¹ at *every* setting — the ceiling does
not touch the ×1.8 linear half of the gap either.

This is now a **second, independent falsification** of R2-capping, on a
different observable from the one ApoPoolDetachment used.

### Step 10 — the lever scan

[`RunLeverScan.m`](RunLeverScan.m) + [`RunLeverScan2.m`](RunLeverScan2.m),
~50 parameters × {×0.6, ×1.7}, slack-only. Needed: `d ln rsK` = **−0.995**.
Ranked and plotted by [`AnalyzeLevers.m`](AnalyzeLevers.m).

The finding is a clean separation into three classes.

**(a) Big movers, all welded to force or `ktr`.** `k2`×0.6 is the largest
kinetic lever (−0.630) but drives `steady` 76 → 93 and `ktr` 52 → 39.
`kstiff2`×1.7 gives −0.479 with `steady` at 119. `ekSE`×0.6 gives −2.18 —
more than needed — with `ktr` at 329 and the feature cost at +546. These are
the "iron law" couplings already documented elsewhere in this repo.

**(b) SRX does nothing, again.** Every SRX rate and gate (`ksr0`, `kmsrd`,
`ksr2srd`, `ksrd2sr`, `sigma1`, `sigma2`, `sigma_srd1`, `sigma_srd2`) moves
`rsK` by |d ln| ≤ 0.09. This independently reproduces
[RestretchVsKtrRecovery](../RestretchVsKtrRecovery/conclusions.md) Part 4 at a
different parameter set: the post-restretch window never engages the trap
because the trap is gated on force and the window never loses force.

**(c) The one orthogonal lever: the series viscoelastic element.**

| lever | `d ln rsK` | ΔfeatCost | `ktr` | `steady` |
|---|---|---|---|---|
| `kSE_M`×0.6 | **−0.247** | **+6.55** | 52.0 *(unchanged)* | 76.2 *(unchanged)* |
| `kSE_M`×1.7 | +0.460 | +11.7 | 51.9 | 76.2 |
| `kSE`×0.6 | −0.247 | +13.1 | 40.9 | 74.8 |
| `mu_neg`×0.6 | −0.110 | +0.31 | 52.3 | 76.2 |
| `kstiff1_n`×1.7 | −0.103 | **+0.07** | 51.8 | 76.2 |

`kSE_M` — the Maxwell element's spring constant — is **monotone in `rsK` and
leaves `ktr` and isometric force completely untouched**. That is the first
lever found in this repo that moves the post-restretch rate without paying the
cross-bridge iron law, and it makes physical sense: it slows the *transmission*
of force rather than the cycling of bridges.

**Caveat that must be checked before believing it:** the scan is slack-only, so
it cannot see the FV or the passive PNB/Mava cost — and the passive experiment
is precisely where a series-viscoelastic change should send its bill.
[`RunTopLeversFull.m`](RunTopLeversFull.m) re-runs the top candidates on the
full battery to expose that.

### Step 11 — I predicted `eta_M` would shift the level, not the slope. Wrong.

I wrote in step 10 that the Maxwell class "cannot fix the amplitude-slope sign,
because a linear series filter slows all amplitudes by the same factor, so it
addresses the ×1.8 half only." [`RunSlopeUnderLevers.m`](RunSlopeUnderLevers.m)
tests that instead of assuming it, and it is **false**:

| variant | slope | `rsK`@0.080 | `rsK`@0.010 |
|---|---|---|---|
| baseline | +714 | 137.2 | 79.8 |
| `eta_M`×0.6 | +218 | 93.0 | 75.1 |
| **`eta_M`×0.3** | **+83** | 80.6 | 70.3 |
| `kSE_M`×0.6 | +300 | 102.2 | 78.4 |
| `k2`×0.6 | +243 | 65.4 | **44.7** |

`eta_M` is not acting as a simple output filter — it flattens the slope 8.6-fold
while barely moving the small-signal level. `k2` does the opposite: it drops the
level onto the data (44.7 vs 43.7) and leaves a slope of +243. **Two levers, two
halves.** So the obvious thing to try is both at once:

| | slope | `rsK`@0.080 | `rsK`@0.010 |
|---|---|---|---|
| **`k2`.6 + `eta_M`.3** | **+21** | **48.7** | **43.7** |
| `k2`.7 + `eta_M`.3 | **−32** | 59.7 | 62.0 |
| `k2`.6 + `eta_M`.3 + `kSE_M`.6 | +23 | 49.0 | 44.5 |
| *data* | *−126* | *~52* | — |

`rsK` ≈ 49 s⁻¹ at **every** amplitude. On the real protocol with the full battery
this is `rsK` = 46.8, **×1.07 of data**, down from ×2.70.

### Step 12 — the bill, and the right compensation

[`RunComboCost.m`](RunComboCost.m), real protocol, full battery:

| variant | `rsK` ratio | `ktr` ratio | steady | feature total |
|---|---|---|---|---|
| baseline | ×2.70 | ×1.05 | 76.2 | **28.2** |
| `eta_M`×0.3 alone | ×1.58 | ×1.05 | 76.2 | **38.0** |
| `k2`×0.6 alone | ×1.44 | ×0.79 | 92.5 | 178.4 |
| `k2`.6 + `etaM`.3 | **×1.07** | ×0.79 | 92.5 | 175.9 |
| + `ka`×1.4 | ×1.05 | ×0.78 | 105.4 | **367.6** |
| + **`kstiff2`×0.85** | ×1.09 | ×0.77 | 81.0 | **89.9** |
| `k2`.7 + `etaM`.3 | ×1.26 | ×0.86 | 87.4 | 124.7 |

Two things worth carrying forward:

1. **`eta_M` alone is the cheap win** — ×2.70 → ×1.58 for +9.8 total, localised in
   `ovrsht_dy` (+3.6), `vall2_dy` (+2.2), `vall_y` (+2.0) and `PS_restretchPeak`
   (+1.2). It *does* send part of its bill to the passive experiment, which is
   what a series-element change should do and is why the slack-only scan
   understated it.
2. **The compensation for `k2`↓ is `kstiff2`↓, not `ka`↑.** `ka`×1.4 doubles the
   damage (force 105 kPa); `kstiff2`×0.85 halves it. This is the opposite of the
   `A0`↓/`ka`↑ pairing that worked for the FV shoulder in
   [RestretchFeatureFit](../RestretchFeatureFit/conclusions.md) — there the goal
   was to add force without inflating peaks; here the goal is to *remove* force
   that `k2`↓ added.

`eta_M` was therefore added to the optimiser's compulsory set alongside
`k2`/`kstiff2`/`ka`.

### Step 13 — the overnight refit, as launched

Two arms, both seeded from `params/optfull7_opt_mov.m`, 18 h budget, 4 workers
each, `MaxRunTime = 600` (≈13× the healthy evaluation, so contention cannot
masquerade as a bad parameter point again):

| tag | `rsK` weight | share of objective | baseline cost | question it answers |
|---|---|---|---|---|
| `rsk_w025` | 0.25 | 48 % | **7.876** | can we get a *usable* snapshot — better `rsK` with the rest of the fit defended? |
| `rsk_w10` | 1.0 | 79 % | **19.178** | how far can `rsK` be pushed at all, and what is given up to do it? |

Both baselines match the values predicted from the saved feature vectors
(7.876 / 19.178) — confirming `rsK` is genuinely inside the objective and not
being scored as a missing feature.

Outputs: `params/rsk_w025_opt.m`, `params/rsk_w10_opt.m` (written only on
improvement), `params/rsk_*_state.mat`, and the per-round `optIter` gallery in
`params/rsk_*_iter/`. The gallery is the frontier: every round records its
converged params *and* the per-feature cost decomposition, so
`analyzeOptIterUniqueness('rsk_w025')` will show whether the near-best solutions
reconverge on the same parameters or are degenerate.

**What to look at first in the morning:** whether either arm got `rsK` below ~×1.5
*without* `steady` leaving 74–79 kPa or `ktr` leaving ×0.95–1.15 — that is the
combination the manual test could not achieve and the joint refit might.

---
