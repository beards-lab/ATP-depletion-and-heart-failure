# The ktr protocol against the slack redevelopment — the model splits them 1.5×, the data does not (2026-09-17)

8 mM ATP. Baseline `params/estiff125_opt.m` (the best 8 mM set when this started);
**outcome `params/fastgate_best.m`, which supersedes it** — see [Adopted](#adopted-paramsfastgate_bestm-is-the-new-8-mm-best).

Reproduce, in order:
`RunKtrVsSlack.m` → `StateDiag.m` → `RunSwaps.m` → `RunVariants.m` → `RunTauScan.m` →
`RunVerify.m` → `RunEstimatorTable.m` → `RunNoiseAndProbe.m` → `RunFeatureEstimatorBias.m` →
`MakeFigure.m` → `RunFVFineGrid.m` → `RunFVLowV.m` → `MakeFVFigure.m` → `TestKtrFeature.m` →
`BuildFastGateBest.m` → `PlotFeatureComparison.m` → `RunFinalCheck.m`.
Held-out check on any paramset: `validateKtrProtocol('<tag>')`.

Figures: `figures/ktr_vs_slack.png` (the defect and the fix), `figures/fv_finegrid.png` (§7),
`figures/features_all.png` + `figures/features_bars.png` (every feature, both variants).

---

## Answer in one line

The model's ktr-protocol redevelopment is **1.9× too fast** while its post-slack redevelopment at
the same length is nearly right. It is **not a missing mechanism and not a contradiction with the
slack features** — it is **one unidentified parameter**, `tau_reg`, which the objective cannot see
because the experiment that constrains it is switched off. Setting it right costs nothing: the
fitted objective goes **0.956 → 0.939**, the staircase is unchanged, and the held-out ktr trace
error drops **73 %** on all three protocol days.

---

## 1. The discrepancy

Both protocols end with force redeveloping isometrically at **SL 2.00 µm** from ~zero — the slack
protocol's cycle 2 releases 1.0999 → 1.0003 ML, and the ktr manoeuvre returns to 1.000 ML. Same
length, same starting force, same endpoint. In the data the two rates coincide:

| `k63` (s⁻¹) | 03/27 | 04/03 | 04/10 | mean |
|---|---|---|---|---|
| post-slack #2 | 50.0 | 44.4 | 39.2 | **44.6** |
| ktr protocol | 49.4 | 42.1 | 38.1 | **43.2** |
| **ratio** | 0.99 | 0.95 | 0.97 | **0.97** |

The model does not. And this is **not an estimator artefact** — all four estimators agree
(03/27; `kfree58` = free amplitude on a 58.5 ms span, the ktr-equivalent control):

| | k63 | kfixA | kfree | kfree58 |
|---|---|---|---|---|
| DATA post-slack #2 | 50.0 | 49.5 | 48.1 | 43.8 |
| DATA ktr protocol | 49.4 | 47.0 | 44.9 | 44.9 |
| **DATA ratio** | **0.99** | **0.95** | **0.93** | **1.02** |
| incumbent post-slack #2 | 55.2 | 56.6 | 59.2 | 69.1 |
| incumbent ktr protocol | 83.2 | 80.9 | 96.6 | 96.6 |
| **incumbent ratio** | **1.51** | **1.43** | **1.63** | **1.40** |

Model-free version, with **no rate estimator at all** — RMS of the normalised trace against the
data, next to the data's own RMS residual about its best single exponential (its noise floor):

| | post-slack #2 | ktr protocol |
|---|---|---|
| data noise floor | 0.0223 | 0.0222 |
| incumbent | 0.0247 = **1.11× noise** | 0.1195 = **5.38× noise** |

The slack window is fitted to within the measurement noise. The ktr window is **5.4× outside it**.
Panel A of the figure shows it: the incumbent's redevelopment is finished while the data is still
at 60 % of F_iso.

Note in passing that the incumbent is *also* fast on the slack — 55.2/56.6/59.2 against a data
50.0/49.5/48.1, i.e. 10–23 % depending on estimator. It is simply far less wrong there. §5 explains
why the objective did not notice.

---

## 2. Why: a 27 ms attachment gate that remembers how long the protocol shortened

`out` at the first sample of each redevelopment (`StateDiag.m`):

| | post-slack #2 | ktr | ratio |
|---|---|---|---|
| attached p1+p2 | 0.1688 | 0.0883 | 0.52 |
| PD (primed detached) | 0.5650 | 0.5961 | 1.06 |
| SR + SRD | 0.110 | 0.248 | 2.25 |
| **A_reg** | **0.2745** | **0.3434** | **1.25** |
| **ka_eff** = ka·A_reg·OV·f_lat·f_sat | **60.8** | **71.1** | **1.17** |

`A_reg` is the registration-availability gate, a first-order state
([dPUdT_CombinedTransitions.m:887](../../Model/dPUdT_CombinedTransitions.m)):

```
A_inf  = A0 + (1-A0)·v/(v + v_ref_reg)          A0 = 0.211,  v_ref_reg = 0.630 µm/s
dA_reg = (A_inf - A_reg)/tau_reg                tau_reg = 26.65 ms
ka_eff = ka · A_reg                             (RegAvailShorteningOnly: v = max(0,-vel))
```

`v_ref_reg` is 0.63 µm/s and every ramp here runs at 60–215 µm/s, so `A_inf ≈ 0.994` during *any*
of them. The gate is therefore a pure **stopwatch on shortening duration**, and the two protocols
shorten for different lengths of time:

| | shortening | predicted A_reg | measured |
|---|---|---|---|
| slack release | 2.2 ms | 0.211 + 0.783(1−e^(−2.2/26.6)) = 0.273 | 0.2745 |
| ktr, first release | 4.8 ms | 0.211 + 0.783(1−e^(−4.8/26.6)) = 0.341 | 0.3401 |
| ktr, + final release | +1.7 ms | — | **0.3434** |

The ktr manoeuvre shortens for **6.5 ms** (4.8 ms down to 0.80 ML plus the 1.7 ms final release);
the slack for **2.2 ms**. With a 27 ms memory the gate is still 25 % higher when redevelopment
starts, and it decays over exactly the ~22 ms the redevelopment takes. Panel C.

### State transplants confirm it carries the effect

Both windows are pure holds, so each is fully determined by its state vector. `RunSwaps.m` replays
each from its captured `PU0` (`seedWindow.m`) and transplants one pool at a time from the post-slack
donor into the ktr recipient. `PT` is the slack variable
([dPUdT_CombinedTransitions.m:222](../../Model/dPUdT_CombinedTransitions.m)) so every single-pool
swap conserves probability automatically.

Replay control: post-slack 57.0, ktr 83.2 — gap **26.2 s⁻¹**.

| transplant into the ktr window | k63 | vs control | share of the gap |
|---|---|---|---|
| **A_reg ← slack** | 66.5 | **−16.7** | **64 %** |
| PD ← slack | 76.9 | −6.3 | 24 % |
| SR, SRD ← slack | 86.3 | +3.1 | — (wrong sign) |
| A_def ← slack | 84.6 | +1.4 | — (wrong sign) |
| A_reg + A_def | 70.3 | −12.9 | |
| A_reg + A_def + SRX | 69.4 | −13.8 | |

`A_reg` alone carries **64 %**; `PD` — itself a consequence of the longer manoeuvre — another 24 %.
**The SRX pool contributes nothing**, reproducing the negative result of
[`../RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) Part 4 (H3) on a different
paramset. `A_def` (compliant realignment, `UseLoadRealign`, τ 21.3 ms — the *other* slow attachment
gate, added since that analysis) also contributes nothing: it is charged by bound mass stripped
**under load**, and both of these windows strip it at F ≈ 0.

> **Invalid swaps, recorded so nobody retries them.** Copying `p1,p2` (or the whole vector) blows up
> to 5×10⁵ s⁻¹. The donor's strain grid is registered to `LXBpivot = 2.20 µm` and the recipient's to
> `2.00 µm`, so a raw element-wise copy is a different strain distribution, not the same one.
> Same caveat as `RestretchVsKtrRecovery`'s H1raw.

---

## 3. Why nobody noticed: the ktr protocol is not in the objective

`params/estiff125_opt.m` has **`params0.RunKtr = 0`** — the protocol is not even simulated. And it
would not have mattered if it were: `runKtrExperiment` returned only a **trace MSE**, and
`params0.OptimizeOn = 'Feats'` reads **only** features. Every `ktr` entry in `params0.fn` comes from
`extractSlackAttributes`, i.e. from the slack protocol.

So `tau_reg` is unidentified. `RunTauScan.m`, three decades:

| tau_reg (ms) | 26.65 | 10 | 5 | 2 | 1 | 0.5 | 0.2 | 0.1 |
|---|---|---|---|---|---|---|---|---|
| **E_t (the fitted objective)** | 0.956 | 1.066 | 1.005 | 0.976 | 1.002 | 0.946 | 0.939 | **0.936** |
| ktr protocol trace MSE | 0.175 | 0.178 | 0.131 | 0.065 | 0.044 | 0.037 | 0.037 | **0.037** |
| k_ktr / k_postSlack2 | 1.51 | 1.50 | 1.36 | 1.20 | 1.15 | 1.07 | **1.04** | 1.05 |

**E_t moves by 0.13 across the whole scan and is minimised at the fast end. The held-out ktr error
moves 5.1× and is monotone.** Panel D. The incumbent's 26.65 ms is not where the fit put `tau_reg`;
it is where the search left it, because nothing was pulling.

---

## 4. The fix, and what it costs

**`tau_reg` ≲ 0.5 ms** — make the availability gate quasi-instantaneous, so it is a function of the
*present* velocity and nothing is remembered once the isometric hold starts. Redevelopment then
becomes protocol-independent by construction: both windows start at `A_reg = A0`.

This is also the more defensible reading of what `A_reg` is. Site availability is a geometric
property of the instantaneous relative sliding; a 27 ms memory would mean the filament remembers
that it was sliding 27 ms ago, and there is no structure that stores that.

Priced against everything (`RunVerify.m`), `tau_reg` 26.65 → 0.2 ms:

| | incumbent | fast gate |
|---|---|---|
| **E_t (fitted objective)** | 0.9559 | **0.9390** |
| ktr k63, 03/27 / 04/03 / 04/10 (data 49.4 / 42.1 / 38.1) | 83.2 / 84.3 / 81.5 | **48.6 / 49.4 / 49.0** |
| ktr trace MSE, three days | 0.175 / 0.201 / 0.243 | **0.037 / 0.043 / 0.053** (−73 / −79 / −78 %) |
| stairs E (held out, `RunStairs = 0`) | 0.0229 | 0.0229 (unchanged) |
| trace RMSE vs noise, ktr | 5.38× | **1.42×** |
| trace RMSE vs noise, post-slack | 1.11× | 1.31× |

Inside the objective it is a **swap, not a loss**. `PlotFeatureComparison.m` →
`figures/features_all.png` (every scored feature, data vs both variants) and
`figures/features_bars.png` (paired per-feature cost + what the fix trades):

| feature | incumbent | fast | Δ |
|---|---|---|---|
| FV_normpowerAvg | 0.2656 | 0.1468 | **−0.119** |
| ktr (slack) | 0.0069 | 0.0906 | **+0.084** |
| rsK | 0.0232 | 0.0499 | +0.027 |
| ktr_rmse | 0.0172 | 0.0002 | −0.017 |
| everything else (20 features) | | | < 0.005 each |
| **total** | **0.9559** | **0.9390** | **−0.017** |

The post-slack window degrades from 1.11× to 1.31× the data's own noise — a change *inside* the
measurement noise — and the ktr window goes from 5.4× outside it to 1.4×.

---

## 5. Is the residual a real trade against the slack? Mostly it is the estimator

The one term that clearly worsens is the objective's slack `ktr` feature: 50.7 → 42.7 against a
target 48.9. But `recoveryWindows` on the *same model trace* reads 47.9. Both comparisons are
like-for-like, and they disagree by 12 % — on the data the same two estimators agree to 2 %.

`RunFeatureEstimatorBias.m` isolates it. `extractSlackAttributes` drops every sample below
`FORCE_THRESHOLD_KPA = 10` (≈13 % of the plateau) and then fits with `F0` pinned at 0 and `t0` free:

| post-slack #2 trace | fit with thr = 0 | fit with thr = 10 kPa | change |
|---|---|---|---|
| DATA | 48.1 | 48.9 | **+1.6 %** |
| incumbent | 57.0 | 50.7 | **−11.1 %** |
| fast gate | 47.0 | 42.7 | **−9.2 %** |

**The cut is harmless on the data and worth −9 to −11 % on the model.** The model's rise has a fast
early foot that the data does not have; removing the first 13 % of the amplitude while holding
`F0 = 0` makes the remainder fit a slower exponential. The bias is **one-sided — model traces only —
so the feature systematically rewards models that are genuinely too fast**, which is exactly the
incumbent's situation (slack `ktr` feature 50.7 vs a target 48.9, while its true rate is 55–59).

Correcting for it, the fast gate's post-slack rate is 47.0 (threshold-free) / 46.6 (k63) /
46.2 (kfixA) against a data 48.1 / 50.0 / 49.5 — **within 4–7 %**, not 13 %.

A real residual of a few percent remains. It is **not** reachable by hand:

- `ka` alone: `E_t` 0.939 → 2.46 at ×1.10, 12.3 at ×1.20. The breakage is `coolDownLS`
  (0.099 → 1.12) and `FV_normpowerAvg` (+0.24) — `ka` is load-bearing for the ramp features.
  And the slack `ktr` feature *falls* (42.7 → 42.3 → 41.4) rather than rising.
- `(ka, kd)` grid: `kd` ×1.3 does lift the rate (slack 42.7 → 45.1, ktr 48.6 → 49.9) but at
  `E_t` 2.85.

**The important result from that grid is the ratio: 1.00–1.06 at every one of the nine points.**
Once the gate is fast, protocol-independence is *structural* — it does not depend on the tuning.
Restoring the absolute rate is an ordinary refit, not a search for a mechanism.

---

## 6. Wiring change: the objective can now see this protocol

Small, additive, verified inert.

| file | change |
|---|---|
| `Model/extractKtrAttributes.m` | **new.** `ktrP_k63` (crossing rate against the known plateau P = F_iso) and `ktrP_F20`. The crossing rule is byte-identical to `recoveryWindows.m` — verified to 0.000 in `TestKtrFeature.m`. |
| `Model/experiments/runKtrExperiment.m` | returns `features_model`, `features_data`. Gated on the protocol actually **being** a ktr manoeuvre (last ramp a release into a ≥30 ms hold) so the FV2 files that share this runner emit nothing. |
| `Model/RunBakersExp.m` | merges them into the feature structs. |
| `Analyses/RestretchVsKtrRecovery/seedWindow.m` | **bug fix.** Its `PU0` length check predated `UseLoadRealign`, so it rejected every paramset with that flag on — which the current best set has. |
| `Auxiliary/plotMultipleFeatures.m` | **bug fix.** It mistook `evalFeatureCost`'s later `w:<varfield>` and `fl:<x>` tokens for an X-axis field name, so `FV_normpowerAvg|w:FV_normpowerVar|…` and any floored feature rendered as an empty "(no data)" tile. Both are scored features, so two of the objective's terms — including its single largest — were silently unplottable. |

Why `k63` and not a free-amplitude fit: the window is 58.5 ms and the data has only reached
0.93 F_iso by its end. Free `A` on a span that short returns ~86 % of the true rate
([`../RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md), "Estimator caveat").
The muscle is back at the pre-release length, so the plateau is *known*.

Regression checks (`TestKtrFeature.m`, `RunFinalCheck.m`), all pass:

- `estiff125_opt` as shipped: `E_t` = 0.9559, identical to before the change.
- `RunKtr = 1` with `fn` unchanged: 0.9559. Features present, unscored, free.
- FV2 file through the runner with `RunKtr = 1`: runs clean, emits nothing.
- `ktrP_k63` scored but unobtainable: +81, the `MISSING_FEATURE_COST` — loud, not silent.

### The ktr protocol is deliberately NOT scored

`ktrP_k63` exists so the protocol can be **measured**, not so it can be fitted. The adopted
paramset keeps `RunKtr = 0` and no `ktrP_*` in `fn`, which leaves the ktr protocol a genuine
out-of-sample test of every future refit. Scoring it would buy a number and spend the only
independent check this model has on a redevelopment rate.

For the record, had it been scored at the weight the repo's own rule gives (`w = λ/s²`, λ anchored
on the slack `ktr` weight; `s` = 0.132, the between-day spread of 49.4/42.1/38.1; floor = one
preparation-width, `fl:0.0175`), the incumbent would pay **0.367** and the fast gate **0.000** — at
0.367 it would have been the largest single term in the objective, larger than F-V's 0.266. That
is a statement about how badly the incumbent misses, not an argument for scoring it.

`validateKtrProtocol('<tag>')` runs the held-out check on all three days for any paramset.

---

## 7. Does the fix break force-velocity? No — but the fine grid finds a different defect

`RunFVFineGrid.m` / `RunFVLowV.m` / `MakeFVFigure.m` → `figures/fv_finegrid.png`.
Velocity grid 0.25 ML/s up to 4, then 1 ML/s to the end of the protocol (6).
**The scored FV uses only five velocities: `FV_velocities = [0 -0.5 -1 -2 -4]`.**

**The fast gate does not break FV.** `F_iso` is 65.46 kPa for both, to the digit.

| v (ML/s) | 0.5 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| incumbent F/F₀ | 0.918 | 0.627 | 0.297 | 0.169 | 0.111 | 0.078 | 0.055 |
| **fast gate F/F₀** | 0.918 | 0.630 | 0.302 | 0.167 | 0.106 | 0.072 | 0.050 |
| data Baker2022 | 0.919 | 0.664 | 0.316 | 0.203 | 0.111 | 0.058 | 0.039 |
| data Isovel2021 | 0.905 | 0.634 | 0.280 | 0.187 | 0.096 | 0.052 | 0.028 |

The two gates differ by ≤ 0.006 in F/F₀ anywhere, always inside the spread between the two
reference datasets. The fast gate is marginally *higher* below 2.5 ML/s and marginally *lower*
above — the incumbent's `A_reg` cannot fully charge during the shortest ramps (at 6 ML/s the ramp
is 16.7 ms against a 26.6 ms time constant), so the slow gate slightly under-attaches at low
velocity and over-attaches at high. Both effects are sub-1 %.

| | incumbent | fast gate | data |
|---|---|---|---|
| peak normalised power | 0.647 @ 1.25 ML/s | **0.654** @ 1.25 ML/s | 0.664 @ 1.00 / 0.634 @ 1.00 |
| `E_fv` (5 scored points) | 0.01895 | **0.01714** | — |
| `FV_normpowerAvg` cost | 0.2656 | **0.1468** | — |

(The data cannot locate the power peak better than ±0.5 ML/s — its grid is 1 ML/s spaced.)

### The defect the 5-point grid hides: the FV curve is non-monotonic

Between the two lowest scored velocities, 0 and 0.5 ML/s, **force rises 21 % above isometric**:

| v (ML/s) | 0 | 0.01 | 0.04 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|
| F/F₀ | 1.000 | 1.029 | 1.156 | **1.211** | 1.162 | 1.079 | 1.035 | 0.918 |
| `A_inf` | 0.211 | 0.223 | 0.300 | 0.401 | 0.517 | 0.596 | 0.652 | 0.695 |

It peaks at **+21.1 % at 0.10 ML/s** and crosses back through F₀ at ~0.40 ML/s, landing exactly on
the data at 0.5. The v→0 limit converges correctly (1.029 at 0.005 ML/s), so this is a genuine,
smooth, well-resolved feature of the model and not a numerical artefact.

**Cause: `A_reg`'s velocity dependence, not its time constant.** The two gates are byte-identical
here — at these velocities the ramps last 0.8–20 s, so `A_reg = A_inf` either way. The `A0 = 1`
control (gate velocity-blind) removes the hump completely: 1.000 → 0.995 → 0.959 → 0.884, monotone.

The reason is that **`v_ref_reg` = 0.63 µm/s is 4 % of Vmax** (~15 µm/s). The gate is therefore not a
graded velocity dependence but very nearly a *step*: isometric availability 0.211, any shortening at
all ≈ full. At 0.10 ML/s availability has already doubled while sliding has barely begun to strip
the bound population, so attachment wins and force goes **up**. A monotonically falling
force-velocity relation is about the most robust fact in muscle mechanics, and shortening produces
force *depression*, not enhancement — so the sign here is wrong, and it is wrong for both gates.

This is the price the model pays for the FV shoulder, and it has never been visible because the
scored grid steps straight from 0 to 0.5 ML/s. It is **independent of this analysis's fix** —
adopt the fast gate on the ktr evidence, and treat the hump as its own item.

---

## Adopted: `params/fastgate_best.m` is the new 8 mM best

`BuildFastGateBest.m` → `params/fastgate_best.m` = `params/estiff125_opt.m` with **one field
changed**, verified by field-by-field diff: `tau_reg` 26.6487 ms → 0.2 ms. The objective is
byte-identical (same 24 `fn` entries, `RunKtr = 0`), so the comparison below is like-for-like.

| | estiff125_opt | **fastgate_best** |
|---|---|---|
| **E_t** (24 scored features) | 0.9559 | **0.9390** (−1.8 %) |
| staircase E (held out, `RunStairs = 0`) | 0.0229 | 0.0229 |
| **ktr protocol k63**, 03/27 / 04/03 / 04/10 (held out; data 49.4 / 42.1 / 38.1) | 83.2 / 84.3 / 81.5 | **48.6 / 49.4 / 49.0** |
| ktr trace error ÷ the data's noise floor | 5.38× / 5.74× / 6.14× | **1.42× / 1.59× / 1.97×** |
| k_ktr / k_postSlack2 (data 0.99) | 1.50 | **1.05** |

Regression: `RunFinalCheck.m` asserts `RunKtr = 0`, no `ktrP_*` in `fn`, `fn` identical to the
incumbent's, and `tau_reg` the only differing field.

## What to do next

1. **Refit from `fastgate_best`** to recover the few percent on the absolute rate. The `(ka, kd)`
   probe says the direction exists and that the ktr/slack ratio holds across it; the naive levers
   are expensive, so this needs the optimizer with compensations, not a hand tune. Re-run
   `validateKtrProtocol` on whatever comes out — that is the whole point of leaving it unscored.
3. **Audit `extractSlackAttributes`' `ktr`.** The 10 kPa threshold costs the feature 9–11 % on model
   traces and 1.6 % on data traces. Either drop the threshold (`F0 = 0` is already anchored, so the
   low-force samples are informative, not harmful) or apply the same cut to both — as it stands the
   feature has a one-sided bias that rewards too-fast models. Same class of defect as
   [`../RestretchMetricAudit`](../RestretchMetricAudit/conclusions.md).
4. **Re-run the 2 mM set** (`ATP-wt-lowatp` worktree) once the gate is fast. The ATP effect on `ktr`
   is measured from the *ktr protocol*
   ([`../ATPEffectReconciliation`](../ATPEffectReconciliation/conclusions.md): ×0.72, CV 4 %), which
   is the protocol this analysis shows the model gets wrong at 8 mM — so the low-ATP `ktr`
   conclusions rest on a window the model was not fitted to.
5. **Not needed:** a new mechanism. The SRX pool, `A_def`, and the bound-pool composition were all
   tested and none of them differentiates these two windows.
6. **Separate item — put a velocity between 0 and 0.5 ML/s in the FV grid.** §7: the model's FV
   curve rises **21 % above isometric at 0.10 ML/s**, unconstrained because `FV_velocities` steps
   0 → 0.5. `v_ref_reg` = 0.63 µm/s is 4 % of Vmax, which makes `A_reg` a step rather than a graded
   velocity dependence. Adding one scored point at ~0.25 ML/s would make the hump pay; whether it
   can be removed while keeping the shoulder is the real question, and it needs a `v_ref_reg` scan.
   There is no such measurement in the current data, so this is also an acquisition request.

---

## Relation to the earlier analysis

[`../RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) (2026-07, on
`ThursdayNightFever.m`) found the same `tau_reg` defect and recommended exactly this fix — its
"Defect 2", with a predicted landing zone of 1–5 ms. Its Defect 1 (`UseMaxwellTensionOnly`) was
adopted; **Defect 2 was not**, and `tau_reg` is back at 26.6 ms two paramsets later. This analysis
re-establishes it on the current best set and differs from it in three ways:

- it prices the fix against the **current feature objective** rather than the old trace cost, and
  finds **no trade-off at all** (0.956 → 0.939) where the old analysis predicted one
  (`E_slack` 173 → 180). The counterweight is F-V, which improves by more than the slack loses;
- it shows the landing zone is **below** 1 ms, not 1–5 ms — the ramps are 1.7–4.8 ms long, so
  `tau_reg` has to be well under that for the gate to return to `A0` before the hold begins;
- it removes the reason the fix keeps getting lost, by putting the protocol in the objective.
