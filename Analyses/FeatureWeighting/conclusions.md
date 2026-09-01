# FeatureWeighting — the objective's weights, taken from the data instead of by hand

**Headline.** Every data-fit weight in `params0.fn` is now `w = λ/s²`, where `s` is the
**measured between-preparation spread** of that feature at 8 mM. One global `λ` makes the
change **cost-neutral at `params/realign4_opt.m` (3.8513)**, so a refit starts from the
same number. Two shares are pinned by hand and flagged as such: **FV at 70 %** (§4) and a
30 % cap on any other term. ATP relevance is measured and *drawn* but does **not** enter
the weights (§4).

> **The diagnostic that fell out of it: in units of its own reproducibility, the model
> already fits everything except the re-stretch recovery.** Sixteen data terms, and
> fifteen sit at χ² = 0.3–3.7. `rsK` sits at **33.4** and the unscored `rsT0` at **36**.
> There is one defect left, and it is the post-re-stretch recovery — the same one
> [RestretchRecoveryFit](../RestretchRecoveryFit/conclusions.md) named structurally.

`RunFeatureSpread.m` → `results/feature_spread.mat`
`BuildWeightedObjective.m` → `results/weight_table.txt`, `results/feature_weighting.png`,
`params/spreadw_seed.m`
`RunOptimSpreadW.m` → the large-scale refit.

---

## 1. Why the old weights had to go

The weights in `realign4_opt`'s `fn` span 0.05 to 62.3 and were each set for a local
reason at the time — a protection weight here, a cost-neutrality constant there. Nothing
made them comparable, so the objective's *composition* was an accident:

| term | share of 3.8513 |
|---|---|
| `FV_normpowerAvg` | **34.8 %** |
| `rsK` | 9.9 % |
| `peak1_dSL` | **8.5 %** |
| `t0_crossing` | 6.8 % |

`peak1_dSL` at 8.5 % is the clearest symptom: it is the *least reproducible* feature in
the whole set (**76 %** between-prep spread) and it was being fitted harder than `ktr`.

## 2. The measurement: how well is each feature determined?

`RunFeatureSpread.m`, three 2026 protocol days, features re-extracted fresh through
`extractSlackAttributes` so every number comes from one code path.

**The Baker slack recordings are excluded on purpose** — different protocol, recovery
windows 2–6× too short, amplitudes ATP-dependently truncated
([SlackDataAnalysis](../SlackDataAnalysis/conclusions.md) §1–2). *(The FV curves are
still Baker-sourced; there is no 2026 replacement for a full F–V curve, and FV's spread
is taken between its own two source datasets instead.)*

Force-dimension features are divided by their own run's mean `steady` first — absolute
force is not transferable between preps, shape is (§4 of SlackDataAnalysis). `steady`
itself keeps its raw CV and is the single absolute anchor.

| | features | `s` (between-prep relative SD) |
|---|---|---|
| **tight** | `peak2` 6.4 %, `Am` 7.9 %, `A` 8.4 %, `vall_y` 8.4 %, `peak1_y` 9.7 % | the trace *shape* is highly reproducible |
| **usable** | `ktr` 11.4 %, `rsR2` 12 %, `rsK` 21.4 %, `ktr_rmse` 23 %, `steady` (absolute) 23.6 %, `rsT0` 23.9 %, `vall2_dy` 24.3 %, `rsA` 24.6 % | |
| **loose** | `restretchSlopeStart` 30 %, `vall2_t` 32.7 %, `t0` 37.7 %, passive `PS_*` 46–48 % | de-weight |
| **barely measured** | `peak1_t` 63.6 %, `peak1_dSL` 76 %, `vall_t` 83.5 % | do not chase |

## 3. The measurement: which features carry the ATP effect?

Ratio 2 mM / 8 mM per prep — rundown-corrected for force-dimension features, raw for
rates ([ATPEffectReconciliation](../ATPEffectReconciliation/conclusions.md)) — and its
discriminability `z = |r̄−1| / √(s² + s_r²)`: how many noise widths ATP moves the feature.

| feature | ratio | z |
|---|---|---|
| `ktr` | **0.540** | **3.69** |
| `rsK` | **0.345** | **2.95** |
| `rsT0` | **0.428** | **2.19** |
| `vall2_t` | 1.924 | 1.66 |
| `steady` (absolute) | **1.346** | 1.34 |
| `rsR2` | 1.187 | 1.20 |
| `t0` | 1.439 | 1.04 |
| `restretchSlopeStart` | 1.265 | 0.61 |
| every force *shape* feature — `A`, `Am`, `peak1_y`, `vall_y`, `peak2`, `vall2_dy` | 0.99–1.04 | **0.07–0.21** |

> **Low ATP scales the trace and slows it. It does not reshape it.**
> The whole ATP signal lives in one scale (`steady` ×1.35) and in the rates
> (`ktr` ×0.54, `rsK` ×0.35, `rsT0` ×0.43). Every force-shape feature is flat to within
> a fifth of a noise width. This is a *falsifiable constraint on the mechanism*: whatever
> ATP gate is added must be nearly shape-neutral at 8→2 mM.

## 4. The rule

```
w_f = λ / s_f²                                    (κ = 0, so R ≡ 1)
```

`evalFeatureCost` forms `Σᵢ |fdᵢ − fsᵢ|² / mean(fd)²` — a **sum of squared relative
residuals** — so dividing by `s_f²` turns every term into a **χ²**: *1.0 means the model
is one preparation-width from the data.* Terms become directly comparable, which is
exactly what the hand-set weights were not.

`λ` is one global constant fixing the reweighted data terms to the subtotal they carry
today (3.4303). Guardrail/boundary terms (`XTOR`, `k2`, `doublePeak`, `coolDownLS`,
`ktr_rmse`, `rsR2`, the passive RMSEs) are **not data** and are untouched, so the grand
total is preserved to 1e-6.

### Why ATP relevance is *not* in the weights

An earlier draft carried `R = 1 + κ·min(z,4)`, κ = 0.5 — a bonus for ATP-responsive
features. **`κ` is now 0, deliberately.** Three reasons, in order of weight:

1. **It breaks the χ².** A weight that also encodes "relevance" is no longer a variance,
   so `cost/s²` stops being a goodness-of-fit and the table stops being readable as one.
   That diagnostic (§6) is the most useful thing this analysis produced; it is not worth
   trading for a tilt.
2. **It is a proxy for a missing term.** The honest way to make ATP matter is a term that
   *scores* the 2 mM condition (§8). A proxy living inside the weights is hard to
   distinguish later from a genuine ATP mechanism — which is exactly the confound the
   tilt was meant to avoid.
3. **It bought little.** `R` spans 1.0–3.0; with FV pinned at 70 % the whole non-FV data
   budget is ~7 %, so the tilt was second-order anyway. And the data already points the
   optimiser at the ATP features on its own: `rsK` is the largest χ² in the set.

ATP relevance is instead **measured, reported and drawn** — red bars in every panel of
`results/feature_weighting.png`, and the `z` column of `results/weight_table.txt`. Set
`KAPPA = 0.5` in `BuildWeightedObjective.m` to put it back into the weights.

### Two conditioning guards — engineering, and labelled as such

* **`FV_SHARE = 0.70` — *this number is arbitrary*.** Force–velocity is a separate,
  independently acquired experiment and the whole muscle work profile is in it, so it is
  held at a fixed share by fiat. Its *measured* spread (11.3 %, and that between two
  **labs** — isometric 56.4 vs 67.4 kPa) would put it at **1.4 %**: the model already
  fits FV to well inside the disagreement between its own two source curves, so `1/s²`
  sees nothing left to buy. That is a statement about the FV data's precision, not about
  FV's importance, and this pin is where the difference is asserted.
* **`SHARE_CAP = 0.30`** on every *other* term. Pure `1/s²` hands `rsK` 70 %, because it
  is 5.8 preparation-widths off — a model *rejection*, not a fitting signal, and an
  objective one feature can buy is how realign2 spent FV monotonically. Applied **after**
  the FV pin, so at `FV_SHARE ≥ ~0.55` it is inert.

### What FV_SHARE = 0.70 costs, in numbers

| | |
|---|---|
| FV | **70.0 %** |
| `rsK` | 12.5 % |
| guardrails (unchanged) | 10.9 % |
| **all fourteen other data terms, combined** | **6.6 %** |

The slack shape is largely along for the ride. **`FV_SHARE = 0.40–0.50` keeps FV dominant
while the slack still has a vote** — one number in `BuildWeightedObjective.m`.

> **⚠ The pin has a side effect worth knowing about.**
> `extractForceVelocityAttributes` normalises the model's power by the **model's own**
> isometric force (`FV_fnorm = FV_f/FV_f(1)`), so **FV carries no information about how
> strong the fibre is**. At 70 %, most of the objective is scale-blind by construction,
> and the weight backing the absolute level (`A`+`peak2`+`peak1_y`+`vall_y`+`steady`)
> falls **125 → 9.4, i.e. 13×**. Not a free ride — a uniform +20 % force drift still
> costs ≈1.9, about half the objective — but 13× cheaper than it was. `ANCHOR_GUARD =
> true` appends one boundary entry `steady[1]|50.6-88.7|5` (the range the three preps
> actually span at 8 mM). It costs exactly 0 inside the range, so cost-neutrality is
> untouched. Default **off**.

## 5. What moved

| term | share before | share after | why |
|---|---|---|---|
| `FV_normpowerAvg` | 34.8 % | **70.0 %** ¹ | pinned. `1/s²` alone would give 1.4 % (χ² 1.67 over 7 points) |
| `rsK` | 9.9 % | **12.5 %** | χ² 33.4. The one broken feature |
| `peak2` | 1.5 % | 1.0 % | the single most reproducible feature (6.4 %) — least de-rated of the slack terms |
| `ktr` | 3.7 % | 0.4 % | χ² 1.11 — already fitted to within a preparation-width |
| `peak1_dSL` | **8.5 %** | **0.2 %** | 76 % spread — was being fitted 5× harder than the data can justify |
| `t0_crossing` | 6.8 % | 0.3 % | 38 % spread |
| `A` | 2.7 % | 0.1 % | χ² 0.29 — fitted 3× better than reproducible |
| `steady` | 0.3 % | 0.003 % | fitted essentially exactly; see the anchor warning above |

¹ pinned by fiat, not by `1/s²`.

![Spread, ATP relevance, χ² and the reweighting](results/feature_weighting.png)

*Red = carries the ATP effect (z ≥ 1). Top left: reproducibility. Top right: where the ATP
signal is — every force-shape feature is flat. Bottom left: the model against the noise
floor. Bottom right: composition before/after, log scale (FV's pin makes a linear axis
useless).*

## 6. The χ² diagnostic — where the model actually stands

Residual in units of the measured between-preparation spread. **χ² ≈ 1 per element means
the model is as close as two preparations are to each other**; asking for better is
asking it to fit one preparation's idiosyncrasy.

| feature | χ² | per element |
|---|---|---|
| **`rsK`** | **33.40** | **6.68** |
| `PS_steady20` | 3.67 | 3.67 |
| `PS_steady22` | 2.84 | 2.84 |
| `peak2` | 2.76 | 0.55 |
| `FV_normpowerAvg` | 1.67 | 0.24 |
| `vall2_dy` | 1.43 | 0.29 |
| `restretchSlopeStart` | 1.28 | 0.26 |
| `peak1_y` | 1.22 | 0.24 |
| `ktr` | 1.11 | 0.22 |
| `t0_crossing` | 0.92 | 0.18 |
| `vall_t` / `peak1_dSL` / `vall_y` | 0.54–0.69 | 0.11–0.14 |
| `PS_restretchPeak` / `A` | 0.29–0.31 | 0.06 |
| `steady` | 0.00 | 0.00 |

**Read this before spending more optimiser time.** Fourteen of sixteen terms are at or
below the noise floor. The gradient that is left is `rsK`, the two passive steady levels
(one absolute offset, not a shape problem), and nothing else.

## 7. Candidates — measured features the objective does not score

| feature | `s` | z | cost if added at its χ² weight | share |
|---|---|---|---|---|
| **`rsT0`** | 23.9 % | **2.19** | **2.47** | **64 %** |
| `vall2_t` | 32.7 % | 1.66 | 0.166 | 4.3 % |
| `rsA` | 24.6 % | 0.29 | 0.052 | 1.3 % |
| `Am` | 7.9 % | 0.11 | 0.024 | 0.6 % |
| `peak1_t` | 63.6 % | 0.31 | 0.020 | 0.5 % |

`rsT0` — the *delay* before the post-re-stretch recovery begins — is **≈6 preparation-widths
off** and is a strong ATP feature (×0.43). It is the same defect as `rsK`, seen on the
time axis rather than the rate axis. It is deliberately **not** added to the default
objective: it would take 64 % on its own and break cost-neutrality. Add it only together
with a mechanism that can move it, or the optimiser will simply pay for it elsewhere.

## 8. Limitation — the honest one

The 2 mM data enters here only as a **weight**. The objective still scores 8 mM alone:
`RunBakersExp` on `main` has no 2 mM branch, so nothing in the cost function actually
*measures* an ATP ratio. The `R` multiplier conditions the baseline for the ATP fit; it
does not perform it.

The complete treatment needs a second slack run at 2 mM inside the same evaluation and
`fn` entries scoring the **ratio** `f₂/f₈` against the measured ratio, weighted by
`1/s_r²` from §3 — i.e. the three ratio features that survive the z-screen
(`ktr` 0.540, `rsK` 0.345, `rsT0` 0.428) plus the scale (`steady` 1.346), and explicitly
**no force-shape ratios**, since those are 1.0 to within noise and would only add
constants. That branch does not exist yet.

## 9. Reproduce

```matlab
run('Analyses/FeatureWeighting/RunFeatureSpread.m')        % data only, ~1 min
run('Analyses/FeatureWeighting/BuildWeightedObjective.m')  % one model eval, ~2 min
DRYRUN = true; run('Analyses/FeatureWeighting/RunOptimSpreadW.m')   % pre-flight
```
