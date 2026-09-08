# DataPooling — how to normalise and pool three preparations, and what that is worth

**Headline.** Between-preparation force variation **is one common scale** (88 % of the
log-variance in PC1, loading flat to 0.955), so one divisor per preparation is the right
model and not just a convenience. Dividing by the pre-slack isometric force at SL 2.2
takes the spread of the force features from **29 % to 15 %** (RMS) / **27 % to 10 %**
(median) and the raw trace from **25 % to 8 %**. **Which** divisor barely matters —
five of the six candidates land within 1.8 pp of each other and a paired bootstrap
cannot separate them. What matters, and is not a matter of taste, is **scope**: the
divisor must be per *preparation*, not per *run*, or the ATP effect is defined away.

> **The number that reframes the exercise.** With *n* = 3, every spread estimate has
> **two degrees of freedom**: a feature measured at *s* = 10 % is only known to lie
> between **5 % and 63 %**. Choosing between two normalisations that differ by 1.3 pp is
> not a decision this dataset can support. Choosing between normalising and not (29 % vs
> 15 %) is.

| script | what it answers | output |
|---|---|---|
| `poolLoadRuns.m` | one code path: traces, anchors, fresh features for all nine runs | `results/pool_raw.mat` |
| `RunNormalizationLadder.m` | which normalisation, at which scope, with which corrections | `results/normalisation_ladder.mat` |
| `RunPoolingComparison.m` | timecourse spread; extract-then-average vs average-then-extract; noise vs biology | `results/pooling_traces.png` |
| `RunFVConsistency.m` | does the in-preparation F–V point rescue the borrowed F–V curve | `results/fv_consistency.png` |
| `buildPooledTarget.m` | **the deliverable**: pooled target + honest error bar per feature | `results/pooled_target.mat` |
| `RunSummaryFigure.m` | the argument on one page | `results/datapooling_summary.png` |
| `WritePooledTarget.m` | freeze the target into version control, next to what it replaces | `../../Model/pooledTarget8mM.m`, `results/pooled_target_report.txt` |
| `BuildPooledSeed.m` | put that target *into the objective* and reweight from its spreads | `../../params/pooled8_seed.m` |
| `RunOptimPooled.m` | refit against it (engine: `../FeatureWeighting/RunOptimSpreadW.m`) | `params/pooled8_opt.m`, `params/pooled8_iter/` |

**[methods.md](methods.md) is the paper-ready paragraph** — 167 words, and it is meant
to be the entire data-manipulation section.

![summary](results/datapooling_summary.png)

---

## 1. The nuisance, named

Tension is force ÷ *estimated* cross-sectional area. Across the three 2026 days the
pre-slack isometric force at SL 2.2 is **79.8 / 70.0 / 50.9 kPa** — a CV of **22 %**.
The generative model this implies is

```
    y[p,c,k]  =  alpha_p^{d_k} * mu[c,k] * exp(eps)
```

`p` preparation, `c` condition, `k` feature; `d_k` = 1 for force-dimension features and
0 for rates, times and ratios. `alpha_p` is a property of the *preparation* — the same
fibre, the same area estimate, both ATP runs.

**Two facts make this testable rather than assumed.**

*One scale is enough.* Centre the log force features by feature and take the SVD:
**PC1 holds 88.3 %** of the between-preparation variance and its loading is flat
(|cos| with the all-ones direction = **0.955**). A pure common scale is exactly a
rank-1 matrix with a flat loading. Both conditions hold.

*Scale-free features are the control group.* `ktr`, `t0`, `rsK`, `rsR2` carry `d_k` = 0,
so no normalisation can touch them. Their spread (11–38 %) is the floor: after
normalisation the force features sit at **4–11 %**, i.e. *better* than the best rate.
Absolute tension goes from the worst-determined quantity in the experiment to among the
best.

## 2. The ladder — six anchors × two scopes × rundown × passive

Every scheme is scored **leave-one-feature-out**: the divisor is re-estimated with
feature *k* held out before *k*'s spread is measured, so a divisor cannot shrink the
spread of the feature it was drawn from.

| scheme | sBar (RMS) | sMed | var. explained | ATP ratio ± sd | z |
|---|---|---|---|---|---|
| **`Fss22` / prep + rundown** | **14.8 %** | 10.5 % | 0.74 | **1.308 ± 0.119** | **1.63** |
| `Fsteady` / prep + rundown | 15.5 % | 9.0 % | 0.74 | 1.308 ± 0.119 | 1.58 |
| `MoR` / prep + rundown | 16.0 % | 7.7 % | 0.69 | 1.308 ± 0.119 | 1.55 |
| `Fiso20a` / prep + rundown | 16.1 % | **5.4 %** | 0.69 | 1.308 ± 0.119 | 1.54 |
| `Fiso20b` / prep + rundown | 16.6 % | 6.6 % | 0.67 | 1.308 ± 0.119 | 1.51 |
| `Fpeak` / prep + rundown | 19.0 % | 10.1 % | 0.57 | 1.308 ± 0.119 | 1.38 |
| none (raw kPa) + rundown | 23.7 % | 21.8 % | 0.33 | 1.308 ± 0.119 | 1.17 |
| none (raw kPa) | 28.9 % | 26.9 % | 0.00 | 1.340 ± 0.274 | 0.85 |

`Fss22` = pre-slack isometric at SL 2.2 (the "8 mM steady state before the slack").
`Fiso20a/b` = the quiescent ML 1.00 plateaus bracketing the slack train. `Fsteady` =
mean of the five slack-segment steady levels. `MoR` = median-of-ratios over all force
features — the reference-free estimator (DESeq size factors, and the same idea in
proteomics). Full 54-row table in the script output.

### 2.1 The anchor choice is not resolvable, and does not need to be

Paired bootstrap over the feature set, 4000 draws, all schemes scored on the same nine
features:

| scheme | sBar | P(best) |
|---|---|---|
| `Fss22`/prep+rd | 14.8 % | **64 %** |
| `Fsteady`/prep+rd | 15.5 % | 18 % |
| `Fiso20a`/prep+rd | 16.1 % | 15 % |
| `Fiso20b`/prep+rd | 16.6 % | 2 % |
| `Fpeak`/prep+rd | 19.0 % | 0.8 % |
| **none** | **23.7 %** | **0.0 %** |

*Normalise* — decisive. *Which anchor* — `Fss22` leads but cannot be separated from
`Fsteady` or `Fiso20a`. So pick it on grounds the bootstrap cannot see: `Fss22` is
measured at the **same length as the reference state**, in the **same activation**, over
~40 ms of clean plateau in the high-rate burst file, and it is **not itself a scored
feature** — the only candidate for which the leave-one-out question does not even arise.

> **`Fpeak` is the one to avoid** (19.0 %, P(best) 0.8 %). The activation peak is reached
> at different times in the two ATP conditions and the fibre is already running down
> through it (RundownCorrection §M.5b), so the peak is contaminated by exactly what the
> anchor is supposed to be free of.

### 2.2 Scope is the decision that actually matters

| scope | sBar | ATP ratio |
|---|---|---|
| per **run** (each condition ÷ its own level) | 17.3 % | *not measurable* |
| per **preparation** (both conditions ÷ the 8 mM level) | 18.6 % | 1.31 ± 0.12 |

A per-run divisor buys 1.3 pp of spread and costs the entire ATP force effect: it
divides each condition by its own level, so whatever it then reports about the level
ratio is the gauge, not the data. **`alpha` is a property of the fibre, so it is shared
by both runs.** This is the one place where "lowest spread wins" gives the wrong answer,
and it is why the criterion in this analysis is spread *and* preserved discriminability.

### 2.3 Rundown is invisible within a condition and decisive across conditions

The lift applied to the run recorded second: 03/27 **+3.4 %**, 04/03 **+6.4 %**,
04/10 **+14.6 %** (2 mM first — it decays ~2× faster and degrades the weaker 8 mM run).

Under per-preparation scope the lift multiplies both the anchor and the features, so it
**cancels exactly** in the 8 mM spread (14.8 % either way). Across conditions it does
not: the ATP ratio's scatter falls **0.274 → 0.119** and z goes **1.09 → 1.63**. And it
explains **50 % of the apparent between-preparation force variance**: the anchor CV falls
**22.0 % → 15.5 %** once each session's second run is lifted.

### 2.4 Subtracting the measured passive does not help overall

The PNB + mavacamten run of each day gives the passive contribution to every feature, at
the right length and time, measured rather than modelled — **6.3 / 12.1 / 8.6 %** of that
day's 8 mM `Fss22`. Subtracting it splits the field:

| helps | hurts |
|---|---|
| `rsA` 23.2 → **8.8 %**, `A` 10.6 → 9.8 %, `vall_y` 10.2 → 9.2 % | `peak2` 7.5 → **11.9 %**, `restretchSlopeStart` 20.1 → **33.5 %**, `vall2_dy` 22.9 → 28.3 % |

Net sBar 14.8 → 16.9 %. The pattern is physical: subtraction helps where passive is a
small quasi-static offset, and hurts on the **re-stretch transient**, where the elastic
and cross-bridge contributions are not additive. **Do not subtract passive as a blanket
correction**; the model has a passive branch and should carry it explicitly.

## 3. What is left after the scale is a LENGTH effect

The trace panel refuses to be a scale story. After dividing by `Fss22` the three
preparations agree to **2.2 %** during the post-re-stretch recovery (at the reference
length ML 1.10) and still differ by **11.0 %** during the slack **hold** (ML 0.94–1.02).
A scale is the same at every length; this is not.

Fit **one length offset per preparation on `A` alone**, gauge Σd = 0, then apply it
unchanged to four other features:

| feature | scale only | + length offset | |
|---|---|---|---|
| `A` | 10.3 % | **0.4 %** | fitted |
| `Am` | 10.1 % | 0.8 % | held out *(nearly the same measurement as `A`)* |
| `vall_y` | 9.9 % | **5.7 %** | held out |
| `peak1_y` | 9.0 % | **4.4 %** | held out |
| `peak2` | 6.8 % | 7.5 % | held out — **at the reference length, where it must do nothing** |

Offsets **+0.048 / −0.060 / +0.011 µm SL**. Two things make this more than curve-fitting:
the genuinely independent held-out features (`vall_y`, `peak1_y`, measured in a different
phase of the trace) improve by ~40 %, and the one measured *at the reference length*
does not — which is the sign convention's own control.

The size matches RundownCorrection's series-elastic creep shift (0.088–0.120 µm), and
`SLslack` is **derived from commanded motor length, not measured**
(PassiveRundownATP). So the offset is either L₀-calibration or accumulated compliance,
and *n* = 3 cannot separate them.

> **Recommendation: report it, do not fit it — yet.** One more per-preparation parameter
> on three preparations spends the dataset. It is quantified here (~5 % residual, the
> largest single term left) and it is a design fix, not an analysis fix: record the L₀
> criterion, or measure SL.

## 4. Pool the FEATURES, never the traces

> **CORRECTED 2026-09-04 by [RestretchMetricAudit](../RestretchMetricAudit/conclusions.md).**
> The −24 % below is real as a number but wrong as a diagnosis. It is **one
> segment**, not a smearing bias: on the pooled trace segments 1–4 give
> `rsK` = 48.8 / 45.3 / 50.9 / 48.3 against a mean-of-features of 53.5 (a few per
> cent), and segment 5 gives **9.9**. Segment 5's feature window is 844 ms, not
> 279 ms — `velocity_segment(5)` there is the end of the protocol at t = 78 s —
> so its plateau, which sets the fit amplitude, lies outside the 570 ms splice
> grid used to build the grand average and is one preparation, rescaled, with a
> 1.5 kPa step at the seam. Averaging the re-stretch window *directly* biases
> `rsK` by −2.0 %.
> The recommendation stands (features, not traces — extraction is non-linear and
> the trap above is easy to re-create), but the *cost* of trace pooling on the
> re-stretch is a few per cent, not a preparation width. And the audited
> `rsK80` supersedes `rsK` throughout: between-preparation SD 17.2 % → 6.1 %,
> ATP `z` 3.2 → 6.8.

Feature extraction is non-linear — exponential fits, peak picking, argmax. Averaging
traces first is a different estimator, and for the re-stretch window a badly biased one:

| feature | mean-of-features | features-of-mean-trace | gap | gap / s |
|---|---|---|---|---|
| **`rsK`** | 53.5 | 40.6 | **−24.0 %** | **1.19** |
| **`rsT0`** | 0.00544 | 0.00420 | **−22.8 %** | **1.11** |
| `rsR2` | 0.734 | 0.793 | +8.0 % | 0.70 |
| `peak1_y` | 0.970 | 1.040 | +7.2 % | 0.70 |
| `A`, `Am`, `ktr`, `t0`, `peak2`, `steady` | | | 0.0–2.3 % | ≤ 0.44 |

The slack-hold recovery survives averaging (`ktr` −1.5 %, `A` 0.0 %) because it is slow
and starts from a hard zero. The **re-stretch** recovery does not: its onset jitters
between preparations, averaging smears the transient, and the fitted rate drops by more
than a preparation width.

> That is the feature the objective is currently failing on and the one `rsK` weights
> most heavily ([FeatureWeighting](../FeatureWeighting/conclusions.md) §6). Pooling
> traces would have manufactured a 24 % error in exactly the wrong place.

*(A practical trap found on the way: `extractSlackAttributes` carries absolute
thresholds in kPa — `FORCE_THRESHOLD_KPA = 10` selects the active window — so a
normalised trace silently returns NaN for every fitted quantity. Normalise the
**features**, or put the pooled trace back on a kPa axis before extracting.)*

## 5. Noise or biology? Mostly biology — so it is preparations we need

Within a run the five slack segments differ only by `SLslack`, so the residual about a
linear fit on `SLslack` is that run's own repeatability. Only `s_meas`/√5 of it reaches
the between-preparation SD.

| feature | s_total | s_meas | contributes | biology share |
|---|---|---|---|---|
| `A`, `t0` | 10.3 %, 37.4 % | 0.4 %, 2.0 % | 0.2 %, 0.9 % | **100 %** |
| `ktr`, `Am` | 11.3 %, 10.0 % | 2.0 %, 2.4 % | 0.9 %, 1.1 % | 99 % |
| `rsK`, `vall_y` | 20.2 %, 9.8 % | 10.1 %, 6.0 % | 4.5 %, 2.7 % | 95 %, 93 % |
| `rsA`, `vall2_dy` | 22.9 %, 22.7 % | 18.6 %, 18.8 % | 8.3 %, 8.4 % | 87 %, 86 % |
| `peak2`, `peak1_y` | 7.0 %, 10.3 % | 8.5 %, 14.2 % | 3.8 %, 6.4 % | 71 %, 62 % |
| `rsT0` | 20.5 % | 29.0 % | 13.0 % | 60 % |

**The features the objective leans on are 95–100 % biological.** No amount of cleaner
extraction moves them; only more preparations do. The exceptions — `peak1_y`, `rsT0`,
`peak2` — carry 30–40 % extraction scatter and *would* respond to a better estimator or
more segments.

### How many measurements do we have, really?

Three preparations. Two degrees of freedom on every spread estimate. **95 % CI on any
*s*: [0.52 s, 6.28 s].** Hence:

* **Moderate the variances.** `s²_mod = (d₀s₀² + d s²)/(d₀+d)` with the prior estimated
  from the observed dispersion of *s²* (Smyth 2004 / limma `fitFDist`), not fixed. The
  estimator returns **d₀ = 2.38** here — prior worth slightly more than the data — and
  it *self-disables* (d₀ → 0) when the features are genuinely not exchangeable, so it
  cannot manufacture agreement. It collapses the weight dynamic range **446× → 15×**.
* **A fixed prior would be wrong.** Shrinking toward a class median pulls `ktr` (11 %)
  toward a class that also holds `vall_t` (84 %). The moment-matched prior is the
  difference between a defensible step and a fudge.

### Is fitting one preparation the same as fitting the pool?

χ² of each preparation against the pooled mean, per feature at 8 mM (the expectation for
three preparations scattering around their own mean is (n−1)/n = **0.67**):

| preparation | χ² | worst feature |
|---|---|---|
| **03/27** — the one the objective currently fits | **0.39** | `vall_t` (1.1) |
| 04/03 | 0.80 | `peak1_dSL` (2.2) |
| 04/10 | 0.75 | `vall2_t` (1.8) |

**03/27 is the most typical of the three**, not an outlier — so the current single-prep
target is not importing a bias, and a refit against the pool should not move the optimum
much. What pooling adds is not a different target; it is **an error bar on every number
in it**, and with it the right to say a residual is or is not significant.

## 6. The in-preparation force–velocity point

`DataCuration/CreateFVRampProtocolVelocityTable.m` built FV2 for 04/10 only. **04/03
runs the identical manoeuvre** — −2.005 ML/s, 34.6 ms, 1.060 → 0.990 ML — in all three
of its conditions. That is **six FV2 measurements, not three**; the configs are now
extended and the 04/03 files are built.

### It validates the borrowed 8 mM curve

| | Baker (2 source sets) | FV2 (2 preparations) | separation |
|---|---|---|---|
| 8 mM, raw | 0.298 ± 0.025 | 0.284 ± 0.062 | **0.21 sd** |
| 8 mM, passive-corrected | 0.298 ± 0.025 | 0.250 ± 0.078 | 0.58 sd |
| 2 mM, raw | 0.145 ± 0.031 | 0.202 ± 0.041 | **1.10 sd** |

P/P₀ at 2 ML/s. **At 8 mM the borrowed curve and the in-preparation point agree to a
fifth of a noise width** — the F–V data is transferable, which until now was assumed.
The ramp is settled: force drifts −1.7 % (04/10) and −3.7 % (04/03) across the 34 ms
window at 8 mM.

### It puts F–V in the same normalisation as everything else

FV2's isometric reference *is* `Fss22` — the same anchor the slack features use. So
P/P₀ and `steady`, `A`, `peak2` now share one denominator, which is precisely what the
scale-blindness warning in FeatureWeighting §4 asks for.

### It corrects two numbers

* **F–V is less reproducible than assumed.** FeatureWeighting takes F–V's spread
  *between the two Baker source sets* (8.5 % at 2 ML/s). Between two actual
  **preparations** it is **21.9 %**; all four together, **13.6 %**. The number that
  applies when the model must predict a new fibre is ~14 %, not ~8 %.
* **The borrowed curve carries a much stronger ATP effect on F–V shape than our own
  preparations do.** In-preparation, P/P₀(2 mM) / P/P₀(8 mM) = **0.719 and 0.705** —
  remarkably tight for two preparations. The Baker sets give **0.528 and 0.441**. As a
  fractional depression of P/P₀ that is **−29 %** measured here against **−51 %**
  borrowed: nearly twice the effect. At FV_SHARE = 0.70 most of the objective is being
  steered by the larger of the two.

### Passive is not negligible at 2 ML/s

The PNB + mavacamten run drives the same ramp and bears **56–58 % of its own isometric
force** — 2.4–4.8 kPa out of a 14.5–19.1 kPa total, i.e. **14–33 %** of the ramp force.
Passive-correcting moves 8 mM P/P₀ from 0.284 to 0.250. Whichever convention is chosen,
it has to be stated; it is not a rounding effect.

## 7. What to do

1. **Adopt `Fss22` / per-preparation / rundown-corrected**, features pooled after
   extraction, spreads moderated. `buildPooledTarget.m` emits it; the methods paragraph
   is [methods.md](methods.md).
2. **Keep the level as one term.** `F_ref` = 69.3 kPa, CV 15.5 %. Without it the
   objective is scale-blind by construction.
3. **Score ATP as ratios, and only five of them**: level ×1.320 ± 8.6 %,
   **`rsK80` ×0.381 (z = 3.42)**, `vall2_t` ×1.924 (2.58), **`rsT20` ×1.752
   (2.36)**, `ktr` ×0.532 (2.32). Every force *shape* ratio is 1.00 ± 0.05 —
   adding those nine terms adds nine constants.
   (`rsK80`/`rsT20` replace `rsK`/`rsT0`: the audited split turns one resolved
   effect into two, because `rsT0`'s 77 % spread resolved nothing while the
   measured 20 % crossing time has 23 % and a real ×1.75 —
   [RestretchMetricAudit §6](../RestretchMetricAudit/conclusions.md).)
4. **Do not subtract passive globally**; let the model's passive branch carry it.
5. **Revisit FV's weight.** Its honest between-preparation spread is 13.6 %, not 8.5 %,
   and its 2 mM ATP effect disagrees with our own preparations by 1.1 sd.
6. **Next dataset:** the binding constraint is *n*, not measurement quality — 95–100 %
   of the spread in the features that matter is biological. Record the L₀ criterion (or
   SL) to close the ~5 % length-offset residual, and keep the FV2 manoeuvre: it is free
   and it is the only F–V measurement tied to these preparations.

## 8. The target, actually in the objective (2026-09-08)

§7 said what to adopt. This section is what adopting it did. Three files, one switch:

* **`WritePooledTarget.m` → `Model/pooledTarget8mM.m`.** `buildPooledTarget` needs the
  six curated protocol files and `data/` is not in git, so the pooled numbers are
  frozen once into a generated function. Regenerate it when a preparation lands — that
  is the point: the target moves as data comes in, and the file header records which
  preparations and which options produced it.
* **`params0.UsePooledTarget`** (default `false`) routes the slack block through
  `Model/pooledTargetFeatures`, which *replaces* `features_data` rather than adding a
  penalty. A pooled target is a different value for the same quantity, so every
  existing `fn` entry, weight and diagnostic keeps working. Modes: `true`/`'all'`,
  `'shape'` (pooled pattern on the incoming force level), `'level'` (force block only).
* **`BuildPooledSeed.m` → `params/pooled8_seed.m`** takes `s_f` from the same pooled
  analysis (`SPREAD_SOURCE = 'pooled'` in `BuildWeightedObjective`) and renames the
  `rsK` term to `rsK80`, so the audited target is scored with the audited estimator.

> **The switch is sticky, and that is deliberate.** `writeParamsToMFile` serialises
> every field, so `params/pooled8_*.m` carry `UsePooledTarget = 1` and re-score
> themselves against the pool wherever they are loaded. A snapshot that names its own
> target is the right default — but it means **a cost quoted for a `pooled8` snapshot is
> not comparable to one quoted for `estiff125_opt`** unless you say which target. To
> compare on the old yardstick, set `params0.UsePooledTarget = false` explicitly after
> `loadParams`. (`getParams` defaults it to `false`, so snapshots predating this change
> keep scoring the single preparation — nothing silently moved.)

### What the level step costs, and what it hides

`estiff125_opt` — the best fit to 03/27 alone — scores **0.955 against 03/27** and
**4.431 against the pool**. The tempting reading is "that's the −13 % level step"
(`F_ref` 69.35 vs 03/27's `Fss22` 79.76). It is not. In χ² units — cost over the
measured between-preparation spread, **per slack segment**, so 1.0 = one
preparation-width:

| term | χ²/segment | term | χ²/segment |
|---|---|---|---|
| **`peak1_y`** | **6.6** | `peak2` | 1.3 |
| **`peak1_dSL`** | **5.2** | `vall2_dy` | 1.3 |
| **`vall_t`** | **4.4** | `steady` | 1.3 |
| `A` | 2.3 | `restretchSlopeStart` | 0.4 |
| `vall_y` | 1.7 | `t0` 0.3 · `ktr` 0.15 · `rsK80` 0.06 · FV 0.02 | — |

Eight of sixteen terms are already at or below one preparation-width, so there is no
broad gain to be had. The residual is **the shape and timing of the re-stretch
transient**, and `peak1_y` is most of it.

**Why 03/27 alone could not see this.** The pool's first re-stretch peak sits **3 %
below** its own isometric; 03/27's sits **12 % above** it. Fitting 03/27 was quietly
*rewarding* the model's lengthening overshoot — the defect flagged at the end of the
spreadw campaign and never scored. Note this is not 03/27 being a bad preparation: its
χ² against the pooled mean is 0.40 where 0.67 is expected (§5). It is one feature where
a single fibre's idiosyncrasy happened to point the same way as a model artefact.

**`rsK80` = 0.06 and `ktr` = 0.15 are essentially exact**, which is the other half of
the reading: the post-re-stretch *rate* is solved (spreadw's win, confirmed against the
pooled mean and with the audited estimator). Any refit that improves `peak1_y` while
those climb is trading a solved feature for a broken one.

### The FV pin moved too — deliberately, and it confounds

`FV_SHARE` 0.70 → **0.40**. At 0.70 the fourteen non-FV data terms held 6.6 % of the
objective between them, so the pooled target's entire contribution would have been
unbuyable. FV's own residual is χ² 0.02/point — inside the disagreement between its two
source curves. The consequence for attribution is real and stated here rather than
buried: this campaign changes the target *and* the pin. §5's `chi2Prep` is the control
that says the target change alone should not have moved the optimum much.

### Held out on purpose

`rsT20` is **not** added. It is the audited re-stretch delay and the model's worst
single feature, and the weight table prices it at 2.7 % — but `rsT0` and `peak1_t`
price at **35 %** each. Adding a term while the target moves makes the result
unreadable. Next campaign: `ADD_FEATS = {'rsT20'}`.

## 9. Limitations

* *n* = 3. Every spread here has a 12× confidence interval. The **rankings** in §2 are
  not resolvable; the **presence or absence** of normalisation is.
* `steady` normalised by `Fss22` is close to a self-normalisation — both are isometric
  force at ML 1.10, seconds apart. Its raw 4 % is a within-run repeatability, not a
  between-preparation spread. Variance moderation happens to fix this (4 % → 15 %); an
  unmoderated objective would hand it the largest weight in the set.
* The FV2 ramp is 34 ms entered through a step, not a sustained isovelocity release.
  Drift across the window is small at 8 mM (−1.7 to −3.7 %) but reaches −14 % for
  04/03's 2 mM run, so the 2 mM FV2 numbers are the softer ones.
* The rundown lift uses φ = 0.452 from a single bracket on 03/27. The correction's
  *direction* is not in doubt and its size is bounded (RundownCorrection §3 finds the
  answer robust over φ = 0.45–0.62), but φ has one calibration.
* The Baker slack recordings remain excluded on
  [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md)'s grounds; the Baker **F–V**
  curves are kept and are now validated at 8 mM by §6.
