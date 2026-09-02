# OptimResultReport — the spread-weighted run (`spreadw`), read against three preparations

**Headline.** `3.8513 → 0.9786`, **−74.6 %** on the spread-weighted objective. But
**87.6 % of that gain is force–velocity alone**, which is what pinning FV at 70 % told
the optimiser to buy. It paid in exactly the places the pin left cheap: the re-stretch
peak and force-redevelopment speed. Two protocols that were never scored — the `ktr`
protocol and the staircase — say the same thing, and the staircase says it *louder* than
before the refit.

> **The one real win:** `rsK`, the post-re-stretch recovery rate that has been the open
> structural defect for months, now sits **on the three-preparation mean** — z = −0.2
> against it, and z = +0.7 against 03/27 itself, which is the lowest-`rsK` prep of the
> three and therefore the harder test. The previous fit sat at ×1.26 of the envelope mean.

Scripts: [`RunDataEnvelope.m`](RunDataEnvelope.m) (data only) → [`RunFitReport.m`](RunFitReport.m)
→ [`MakeSlideFigures.m`](MakeSlideFigures.m) → [`MakeSlide.m`](MakeSlide.m) (a `.pptx`)
→ [`BuildReport.m`](BuildReport.m) (single-file HTML). All TAG-driven; re-run for any
future tag with `TAG_NEW`/`TAG_BASE`.
**HTML report: [`results/report_spreadw_opt.html`](results/report_spreadw_opt.html)** —
1.6 MB, all figures inlined, nothing external to lose.
**Slide: [`results/slide_spreadw_opt.pptx`](results/slide_spreadw_opt.pptx)**.

---

## 0. Provenance, and the one iteration that does not reproduce

The run was executed off this machine; only `params/spreadw_iter_/` came back — **no
`spreadw_state.mat`**, so there is no `optIter`/`featCost` history and no `_opt.m`. The
gallery was copied to `params/spreadw_iter/` and **iteration 18** written to
`params/spreadw_opt.m` with its provenance in the header.

Re-evaluating each iteration here against its recorded header cost:

| iter | header | local | Δ |
|---|---|---|---|
| 1 | 3.4470 | 3.4520 | +0.1 % |
| 4 | 2.0920 | 2.1153 | +1.1 % |
| 7 | 1.7811 | 1.7883 | +0.4 % |
| 10 | 1.6630 | 1.6686 | +0.3 % |
| 12 | 1.6231 | 1.6333 | +0.6 % |
| **13** | **1.4769** | **1.6615** | **+12.5 %** |
| 14 | 1.3344 | 1.3362 | +0.1 % |
| 15 | 1.0626 | 1.0674 | +0.5 % |
| 16 | 1.0231 | 1.0235 | +0.0 % |
| 17 | 1.0010 | 1.0035 | +0.3 % |
| **18** | **0.9715** | **0.9786** | **+0.7 %** |

Seventeen of eighteen reproduce inside 1.1 %. **Iteration 13 does not** — and the model is
deterministic here (evaluated twice, identical to the digit). Only **9 scalars move**
between iter 12 and 13, all by 0.4–1.6 % (`eta_M` +1.6 %, `v_ref_reg` −1.6 %, `kSE`
−1.3 %, `ka` +0.4 %). A 1 % parameter step cannot honestly produce a 9 % cost drop.

> **That step was a numerical wobble the run machine pocketed, not a gain.** It does not
> contaminate the result — the search climbed past it and round 14 resumed from a genuine
> point — but `optimizeFeatures` has **no repeat-evaluation guard on accept**, and this is
> the second time the objective's roughness has bought a fake step.

**Not converged.** A 20-round stall (rounds 14–33) broke, then five straight improvements
to round 38, the last worth −2.9 %, still descending when the budget ran out.

## 1. Where the −2.87 went

Both snapshots scored on one common `fn` (realign4's own list still carries the
duplicated `rsR2`; it costs 0, so its headline is unchanged).

| feature | realign4 | spreadw | Δ | share of the gain |
|---|---|---|---|---|
| **`FV_normpowerAvg`** | 2.6959 | 0.1803 | **−2.5156** | **87.6 %** |
| **`rsK`** | 0.4807 | 0.0635 | **−0.4171** | **14.5 %** |
| `coolDownLS` | 0.1347 | 0.0797 | −0.0550 | 1.9 % |
| `doublePeak`, `peak2`, `ktr_rmse`, `PS_*` | | | −0.03 | 1.2 % |
| **`peak1_dSL`** | 0.0082 | 0.0448 | **+0.0366** | −1.3 % |
| **`restretchSlopeStart`** | 0.0184 | 0.0492 | **+0.0308** | −1.1 % |
| **`peak1_y`** | 0.0175 | 0.0399 | **+0.0224** | −0.8 % |
| `ktr` | 0.0159 | 0.0349 | +0.0190 | −0.7 % |
| `vall2_dy`, `PS_holdDecayRMSE`, `vall_t`, … | | | +0.03 | −1.3 % |
| **total** | **3.8513** | **0.9786** | **−2.8727** | |

![per-feature cost shift](results/fit_costshift.png)

**The FV gain is real, not bookkeeping.** `realign4` fell *below* the band spanned by the
two FV source curves at 2 ML/s (0.502 vs data 0.596); `spreadw` now tracks the data mean
across the whole power curve, at χ² well inside the two sources' own disagreement. The
long-standing FV shoulder is fixed.

![force-velocity](results/fit_fv.png)

## 2. The model against three preparations, not one

`RunDataEnvelope.m` builds a proper reference. All three 2026 days share 4 kHz sampling,
5 cycles at 0.600 s, the same release velocities and the same 4/5/6/7/3 ML/s re-stretch
ladder — only the session-time offset differs — so aligning on each release is a pure
index shift. Every trace is divided by its **own** mean `steady` (77.3 / 64.6 / 47.9 kPa)
because absolute force is not transferable but shape is. The model is normalised the same
way, so nothing here can be bought by moving the force level.

**The envelope is tight: SD is 3.2 % of the mean (median), 6.2 % at p90.** It is a real
yardstick, not a smear.

![model on the timecourse envelope](results/fit_envelope_overlay.png)

Distance from the envelope, in its own SDs, averaged over the five cycles. **03/27 is
listed separately** — it is the only preparation the model was ever shown, so a miss
against the *mean* is a defect only if it is also a miss against *that* trace.

| feature | 3-prep mean | 03/27 | realign4 | spreadw | z vs mean | z vs 03/27 |
|---|---|---|---|---|---|---|
| `rsT0` | 0.00544 | 0.00494 | 0.00206 | 0.00213 | **−4.1** | **−3.9** |
| `vall_t` | 0.00665 | 0.00800 | 0.01317 | 0.01443 | **+2.9** | **+2.7** |
| `peak1_dSL` | 0.0178 | 0.0176 | 0.0307 | 0.0400 | **+2.8** | **+2.8** |
| **`peak1_y`** | 1.027 | 1.065 | 1.211 | **1.233** | **+2.2** | +1.7 |
| `ktr` | 47.85 | 49.21 | 50.49 | 52.93 | +1.0 | +0.7 |
| `steady` | 63.3 | 77.4 | 77.0 | 77.2 | +0.9 | **−0.0** |
| `restretchSlopeStart` | 1324 | 1588 | 1435 | 1235 | −0.8 | −1.5 |
| `rsA` | 0.176 | 0.164 | 0.180 | 0.188 | +0.7 | +0.8 |
| `vall2_dy` | −0.175 | −0.163 | −0.180 | −0.188 | −0.7 | −0.8 |
| `Am` | 0.754 | 0.800 | 0.807 | 0.789 | +0.6 | **−0.2** |
| `A` | 0.760 | 0.804 | 0.809 | 0.789 | +0.4 | **−0.3** |
| `vall_y` | 0.868 | 0.902 | 0.903 | 0.890 | +0.4 | **−0.1** |
| **`rsK`** | 53.47 | **43.74** | 67.11 | **51.32** | **−0.2** | **+0.7** |
| `peak2` | 1.019 | 1.014 | 1.041 | 1.016 | −0.1 | +0.1 |

> **The last two columns are why marking 03/27 was worth doing.**
>
> Four features that look mildly off against the mean are essentially **exact** against
> the prep the model was fitted to: `steady` +0.9 → **−0.0**, `Am` +0.6 → **−0.2**,
> `A` +0.4 → **−0.3**, `vall_y` +0.4 → **−0.1**. Those are not defects — 03/27 is simply
> the strongest and stiffest of the three, and the model is sitting on it.
>
> **`rsK` is now inside on both readings** — z = −0.2 against the mean and **+0.7 against
> 03/27**, which is the *lowest*-`rsK` prep of the three (43.7 vs a mean of 53.5) and so
> the harder test. realign4 was at 67.1. This is the first fit for which the
> post-re-stretch recovery rate is not a defect.

**Ten of fourteen features are now inside 1 SD of the envelope.** What is left is a tight
cluster: `rsT0` (−4.1), `vall_t` (+2.9), `peak1_dSL` (+2.8), `peak1_y` (+2.2). `vall_t` and
`peak1_dSL` carry 68 % and 64 % between-prep spread — barely measured, which is why the
weighting de-rated them, and both got *worse* here because nothing defends them. `rsT0` and
`peak1_y` are the real ones, and §3 shows they are the same defect.

![features against the envelope](results/fit_features.png)

## 3. The predictions — and they agree with each other

Neither protocol was scored during the fit (`RunKtr = 0`, `RunStairs = 0`), so these are
genuine out-of-sample tests.

![ktr prediction](results/fit_ktr_pred.png)
![staircase prediction](results/fit_stairs_pred.png)

| observable | data | realign4 | spreadw |
|---|---|---|---|
| `ktr` protocol, re-stretch overshoot | 0.724 | 1.000 (+38 %) | 1.006 (**+39 %**) |
| `ktr` protocol, t to 90 % of final | 31.5 ms | 12.6 ms (**2.5× fast**) | 12.4 ms (**2.5× fast**) |
| staircase, mean step-peak overshoot | — | +7.3 % | **+9.3 %** |

> ### One defect, seen three ways — and the refit did not touch it
>
> **The model's transient response to imposed LENGTHENING is systematically too large:**
> slack re-stretch `peak1_y` **+20 %** over the envelope (z = +2.2), `ktr`-protocol
> re-stretch overshoot **+39 %**, staircase step peaks **+9.3 %**. Three independent
> protocols, one sign, one mechanism.
>
> **And force redevelopment after a length change is systematically too fast:** `ktr`
> protocol **2.5×** too fast, and `rsT0` — the delay before recovery starts — is
> **61 % too short** and the single worst feature in the set at z = −4.1.
>
> **The staircase got worse** (+7.3 % → +9.3 %), and the `ktr` overshoot marginally so.
> That is the honest cost of this run: nothing in the objective defends them. `peak1_y`
> holds ~4 % of it, `ktr` ~4 %, and the two protocols above 0 %.

## 4. What this says about the objective

The run did exactly what it was told, and the pin worked: FV was genuinely broken and is
now genuinely fixed. But **the yardstick has now served its purpose.** With FV's cost down
to 0.18 it still holds 18 % of the total while sitting well inside its data's own
disagreement — the objective is now spending its largest single weight defending a solved
problem, and `rsK`, the other thing it bought, is also done.

Everything that is left is either unmeasurable (`vall_t`, `peak1_dSL` — 68 % and 64 %
spread) or unscored (`rsT0`, the two prediction protocols). **That is the definition of an
objective that has run out of information.**

Three concrete options, in the order I would take them:

1. **Re-pin FV lower and rebuild on `spreadw_opt`** — `FV_SHARE = 0.40` in
   `BuildWeightedObjective.m`. Holding 70 % on a solved feature is the main reason the
   last 20 rounds stalled.
2. **Score `rsT0`.** The worst feature in the set (z = −4.1), data-backed, and a strong
   ATP feature (×0.43 — see [FeatureWeighting §7](../FeatureWeighting/conclusions.md)).
   It is the time-axis twin of the `rsK` defect this run just closed, which is direct
   evidence the mechanism *can* move: `rsK` went 67.1 → 51.3 without `rsT0` following, so
   the two are not welded and scoring `rsT0` should be able to separate them.
3. **Put the lengthening overshoot in the objective.** The one defect three independent
   protocols agree on, and it got *worse* here because nothing scores it but `peak1_y` at
   4 %. Turning on `RunKtr`/`RunStairs` makes them constraints — at the cost of losing the
   only out-of-sample tests this project has.

## 5. Reproduce

```matlab
run('Analyses/OptimResultReport/RunDataEnvelope.m')   % data only, ~40 s
run('Analyses/OptimResultReport/RunFitReport.m')      % two full runs, ~4 min
run('Analyses/OptimResultReport/MakeSlideFigures.m')  % 4 slide-grade PNGs, ~5 s
run('Analyses/OptimResultReport/MakeSlide.m')         % -> results/slide_<tag>.pptx
run('Analyses/OptimResultReport/BuildReport.m')       % single-file HTML, ~3 s
```

`REUSE = true` before `RunFitReport.m` replays the saved run instead of re-simulating —
for figure and axis tweaks, which otherwise cost four minutes of ODE solving each time.

Two generic pieces came out of this, both in `Auxiliary/` and reusable anywhere in the
project:

* **`buildHtmlReport.m`** — fills `{{PLACEHOLDER}}`s in an HTML template, inlines every
  value that is a path to an image as a base64 data URI, and **errors** on any placeholder
  left unfilled, so a template typo fails loudly instead of leaving a hole in the page.
* **`buildPptxSlide.m`** — writes a single-slide `.pptx` from images and text boxes by
  emitting the OOXML parts and zipping them. No Report Generator toolbox (not licensed
  here), no Office automation, no Python. Verified: 11 XML parts well-formed, every
  relationship id resolved, every target present, every part declared in
  `[Content_Types].xml`.
