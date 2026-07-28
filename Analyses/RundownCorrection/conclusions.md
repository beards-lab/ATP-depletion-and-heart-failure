# RundownCorrection — how the preparation degrades, and how to correct for it

**Headline.** Rundown is **linear in effective activated time**, not wall-clock time;
only **φ ≈ 0.45** of each recording's within-run force decline is permanent; and
mechanically it behaves as **creep in the series elastic** — the element gets
~35 % softer and ~0.1 µm longer — which alone reproduces the force loss, the
`ktr` loss *and* the bend in the length–tension curve. A correction built on this
collapses the cross-prep spread of the ATP force effect from **CV 20 % → 3–6 %**.

> **[methods.md](methods.md)** is the paper-ready appendix: notation, equations,
> calibration, validation, limitations and a worked example. Read that to *apply*
> the correction; read this page for *what it found*.

| script | what it does | figure |
|---|---|---|
| `RunBracketExplainer.m` | what "the bracket" is, from the raw timecourses | `results/bracket_explained.png` |
| `RunRundownCorrection.m` | builds + calibrates the correction on effective activated time | `results/rundown_correction.png` |
| `RunRundownMechanisms.m` | data-side test: length–tension shape and sampling-matched `ktr` | `results/rundown_mechanism.png` |
| `RunMechanismSimulation.m` | imposes each candidate lesion **on the model** and compares | `results/mechanism_simulation.png` |

The whole analysis rests on one design feature: **03/27 recorded 8 mM twice**,
fresh (t = 0) and again at the end of the session (t = +694 s, after the 2 mM and
PNB+Mava runs). That bracket measures rundown with no modelling assumptions.
*(Its file is named `04_..._repeat` but was recorded last, at 15:29:40 — the file
numbering is not chronological.)*

**"Bracket" means bracketed in TIME, not in force.** Force never rises: the 8 mM
condition falls monotonically 68.43 → 56.45 kPa. The 2 mM run reads higher than
both only because low ATP increases force. See `results/bracket_explained.png`.

---

## 1. The clock is activated time, and the decay is linear

A skinned fibre degrades while it generates force, not while it sits in relaxing
solution — so the natural clock is accumulated **activated** time and the natural
rate is each recording's own within-run force slope (kPa/s, measured on the
isometric plateaus that bracket the slack train).

| | fresh (t=0) | repeat (t=+694 s) |
|---|---|---|
| force at SL 2.0 | 68.43 kPa | 56.45 kPa (−17.5 %) |
| within-run slope | **−0.4592 kPa/s** | **−0.4618 kPa/s** |

The slopes are identical in **absolute** terms even though the second run is
17.5 % weaker. Exponential (force-proportional) rundown would have predicted
−0.379. **Rundown is a constant kPa/s while activated.** The gap between two runs
is therefore naturally an *effective time*, `τ_eff = ΔF/|slope| = 26.1 s` — versus
694 s of wall clock. This is the `T0` quantity `DataCuration/FitRundownCorrection.m`
already computes.

## 2. Most of the within-run decline is reversible

Activation lasts a very consistent **15.3 ± 0.4 s** in every force-generating run
(PNB+Mava is blebbistatin-suppressed, slope ≈ −0.016 kPa/s, and is charged zero
damage). Integrating each run's own slope over its activation predicts **26.5 kPa**
of loss across the bracket; the bracket measures **11.98 kPa**.

> **φ = 0.45** of the within-run decline is permanent damage; **~55 % recovers**
> between activations (product accumulation washing out).

Damage separating two runs:  `loss = φ · Σ |slope_i| · T_act_i`  over the
intervening force-generating activations.

## 3. Two independent calibrations of φ agree

| source of φ | value |
|---|---|
| the 03/27 bracket (mechanical, no ATP assumption) | **0.452** (SL 2.0) / 0.342 (SL 2.2) |
| minimising cross-prep CV of the ATP force ratio | **0.625** |

Nothing links these two routes, yet they land in the same region — and across the
whole range the corrected ATP force effect is **≈ ×1.35**:

| φ | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| 0.00 (raw) | 1.222 | 1.194 | 1.685 | 1.367 | **20.2 %** |
| 0.45 (bracket) | 1.281 | 1.325 | 1.441 | 1.349 | **6.1 %** |
| 0.62 (optimum) | 1.305 | 1.384 | 1.365 | 1.351 | **3.0 %** |
| 1.20 | 1.393 | 1.621 | 1.162 | 1.392 | 16.5 % |

An alternative damage model — one universal kPa/s rate rather than
slope-proportional — gives 1.335 / 1.341 / 1.486 (CV 6.2 %). **The answer is
robust to which damage model is assumed.**

> **This supersedes the earlier wall-clock exponential treatment**, which
> under-corrected (it ignored that damage accrues during activation) and therefore
> wrongly labelled 04/10 an irreconcilable outlier.

## 4. Mechanically, rundown is series-elastic creep

Each recording contains a **6-point length–tension curve** (quasi-isometric force
at the end of each of the five slack holds, ML 0.94–1.02, plus the pre-slack
plateau at ML 1.10). That curve is the discriminator, because the candidate
lesions deform it differently.

### Data side (`RunRundownMechanisms.m`)

| model | parameter | RMSE | ΔAIC |
|---|---|---|---|
| pure vertical **scale** (fewer heads) | s = 0.832 | 2.11 kPa | +9.0 |
| pure horizontal **shift** (lost length) | d = 0.060 ML = 0.120 µm SL | 1.24 kPa | +2.6 |
| **both** | s = 0.944, d = 0.044 ML | **0.84 kPa** | **0** |

A pure scale is refuted by the ML 1.10 point alone: it needs **+4.4 kPa** there
while every other point sits ~1.2 kPa low. **A scale cannot bend, and the curve
bends.**

`ktr` (sampling-matched — the repeat is 10 ms while τ ≈ 20 ms, so the fresh trace
was resampled onto the repeat's *own* sample times and both fitted identically):
**×0.884 ± 0.021**, all five segments down, *t* = 5.5. Sampling alone accounts for
only ×0.956.

### Model side (`RunMechanismSimulation.m`)

Each lesion imposed on the model, then interpolated to the level that reproduces
the observed force loss exactly, so `ktr` and the L–T shape become free
discriminators:

| id | lesion | lever | ktr at matched force | L–T shape |
|---|---|---|---|---|
| M1 | fewer attached heads | `kstiff1,2` ×0.84 | 0.970 (obs 0.884) | **scale** ✗ |
| M2 | serial creep | `kSE` ×0.4 → force barely moves | huge ktr lever (×0.62) | — |
| M3 | lost SL at fixed ML | `SL0` −0.098 µm | 0.982 ✗ | **SHIFT** ✓ |
| M4 | uniform kinetic slowdown | `xrate` → force *unchanged* | pure ktr lever | — |
| M5 | reduced attachment ("less power") | `ka` ×0.73 | **1.045 — wrong sign** ✗ | scale ✗ |

**No single lesion works.** Force + the bend need M3; `ktr` needs M2 or M4. But
**M2 and M3 are the same physical lesion**: a series elastic element that creeps
longer *and* softer simultaneously lets the sarcomeres sit shorter at a given ML
(M3) and adds compliance (M2). Testing that pairing:

| combination | force | ktr | L–T shape |
|---|---|---|---|
| **C1 series-elastic damage** (`kSE`×0.65, SL −0.098 µm) | **0.827** | **0.877** | **SHIFT (bends)** |
| C2 fewer heads + slower motors (`kstiff`×0.84, `xrate`×0.85) | 0.841 | 0.914 | scale ✗ |
| **OBSERVED** | **0.829** | **0.884** | **SHIFT (bends)** |

Both combinations can reach the observed *(force, ktr)* point — that pair of
numbers alone does **not** identify the lesion. **The length–tension shape does**,
and it selects series-elastic damage.

### What this rules out

* **Fewer attached heads / myosin damage (M1)** — cannot bend the curve, and
  leaves `ktr` far too high.
* **A dead region in PARALLEL** (torn myofibrils still attached at both ends) —
  on the active side it is indistinguishable from M1, and it differs only by
  adding parallel passive force, which is ~5 kPa here (the PNB+Mava level) and
  cannot produce a 4.4 kPa *bend* in the right direction.
* **"Less power" per sarcomere (M5, reduced attachment)** — moves `ktr` the wrong
  way in this model.
* **Uniform kinetic slowdown (M4)** — leaves force untouched (duty ratio is
  invariant to a uniform rate scaling), so it cannot be the primary lesion,
  though it remains an admissible minor add-on.

## 5. Do not correct `ktr`

Tempting, but it must be judged on whether it improves cross-prep consistency:

| | 03/27 | 04/03 | 04/10 | Baker | CV |
|---|---|---|---|---|---|
| raw `ktr` ratio | 0.537 | 0.593 | 0.492 | 0.566 | **8 %** |
| rundown-corrected | 0.550 | 0.607 | 0.465 | — | **13–16 %** |

The correction makes it **worse**, under both the wall-clock and the
effective-time model. The bracket's `ktr` loss does not transfer to the short
gaps. Leave `ktr` raw — over a 132–138 s gap the bias is ~2 % — and never compare
`ktr` across long session gaps.

---

## How to use this

1. **Correct force with** `loss = φ · |slope_earlier| · T_act_earlier`, φ = 0.45,
   referenced to the **later** run's own force (the later run is the degraded one).
2. **Do not correct `ktr`.**
3. **Prefer a length offset over a force scale** when representing rundown in the
   model: the dominant term is a ~0.1 µm SL loss plus ~35 % extra series
   compliance, not a loss of heads. A consequence is that the nominal `SLslack`
   axis (2.04 … 1.88 µm) is systematically too long for later activations.
4. `DataCuration/FitRundownCorrection.m`'s `r(F,SL)` surface remains the right
   *force* correction, and §4 explains **why it has to be SL-dependent**: most of
   the loss is a length shift, and the surface absorbs it in force coordinates.
   Its hand-tuned `r0 = 1.214, k = −0.6` reproduces the bracket to within 2 %.

## Next

* **Record a hi-res `*_slack.txt` on the end-of-session repeat.** The `ktr` and
  length–tension parts are limited by the repeat being log-only (10 ms). Note
  `mergeLogsAndBursts` currently *excludes* filenames containing `repeat` from
  burst matching — a one-line change.
* **Measure stiffness directly.** `restretchSlopeStart / v_restretch` separates
  "fewer heads" (stiffness ∝ force) from "series creep" (stiffness falls faster
  than force) in one number, but needs hi-res on both runs.
* Repeat the bracket on a second prep to test whether φ ≈ 0.45 is universal.

Downstream: [ATPEffectReconciliation](../ATPEffectReconciliation/conclusions.md)
applies this correction; [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md)
is the cross-dataset survey it came from.
