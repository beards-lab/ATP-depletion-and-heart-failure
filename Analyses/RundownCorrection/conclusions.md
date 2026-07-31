# RundownCorrection — how the preparation degrades, and how to correct for it

**Headline.** Rundown is **linear in effective activated time**, not wall-clock time;
only **φ ≈ 0.45** of each recording's within-run force decline is permanent; and
mechanically it behaves as **creep in the series elastic** — the element gets
~35 % softer and ~0.1 µm longer — which alone reproduces the force loss, the
`ktr` loss *and* the bend in the length–tension curve. A correction built on this
collapses the cross-prep spread of the ATP force effect from **CV 20 % → 3–6 %**.

> **[strategy.md](strategy.md)** — the decision document: what to correct, which
> running order to request, where to put repeats, how to pool. **Start there** if you
> are designing or analysing a dataset. It also corrects two recommendations made
> earlier on this page.
>
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
both only because low ATP increases force.

![The bracket, from the raw timecourses](results/bracket_explained.png)

*(a) the four raw 03/27 timecourses — 2 mM (orange) highest, 8 mM fresh (light blue),
8 mM repeat (dark blue) lowest. (c) the damage staircase: force drops only during the
shaded activation windows and lands on the measured repeat. The orange bar is the ATP
effect — the measured 2 mM force against where 8 mM would have been at that instant.*

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

![Calibration and validation](results/rundown_correction.png)

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
numbers alone does **not** identify the lesion.

> **⚠ Revised by §4b.** The length–tension verdict above rests on the single
> **ML 1.10** point. Testing against the *five-segment* profile instead reverses the
> force half of the conclusion: the force loss looks like **fewer heads**, not lost
> length. The `ktr` half (added series compliance) survives unchanged.

![Model perturbation study](results/mechanism_simulation.png)

*Each lesion is tuned to reproduce the observed force loss exactly, so `ktr` (top
right) and the length–tension shape (bottom middle) become free discriminators. Only
the lesions whose "shift" bar beats their "scale" bar bend the curve as the data does.*

![Data-side evidence: length-tension and ktr](results/rundown_mechanism.png)

## 4b. Revision — the five-segment profile overturns the force half

`RunModelRundownFit.m` → `results/model_rundown_fit.png`. The bracket gives not one
number but a **profile**: `ktr` and amplitude at each of the five slack segments, for
fresh and repeat. A lesion tuned to the mean is fitting one number; a lesion that
reproduces the profile is being tested against five.

![Rundown as a decaying model parameter](results/model_rundown_fit.png)

**Observed** (repeat / fresh, sampling-matched):

| | s1 (2.04 µm) | s2 | s3 | s4 | s5 (1.88 µm) | mean |
|---|---|---|---|---|---|---|
| `ktr` ratio | 0.947 | 0.906 | 0.847 | 0.829 | 0.893 | 0.884 |
| amplitude ratio | 0.823 | 0.805 | 0.808 | 0.802 | 0.792 | **0.806** |

**The amplitude profile is essentially FLAT** (gradient −0.031 across the ladder).
That is the discriminator:

| candidate | RMSE `ktr` | RMSE amplitude | **total** |
|---|---|---|---|
| SL0 only (−0.098 µm) | 0.116 | 0.058 | 0.130 |
| kSE only (×0.65) | 0.069 | 0.187 | 0.199 |
| SL0 + kSE | 0.066 | 0.052 | 0.084 |
| kstiff only (×0.84) | 0.108 | **0.023** | 0.111 |
| **kstiff + kSE (×0.84, ×0.65)** | **0.058** | **0.017** | **0.060** |

A length shift of the required size predicts a **steep** amplitude gradient
(0.910 → 0.726, −0.184) because the deep segments lose more length–tension. The data
show −0.031. A force scale predicts −0.001. **On five points the force loss is a
scale, not a shift.**

**So why did §4 say shift?** Because the 6-point length–tension curve includes the
**ML 1.10** pre-slack plateau (ratio 0.887), which sits *outside* the slack ladder and
is the only point that breaks the flat pattern — the other five run 0.79–0.83. One
point drove that verdict. Partial explanations for it: passive force is largest there
(~5 kPa of 80, and passive does not run down), which moves the expected ratio from
0.81 to ~0.82 — some of the way to 0.887, not all. **Treat the length-shift evidence
as provisional and resting on a single measurement.**

**Revised mechanism.** Rundown = **loss of force-generating material in parallel
(≈16 %) plus added series compliance (≈35 %)**. That is your "tearing the cell"
hypothesis: myofibrils rupturing near the attachments lose parallel force-generating
material *and* leave a compliant damaged end region in series. One lesion, two
consequences, and it fits both profiles (total RMSE 0.060 versus 0.084 for the
series-creep picture).

**What is unchanged:** the `ktr` loss still requires **added series compliance**. No
alternative reaches it — `kstiff`, `SL0` and `ka` leave `ktr` at 0.98–1.05, and while
`xrate` moves `ktr` it leaves force untouched.

### Using rundown as a model nuisance parameter

This is the practical payoff. Rather than correcting data, give each run a **single
rundown coordinate** mapped onto (`kstiff`, `kSE`) jointly, fit it alongside, and the
ATP parameters are estimated on a common footing automatically. One nuisance
dimension per run, not two free knobs — the two model parameters move together
because they are two consequences of one lesion.

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

## 5. Decay geometry: when rundown acts, and why it biases low ATP

`RunDecayGeometry.m` → `results/decay_geometry.png`. Three isometric ML = 1.0
windows are used: **W1** 68.20–69.20 s (before any perturbation — the only clean
one), **W2** 71.20–72.35 s (after the ktr burst), **W3** 77.70–78.30 s (after the
last slack). W2/W3 are post-perturbation recoveries; they sit at identical protocol
times in every run so the bias cancels in between-run ratios, but they are not
absolute isometric values.

![Decay geometry](results/decay_geometry.png)

### 5a. The activation rise overlaps the decay — and it is ATP-dependent

Force takes **3.3–5.8 s to peak** after activation begins (calcium and substrate
diffusing into a skinned preparation), and rundown is already running during that
rise. So the observed peak is *not* the undamaged force. More importantly the two
conditions do not develop force at the same rate:

| | time to peak after onset | at W1 (68.7 s) |
|---|---|---|
| 8 mM (03/27, 03/27 rpt, 04/03) | 5.3 – 5.8 s | **still rising** |
| 2 mM (all three) | 3.6 – 4.6 s | **already falling** |

*(04/10's 8 mM peaks at 3.3 s — the same run that is anomalous everywhere else.)*

Consequence: at the W2 window the 2 mM run has already been decaying ~3.0 s while
the 8 mM run has been decaying ~0.7 s, **and** its slope is ~2× steeper. On 03/27
that is 3.65 kPa lost from the 2 mM peak against 0.34 kPa from the 8 mM peak — a
systematic understatement of the ATP effect by ~4 % that no φ-based correction
touches, because it happens *inside* each run.

### 5b. Does "no decay between activations" hold? — very nearly

Extrapolating each run's decay line to its own peak and its own end gives a
continuity chain on 03/27. With `a` = ATP effect and `r` = between-activation
recovery (`r > 1` means force partially recovers while relaxed):

```
  8 mM #1 : peak 68.64 -> end 64.27      (loses 4.36 kPa while activated)
  2 mM #2 : peak 87.08 -> end 71.51      (loses 15.57 kPa)
  8 mM #4 : peak 57.05 -> end 52.46      (loses 4.59 kPa)

  8mM#1_end -> 2mM#2_peak :  x1.355  =  r * a
  2mM#2_end -> 8mM#4_peak :  x0.798  =  r^n / a
```

| rest intervals assumed | ATP effect `a` | recovery `r` |
|---|---|---|
| one each | **1.303** | +4.0 % per rest |
| two after the 2 mM run | **1.320** | +2.6 % per rest |

> **The assumption essentially holds.** The residual is a small **recovery** of
> 2.6–4.0 % per rest interval — not further decay. And this route yields an ATP
> effect of **×1.30–1.32 from 03/27 alone, with no φ and no damage model assumed** —
> an independent confirmation of the ×1.36 consensus reached the other way.

This is the strongest validation in the analysis: it uses only the within-run decay
lines and the assumption of continuity, and it lands on the same number as the
φ-calibrated route.

### 5c. Confirmed: 2 mM runs down ~2× faster

| | 8 mM | 2 mM | ratio |
|---|---|---|---|
| absolute (kPa/s) | 0.560 ± 0.142 | 1.092 ± 0.113 | **×1.95** |
| relative (%/s) | 1.06 ± 0.37 | 1.50 ± 0.12 | ×1.42 |
| paired within prep | | | ×2.66, ×1.37, ×1.77 (**all > 1**) |

Two things worth noting. The **2 mM relative slope is remarkably tight** across
preps (1.41–1.63 %/s, CV 8 %) while the 8 mM is scattered (0.67–1.42 %/s, CV 35 %) —
low ATP appears to impose a stereotyped decay rate. And **03/27's 8 mM has the
shallowest relative slope of any run (0.67 %/s)**: the preparation the entire
correction is calibrated on is the most robust one in the set, which is a reason to
want a second bracket on a more typical fibre.

### 5d. The improved correction

Referencing each run to **its own peak** rather than to a fixed protocol window:

| referencing | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| fixed window W2 | 1.225 | 1.196 | 1.678 | 1.366 | **20 %** |
| own-peak | 1.269 | 1.233 | 1.649 | 1.384 | **17 %** |
| own-peak + damage over its own decay phase | 1.306 | 1.313 | 1.486 | **1.369** | **7 %** |

The own-peak step is **free** — no φ, no damage model — and it removes the intra-run
bias of 5a. Adding the damage term over each run's own decay phase (rather than a
common activation duration) reaches CV 7 % at ×1.37, converging with the ×1.36 from
the independent φ route.

## 5e. Can the refinements rescue a `ktr` correction? No — and now we know why

`RunJointCorrection.m` → `results/joint_correction.png`. Re-derives φ on the
own-peak / `T_dec` footing (**φ = 0.581**, versus 0.452 on the `T_act` footing —
the *product* φ·|s|·T is what the bracket constrains, so this is the same model
re-expressed) and then asks the only question that matters for each observable:
**does the correction reduce between-preparation scatter?**

![Joint correction](results/joint_correction.png)

**Force — yes.** Own-peak referencing plus damage over each run's own decay phase:

| | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| damage over the gap | 2.9 % | 6.4 % | 14.1 % | | |
| raw (own-peak) | 1.269 | 1.233 | 1.649 | 1.384 | **17 %** |
| corrected | 1.306 | 1.311 | 1.445 | **1.354** | **6 %** |

Two different footings (`T_act` with φ=0.452, `T_dec` with φ=0.581) give ×1.37
and ×1.35 — the framework is self-consistent.

**`ktr` — no.** The coupling to the force correction was left free and scanned:

```
ktr_frac(gap)  =  gamma * force_frac(gap)
```

| γ | meaning | ktr ratios | CV |
|---|---|---|---|
| 0 | no correction | 0.537 / 0.593 / 0.492 | **9.4 %** |
| +0.68 | implied by the bracket (force ×0.829, `ktr` ×0.884) | 0.548 / 0.618 / 0.449 | 15.8 % |
| +0.71 | predicted by the series-creep lesion | 0.548 / 0.620 / 0.447 | 16.1 % |
| −0.50 | scan edge | 0.529 / 0.574 / 0.529 | 4.7 % |

`CV(γ)` is **monotone across the whole range** — there is no interior optimum. The
apparent improvement at negative γ is not a mechanism (a negative coupling would
mean rundown makes `ktr` *faster*); it works only by lifting 04/10, the lowest
ratio, toward the other two — fitting one preparation's idiosyncrasy with one free
parameter on three points.

And there is no dose-response to follow: `ktr` ratio versus force-damage dose gives
**r = −0.63, p = 0.57 (n = 3)**, and the pattern is not even monotone (0.537 at
2.9 %, 0.593 at 6.4 %, 0.492 at 14.1 %).

**Why the two observables cannot share a correction.** The lesion is series-elastic
creep (§4). It reaches **force** through the length–tension relation — a strong
channel — and **`ktr`** through added compliance — a weak one. Over the bracket that
is force ×0.829 against `ktr` ×0.884, and the series-creep lesion predicts ×0.827 /
×0.877. A scaled correction assumes exactly the proportionality that the mechanism
denies.

**Could a within-run `ktr` decay be measured instead?** Attempted, and it is below
the noise floor. Comparing `ktr` from the ktr burst (t ≈ 71.6 s) with slack segment 2
(t ≈ 75.5 s, same SL 2.0, same activation, ~4 s apart) gives ratios of
1.09 / 0.86 / 1.45 / 1.22 / 0.94 / 0.73 across the six runs — scatter far larger than
any plausible 4-second decay. The two are mechanically different perturbations and
the burst fit window is only ~2.8 time constants. **This measurement needs a hi-res
repeat run to be made properly.**

## 5f. Why the slack `ktr` falls with rundown — and why the two ktr protocols disagree about ATP

`RunKtrProtocolSensitivity.m` → `results/ktr_protocol_sensitivity.png`.

![ktr protocol sensitivity](results/ktr_protocol_sensitivity.png)

Two observations needed one explanation:

1. **The rundown `ktr` drop is real and is not biological variability.** Between
   03/27's 8 mM run and its 8 mM *repeat* — same fibre, same condition — the slack
   `ktr` falls ×0.884 in **every one of the five segments** (t = 5.5).
2. **The ATP effect on `ktr` is protocol-dependent.** Fits are well conditioned
   (R² 0.974–0.995, tight CIs):

| prep | ktr protocol | slack (seg 2) |
|---|---|---|
| 03/27 | 0.746 | 0.591 |
| 04/03 | 0.732 | 0.616 |
| 04/10 | 0.691 | 0.537 |
| **mean** | **0.723** (CV 4 %) | **0.581** (CV 7 %) |

Both are internally tight, yet **24 % apart**. And they diverge only at low ATP —
slack/protocol is 1.09 / 1.45 / 0.94 at 8 mM but 0.86 / 1.22 / 0.73 at 2 mM.

### The explanation: the slack `ktr` is compliance-contaminated

After an unloaded shortening the contractile element must shorten *again* to
re-stretch the series element before force appears. The ktr protocol restretches to
1.05 ML and returns to 1.0, starting much closer to its force-bearing configuration.
Model test:

| kSE | slack `ktr` | protocol `ktr` |
|---|---|---|
| ×1.00 | 44.88 (—) | 69.01 (—) |
| ×0.80 | 41.59 (×0.927) | 67.00 (×0.971) |
| **×0.65** | 37.26 (**×0.830**) | 66.63 (**×0.965**) |
| ×0.50 | 32.46 (×0.723) | 65.27 (×0.946) |

**The slack measurement is ~5× more compliance-sensitive.** That single fact explains
both observations:

* **The rundown drop** — rundown adds ~35 % series compliance (§4b), which costs the
  slack `ktr` 17 % in the model against the −12 % observed. Same sign, right order.
* **The protocol disagreement** — a measurement 5× more compliance-sensitive reports a
  different rate whenever compliance differs between the conditions compared.

> **Consequence.** The slack `ktr` is the right probe for **rundown** (a compliance
> lesion) and the **wrong** probe for cross-bridge kinetics. For the ATP effect prefer
> the ktr-protocol value **×0.72 (CV 4 %)** over the slack value ×0.58 (CV 7 %).

**Caveat.** The model's absolute protocol-`ktr` (69/s) is much faster than its slack
`ktr` (45/s), whereas the data have them within 10 % at 8 mM — the defect documented
in [RestretchVsKtrRecovery](../RestretchVsKtrRecovery/conclusions.md). The
sensitivity *ratio* is a within-model comparison and survives that, but its magnitude
may be overstated.

**What could not be done.** The 03/27 repeat has no hi-res burst file, only the 10 ms
log. Its ktr-protocol redevelopment lasts ~60 ms — 6 samples — so the rundown drop is
measurable only in the slack, whose 298 ms hold gives 30. Recording a hi-res
`*_stiff_ktr.txt` on an end-of-session repeat would test the compliance explanation
directly: it predicts the protocol `ktr` should fall only ~3 % where the slack `ktr`
falls 12 %.

## 6. Do not correct `ktr`

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
