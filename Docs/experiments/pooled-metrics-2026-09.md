# The pooled 8 mM fit target — methods, current numbers, and whether to identify now

**2026-09-10 · F. Jezek · three preparations (03/27 M, 04/03 F, 04/10 M), 8 mM and 2 mM MgATP**

> ## The ask, in five lines
>
> We have replaced a fit target built from **one** preparation with one built from
> **three**, carrying an error bar on every number. The target moved by more than the
> model's own residual, so this is not cosmetic: the old target was quietly rewarding a
> known model artefact.
>
> **The recommendation on identification is: not yet.** Do the structural refit now;
> start hard identification when *n* ≥ 6. The reason is arithmetic, not caution — at
> *n* = 3 each objective weight is known only to a factor of ~18, so a parameter
> confidence interval computed today would mostly be reporting our ignorance of the
> weights, not the data's information about the parameter.

---

## 1. Methods

### 1.1 The nuisance being removed

Tension is force ÷ *estimated* cross-sectional area. Across the three days the pre-slack
isometric force at SL 2.2 µm is **79.8 / 70.0 / 50.9 kPa** — a CV of **22 %**. Almost
none of that is biology of interest; it is the area estimate.

The generative model this implies is

```
    y[p,c,k]  =  alpha_p^(d_k) · mu[c,k] · exp(eps)
```

`p` preparation, `c` ATP condition, `k` feature, `d_k` = 1 for force-dimension features
and 0 for rates, times and ratios. **`alpha_p` is a property of the preparation** — the
same fibre, the same area estimate, both ATP runs.

Two facts make this testable rather than assumed:

- **One scale suffices.** Centre the log force features and take the SVD: **PC1 holds
  88.3 %** of the between-preparation variance and its loading is flat (|cos| with the
  all-ones direction = **0.955**). A pure common scale is exactly a rank-1 matrix with a
  flat loading. Both conditions hold.
- **Scale-free features are the control group.** `ktr`, `t0`, `rsK`, `rsR2` carry
  `d_k` = 0, so no normalisation can touch them; their 11–38 % spread is the floor.
  After normalisation the force features sit at **4–11 %**, i.e. *better than the best
  rate*. Absolute tension goes from the worst-determined quantity in the experiment to
  among the best.

### 1.2 Four steps, in this order

1. **Scale.** Divide every force-dimension feature of preparation *p* by that
   preparation's own `Fss22`, taken from its **8 mM** run. **One divisor per
   preparation, shared by both ATP conditions** — so the 2 mM/8 mM ratio survives
   untouched and stays a *measurement* rather than a definition.
2. **Rundown.** Multiply the run recorded **second** in each session by
   `1 + φ·|slope|·T / F`, with **φ = 0.452**. Without this the ATP ratio depends on
   running order, which differs across the three days by design.
3. **Pool.** Mean of the per-preparation **features** — never features of an averaged
   trace. Averaging traces biases the re-stretch rate by 24 %, more than a preparation
   width, and that bias traces to a single segment whose feature window (844 ms) is not
   covered by the splice grid (570 ms).
4. **Spread.** *s* = between-preparation relative SD, **variance-moderated** toward the
   median of its own class (the standard low-replicate remedy; Smyth 2004). Prior
   **d₀ = 2.58** against *d* = 2 from the data, giving a moderated **4.58 dof**.

**Step 5, added 2026-09.** *s* is measured on **exactly the slack cycles the objective
scores.** For the re-stretch rate that is cycles 1–4 (§4.4). Measuring the spread over a
different set than the fit scores makes the χ² not a χ².

### 1.3 Two choices worth defending explicitly

**Which anchor barely matters; the scope decides everything.** A paired bootstrap over
the feature set (4000 draws) cannot separate the top three anchors — P(best) 64 / 18 /
15 % for `Fss22` / `Fsteady` / `Fiso20a`, all within 1.8 pp of spread. *Normalising
versus not* is decisive (29 % → 15 %); *which divisor* is not a decision this dataset can
support. `Fss22` is chosen on grounds the bootstrap cannot see: it is measured at the
same length and activation as the reference state, over ~40 ms of clean plateau, and it
is **not itself a scored feature**.

**Per preparation, not per run.** A per-run divisor buys 1.3 pp of spread and costs the
entire ATP force effect — it divides each condition by its own level, so whatever it then
reports about the level ratio is the gauge, not the data. This is the one place where
"lowest spread wins" gives the wrong answer.

*Rundown is worth its own line: the lift explains **50 %** of the apparent
between-preparation force variance (anchor CV 22.0 % → 15.5 %), and it takes the ATP
ratio's scatter from 0.274 to 0.119.*

### 1.4 Scope — what is pooled and what is not

| | status |
|---|---|
| slack block, 20 features × 5 cycles, 8 mM **and** 2 mM | **pooled, with error bars** |
| the absolute level (`F_ref`) | **pooled**, kept as its own term |
| force–velocity (borrowed curve + in-preparation FV2 ramps) | not pooled |
| ktr-burst protocol, staircase ramp-up | not pooled — held as **predictions** |

---

## 2. Results — the current pooled target

**Level.** `F_ref` = **69.35 kPa**, CV **15.5 %** (per preparation 79.8 / 70.0 / 58.3
after the rundown lift). This is the one absolute number the normalisation removes, and
it is kept rather than discarded — without it the objective is scale-blind by
construction.

**ATP effect on the level:** force ×**1.320 ± 8.6 %** (per preparation 1.219 / 1.298 /
1.443), i.e. ~4.8 SEM from unity.

### 2.1 The twenty features

Sorted by ATP discriminability *z* = |ratio−1| / √(s² + s_ratio²).

| feature | cls | pooled 8 mM | *s* used | *s* raw | 2/8 ratio | *z* | w = 1/s² |
|---|---|---|---|---|---|---|---|
| **rsK80** | 0 | 48.87 s⁻¹ | 15.8 % | 7.7 % | **0.381** | **3.63** | 40.1 |
| **vall2_t** | 0 | 6.38 ms | 26.3 % | 32.7 % | **1.924** | **2.58** | 14.5 |
| **rsK** | 0 | 53.47 s⁻¹ | 20.5 % | 21.4 % | **0.345** | **2.46** | 23.7 |
| **ktr** | 0 | 47.85 s⁻¹ | 16.7 % | 11.4 % | **0.532** | **2.37** | 35.8 |
| **rsT20** | 0 | 8.54 ms | 21.4 % | 23.2 % | **1.752** | **2.37** | 21.8 |
| steady | F | 0.945 | 15.2 % | 4.0 % | 1.320 | 1.84 | 43.5 |
| vall_y | F | 0.820 | 16.4 % | 10.2 % | 1.384 | 1.75 | 37.2 |
| A | F | 0.719 | 16.5 % | 10.6 % | 1.337 | 1.69 | 36.8 |
| peak2 | F | 0.963 | 15.7 % | 7.5 % | 1.365 | 1.62 | 40.4 |
| Am | F | 0.713 | 16.4 % | 10.2 % | 1.331 | 1.59 | 37.3 |
| peak1_y | F | 0.970 | 16.5 % | 10.5 % | 1.390 | 1.50 | 36.9 |
| t0 | 0 | 10.96 ms | 29.0 % | 37.7 % | 1.434 | 1.29 | 11.9 |
| vall2_dy | F | −0.165 | 21.3 % | 22.9 % | 1.247 | 1.03 | 22.1 |
| rsR2 | 0 | 0.734 | 16.9 % | 12.0 % | 1.193 | 1.02 | 35.0 |
| rsA | F | 0.166 | 21.4 % | 23.2 % | 1.237 | 0.96 | 21.8 |
| restretchSlopeStart | F | 19.8 s⁻¹ | 20.0 % | 20.1 % | 1.237 | 0.89 | 25.1 |
| peak1_dSL | 0 | 17.8 mML | 52.4 % | 76.0 % | 1.558 | 0.78 | 3.6 |
| rsT0 | 0 | 5.44 ms | 21.7 % | 23.9 % | 0.664 | 0.42 | 21.2 |
| peak1_t | 0 | 1.83 ms | 44.6 % | 63.6 % | 1.301 | 0.42 | 5.0 |
| vall_t | 0 | 6.65 ms | 57.1 % | 83.5 % | 1.094 | 0.11 | 3.1 |

Force-dimension ('F') values are in units of `F_ref`; multiply by 69.35 kPa for absolute.

### 2.2 The two readings that matter for the science

**The ATP effect on force is ONE number, not nine.** Divide each force feature's 2/8
ratio by the level and every force *shape* ratio is **1.00 ± 0.05**. Nine force ratio
terms would be nine copies of one measurement. What ATP genuinely and separately moves
is **kinetic**: the re-stretch rate (×0.381), `ktr` (×0.532), and two timing features
(`rsT20` ×1.752, `vall2_t` ×1.924).

**Variance moderation is doing real work, and it must.** Compare *s* raw with *s* used.
`steady`'s raw 4.0 % is a *within-run repeatability* masquerading as a between-preparation
spread — it is the post-re-stretch isometric divided by the pre-slack isometric at the
same length, seconds apart, so it is close to a self-normalisation. Unmoderated, it would
take the largest weight in the set. Moderation lifts it to 15.2 % and compresses the
whole 1/s² dynamic range from **446× to 13×**.

---

## 3. What changed against the previous target

The previous target was **03/27 alone**, in that fibre's own kPa.

### 3.1 The target values

| | change, old → new |
|---|---|
| whole force block | **−13.1 % uniform** (`F_ref`/`Fss22`(03/27) = 69.35/79.76 = 0.869) |
| `peak1_y` | **−24.7 %** — i.e. ~−11.6 % of *shape* on top of the level |
| `A`, `Am` | −19.9 %, −20.0 % |
| `peak2`, `steady` | −14.9 %, −15.3 % (level only) |
| `peak1_t`, `vall_t`, `peak1_dSL` | −36.8 %, −36.1 %, −30.3 % (timing) |
| `rsK` | **+22.3 %** — 03/27 is the lowest-rate preparation of the three |

### 3.2 Is fitting one preparation the same as fitting the pool?

**Almost.** χ² of each preparation against the pooled mean is **0.40 / 0.77 / 0.77**
against an expectation of 0.67 for *n* = 3. 03/27 is the *most* typical of the three —
pooling adds error bars, it does not relocate the target wholesale.

**With one exception, and it is the one that mattered.** The pool's first re-stretch peak
sits **3 % below** its own isometric; 03/27's sits **12 % above** it. Fitting 03/27 alone
was *rewarding* the model's known lengthening overshoot. This is not 03/27 being a bad
preparation — it is one feature where a single fibre's idiosyncrasy happened to point the
same way as a model artefact.

### 3.3 Two metric corrections folded in

**`rsK` → `rsK80`.** The old re-stretch rate took its amplitude from the tail of its own
fit window, so a known ×0.6 time dilation returned 24 % wrong and a 60 ms window shift
moved it **+77 %**. `rsK80` = ln4/(t₈₀−t₂₀) takes the amplitude from `steady` and reads
the rate off the 20–80 % band. Spread 20.5 % → **15.8 %**, *z* 2.46 → **3.63**.

**`rsK80` measured on cycles 1–4 only.** Injecting each run's own measured transducer
noise into its own trace 200× and re-extracting gives each feature's *measurement* noise
with biology held exactly fixed. `rsK80` came back at **32.7 %** — *above* the 16.9 %
between-preparation spread its weight was built on, which is impossible: a feature cannot
be better determined between preparations than one measurement allows. The whole apparent
spread was the transducer, and the objective was weighting it as a joint-largest term.

It is one cycle. Cycle 5 re-stretches to a **different final length** (1.000 ML, not
1.100) by the **smallest** amplitude (0.061 ML against 0.080–0.140) — it is not on the
ladder. Its `rsK80` noise is 29 / 44 / **245 %** against 6–15 % for cycles 1–4, in every
run and both ATP conditions. Dropping it: noise **32.7 % → 5.2 %**, and the weight goes
**up** (1/s² 35.0 → 40.1). Cycle 5's *values* are retained; only its contribution to the
spread is removed.

### 3.4 A separate question answered on the way: is the recovery ripple real?

Every force recovery carries a small oscillation on the exponential. Ninety recovery
windows across the three preparations, at 8 mM, 2 mM and with cross-bridges chemically
inhibited (PNB + mavacamten):

- **Nine tenths of it is plain white transducer noise**, and there is no discrete
  resonance line anywhere in the residual spectra.
- With **cross-bridges switched off**, σ = 1.19 kPa against **1.24 kPa** active.
- The reproducible part's period moves ×0.70 with ATP while every measured rate moves
  ×1.9–2.6 — it does not track the kinetics.
- Its amplitude does not scale with the rise it sits on (r = −0.38).

**It is the rig.** The recoveries are single exponentials, as they appear. Consequence
for the objective: the data's own exponential reaches R² = 0.795, but **a model sitting
exactly on the truth would reach only 0.834**, because the comparison is against a noisy
trace. That gap must never be scored as a defect — `rsR2` stays a *boundary* term and
must not become a fit-to-data term.

### 3.5 What it did to the model

| | old objective (03/27) | pooled objective |
|---|---|---|
| `estiff125_opt` — best fit to 03/27 | **0.9547** | 4.4309 |
| `pooled8_opt` — refit to the pool (4 h) | 1.0876 | **2.6999** |
| | **+13.9 %** | **−39.1 %** |

**This is a relocation to a better-justified target, not a free win** — it gives up 14 %
on one fibre's objective to gain 39 % on the pool. The asymmetry is exactly what §3.2
predicts.

Against the pooled error bars, the refit **halved the resolved misses**:

| | previous best (`spreadw_opt`) | `pooled8_opt` |
|---|---|---|
| inside 1 SD | 13 / 20 | 13 / 20 |
| **beyond 2 SD** | **4** (all timing) | **1** |
| `rsT0` | −2.6 | −2.88 |
| `peak1_t` | +2.6 | **+1.96** |
| `peak1_dSL` | +2.1 | **+1.64** |
| `vall_t` | +2.0 | **+1.52** |

**And the gain generalised to a protocol the fit never saw.** On the ktr burst — not in
the objective — the re-stretch overshoot went 1.079 → **1.061** against a data value of
1.032. A gain in an unscored protocol is the strongest available evidence that the pooled
target found something real rather than fitting noise.

**Honestly, two things got worse**, both unscored: ktr redevelopment is now +41 % too fast
(was +30 %), and the staircase force–length end point +5.9 % (was +0.9 %). The
redevelopment is the same defect the scored slack `ktr` reports, read by a second
protocol — not an independent problem.

---

## 4. What the coming datasets will buy

The planned round is **6 preparations**, `8-4-2-8` alternating with `2-8-2`, a rundown
bracket in every session, sex balanced at every rung.

**Three different quantities improve at three different rates, and conflating them is how
sample-size arguments go wrong:**

| *n* | raw dof | moderated dof | 95 % CI on *s* | …moderated | ⇒ CI on **w = 1/s²** | SEM of a pooled mean |
|---|---|---|---|---|---|---|
| **3 (now)** | 2 | 4.58 | **12.1×** | 4.2× | **~18×** | 9.2 % |
| **6** | 5 | 7.58 | **3.9×** | 2.9× | **~8×** | 6.5 % |
| 9 | 8 | 10.58 | 2.8× | 2.4× | ~6× | 5.3 % |
| 12 | 11 | 13.58 | 2.4× | 2.2× | ~5× | 4.6 % |

- **Improves as 1/√n — the target VALUES.** Every pooled mean's SEM falls 9.2 % → 6.5 %.
  On the primary readout, the 2/8 ratio, SEM goes **4.6 % → 3.3 % → 2.7 %**.
- **Improves with dof — the WEIGHTS, and this is the big one.** At *n* = 3 the raw 95 %
  interval on a spread spans **12×**: a feature measured at *s* = 10 % is only known to
  lie between 5 % and 63 %. Since the weight is 1/s², today's weights are individually
  uncertain by about a factor of **18**. Moderation is what makes them usable at all, and
  going to *n* = 6 collapses the raw interval to **3.9×**. **The step from 3 to 6 is by
  far the most valuable; beyond 6 the returns flatten.**
- **Does NOT improve with *n*.** The between-preparation SD itself — that is the biology,
  and more preparations measure it better without making it smaller. Nor the rig noise
  floor (σ = 0.96–1.58 kPa per run, a transducer property). Nor cycle 5's unsuitability,
  which is protocol geometry.

**Which features gain, and which will not.** From the noise-injection experiment, the
share of each weight that is rig rather than biology:

| gains fully from more preparations (noise ≪ spread) | gains little (noise ≈ spread) |
|---|---|
| every force feature (noise 0.0–3.9 % vs spread 15–17 %), `ktr` 0.5 %, `steady` 0.4 %, `A` 0.0 %, `t0` 1.1 % | `rsT0` (noise 23.1 % vs spread 21.7 % — **never score it**), `rsT20` (43 % of the weight is rig), `peak1_t` / `vall_t` / `peak1_dSL` (18–22 %) |

**Newly possible at *n* ≥ 6, and impossible now: leave-one-preparation-out validation.**
At *n* = 3, holding one out leaves a two-preparation spread with **1 dof** — the held-out
test is uninterpretable. This is the single largest methodological gain of the round, and
it is what turns "the model fits the pool" into "the model predicts a preparation it never
saw."

**What no affordable *n* will buy:** the *shape* of the ATP dose-response.
Discriminating hyperbolic saturation from log-linearity needs ~42 preparations per dose at
the realistic CV. If the purpose of the 4 mM arm is "what is the functional form", this
preparation's scatter forecloses it. If the purpose is "give the model a third point it
was not fitted to", 3 sessions is right and it is a genuinely strong test — 89 % power
that 4 mM differs from 8 mM, 77 % to catch a 15 % model error.

---

## 5. Hard identification: now, or wait?

**Recommendation: do the structural refit now; begin hard identification when *n* ≥ 6
*and* the overshoot residual is at or below ~1 preparation-width.**

This is not a counsel of perfection. There are four specific reasons, in descending
order of force.

**1. The weights are not yet known well enough for a parameter CI to mean anything.**
Hard identification propagates target uncertainty into parameter confidence intervals. At
*n* = 3 each target value carries ~9 % SEM *and* each weight is uncertain by a factor
~18 (§4). An interval computed against that would largely be reporting our ignorance of
the objective, not the data's information about the parameter — and it would look
authoritative while doing so.

**2. There is a known structural defect, and identification against an unreachable
target is biased.** The entire remaining residual is one coupled object — the re-stretch
transient — at **4.2 / 3.8 / 2.8** preparation-widths on `peak1_y` / `peak1_dSL` /
`vall_t`, against nine of seventeen terms already at or below one width. It is not a
parameter problem: it survives into the unscored ktr protocol and the staircase. Fitting
confidence intervals around a point the model cannot structurally reach yields precise
answers to the wrong question.

**3. We cannot currently say whether the model is identifiable at all, and the evidence
we have is weaker than it looks.** The near-optimal cloud from the 4 h refit contains
**2** solutions, produced by a monotonically descending 9-round chain. Those are
correlated iterates, not a sample of the near-optimal manifold, so the "unique, 1.4 %
spread" verdict it reports is not a result. The remedy is already built: a **30 h
multi-start campaign in four profiles** (control / wider draw / hard kick / global
surrogate) that produces *independent* near-optimal points — which is precisely the input
an identifiability analysis requires.

**4. Cross-validation is impossible at *n* = 3** (§4). Any identification done now could
not be tested out-of-sample.

### What to do in the meantime — all of it prerequisite, none of it wasted

| now | why it does not depend on more data |
|---|---|
| **Run the 30 h multi-start campaign** | buys the descent that was still ongoing when the 4 h budget expired (9 improvements in 9 rounds, zero stalls) **and** the independent solutions identifiability needs |
| **Practical identifiability *screening*** — sensitivity + spread across the four profiles | tells us which parameters this protocol can constrain *even in principle*. Worth knowing **before** the data arrives, because it can still change the protocol |
| **Close two objective inconsistencies** | `coolDownLS` is a least-squares residual against 03/27's *raw trace in absolute kPa*, so it pulls toward 79.76 kPa while the pooled features pull toward 69.35. And the force–length end point is scored only through slack cycle 5, while the staircase measures it directly and **disagrees by 10 %** (slack 1.302 vs staircase 1.187) |
| **Decide whether to score `rsT20`** | the recovery is now known to be a pure exponential, i.e. three numbers — amplitude, rate, delay. We score the amplitude and the rate; the **delay is unconstrained**. `rsT20` is the audited delay, ATP effect resolved (*z* = 2.37), noise 14 % inside its 21 % spread |

---

## 6. Limitations, stated plainly

- ***n* = 3.** Every spread here has a ~12× confidence interval. The **rankings** among
  normalisation schemes are not resolvable; the **presence or absence** of normalisation
  is. Everything downstream is planned at CV 8 %, bracketed by 4 % and 11 %.
- **φ has one calibration.** The rundown lift uses φ = 0.452 from a single bracket on
  03/27. The correction's *direction* is not in doubt and its size is robust over
  φ = 0.45–0.62, but a second bracket is the single most valuable next measurement — and
  the `2-8-2` bracket additionally tests whether φ is ATP-dependent, which the 2 mM run's
  ×1.95 faster decay gives direct reason to suspect.
- **`steady` is nearly a self-normalisation** against `Fss22` (both isometric at ML 1.10,
  seconds apart). Moderation happens to correct this; an unmoderated objective would hand
  it the largest weight in the set.
- **The 2 mM side is not yet in the objective.** The ATP block above is measured and
  reported but the fit runs at 8 mM only; the low-ATP gate lives in a separate worktree.
- **FV, ktr and the staircase are not pooled.** They are held as predictions, which is
  their most valuable use, but it means their error bars are not comparable with the
  slack block's.

---

### Provenance

Target frozen in `Model/pooledTarget8mM.m` (regenerate with
`Analyses/DataPooling/WritePooledTarget.m` when a preparation lands). Analyses:
`Analyses/DataPooling` (normalisation, pooling, the target),
`Analyses/RestretchMetricAudit` (`rsK80`, cycle 5),
`Analyses/RecoveryOscillation` (the ripple, the noise budget, the weight correction),
`Analyses/ProtocolDesign` (sizing, order, sex balance),
`Analyses/RundownCorrection` (φ). Campaign: `Analyses/RecoveryOscillation/RunOptimRigAware.m`.
Every snapshot cost quoted here was **re-evaluated from disk**, not read from a file
header — a header has been wrong before.
