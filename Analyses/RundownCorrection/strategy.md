# Strategy — rundown, running order, and pooling for model validation

*A decision document. What to do, what not to do, and why. Supersedes the
per-feature recommendations scattered through `conclusions.md`.*

---

## 0. First, two corrections to what I previously said

**(a) "The ktr protocol is rundown-proof" is not established.** What is established
is a *model prediction*: the protocol `ktr` is ~5× less sensitive to added series
compliance than the slack `ktr`. It has **never been measured** on a run-down fibre,
because 03/27's repeat has no hi-res burst file. Treat it as a strong hypothesis with
a sharp prediction attached (§4), not as a fact.

**(b) "Correct force, do not correct `ktr`" was wrong.** It produces an internally
inconsistent dataset — a force belonging to a fresh fibre paired with a `ktr`
belonging to a degraded one, a state that never existed. The model is then asked to
fit an impossible fibre, and will distort real parameters trying.

The error was in the criterion. I rejected the `ktr` correction because it *increased
between-preparation CV*. But CV conflates two different things: whether the correction
is right, and whether preparations are biologically identical. For force, correcting
reduced CV — good evidence. For `ktr`, it increased CV — which is equally consistent
with the correction being right and the preparations genuinely differing. **Internal
consistency is the stronger requirement.**

---

## 1. Intra-record rundown — correct it, but only in force

Within one activation, force decays 1.0–1.6 %/s. The five slack segments span 2.7 s,
so the last is measured on a baseline 1.8–4.4 % below the first — and the 2 mM records
decay ~2× faster, so this does **not** cancel in an ATP ratio. It accounts for 8–16 %
of the apparent `A(1→5)` decline, which is a fit target.

**Do:** subtract each record's own linear slope (fitted on the two isometric ML = 1.0
plateaus that bracket the slack train) from the force trace before extracting features.
It is free — the slope is already measured — and it removes a systematic,
ATP-dependent distortion of the length dependence.

**Do not** attempt an intra-record `ktr` correction. The drift over the 2.7 s train is
~2 %, and segment-to-segment `ktr` scatter is ~5 %; it is below the noise floor. The
attempt to measure it directly (ktr burst at t≈71.6 s vs slack segment 2 at t≈75.5 s,
same SL, same activation) gave ratios of 1.09/0.86/1.45/1.22/0.94/0.73 — pure noise.
The residual inconsistency this leaves is ~2 % in `ktr`, an order below the effects
being studied.

**Also do:** reference each record to **its own activation peak**, not to a fixed
protocol time. 2 mM peaks ~2 s earlier than 8 mM, so a fixed window compares a 2 mM
record that has been decaying ~3.0 s with an 8 mM record that has been decaying ~0.7 s.
This costs nothing and removes ~4 % of bias that sits *inside* each record, where no
between-record correction can reach.

---

## 2. Sequence rundown — where the bias comes from, in full

### 2.1 What is biased against what

**BIAS = (the measured 2 mM / 8 mM force ratio) ÷ (the true ATP effect at a fixed
fibre state).** The true effect is what you would get if both conditions could be
measured on the *same* fibre at the *same* instant. They cannot, so one of them is
always measured on a more degraded fibre.

Notation:

* `a` — the true ATP effect, `F_2mM / F_8mM` at one fibre state.
* `L` — the **absolute** damage in kPa accrued during the **earlier** run's
  activation, `L = φ · |slope_earlier| · T_earlier`. Rundown is linear in kPa/s
  (§M.1), so the damage is an absolute subtraction, **not** a fractional one. *This
  is the root of the asymmetry.*
* The later run is the degraded one: `F_later_measured = F_later_true − L`.
* `f = L / F_later_measured` — the damage as a fraction of the run it actually hit.

### 2.2 The two orders, derived

**8 → 2** — the 8 mM runs first and does the damage; the **2 mM** run is degraded,
and it is the **numerator**:

```
    R_meas = F2_meas / F8            F2_true = F2_meas + L
    a      = (F2_meas + L)/F8  =  R_meas · (1+f)
    BIAS   = R_meas / a        =  1/(1+f)   < 1     the ratio UNDER-states a
```

**2 → 8** — the 2 mM runs first and does the damage; the **8 mM** run is degraded,
and it is the **denominator**:

```
    R_meas = F2 / F8_meas           F8_true = F8_meas + L
    a      = F2/(F8_meas + L)  =  R_meas / (1+f)
    BIAS   = R_meas / a        =  (1+f)     > 1     the ratio OVER-states a
```

Same `f` produces a **larger** bias when it lands in the denominator: `1/(1+f)`
deviates from 1 by `f/(1+f)` while `(1+f)` deviates by `f`.

### 2.3 The numbers, and why they are asymmetric

| | `L` (kPa) | degraded run | `f` | position | **BIAS** |
|---|---|---|---|---|---|
| **8 → 2** | 0.581 × 0.560 × 10.4 = **3.38** | 2 mM, F ≈ 70 kPa | **0.048** | numerator | **0.954** (−4.6 %) |
| **2 → 8** | 0.581 × 1.092 × 11.6 = **7.36** | 8 mM, F ≈ 50 kPa | **0.147** | denominator | **1.147** (+14.7 %) |

Three compounding factors, all pushing the same way:

1. **The 2 mM run decays ×1.95 faster** (1.092 vs 0.560 kPa/s) and runs ~1 s longer,
   so it inflicts **×2.17** the damage.
2. **It degrades the weaker run.** The 8 mM plateau is ~50 kPa against the 2 mM's
   ~70, so the same kPa is a **×1.40** larger fraction. Net: `f` is **×3.0** larger.
3. **That larger `f` then lands in the denominator**, where the same fraction
   produces a bigger bias.

Result: **the 2 → 8 order carries ≈3.2× the bias of 8 → 2**, and in the opposite
direction.

*(F ≈ 50 / 70 kPa are representative plateau forces. The real bias is prep-specific
and is computed per prep in [ATPEffectReconciliation](../ATPEffectReconciliation/conclusions.md);
04/10's actual `f` was 0.141.)*

### 2.4 Consequence: balancing does *not* cancel the bias

| design (n = 6) | bias, uncorrected | bias, corrected | can validate? |
|---|---|---|---|
| 6× 8→2 | **−4.6 %** | −1.2 % | no |
| **5× 8→2, 1× 2→8** | **−1.4 %** | **−0.3 %** | **yes** |
| 4× 8→2, 2× 2→8 | +1.8 % | +0.5 % | yes |
| 3× 8→2, 3× 2→8 | +5.1 % | +1.3 % | yes |
| 6× 2→8 | +14.7 % | +3.7 % | no |

A "balanced" 3+3 design leaves **+5.1 %** — worse than all-8→2 (−4.6 %), because it
averages a small negative bias against a large positive one.

**So correction is required regardless of design.** The role of the 2→8 records is
not to cancel the bias by averaging — it is to **detect a broken correction**, which
shows up as a systematic difference between the two order groups.

---

## 3. Should rundown be considered at all? — Yes, and here is the decision

**Yes.** At 5–15 % of the measured ATP effect (which is ~34 %), it is not negligible;
uncorrected, the running order alone moves the answer from ×1.23 to ×1.65.

**How** depends on what the numbers are for:

### Building ONE fitting target: six strategies compared

Identification is expensive, so the practical need is **one pooled 8 mM target and one
pooled 2 mM target** whose difference is the ATP effect and not a rundown artefact.
`../ATPEffectReconciliation/RunPoolingStrategies.m` evaluates six ways, by simulation
against a known truth (3 preps per order):

| strategy | bias | SEM | RMSE | assumes |
|---|---|---|---|---|
| P1 naive (ignore order) | +5.1 % | 0.025 | 0.073 | nothing |
| P2 balanced-count, arithmetic | +5.0 % | 0.025 | 0.072 | equal n both orders |
| P3 balanced-count, geometric | +4.6 % | 0.025 | 0.067 | equal n both orders |
| P4 first activations only | **0 %** | **0.219** | 0.219 | unbiased but *unpaired* |
| P5 correct with a transferred φ | −0.0 %\* | 0.022 | 0.022 | **φ transfers** |
| **P6 two-order solve (a *and* φ)** | **−0.1 %** | 0.024 | **0.024** | both orders present |

\* *P5's zero bias is an artefact of the simulation handing it the true φ.*

**Why the balanced-count proposal (P2/P3) leaves ~5 %.** Averaging three 8 mM records
from each order against three 2 mM records from each order *would* cancel if both
orders inflicted equal damage. They do not: `f ≈ 0.05` for 8→2 against `f ≈ 0.15` for
2→8 (§2.3). Averaging a 5 % deficit against a 15 % deficit leaves ~5 %.

**Why P4 fails despite being unbiased.** Every prep's *first* activation is undamaged,
so taking 8 mM from the 8→2 preps and 2 mM from the 2→8 preps is bias-free. But it is
**unpaired** — the two conditions come from different animals — so the ~20 %
between-prep force scale no longer cancels. Its SEM is **9× worse**. Pairing within
preparation is worth more than removing the bias.

**Why P6 wins.** With `f_i = φ·g_i` and `g_i = |slope_earlier|·T_dec_earlier / F_later`
measured directly per prep, the two orders give

```
    8 -> 2 :   R_i  =  a / (1 + phi*g_i)
    2 -> 8 :   R_j  =  a * (1 + phi*g_j)
```

— two unknowns in n equations. **φ is fitted, not transferred.** Robustness check:

| true φ | P5 bias | P5 RMSE | P6 bias | P6 RMSE |
|---|---|---|---|---|
| 0.40 | −1.3 % | 0.028 | −0.0 % | 0.025 |
| 0.58 *(= the transferred value)* | 0.0 % | 0.022 | −0.1 % | 0.025 |
| 0.80 | +1.5 % | 0.030 | −0.1 % | 0.024 |
| 1.00 | +3.0 % | 0.047 | −0.0 % | 0.024 |

P5 is unbiased **only** when the borrowed φ happens to be right; P6 is unbiased at
every φ, at a negligible cost in SEM. Since φ transferring between preparations is the
single largest untested assumption in the whole correction, P6 removes it.

**On the current data P6 returns `a = 1.344`** — identical to the independent
model-lesion answer — with φ = 0.889 and residuals of 0.016 / −0.012 / −0.004.
With 3 preps and 2 parameters there is only 1 degree of freedom, so that agreement is
encouraging but not yet a test. **3 preps per order gives SEM 0.024; that is the design
target.**

### For MODEL FITTING — correct the second record of each prep, and fit one prep at a time

**This replaces the λ recommendation I gave earlier, which was over-engineered.**

The λ idea was: fit a per-record rundown coordinate mapped onto
`kstiff × (1 − 0.16λ)`, `kSE × (1 − 0.35λ)`. Three objections, all fair:

* **λ is not a timecourse.** It is *one constant per record* — the fibre's state
  during that recording — not a within-record time dependence. I wrote it in a way
  that invited the opposite reading.
* **It does not match the objective.** You are fitting two *states* (8 mM, 2 mM), not
  a degradation timecourse. λ adds a dimension that is not the physics of interest.
* **It needs infrastructure you do not have.** Fitting λ requires a multi-record cost
  function with shared physiology and per-record nuisances. The current workflow is
  one `params0`, one `velocitytableonfile`, one `RunBakersExp`.

**The practical route instead.** Within one preparation the rundown problem is
exactly one number: the second record is degraded relative to the first. So:

1. **Fit one preparation at a time**, its 8 mM and 2 mM records together.
2. **Apply a single correction to that preparation's SECOND record only**, scaling
   *all* its force features and its slack `ktr` by that prep's own factors
   (computed from its own slopes; measured directly if it has a repeat).
3. The model then sees **two ATP states of one fibre at one state of damage** —
   which is exactly the objective, with no extra parameters.
4. **Repeat per preparation.** The physiology parameters should agree across preps;
   that agreement *is* the validation. Prep-to-prep force scale stays a free
   parameter, as it must (absolute force CV is 20 %).

This satisfies the consistency requirement — force and `ktr` are corrected by the
same lesion, so the corrected record describes a real fibre state — without any new
fitting machinery. `features_data` is already the fit target, and correcting stored
features is a small, auditable step.

**When λ would be worth building:** only if you later want to fit *all* preparations
simultaneously with shared physiology. Then λ is the right structure — one nuisance
dimension per record, both parameters moving together because they are two
consequences of one lesion. It is a fallback, not the recommendation.

### For REPORTING a summary ATP effect — correct, consistently

Use the model-lesion correction (`../ATPEffectReconciliation`), which corrects force
**and** `ktr` from the same λ. Do not correct one and not the other.

Current best estimates:

| quantity | value | basis |
|---|---|---|
| ATP effect on **force** | **×1.34** (CV 4 %) | 3 preps, model-lesion corrected |
| ATP effect on **`ktr`** | **×0.72** (CV 4 %) | ktr protocol, 3 preps, *uncorrected* |

The `ktr` protocol value is quoted uncorrected because its rundown sensitivity is
predicted to be ~3 %, below its 4 % CV — and because that prediction is untested (§0a).

### Which `ktr` to use where

| purpose | measurement | why |
|---|---|---|
| cross-bridge **kinetics** (ATP effect, model fitting) | **ktr protocol** | ~5× less compliance-contaminated |
| **rundown / compliance** state | **slack `ktr`** | it is the compliance-sensitive probe — that is a feature here |

Do not average them, and do not fit them with a shared rate parameter until
[RestretchVsKtrRecovery](../RestretchVsKtrRecovery/conclusions.md) is resolved.

---

## 4. Design for the three incoming datasets

**Recommended: 2× (8 → 2) and 1× (2 → 8), with repeats on two of them.**

| dataset | order | repeat | what it buys |
|---|---|---|---|
| A | 8 → 2 | **8-2-8** | measures φ for that prep instead of assuming it; directly gives the correction the 8→2 comparison needs |
| B | 8 → 2 | — | replication at the low-bias order |
| C | 2 → 8 | **2-8-2** if possible | validates the correction (opposite bias sign) *and* tests whether φ is ATP-dependent |

Rationale:

* **Majority 8→2** because its bias is 3.6× smaller — less to correct, less correction
  error.
* **At least one 2→8** because without it a systematic error in the correction is
  invisible. This is the single most valuable design feature after the repeat.
* **8-2-8 is the most useful repeat**: it brackets the 2 mM run in 8 mM units, which
  is exactly the correction the majority design needs.
* **2-8-2 tests the largest remaining assumption** — that φ transfers between ATP
  conditions. Given the 2 mM record decays ×1.95 faster, φ plausibly differs, and a
  single 03/27 bracket cannot tell.

**Also request, on at least one repeat: a hi-res `*_stiff_ktr.txt` and `*_slack.txt`.**
This is the decisive test of §0a and costs one extra burst recording. It predicts:
**the protocol `ktr` should fall ~3 % where the slack `ktr` falls ~12 %.** If they
fall equally, the compliance explanation is wrong and the slowing is genuinely
kinetic — which would change the ATP `ktr` number back toward ×0.58.

*Note for the analyst:* `mergeLogsAndBursts` currently **excludes** filenames
containing `repeat` from burst matching. That is a one-line change, needed before any
repeat's hi-res data will merge.

---

## 5. Pooling for model validation

1. **Fit ratios, never absolute levels.** Absolute 8 mM force spans 37.5–62.2 kPa
   across preparations (CV 20 %); normalised shape collapses to CV 4–7 %. Use one free
   force scale per preparation.
2. **Pair within preparation.** The ATP effect is a within-prep quantity; every
   prep-specific nuisance cancels in the ratio and none of it cancels in the levels.
3. **Weight by design quality**, not equally: preparations with a measured repeat
   carry a measured correction; those without carry an assumed one.
4. **Report the two order groups separately as well as pooled.** Agreement between
   them is the validation; a discrepancy is a broken correction, not noise.
5. **Precision:** residual between-prep SD after correction is 0.053, so SEM is 0.031
   at n = 3, 0.022 at n = 6, 0.018 at n = 9. Going from 3 to 6 preparations buys a
   third off the random error — but the systematic order bias does **not** shrink with
   n, which is why §2 matters more than sample size.
6. **Exclude Baker amplitudes** entirely (different protocol, truncated recovery); its
   `ktr` is usable.

---

## 6. Summary — the short version

| question | answer |
|---|---|
| Correct intra-record? | **Yes, force only** — detrend by the record's own slope, and reference to its own activation peak. `ktr` drift over the train is ~2 %, below the ~5 % segment noise. |
| Correct for sequence? | **Yes.** Bias is −4.6 % (8→2) or +14.7 % (2→8) — see §2 for the derivation. |
| Balance the design? | **No.** The biases are unequal, so 3+3 leaves +5.1 %, worse than all-8→2. Prefer 8→2; keep one 2→8 to *validate*, not to cancel. |
| How to build ONE pooled target? | **P6 — solve `a` and `φ` jointly from both orders** (§3). Unbiased without assuming φ transfers; SEM 0.024 at 3 preps per order. |
| Does balanced averaging cancel the damage? | **No, it leaves ~5 %** — the two orders inflict unequal damage (f ≈ 0.05 vs 0.15). Useful as a sanity check, not as the estimator. |
| How to correct for fitting? | **Correct the second record of each prep** (force *and* slack `ktr`, same lesion), fit one prep at a time. **Not** a per-record λ — that was over-engineered and needs machinery you do not have. |
| Correct for reporting? | **Yes, force and `ktr` together.** Never one without the other — that describes a fibre that never existed. |
| Which `ktr`? | **Protocol** for kinetics (×0.72); **slack** for rundown/compliance. |
| Biggest open risk | φ is calibrated on **one** bracket and assumed transferable. A repeat of the opposite order closes it. |

## 7. The recipe, end to end

1. **Per prep**, extract features from both records (detrended, own-peak referenced).
2. **Per prep**, form the 2 mM / 8 mM ratio for every feature — this cancels the ~20 %
   prep-to-prep force scale, which is larger than the rundown effect.
3. **Compute `g_i = |slope_earlier| · T_dec_earlier / F_later`** for each prep — all four
   quantities are directly measured.
4. **Fit `a` and `φ` jointly** across all preps of both orders (P6). If any prep has a
   sandwich repeat, add its directly measured damage as an extra constraint on `φ`.
5. **Build the two targets**: undo `φ·g_i` on each prep's degraded record (force *and*
   slack `ktr`, same lesion), normalise each prep, then average. That gives one 8 mM and
   one 2 mM target at a common, undamaged fibre state.
6. **Fit the model to those two targets** with a single free force scale.
7. **Validate** by checking the two order groups agree, and that `φ` from the fit matches
   `φ` from any sandwich prep.

Reproduce: `RunDesignStrategy.m` → `results/design_strategy.png`;
`../ATPEffectReconciliation/RunPoolingStrategies.m` → `pooling_strategies.png`.

![Design strategy](results/design_strategy.png)

![Pooling strategies](../ATPEffectReconciliation/results/pooling_strategies.png)
