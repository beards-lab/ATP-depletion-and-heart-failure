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

## 2. Sequence rundown — the order matters, and it is *not* symmetric

This is the finding that should drive the experimental design.

| | damage | on which run | as a fraction | bias on the ATP ratio |
|---|---|---|---|---|
| **8 → 2** | 3.4 kPa | the 2 mM run (F ≈ 70 kPa, the **strong** one) | 4.8 % | **×0.952** |
| **2 → 8** | 7.4 kPa | the 8 mM run (F ≈ 50 kPa, the **weak** one) | 14.7 % | **×1.173** |

Two effects compound: the 2 mM record decays **×1.95 faster**, and the record it
degrades is **×1.40 weaker**, so the same absolute damage is a much larger fraction.
**The 2 → 8 order carries 3.6× the bias.**

### Consequence: balancing the design does *not* cancel the bias

| design (n = 6) | bias, uncorrected | bias, corrected | can validate? |
|---|---|---|---|
| 6× 8→2 | **−4.8 %** | −1.2 % | no |
| **5× 8→2, 1× 2→8** | **−1.2 %** | **−0.3 %** | **yes** |
| 4× 8→2, 2× 2→8 | +2.5 % | +0.6 % | yes |
| 3× 8→2, 3× 2→8 | +6.2 % | +1.6 % | yes |
| 6× 2→8 | +17.3 % | +4.3 % | no |

A "balanced" 3+3 design leaves **+6.2 %** — worse than all-8→2. Averaging a small
negative bias against a large positive one does not cancel.

**So: correction is required regardless of design.** The role of the 2→8 records is
not to cancel the bias by averaging — it is to **detect a broken correction**, which
shows up as a systematic difference between the two order groups.

---

## 3. Should rundown be considered at all? — Yes, and here is the decision

**Yes.** At 5–15 % of the measured ATP effect (which is ~34 %), it is not negligible;
uncorrected, the running order alone moves the answer from ×1.23 to ×1.65.

**How** depends on what the numbers are for:

### For MODEL FITTING — do not correct the data. Fit the fibre state.

Give each record **one rundown coordinate** λ, mapped jointly onto
`kstiff × (1 − 0.16λ)` and `kSE × (1 − 0.35λ)` (the lesion identified from the
5-segment profile). Fit λ per record alongside the physiology.

This is the answer to the consistency objection: nothing is corrected, so nothing can
be inconsistent. The model reproduces what was actually measured, on a fibre whose
state is an explicit parameter. It is **one nuisance dimension per record**, not two
free knobs — the two model parameters move together because they are two consequences
of one lesion.

It also makes φ, `T_dec` and own-peak referencing unnecessary for fitting: the model
carries the state.

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
| Correct intra-record? | **Yes, force only** (detrend + own-peak reference). `ktr` drift is below noise. |
| Correct for sequence? | **Yes** — the order bias is 5–15 %, and it is *asymmetric*: 2→8 is 3.6× worse. |
| Balance the design? | **No** — unequal biases do not cancel. Prefer 8→2, keep one 2→8 to *validate*. |
| Correct the data for fitting? | **No** — fit a per-record rundown coordinate λ instead. Avoids the consistency problem entirely. |
| Correct for reporting? | **Yes, force and `ktr` together from the same λ.** Never one without the other. |
| Which `ktr`? | **Protocol** for kinetics (×0.72); **slack** for rundown/compliance. |
| Biggest open risk | φ is calibrated on **one** bracket and assumed transferable. A second repeat, ideally of the opposite order, closes it. |

Reproduce: `RunDesignStrategy.m` → `results/design_strategy.png`.

![Design strategy](results/design_strategy.png)
