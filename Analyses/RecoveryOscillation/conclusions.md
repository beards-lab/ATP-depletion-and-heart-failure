# RecoveryOscillation — the ripple on every force recovery is the rig, and the objective was already nearly immune to it

*Asked: the recoveries look like exponentials with a small oscillation on top;
if that is the measurement rig, does the objective need reweighting?
Answer: yes it is the rig — 90 % of it is plain transducer noise, not an
oscillation — and the objective survives that almost untouched. **One** term
did not, and it was the largest weight in the set.*

Scripts, in order: [RunOscillationID.m](RunOscillationID.m) ·
[RunNoiseBudget.m](RunNoiseBudget.m) · [RunObjectiveFloor.m](RunObjectiveFloor.m) ·
[RunFeatureNoise.m](RunFeatureNoise.m) → [BuildRigAwareSeed.m](BuildRigAwareSeed.m)
→ [RunOptimRigAware.m](RunOptimRigAware.m) → [RunRigAwareReport.m](RunRigAwareReport.m).

---

## 0. The short version

| claim | verdict |
|---|---|
| the recoveries are single exponentials | **yes.** What sits on top is 90 % white noise, 6 % a slow ramp-locked wander, 4 % unexplained |
| the ripple is a rig property, not biology | **yes.** σ is 1.19 kPa with the cross-bridges chemically switched off against 1.24 kPa with them working, and the locked part's period does not move with ATP |
| the objective is chasing it | **almost nowhere.** The three trace-least-squares terms carry 12.8 % of the objective and only **0.7 %** of it is below the noise floor |
| so the weights need adapting | **one term.** `rsK80`'s apparent between-preparation spread (16.9 %) is *smaller than its own measurement noise* (32.7 %), which is impossible — the feature was carrying no resolvable between-preparation information at all |
| the fix is to de-weight it | **no** — it is to stop scoring slack cycle 5. Noise 32.7 % → **5.2 %**, and the weight goes *up* |

---

## 1. The residual of a one-exponential fit, taken apart

[RunNoiseBudget.m](RunNoiseBudget.m) fits the same estimator the objective uses
(baseline = the minimum in the first 80 ms, amplitude fixed from the plateau,
only `k` and `t0` free) to all 90 recovery windows — post-slack and
post-re-stretch, three preparations, 8 mM / 2 mM / PNB — and splits what is
left over three ways, with no model of the ripple:

* **σ_w**, white transducer noise, from the flat plateau of the *same* window
  above 60 Hz. Two independent estimators (the high-passed plateau and the
  first difference) agree to within a few per cent across all 18 run × window
  types — the middle-bottom panel of `results/noise_budget.png`.
* **σ_c**, the part **phase-locked to the ramp**, by coherently averaging the
  five slack cycles of a run on their common post-ramp clock. Anything locked
  survives; white noise falls as 1/√5. The excess of the average's variance
  over σ_w²/5 is the locked part.
* **σ_s**, whatever is left.

```
post-re-stretch, active runs        post-slack, active runs
  residual RMS  1.31 kPa              residual RMS  1.50 kPa
    = 10.6 % of the rise                = 2.9 % of the rise
  white        90 %  (σ 1.24)         white        69 %  (σ 1.25)
  locked        6 %  (σ 0.28)         locked       21 %  (σ 0.68)
  unexplained   4 %  (σ 0.25)         unexplained  10 %  (σ 0.45)
```

**The "oscillation" is mostly not an oscillation.** Nine tenths of it is plain
broadband noise. The genuinely reproducible part is small and *slow* — its
dominant period is 190–260 ms, comparable to the window itself, so over one
recovery it is a wander, not a ring. There is no discrete resonance line
anywhere in the residual spectra ([RunOscillationID.m](RunOscillationID.m),
`results/oscillation_id.png`): the power is concentrated below 8 Hz and flat
above it.

### Three controls, all pointing the same way

* **Cross-bridges off.** In the PNB + mavacamten runs — same rig, same ramps,
  no active force — σ_w is **1.19 kPa** against **1.24 kPa** in the active
  runs. The noise is the transducer, full stop.
* **Nucleotide.** The locked part's period moves ×0.70 from 8 mM to 2 mM while
  every measured rate moves ×1.9–2.6 (`ktr` ×1.9, `rsK80` ×2.6). It does not
  track the kinetics.
* **Amplitude.** The ripple does not scale with the rise it sits on
  (r = −0.38), which a cross-bridge process producing it would.

### What that costs a model that has no noise

The data's own exponential reaches R² = 0.795 on the post-re-stretch window.
A model sitting *exactly on the truth* would reach **0.834**, because the
comparison is against a noisy trace. The model reaches 0.999. That gap is not
a defect and must never be scored as one.

---

## 2. How much of the objective is unreachable

[RunObjectiveFloor.m](RunObjectiveFloor.m). Three `fn` terms score a model
trace against a data trace by least squares with the data value pinned at 0.
Zero is not attainable — the floor is σ.

| term | σ in that window | floor | cost at `pooled8_opt` | irreducible |
|---|---|---|---|---|
| `coolDownLS` | 1.59 kPa | 0.378 | 0.1340 | 5 % |
| `PS_rampupRMSE` | 1.23 kPa | 1.232 | 0.1241 | 3 % |
| `PS_holdDecayRMSE` | 1.44 kPa | 1.444 | 0.0883 | 11 % |
| | | | **0.3464 (12.8 %)** | **0.0198 (0.7 %)** |

Those terms sit 2–6× above their own floor, so they are fitting signal. **No
weight change is warranted, and none was made.** The reason the objective is
this robust is not luck: [RestretchMetricAudit](../RestretchMetricAudit/conclusions.md)
already moved the recovery onto rates and levels and off trace least squares.
This analysis is that decision's receipt.

Two guardrails are also confirmed rather than changed:

* **`rsR2` must stay a boundary** (`0.4-1.0`), where it costs exactly 0. Turning
  it into a fit-to-data term would ask the model to reproduce 0.795 — a number
  that is 90 % noise and that a noise-free model cannot reach.
* **`ktr_rmse`'s `0-0.5` bound** asks the model's ktr rise to be well described
  by *one* exponential. The measurement has now shown that is the right shape,
  so the guardrail is earned rather than assumed.

---

## 3. The one term that was measuring the rig

[RunFeatureNoise.m](RunFeatureNoise.m) adds each run's **own measured noise**
to its **own trace**, 200 times, and re-extracts every feature through the
production extractor. Biology is held exactly fixed, so the spread across
realisations *is* the measurement noise of that feature, `s_noise`. Compare it
with `s_prep`, the between-preparation spread the weights `w = λ/s²` are built
on:

| feature | s_noise | s_prep | verdict |
|---|---|---|---|
| **`rsK80`** | **32.7 %** | **16.9 %** | **impossible — s_prep is below its own noise** |
| `rsT0` | 23.1 % | 22.0 % | same, but not scored — *never score it* |
| `rsT20` | 14.2 % | 21.7 % | usable |
| `peak1_dSL`, `vall_t` | 23 %, 24 % | 52 %, 56 % | loose but honest |
| every force feature | 0.0–3.9 % | 16–22 % | biology |

A feature cannot be better determined *between preparations* than the noise on
a single measurement allows. `rsK80`'s 16.9 % is a small-*n* fluke — exactly the
trap [DataPooling](../DataPooling/conclusions.md) warns about, where *n* = 3
gives 2 dof and a spread measured at 10 % is only known to lie between 5 % and
63 %. With cycle 5 in, **the whole apparent spread of the re-stretch rate is the
transducer**, and the objective was weighting it as its joint-largest term.

### It is one slack cycle, and it was already known to be different

```
per-cycle s_noise of rsK80, %      cycle 1     2     3     4       5
03/27  8 mM                           15.1  12.8  10.5  11.9    28.9
04/03  8 mM                            7.1   7.9   5.8   5.7    44.0
04/10  8 mM                           11.1  13.1  10.4   7.1   245.0
03/27  2 mM                           10.3   9.2   8.9   8.7    35.6
04/03  2 mM                           11.7   9.3   6.6   6.2       -
04/10  2 mM                            9.1  11.5   8.2  11.3    30.8
```

Every run, both nucleotide conditions, no exception. `rsK80` is `ln4/(t₈₀−t₂₀)`,
so its precision is set by how well two crossings are located, i.e. by
σ/amplitude. Cycle 5 re-stretches **to a different final length** (1.000 ML,
not 1.100) by the **smallest** amplitude (0.061 ML against 0.080–0.140) — it is
not on the ladder at all, which
[RestretchMetricAudit §3](../RestretchMetricAudit/conclusions.md) established
and which was never carried through to the feature the objective scores.

Dropping it: **s_noise 32.7 % → 5.2 %**, `s_prep` 10.2 % → 7.7 % raw
(15.8 % moderated). There is signal again.

---

## 4. What changed, and what it is worth

`rsK80` → `rsK80[1:4]`. That is the whole change.

It is a **scoring-window** change, not a de-weighting: the weight follows from
the tighter spread and goes **up** (1/s² 35.0 → 40.1), which is the right
direction — cycle 5 was adding scatter, not information. Cost at `pooled8_opt`
falls 2.6999 → **2.6858**, and that 0.0141 is deterministic, not the objective's
0.076 roughness: it is exactly cycle 5's `rsK80` residual leaving.

Machinery added, all default-off:

* `applyFeatSelector` now takes ranges and lists — `rsK80[1:4]`, `peak2[1,3,5]`.
  Tested in [`Workbench/Tests/TestFeatSelector.m`](../../Workbench/Tests/TestFeatSelector.m),
  including the trap that `str2double('1,3,5')` is **135** (MATLAB reads the
  comma as a thousands separator), so a list must be matched before any numeric
  parse. `hasTopLevelComma` / `splitTopLevelComma` are now bracket-aware so a
  list cannot be mistaken for a grouped `fn` entry.
* `buildPooledTarget` takes `opts.segments`, per-feature, and the frozen
  `Model/pooledTarget8mM.m` records `P.seg` — **s must be measured on exactly
  the cycles the objective scores**, or the χ² is not a χ².

`P.val` is byte-identical: the target values did not move, only the scale
against which residuals are judged.

> **The honest headline is that the premise was right and the consequence was
> small.** The recoveries are exponentials, the ripple is the rig, and the
> objective — because it scores rates and levels rather than traces — was
> already immune everywhere except one feature. That one feature happened to be
> the audited re-stretch rate carrying the largest ATP signal in the set
> (z = 3.6), so it was worth finding.

---

## 5. Where the residual actually is

At the rig-aware seed (2.6858), in χ² units — 1.0 = one preparation-width,
per slack segment in brackets:

```
peak1_y     20.8  [4.2]   <- the lengthening overshoot
peak1_dSL   18.8  [3.8]   <- and where that peak sits after the release
vall_t      13.9  [2.8]
vall2_dy     8.0  [1.6]
A            6.7  [1.3]   steady 6.5 [1.3]   ktr 5.2 [1.0]
peak2        4.5  [0.9]   restretchSlopeStart 2.7 [0.5]
vall_y       2.0  [0.4]   rsK80 1.5 [0.3]    t0 1.5 [0.3]
FV          0.06 [0.01]
```

Nine of seventeen terms are at or below one preparation-width. **The entire
remaining residual is the shape and timing of the re-stretch transient**, and
those four terms are four readings of one overshoot. Nothing in this analysis
touches it — it is the open structural defect, and it is what the 30-hour
campaign is aimed at.

An argued option, deliberately left **off**: if the recovery is a pure
exponential then it *is* three numbers — amplitude, rate, delay. The objective
scores the amplitude (`vall2_dy`, `steady`) and the rate (`rsK80[1:4]`) but not
the delay, so a third of a description now known to be sufficient is
unconstrained. `rsT20` is the audited delay, its ATP effect is resolved
(z = 2.4) and its noise (14 %) is inside its spread (22 %).
`ADD_RST20 = true` in [BuildRigAwareSeed.m](BuildRigAwareSeed.m) turns it on —
but then the campaign tests two changes, not one.

---

## 6. The campaign

[RunOptimRigAware.m](RunOptimRigAware.m), four profiles differing in *how* they
search, not in what they score, run as independent sessions
([launch_rigaware.ps1](launch_rigaware.ps1) / [.sh](launch_rigaware.sh)):

| | draw | simplex | kick | surrogate | why |
|---|---|---|---|---|---|
| **a** | 10 | 60 | 0.10 | – | the pooled8 setting: the control, and it was still descending monotonically when its 4 h ran out |
| **b** | 16 | 90 | 0.10 | – | the four open terms are one coupled overshoot; a 10-parameter draw may not span the direction that moves them together |
| **c** | 6 | 40 | 0.25 | – | many more rounds and the hardest escape from a flat ridge — how `estiff125` stalled ([StairsPrediction](../StairsPrediction/conclusions.md)) |
| **d** | 10 | 60 | 0.10 | 120 | the only member that can leave the basin |

One evaluation is ~38 s at 5 workers, so profile **a** is ~40 min/round and
~45 rounds in 30 h. A session saturates at 5–6 workers because the parallelism
is *inside* one evaluation (over the five slack segments) — extra cores buy
more profiles, not faster rounds.

**Always dry-run first.** `DRYRUN=true` validates the seed, the objective's
coverage and the pool for one evaluation, and is the difference between finding
a broken seed now and finding it in 30 hours.

[RunRigAwareReport.m](RunRigAwareReport.m) consolidates, and refuses two
mistakes this project has already made: it **re-evaluates every snapshot from
disk** rather than trusting its header ([StairsPrediction](../StairsPrediction/conclusions.md):
a header was fake and a campaign was ranked on it), and it reports each
profile's cost at the round count of the *shortest* run, so the ranking is
budget-matched and not just wall-clock-matched.
