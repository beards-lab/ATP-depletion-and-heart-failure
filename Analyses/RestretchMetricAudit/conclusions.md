# RestretchMetricAudit — rsK was a bad estimator, and two of the results built on it were artefacts

*Asked: is the re-stretch recovery rate badly measured? Answer: yes, and fixing it
retires two conclusions that were resting on it.*

Scripts, in order: [RunMetricAudit.m](RunMetricAudit.m) ·
[CheckPoolingArtefact.m](CheckPoolingArtefact.m) ·
[RunSegment5Check.m](RunSegment5Check.m) ·
[RunAmplitudeIsTheWeakLink.m](RunAmplitudeIsTheWeakLink.m) ·
[RunFinalMetricChoice.m](RunFinalMetricChoice.m).
The estimators themselves are in [rsEstimators.m](rsEstimators.m); the chosen one is
now `rsK80` in `Model/extractSlackAttributes.m`.

---

## 0. The short version

| claim before | claim after |
|---|---|
| trace averaging biases `rsK` by **−24 %**, more than a preparation width | it biases segments 1–4 by **−2 %**. The −24 % is **one segment** (5) whose 844 ms feature window is not covered by DataPooling's 570 ms splice grid |
| `dk/dSL` in the data is **−126** and the model's **+714**, so the model's amplitude dependence has the wrong sign | on the four comparable segments the data slope is **−4 to −30**. Segment 5 — a *different manoeuvre* — supplied 76 % of it. The model is still wrong, but it is wrong against a near-flat target |
| `rsK`'s 21 % between-preparation spread is 95 % biology | 21 % → **6.1 %** by changing the estimator. It was mostly the estimator |
| `rsK` is a re-stretch **rate** | it is a rate *plus* the tail of its own fit window. A known ×0.6 time-dilation comes back 24 % wrong; truncating the window to 60 ms moves it +77 % |

The fix, in one line: **take the amplitude from `steady` instead of from the tail of
the fit window, and read the rate off the 20–80 % band of the rise.**

```
A_ext = steady - B                    B = the located post-re-stretch minimum
z(t)  = (y - B) / A_ext               the shape, 0 -> 1
rsK80 = ln(4) / (t80 - t20)           the rate
rsT20 = t20                           the delay, kept separately
```

---

## 1. What was tested, and why those tests

Seven estimators, all on the same shape `z` so none of them can absorb a
force-scale error, differing only in *which part of the rise they look at*:

| | what it uses |
|---|---|
| `k_fit` **(current)** | least squares of `1-exp(-k(t-t0))` over the whole window, `k` and `t0` free |
| `k_t0z` | the same with the clock started at the valley (`t0 = 0`) |
| `k_band` | log-linear regression of `log(1-z)` on `t`, `z ∈ [0.15, 0.90]` — the top of the rise |
| `k_t50` | `ln2 / t50` |
| `k_2080` | `ln4 / (t80 - t20)` — toe- and plateau-free |
| `k_mrt` | `1/MRT`, `MRT = ∫(1-z)dt`, finite-window corrected — the only *linear* functional |
| `k_area` | `1/∫(1-z)dt`, uncorrected |

The suggestion to start the fit at the post-re-stretch valley was already in the
code (`fitRestretchRecovery` anchors on the `vall2` minimum); the suggestion to
rescale completely is what `z` is; the suggestion to fit only the top of the rise
is `k_band` / `k_2080`. All three were tried. **The linear-functional idea, which
looked like the principled answer to a pooling bias, fails outright**: the
post-re-stretch overshoot pushes `z` above 1, the integral cancels against
itself, and on 04/10 `k_mrt` returns 987 s⁻¹ against a true ~50. An estimator that
commutes with averaging is worthless if it does not measure the rate.

Six tests, the first of which is the only one with a ground truth:

* **T1 identity** — time-dilate a real trace by a known `c`; a rate estimator must
  return `c·k`. Run on the actual data, all 30 windows.
* **T2 delay** — inject a pure 5 ms onset delay; a *rate* must not move.
* **T3 span** — truncate the window to 60 / 100 / 279 ms.
* **T4 noise** — add the measured noise (σ = 1.29 kPa, read off the plateaus) 200×.
* **T5 pooling** — mean-of-features vs features-of-mean-trace, two alignments.
* **T6 target quality** — between-preparation SD, within-run repeatability,
  and the 2 mM/8 mM discriminability `z`.

T2 is passed by everything (−0.0 %): the delay is *not* leaking into any of the
rates. That eliminates the obvious suspect straight away.

---

## 2. The −24 % pooling result does not survive re-running

[CheckPoolingArtefact.m](CheckPoolingArtefact.m) reproduces DataPooling's number and
then takes it apart per segment:

| segment | 1 | 2 | 3 | 4 | 5 | mean |
|---|---|---|---|---|---|---|
| `rsK` from the pooled trace | 48.8 | 45.3 | 50.9 | 48.3 | **9.9** | 40.6 |
| mean-of-features, per prep | | | | | | 53.5 |

Segments 1–4 sit at 45–51 against a mean-of-features of 53.5 — a few per cent, the
same as [RunMetricAudit](RunMetricAudit.m)'s own window-local averaging finds
(−2.0 %). One number, 9.9, drags the average down by 24 %.

The reason is mechanical. `velocity_segment(5)` for the fifth segment is not the
end of a slack cycle but the end of the protocol at t = 78 s, so that window is
**0.844 s long against 0.279 s for the others**. DataPooling's grand average is
built on a grid ending at release + 0.570 s and spliced into 03/27's own
datatable, so the last ~500 ms of segment 5's window — including the plateau that
sets the fit amplitude — is one preparation, rescaled, with a 1.5 kPa step at the
seam. The amplitude is wrong and the fixed-`A` fit collapses.

> The conclusion "pool the features, never the traces" is still right, but not for
> the reason given. The trace-averaging penalty on the re-stretch is a few per
> cent, not a preparation width. What actually bit was a **window that leaves the
> splice**, and the same trap will bite anything else that averages on a relative
> clock and re-extracts on an absolute one.

---

## 3. Segment 5 is a different manoeuvre and should never have been in the ladder

Read straight off the velocity table ([RunSegment5Check.m](RunSegment5Check.m)):

| segment | slack to (ML) | re-stretch **back to** (ML) | re-stretch **amplitude** | ramp (ML/s) | window |
|---|---|---|---|---|---|
| 1 | 1.020 | 1.1005 | 0.080 | 3.9 | 0.279 s |
| 2 | 1.000 | 1.1006 | 0.100 | 4.9 | 0.279 s |
| 3 | 0.980 | 1.1006 | 0.121 | 5.9 | 0.280 s |
| 4 | 0.960 | 1.1006 | 0.140 | 6.8 | 0.279 s |
| **5** | 0.940 | **1.0002** | **0.061** | **2.9** | **0.844 s** |

Segments 1–4 are a clean ladder: one final length, amplitude rising 0.080 → 0.140 ML.
Segment 5 returns to a **different final length**, by the **smallest** amount, at the
**slowest** ramp. Plotted against `SLslack` it looks like the far end of the ladder;
it is in fact off the ladder entirely.

Consequences, all measured:

| | segments 1–5 | segments 1–4 |
|---|---|---|
| `dk/dSLslack`, current estimator | **−85** | **−20** |
| `dk/dSLslack`, `rsK80` | −36 | **−4** |
| between-prep SD, current | 20.2 % | 17.2 % |
| between-prep SD, `rsK80` | 8.2 % | **6.1 %** |
| ATP `z`, current | 2.43 | 3.09 |
| ATP `z`, `rsK80` | 6.14 | **6.79** |

Between-preparation spread is concentrated in segment 5 (39.9 % there against
9.4 % in segment 1 for the current estimator).

> **This is the finding with the longest reach.** `[[restretch_recovery_rate]]`
> recorded the model's amplitude slope as the open structural defect: data −126,
> model +714. On the comparable segments the data slope is between −4 and −30
> depending on the estimator, i.e. **the data says the re-stretch recovery rate
> barely depends on the re-stretch amplitude at all.** The model's +714 is still a
> defect, but the diagnosis changes from "wrong sign" to "a strong dependence
> where there is essentially none", and three quarters of the apparent data slope
> was one non-comparable segment.
>
> Segment 5 should still be *simulated and scored* — the model drives the same
> velocity table, so it is a legitimate constraint. What it must not be is the
> fifth rung of an amplitude ladder.

---

## 4. The weak link was the amplitude, not the exponential

Every estimator failed T1 and T3 in the same direction, which is the signature of
a shared cause rather than a formula problem. The shared thing is
`A = median(last 15 % of the window)`. Both failing tests move the *end* of the
window: dilating by `c = 0.6` pulls the plateau sample back to where `z ≈ 0.99`,
truncating to 60 ms pulls it back to `z ≈ 0.85`. `A` comes out small and every rate
comes out fast.

Supplying `A` from `steady` instead — the same number to within 1 % (11.59 vs
11.62 kPa on segment 1), but measured *outside* the fit window, at the reference
length, over the last 20 ms before the segment ends:

| | identity error, worst `c` | span error at 60 ms |
|---|---|---|
| current, `A` from window tail | 24.2 % | **+77.0 %** |
| current, `A` from `steady` | 8.0 % | −6.7 % |
| `k_2080`, `A` from window tail | 17.9 % | +59.1 % |
| `k_2080`, `A` from `steady` | **6.0 %** | **0.0 %** |

The model-free estimators go to *exactly* zero span error once `A` is external,
because they never touch the tail. That is the whole argument for them.

---

## 5. The decision

Ten combinations (5 rate formulas × 2 amplitude sources), segments 1–4, crossing
smoother fixed at 5 ms a priori:

| rate | A from | k (1/s) | s_betw | s_within | ATP r | z | dk/dSL | identity | span60 |
|---|---|---|---|---|---|---|---|---|---|
| **t20–t80** | **steady** | **47.5** | **6.1 %** | 3.6 % | **0.379** | **6.79** | −4 | **6.0 %** | **0.0 %** |
| t50 | steady | 34.5 | 6.8 % | 6.7 % | 0.493 | 4.52 | −30 | 8.9 % | 0.0 % |
| band 15–90 | steady | 49.7 | 15.5 % | 4.6 % | 0.342 | 4.16 | −30 | 8.7 % | 4.0 % |
| t0 = 0 | steady | 41.7 | 12.7 % | 3.4 % | 0.400 | 3.38 | −8 | 6.3 % | 9.2 % |
| current | steady | 51.2 | 15.9 % | 3.6 % | 0.335 | 3.25 | −9 | 8.0 % | 9.4 % |
| current | tail *(status quo)* | 50.6 | 17.2 % | 5.3 % | 0.348 | 3.23 | −20 | 24.2 % | 77.0 % |

**`ln(4)/(t80−t20)` with the amplitude from `steady` wins on every column that
matters** and is not a tuned choice: its *value* is invariant to the crossing
smoother (47.7 / 47.5 / 47.5 at 2 / 5 / 10 ms) where `t50`'s drifts 38.3 → 32.5.
`t50` at its best smoother has a lower spread than `t20–t80` at a bad one, so the
2 ms optimum in `RunAmplitudeIsTheWeakLink` panel (c) is not to be read as a
ranking — with `A` external, `t20–t80` is stable and `t50` is the one that drifts.

Per preparation, segments 1–4: **44.5 / 47.7 / 50.3 s⁻¹** (mean 47.5, SD 6.1 %).

### What it costs

`rsK80` is not a first-order rate constant. If the recovery is not first-order,
`ln4/(t80−t20)` is a *shape-weighted* rate, and comparing it to a model rate
constant is only legitimate because the model is scored through the same
estimator. That is the trade: it is reproducible and it resolves ATP, but it is
defined by its estimator rather than by a mechanism. The exponential's `rsR2` is
kept as the separate statement of "is it even first-order" (data 0.77, model
0.999 — see §7).

---

## 6. Both new features clear the noise floor for ATP

| feature | 2 mM / 8 mM | s | z | |
|---|---|---|---|---|
| **rsK80** | **0.381** | 6.5 % | **3.42** | the rate: ×0.38 at 2 mM |
| vall2_t | 1.924 | 24.4 % | 2.58 | |
| rsK *(old)* | 0.345 | 17.0 % | 2.44 | |
| **rsT20** | **1.752** | 23.4 % | **2.36** | the delay: ×1.75 at 2 mM |
| ktr | 0.532 | 10.4 % | 2.32 | |
| rsT0 *(old)* | 0.664 | 77.0 % | 0.42 | unresolvable |

Splitting rate from delay honestly turns **one** resolved ATP effect into **two**:
`rsT0`, the free delay of the old fit, had a 77 % spread and carried no
information; `rsT20`, the measured 20 % crossing time, has 23 % and a resolved
×1.75. Low ATP both slows the recovery *and* delays its onset, and those are now
two separate constraints on the model instead of one.

---

## 7. What the fixed metric now says about the model

`spreadw_opt` scored against the pooled target
([../DataPooling/RunTargetComparison.m](../DataPooling/RunTargetComparison.m)):

* `rsK80` 49.8 vs 47.5 ± 6.1 % → **z = 0.29**. The rate is right.
* `rsT20` 7.1 ms vs 8.6 ms → z = −0.81. The onset is slightly early but inside.
* `rsR2` **0.999** vs **0.770** → z = +1.71. The model's recovery is a *pure*
  exponential; the data's is not.

That last line is where the residual actually lives, and the timecourse figure
shows it plainly: the model reaches its plateau by ~360 ms while the pooled data
is still climbing at 420 ms. The model matches the fast phase and has no slow
tail. Because `rsK80` takes its amplitude from `steady`, a missing slow tail
barely moves it — which is correct behaviour for a *rate* metric, but it means
**the slow phase needs its own feature.** Candidate: the fraction of the recovery
completed by a fixed time (e.g. `z` at 60 ms), which is what a two-phase recovery
differs in and a one-phase one cannot fake.

---

## 8. Limitations

* n = 3. Every spread here has a 12× confidence interval. The *ranking* of the top
  two estimators is not resolvable; the gap between them and the status quo is.
* The choice of the 20–80 % band, like any band, is a convention. It was fixed
  before the target numbers were computed, and the value's insensitivity to the
  smoother is evidence it is not tuned, but a different band would give a
  different number for the same physiology.
* `rsK80` is defined by its estimator (§5). Any comparison to a literature rate
  constant needs that estimator applied to the literature trace.
* Segment 5 is kept in the target (all per-segment values are stored) but excluded
  from the segment-averaged summaries and from any amplitude trend.
* The measured noise test (T4) adds a full measurement noise on top of data that
  already contains one, so its CVs are upper bounds; the between-preparation SD is
  the observable that matters and it is reported separately.
