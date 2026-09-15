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

---

## 7. The campaign ran, and it bought the wrong thing (2026-09-11)

14 h on Great Lakes, four profiles, killed early (§8). Every snapshot
re-evaluates here within 0.016 of its cluster header, so the numbers are real:

| profile | rounds | cluster header | re-evaluated here | vs seed |
|---|---|---|---|---|
| **a** (control) | 13 | 2.4153 | **2.4088** | **−10.3 %** |
| **b** (wide draw) | **2** | 2.4063 | **2.4183** | −10.0 % |
| c (hard kick) | 65 | 2.6598 | 2.6758 | −0.4 % |
| d (surrogate) | 1 | 2.6748 | 2.6756 | −0.4 % |

The cluster reads the seed at **2.7008** against 2.6858 here — a systematic
+0.015 platform offset (different MATLAB release). It is inside the 0.076 σ
noise and cancels in any within-machine comparison, but absolute costs must not
be quoted across machines without it.

**Then look at where the 10 % came from.**

```
total drop          2.6858 -> 2.4088   = -0.2770
FV_normpowerAvg     1.0743 -> 0.7690   = -0.3053
everything else, net                   = +0.0283
```

**The F–V term alone accounts for 110 % of the gain. Strip it and the fit got
1.8 % WORSE** (non-FV 1.6115 → 1.6398). The terms the campaign existed to move
did not move, or moved backwards:

| | seed → best | |
|---|---|---|
| `peak1_y` | 0.2619 → 0.2473 | −0.0145, ~5 % of the term |
| `peak1_dSL` | 0.2364 → 0.2446 | **worse** |
| `vall_t` | 0.1751 → 0.1801 | **worse** |
| `ktr` | 0.0650 → 0.0737 | **worse** — a solved rate |
| `rsK80` | 0.0188 → 0.0244 | **worse** — a solved rate |
| `coolDownLS` | 0.1340 → 0.1483 | **worse** |

**This is the spreadw campaign again** (3.85 → 0.98, "88 % of it FV"), and the
cause is structural, not bad luck. `FV_SHARE = 0.40` hands F–V 1.07 of a 2.69
objective while its χ² is **0.10 — a third of one preparation-width**. It is the
smoothest, most scale-blind, most buyable term in the set, and it already sits
far inside the disagreement between its own two source curves (11.3 %, between
LABS). Any search spends there first. **Lowering the pin 0.70 → 0.40 was not
enough: 0.40 still defends a solved problem.**

`c` is the control that proves the mechanism. Its draw is the six compulsory
parameters only, it cannot reach F–V's levers, and it spent 65 rounds going
nowhere (−0.4 %, every ±25 % kick returning 2.79–3.34 against 2.6598). `b`
reached the same place as `a` in **2 rounds against 13** — the wider draw does
search better, but what it found faster was still F–V.

**Do not relaunch this objective.** The next campaign needs F–V *floored*, not
merely de-weighted: a term below one preparation-width should stop paying.

## 8. Why three of four died

`a`, `b` and `d` end **mid-round, mid-simplex, with no error, no budget message
and `stall = 0`** — `a` in round 14 at 2.4153, `b` in round 3 at 2.4063, both
still descending. An abrupt end with nothing in the log is a signal, not a
crash: four sessions × (1 client + 5 workers) = 24 MATLAB processes, and the
narrowest profile (`c`) is the one that survived 20 h. Budget ~2 GB per
process, or run two profiles rather than four.

`d` is a separate, fixable defect: **10,380 `Failure at t=0. Unable to meet
integration tolerances`** warnings, each with a full stack trace, produced a
**19.8 MB** log in a single round. `surrogateopt` samples the raw 88-parameter
box globally and most of that box is not integrable. Either give the surrogate
tighter bounds or suppress that warning inside the objective — as it stands,
profile `d` spends its budget writing a log.

## 9. The fix: a floor, not another pin (2026-09-11)

§7's failure is not fixable by re-weighting. Any nonzero weight on a smooth,
nearly-solved term still offers the cheapest descent in the problem, so a fixed
budget will always be spent there — 0.70 → 0.40 was tried and §7 is the result.
What has to change is not how much a term is worth but **whether it pays at all
once it is inside the noise**.

`evalFeatureCost` gains an `fl:X` token: subtract X from the residual, clamp at
zero. The cost is a sum of squared relative residuals, so one preparation-width
over `nEl` elements is `nEl*s²`; `BuildWeightedObjective` gains `FLOOR_CHI2` and
emits `fl: FLOOR_CHI2 * nEl * s²` for every data term. Ten tests in
[`TestFeatFloor.m`](../../Workbench/Tests/TestFeatFloor.m), including the two
properties that matter: it is a **deadband, not a deletion** — a floored term
still costs more as it drifts further out, so it defends itself — and the NaN
penalty is applied *after* the floor, because a NaN is a broken evaluation, not
a small residual.

`BuildFlooredSeed.m` → `params/floored_seed.m`, from `rigaware_a_opt` (the
campaign's best, 2.4088 here) at `FLOOR_CHI2 = 1`:

```
2.4088 -> 0.9081   (1.5007 = 62 % of the objective was residual inside
                    one preparation-width, and now earns nothing)
```

and the composition is the point:

| | before | after |
|---|---|---|
| **FV_normpowerAvg** | **40 %** | **0 %** |
| `ktr` | 3 % | 1.0 % |
| `rsK80` | 1 % | 0 % |
| **re-stretch block** | ~28 % | **47.7 %** |
| `peak1_y` + `peak1_dSL` | ~20 % | **34 %** |

F–V, `ktr` and `rsK80` can no longer be bought. The search has nothing cheap
left, which is the whole intent.

**Watch this instead:** the unfloored trace-least-squares terms are now the
second block — `coolDownLS` 16.3 %, `PS_rampupRMSE` 14.4 %, `PS_holdDecayRMSE`
10.5 %, together **41 %**. They are *correctly* unfloored (§2: they sit 2–6×
above their own noise), but `coolDownLS` is anchored to 03/27's raw trace in
absolute kPa while the feature block targets the pooled level
([DataPooling §8](../DataPooling/conclusions.md)), so it pulls against the
pooled target. If the next campaign's gain lands there, that contradiction is
the reason and the fix is to rescale the trace windows by `Fref/Fss22`.

**Two bugs found on the way, both of which had already propagated:**

* `isData` classified by bare feature name, so the **boundary** entry
  `steady[1]|52.5-87.7|5` counted as a fit term and had its deliberate weight
  of 5 rewritten to 0.46 — a 10× weakening of the anchor guard. Boundary
  entries never read `feats_data`, so `1/s²` is meaningless for them.
* Because of that, `rigaware_a_opt` carries the guard **twice**, at 5 and at
  0.547416. Both cost 0 at every point evaluated so far, so nothing published
  is wrong, but the guard's strength was not the number anyone chose. The
  anchor block now replaces every existing guard on the same feature+range
  rather than appending.

---

## 10. How to run it

One name selects the whole campaign. [`campaignSpec.m`](campaignSpec.m) holds
what each starts from, and the launcher, the SLURM job, the campaign script,
`campaignStatus` and the consolidation report all read it from there — because
consolidating one campaign's snapshots against *another* campaign's seed is
exactly how a "−10 %" gets printed for a search that went nowhere.

| campaign | seed | cost | status |
|---|---|---|---|
| `rigaware` | `params/rigaware_seed.m` | 2.6858 | **ran, superseded** (§7). Kept only so the result reproduces |
| **`floored`** | `params/floored_seed.m` | **0.9081** | the one to run |

```bash
# 1. pre-flight -- one evaluation, no search. Always do this first.
matlab -batch "addpath(genpath('.')); PROFILE='a'; DRYRUN=true; \
    run('Analyses/RecoveryOscillation/RunOptimRigAware.m')"

# 2. launch. CAMPAIGN defaults to `floored`; -Go / GO=1 is the safety catch.
.\Analyses\RecoveryOscillation\launch_rigaware.ps1 -Hours 30 -Go       # Windows
GO=1 HOURS=30 Analyses/RecoveryOscillation/launch_rigaware.sh          # Linux
sbatch --account=<acct> Analyses/RecoveryOscillation/rigaware.sbatch   # Great Lakes

# 3. watch, from the state files rather than the logs
matlab -batch "addpath(genpath('.')); campaignStatus"

# 4. consolidate
matlab -batch "addpath(genpath('.')); CAMPAIGN='floored'; \
    run('Analyses/RecoveryOscillation/RunRigAwareReport.m')"
```

Pre-flight passes at exactly 0.9081, and the composition is what the exercise
was for — F–V is gone from the top ten entirely:

```
peak1_y            0.1558   17.2 %      PS_holdDecayRMSE   0.0954   10.5 %
peak1_dSL          0.1535   16.9 %      PS_steady20        0.0277    3.0 %
coolDownLS         0.1483   16.3 %      vall2_dy           0.0252    2.8 %
PS_rampupRMSE      0.1306   14.4 %      PS_steady22        0.0199    2.2 %
vall_t             0.0990   10.9 %      steady             0.0194    2.1 %
```

**Run two profiles, not four** (§8: four sessions is 24 MATLAB processes on one
node and three of them died). `a` and `b` are the pair worth the cores — `c`
cannot reach the levers that matter and `d` spends its budget writing a log.

**Read the result against §7, not against the total.** The question is not
whether the number fell, it is whether `peak1_y` and `peak1_dSL` fell. The
report flags a protected term (`rsK80`, `ktr`, F–V) that grows by more than
20 % of itself, which is what §7 needed and did not have.

---

## 11. The campaign was resumed, and this time it moved the overshoot (2026-09-15)

`params/rigaware/` is a **continuation** of §7, not a new run: same unfloored
objective (`rsK80[1:4]`, no `fl:` tokens), same seed, resumed from the §7 state.
Rounds 13→22 for `a`, 2→15 for `b`, 65→98 for `c`. Every snapshot re-evaluated
here.

| profile | rounds | §7 best | now | vs seed | stall |
|---|---|---|---|---|---|
| a (control) | 13 → **22** | 2.4088 | 2.2944 | −14.6 % | 1 — **converged** |
| **b (wide draw)** | 2 → **15** | 2.4183 | **2.0271** | **−24.5 %** | 0 — **still descending** |
| c (hard kick) | 65 → **98** | 2.6758 | 2.6758 | −0.4 % | 0 |
| d (surrogate) | — | 2.6756 | 2.6756 | −0.4 % | no state file returned |

### The profile question is now settled, and on equal evaluations

`b` did 15 rounds × 90 evaluations = **1350**; `a` did 22 × 60 = **1320**. Same
budget to within 2 %, and `b` is **0.267 lower**. This is not `b` being further
along — it is `b` searching better, which is what §6 predicted for the right
reason: the four open terms are four readings of one overshoot, and a
10-parameter draw does not span the direction that moves them together.

`a`'s history flattens (2.284 for four rounds, then a failed kick at 2.406);
`b`'s is **monotone over all 15 rounds** and was still falling when the run
ended. **`b` is not converged.**

`c` is finished as an idea. Its 98 rounds are a 49-cycle sawtooth — incumbent
2.66, kick to 2.8–10.7, revert to 2.66, repeat — so the six compulsory
parameters genuinely cannot move this objective, and 33 further rounds bought
**nothing**. Do not give it cores again.

### Is it F–V again? Half the disease, not all of it

```
total          2.6858 -> 2.0271   = -0.6587
F-V            1.0743 -> 0.5087   = -0.5657   86 % of the drop
everything else 1.6115 -> 1.5184  = -0.0930
```

86 % is still F–V, so the **headline number remains mostly unreal** and §9's
diagnosis stands unchanged. But the difference from §7 is the one that matters:
there the non-F–V fit got **1.8 % worse**; here it got **5.8 % better**, and the
gain is in the right place — `peak1_y` −29 %, `peak1_dSL` −21 %, `vall_y` −48 %.

### The honest test: score it where F–V pays nothing

Every candidate re-evaluated under the **floored** objective, so F–V is
deadbanded to zero and only movement outside one preparation-width counts:

| | floored cost | |
|---|---|---|
| `floored_seed` ( = §7's `rigaware_a_opt`) | 0.9081 | the current campaign seed |
| `rigaware/rigaware_a_opt` | 0.8961 | −1.3 % |
| **`rigaware/rigaware_b_opt`** | **0.8072** | **−11.1 %** |

**The gain is real and it is in the defect.** Under the floor, where F–V
contributes exactly nothing:

```
peak1_y     0.1558 -> 0.1033   -34 %     <- the lengthening overshoot
peak1_dSL   0.1535 -> 0.1046   -32 %     <- and where that peak sits
vall_t      0.0990 -> 0.0896    -9 %
A, steady, vall2_dy, PS_steady22          all improved
```

This is the **first movement on the open structural defect** in this whole
line of work. §7's campaign did not manage it; `b` did.

### What it cost, read through the floor

The unfloored table flags `rsK80` at **+91 %** — but under the floor `rsK80`
still pays **exactly zero**: even after that regression it is inside one
preparation-width. The floor is doing precisely its job, telling us the
headline regression is not worth defending. What the floor does *not* forgive:

* **`doublePeak` 0.0000 → 0.0056.** A shape guardrail that was exactly zero at
  every point evaluated so far is now nonzero. The model is starting to produce
  a second peak. **This is the one to watch**, because it is a qualitative
  change in the trace, not a quantitative one.
* `PS_holdDecayRMSE` +17 %, `coolDownLS` +4 % — the unfloored trace-LSQ block
  again, exactly where §9 said to watch it.
* `ktr` +24 % of a small floored term (0.0091 → 0.0114).

### What to do

**Re-seed the floored campaign from `params/rigaware/rigaware_b_opt.m`, not from
`rigaware_a_opt`.** It is 11 % better on the objective the campaign will
actually use, and the improvement is in the two terms that campaign exists to
move. `BuildFlooredSeed.m` needs one line changed — `BASE = 'rigaware/rigaware_b_opt'`.

**Run profile `b` alone, or `b` and `a`.** `c` and `d` have now spent two
campaigns producing nothing between them; §8's advice to run two profiles rather
than four is reinforced by the result as well as by the memory ceiling.

**Before adopting `b`, check it is not degenerate**:
`analyzeOptIterUniqueness('rigaware_b')`. Fifteen monotone rounds of a
16-parameter draw is exactly the shape that produces a sloppy solution, and a
best cost reached along a direction the data cannot see is not a result.
