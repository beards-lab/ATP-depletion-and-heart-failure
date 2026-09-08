# Methods — the paragraph

*This is the whole data-manipulation section. Everything behind it is in
[conclusions.md](conclusions.md); nothing behind it needs to appear in a paper.*

---

## The text

> **Normalisation and pooling.** Force-dimension observables were divided, per
> preparation, by that preparation's isometric force measured at SL 2.2
> immediately before the slack train of its 8 mM run (*F*<sub>ref</sub>). One
> divisor per preparation, shared by both ATP conditions, so the 2 mM / 8 mM
> ratio is independent of it. Because the two conditions are recorded
> sequentially and the running order differed between preparations, the run
> recorded second in each session was corrected for rundown by
> 1 + φ·|*s*<sub>1</sub>|·*T*<sub>1</sub> / *F*<sub>2</sub>, with φ = 0.452 and
> *s*<sub>1</sub>, *T*<sub>1</sub> the first run's force slope and activation
> duration (Rundown­Correction, methods §M.5). Features were extracted from each
> preparation separately and then averaged; traces were never averaged before
> extraction. Each feature's uncertainty is the between-preparation SD of the
> normalised value, moderated by an empirical-Bayes prior estimated from the
> feature set as a whole (Smyth 2004), because *n* = 3 leaves two degrees of
> freedom per estimate. Absolute tension, which the normalisation removes, is
> reported separately: *F*<sub>ref</sub> = 69.3 kPa, CV 15.5 % (*n* = 3,
> rundown-corrected).

167 words. Four sentences of procedure, one of statistics, one of disclosure.

## What each clause is doing, in one line

| clause | why it is there |
|---|---|
| divide by *F*<sub>ref</sub> | cross-sectional area is estimated, not measured; it is a per-preparation multiplicative nuisance worth 22 % |
| **one** divisor per preparation | makes the ATP ratio nuisance-free — a per-*run* divisor would define it away |
| *F*<sub>ref</sub> at SL 2.2, pre-slack, 8 mM | same length, same activation, ~40 ms of clean plateau; nothing else in the record is measured better |
| rundown lift on the second run | the running order is not the same in every preparation, so without it the ATP effect depends on the order |
| extract then average | feature extraction is non-linear (exponential fits, peak picking, argmax), so the two routes are different estimators — and on the re-stretch, a feature window that runs past the averaging grid breaks outright (RestretchMetricAudit §2) |
| moderated SD | a 2-dof SD has a 95 % CI spanning 12×; unmoderated, one lucky feature takes over the objective |
| *F*<sub>ref</sub> reported | the level is real data; normalising it away and not saying so leaves the fit scale-blind |

## The objective this produces

```
chi2  =  SUM_k  || f_k^model - fbar_k ||^2 / s_k^2          shape, normalised
      +  (F_ref^model / 69.3 - 1)^2 / 0.155^2               the one level term
      +  SUM_{k in ATP}  (r_k^model - rbar_k)^2 / sr_k^2    the ATP ratios
```

`fbar_k`, `s_k`, `rbar_k`, `sr_k` come from `buildPooledTarget.m`. Five features
clear *z* ≥ 2 for the ATP block: the force **level** (×1.32), the re-stretch
recovery **rate** `rsK80` (×0.381, *z* = 3.4) and its **delay** `rsT20` (×1.75,
2.4), `vall2_t` (×1.92, 2.6) and `ktr` (×0.53, 2.3). Every force *shape* ratio is
1.00 ± 0.05 and would only add a constant.

`rsK80` and `rsT20` replace the older `rsK`/`rsT0` pair on the strength of
[RestretchMetricAudit](../RestretchMetricAudit/conclusions.md): the amplitude is
taken from `steady` rather than from the tail of the fit window and the rate is
read off the 20–80 % band of the rise, which drops the between-preparation SD
from 17.2 % to 6.1 % and makes the estimate invariant to the window (a known
time-dilation is recovered to 6 % instead of 24 %; truncating the window to 60 ms
moves it by 0 % instead of 77 %). Segment-averaged quantities use segments 1–4:
segment 5 re-stretches to ML 1.00 rather than 1.10, by the smallest amount, and
is not part of the amplitude ladder.

## Reproduce

```matlab
T = buildPooledTarget();                              % the target + its error bars
run('Analyses/DataPooling/RunTargetComparison.m')     % old vs new target vs model
run('Analyses/DataPooling/RunNormalizationLadder.m')  % why this normalisation
run('Analyses/DataPooling/RunPoolingComparison.m')    % why extract-then-average
run('Analyses/DataPooling/RunFVConsistency.m')        % the in-preparation F-V point
run('Analyses/DataPooling/RunSummaryFigure.m')
```
