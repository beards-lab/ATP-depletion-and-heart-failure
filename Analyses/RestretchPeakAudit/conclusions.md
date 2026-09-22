# The re-stretch peak detector was reading noise — 57 % of the objective was chasing it (2026-09-21)

Reproduce: `RunPeakAudit.m` → `RunSpreadReconcile.m` → `RunInstrument0410.m`.
Fix: `Model/extractSlackAttributes.m`. Re-curation: `DataCuration/RecurateSlackFeatures2026.m`.
Figure: `results/instrument_0410.png`.

---

## The bug

`extractSlackAttributes` used a **fixed `PEAK_MIN_PROMINENCE = 0.1 kPa`**. Against these
recordings' own noise that is:

| prep | F_iso | noise σ | 0.1 kPa = |
|---|---|---|---|
| 03/27 | 80.1 | 1.360 | **0.07 σ** |
| 04/03 | 70.0 | 1.324 | **0.08 σ** |
| 04/10 | 50.9 | 0.907 | **0.11 σ** |

So `findpeaks` returned **11–15 "peaks" per re-stretch**, almost all noise ripples — and the code
takes `peak1_t(1)`, **the first in time, not the largest**. The real overshoot has prominence
12–19 kPa; the ripples taken instead had 0.21–2.01 kPa:

```
03/27 seg1: TAKEN +0.30 ms prom  1.34 -> dSL 0.0006   LARGEST +3.05 ms prom 19.44 -> 0.0225
03/27 seg5: TAKEN +1.65 ms prom  2.01 -> dSL 0.0091   LARGEST +4.65 ms prom 15.41 -> 0.0272
04/10 seg1: TAKEN +0.80 ms prom  0.42 -> dSL 0.0033   LARGEST +3.80 ms prom 12.77 -> 0.0281
04/10 seg5: TAKEN +0.45 ms prom  0.63 -> dSL 0.0012   LARGEST +6.45 ms prom  6.07 -> 0.0375
```

**7 of 15 preparation × segment cells** were mis-detected — not one bad preparation, though 04/10
(the weakest, 50.9 kPa) was worst with 4 of 5.

> Taking the *largest* is not the fix either: on 04/03 seg4 and 04/10 seg4 the largest peak sits at
> +20 ms, the end-of-ramp rise. Raising the threshold is what makes *first* the right one.

### Why it mattered so much

`peak1_dSL = SL(t_peak) − SLslack`, and **SL is the commanded motor length ×2, not a measurement** —
identical across all three days to 0.1 sample. So `peak1_dSL` carries no length information at all:
it is peak *timing* re-expressed through a known ramp. During the re-stretch SL advances
**0.00247 µm per 0.25 ms sample**, one sample is **14 %** of the feature's mean, and the pooled
target's entire across-segment range was **9.3 samples**. A detector that slips a few samples spans
the whole feature.

## The fix

```matlab
PEAK_PROM_K_SIGMA  = 3;     % prominence floor in units of the trace's noise sigma
PEAK_PROM_MIN_FRAC = 0.05;  % ... and at least this fraction of the window's force range
peak_prom = max(PEAK_PROM_K_SIGMA*noise_sigma, PEAK_PROM_MIN_FRAC*range(y));
```

`noise_sigma` is a MAD-based estimate of the trace's own high-pass residual, so it degrades sanely
to ~0 on a smooth model trace and the fractional floor takes over. A threshold sweep gives
**identical results at 3, 4 and 6 σ** — a plateau, which is what a correct threshold looks like:

| threshold | between-prep rel SD |
|---|---|
| 0.07 σ (shipped) | 0.641 |
| 1 σ | 0.158 (unstable — 2 σ jumps to 0.284) |
| **3 σ** | **0.171** |
| 4 σ / 6 σ | 0.171 |

## What re-curation did

All three 8 mM files re-extracted through one code path. (`UpdateSlackFeatures`' list did not even
contain **04/03**, and its `featRunIdx` defaults to the 04/10 rows, so the three preparations'
stored features had not all come from the same code.)

| feature | old spread | new | change |
|---|---|---|---|
| **peak1_dSL** | 0.525 | **0.172** | **0.33×** |
| **vall_t** | 0.548 | **0.184** | **0.33×** |
| peak1_y | 0.320 | 0.239 | 0.75× |
| vall_y | 0.262 | 0.244 | 0.93× |
| A / steady / peak2 / ktr | — | unchanged | 1.00× |

04/10's `peak1_dSL` went `[0.0033 0.0045 0.0367 0.0180 0.0012]` → `[0.0281 0.0280 0.0367 0.0353
0.0375]`; its `vall_t` `[0.0010 0.0010 0.0070 0.0050 0.0005]` → `[0.0097 0.0085 0.0070 0.0050
0.0092]`.

Also surfaced: `rsK80` was **empty** in all three stored feature sets and is now populated. It was
only ever supplied through the pooled swap.

## Consequences for the objective

**The "wrong sign of segment-dependence" finding is retracted for three of its four features.**
The pooled `peak1_dSL` target was `[0.0074 0.0175 0.0303 0.0250 0.0090]` — a non-monotonic arc the
model was penalised for not reproducing. It is now `0.0278` mean and monotone. The model always
had the right shape. Same for `vall_t` and `peak1_y`. **Only `ktr`'s segment slope survives as a
genuine sign mismatch** (two of three preps fall across the ladder, the model rises).

Scored on the re-curated target with floors and weights re-derived from the new spreads
(`params/refit4_seed.m`, same parameters as `refit3_opt`):

| | before re-curation | after |
|---|---|---|
| peak1_dSL | 0.1326 | **0.0568** — despite a 5.4× tighter floor and 7× higher weight |
| peak1_y, vall_y, peak2, ktr, steady, rsK80, FV, all PS_* | paying | **all zero** |
| vall_t | 0.1112 | 0.1260 (floor 6.2× tighter, weight 8.3× higher) |

The objective is now `vall_t` 39 %, `peak1_dSL` 18 %, `ktrP_k63` 16 %, `vall2_dy` 11 %,
`restretchSlopeStart` 9 % — re-stretch transient timing plus the ktr rate. Coherent, and pointed
at things that are actually wrong.

## A methodological consequence, recorded

`buildPooledTarget`'s `fitFDist` now returns **d0 = Inf**: with the two wild features fixed, the
remaining twenty spreads (4–19 %) are statistically indistinguishable at n = 3. The formula had no
limit at that end and silently produced `s = NaN` for every feature, which would have poisoned
every floor and weight downstream. Now capped at `D0_CAP = 5` with a warning.

**The honest reading is that the old differential weighting was itself partly an artefact of this
bug** — `peak1_dSL` at 52 % and `vall_t` at 57 % supplied most of the apparent dispersion that
justified weighting features differently. Post-fix the weights are near-uniform (25–34). That is a
judgement call at the cap, and it should be argued with rather than inherited.
