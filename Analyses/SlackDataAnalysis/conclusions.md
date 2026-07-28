# SlackDataAnalysis — survey of every slack dataset

**Headline.** Nine slack recordings audited through one extraction path. The Baker
set is a **different protocol** whose recovery windows are 2–6× too short, so its
amplitude features are truncated and its apparent low-ATP force *loss* is an
artifact. Absolute force is not transferable between preps (CV ~20 %), but the
**shape collapses to CV 4–7 %** once each prep's own strength is divided out.

`RunSlackDataAnalysis.m` → `results/`. Data only, no model.
Nine active datasets (4 protocol days × 8/2 mM, plus Baker 0.2 mM) and three
PNB+Mava passive controls.

> Rundown and the ATP effect were spun out into their own analyses:
> **[RundownCorrection](../RundownCorrection/conclusions.md)** (what rundown is
> and how to correct it) and
> **[ATPEffectReconciliation](../ATPEffectReconciliation/conclusions.md))**
> (the reconciled ATP numbers: force ×1.36, `ktr` ×0.55).

| figure | what it shows |
|---|---|
| `fit_<tag>.png` (12) | feature-extraction QC per dataset — the fits behind every number |
| `overview_traces.png` | all force traces, one tile per day, aligned on the first slack; passive control below |
| `features_all.png` | all datasets on the standard feature panels. Colour = day, **solid = 8 mM, dashed = 2 mM, dash-dot = 0.2 mM** |
| `features_normalised.png` | same, force features divided by each prep's own mean `steady` |

**Why features are re-extracted in-script.** The stored `features_data` were
written at different times by different script revisions (some smoothed, some
not). A cross-dataset comparison is only meaningful if every number comes from one
code path, so this script re-runs `extractSlackAttributes` with identical settings
on every dataset and never writes back.

**Acquisition order.** All days ran 8 mM → 2 mM **except 04/10**, which ran
2 mM → 8 mM. 03/27 additionally repeated the 8 mM at the end of the session.

---

## 1. The Baker dataset is a different experiment

Slack depth (`SLslack` 2.04 → 1.88) and the restretch ladder (4,5,6,7,3 ML/s) are
shared. **Nothing else is.**

| | Baker | 03/27, 04/03, 04/10 |
|---|---|---|
| inter-slack interval | 0.30 / 0.35 / 0.45 / 0.50 s | 0.60 s uniform |
| **slack hold** (= the `ktr` fit window) | **0.05 → 0.15 s** | **0.298 s** |
| release velocity | −89 … −174 ML/s | −36 … −71 ML/s (2.5× slower) |
| restretch velocity | 4, 5, 6, 7, 3 | same ✓ |

The slack hold **is** the window `A`, `Am`, `ktr` and `t0` are measured in.

## 2. …and its recovery is truncated, ATP-dependently

`A` is the fitted asymptote of `F = A(1−e^{−ktr(t−t0)})`; `Am` is the last
*measured* point in the same window. `Am/A` is the fraction of the asymptote
actually reached.

| dataset | ktr × hold (>3 ⇒ complete) | Am/A |
|---|---|---|
| Baker 8 mM | 1.8 – 3.7 | 0.84 – 0.96 |
| Baker 2 mM | **1.2 – 1.8** | **0.70 – 0.83** |
| Baker 0.2 mM | 1.1 – 1.7 | 0.67 – 0.84 (2/5 fits fail) |
| 03/27, 04/03, 04/10 (8 & 2 mM) | 4.0 – 16.7 | 0.94 – 1.05 |

Because `ktr` is ~0.55× at low ATP, the same fixed hold truncates the low-ATP
recovery *harder*. That manufactures an apparent low-ATP force loss in exactly the
truncation-sensitive features and in no others:

| Baker, 2 mM / 8 mM | value | truncation-sensitive? |
|---|---|---|
| `Am` (last measured) | **0.713** | yes — most sensitive |
| `A` (asymptote fitted from 1.2 τ of data) | **0.856** | yes — ill-conditioned |
| `steady`, `peak2` (read long after the restretch) | **1.006, 1.083** | no |

**Do not pool Baker amplitude features with the 2026 preps.** Its `ktr` is fine —
a *rate* is identifiable from 1.2 τ; an *asymptote* is not.

## 3. `ktr` is the most reproducible quantity in the dataset

`ktr(2 mM)/ktr(8 mM)` across **4 preps, 2 ATP orders, 2 eras, 2 protocols**:

| Baker | 03/27 | 04/03 | 04/10 | mean | **CV** |
|---|---|---|---|---|---|
| 0.566 | 0.537 | 0.593 | 0.492 | **0.547** | **8 %** |

## 4. Absolute force is not transferable; shape is

Between-prep CV at fixed ATP, raw and after dividing each prep's force features by
its own mean `steady`:

| feature | raw CV (8 mM) | normalised CV |
|---|---|---|
| `A` / `Am` | 23 % / 22 % | **7 % / 7 %** |
| `peak1_y` | 22 % | **4 %** |
| `peak2` | 19 % | **5 %** |
| `vall_y` | 21 % | **6 %** |
| `vall2_dy` | 28 % | 19 % |
| `restretchSlopeStart` | 23 % | 19 % |

Absolute 8 mM `A` spans 37.5 – 62.2 kPa across preps — biological variation in
cross-section, damage and rundown, not model error. **Fit normalised shape plus
one per-prep force scale**, not absolute kPa. Fitting absolute force to a single
prep is what ThursdayNightFever does: it matches 03/27's `A` to 0.1 % and
therefore misses 04/10 by 1.6–2×, baking one prep's scale into the kinetics.

## 5. Feature reliability ranking

| tier | features | evidence |
|---|---|---|
| **use** | `SLslack` (CV 0 %), `ktr` ratio (8 %), normalised `peak1_y` `peak2` `vall_y` `A` `Am` (4–7 %) | §3, §4 |
| **de-weight** | `t0`, `restretchSlopeStart`, `vall2_dy` (19–37 %) | |
| **do not fit** | `ovrsht_dy` (ratio CV **115 %**, sign flips: +1.35 / +0.34 / +0.37 / −0.01); `ktr2_*` (`ktr2_k` CV 41 %; 03/27 8 mM sits at 46 while every other dataset is 14–22; `ktr2_omega` = 158 for Baker 2 mM) | the 2-component fit is not identifiable on this data |

## 6. Remaining inconsistencies worth a look

* **04/03 `t0` is 1.5–2× everyone else** (8 mM 0.0152 s vs 0.0075–0.0102) under an
  identical protocol. `t0` is the delay before force redevelops = time to take up
  the imposed slack by unloaded shortening; same `SLdiff`, longer delay ⇒ **lower
  unloaded shortening velocity (Vmax) in that prep**. Its `vall2_dy` (−20.9) is
  also the extreme, and its length dependence stays steepest even after
  detrending. Check against that day's FV before using it.
* **Baker 0.2 mM shows a rigor signature**: normalised `peak2` and `vall_y` are the
  *highest* of all datasets (1.3–1.4) while `steady` is the lowest — a large
  elastic restretch response with little active force, i.e. many stiff
  non-force-generating bridges. Physiologically sensible and interesting, but its
  recovery fits fail on 2 of 5 segments (`t0` → 2e−10), so `A`/`ktr` are unusable
  there.
* **Passive (PNB+Mava) level varies 2–3× between days**: 03/27 restretch peaks
  21–24 kPa on a ~0 kPa baseline; 04/03 15–19 kPa on ~10 kPa; 04/10 8–12 kPa on
  2–5 kPa. TNF's default passive reference is the 03/27 file — the stiffest of the
  three, and still built by a stale script revision. A per-day passive reference is
  needed before passive↔active are jointly identified.
* **Segment 5 is protocol-different by design** in every dataset (restretch drops
  to ~2.9 ML/s, final hold is longer), which is why normalised `steady` is flat at
  ~1.05 for segments 1–4 then falls to ~0.8. Reproducible across all datasets — the
  model must reproduce it, but it is not an anomaly.

## 7. Detrend before extracting

Each recording decays during its own activation; the five slacks span ~2.7 s of it,
and the 2 mM recordings decay **1.5–2× faster**, so the drift does not cancel in the
ATP ratio. It accounts for 8–16 % of the apparent `A(1→5)` decline. Removing each
recording's own slope (fitted on the isometric plateaus that bracket the slack
train, so it costs nothing) shrinks the apparent length dependence by 2–4 pp:

| day | cond | A(1→5) raw | detrended |
|---|---|---|---|
| 03/27 | 8 / 2 mM | −21.5 / −24.6 % | −19.8 / **−20.8 %** |
| 04/03 | 8 / 2 mM | −31.2 / −32.4 % | −27.5 / −28.1 % |
| 04/10 | 8 / 2 mM | −23.6 / −24.5 % | −20.1 / −20.7 % |

It barely moves the ATP *ratio* but cleans up the **`A`-vs-`SLslack` slope, which
is a fit target**, and brings 03/27's two conditions into agreement where they were
3 pp apart.

---

## How to fit better

1. Fit **ratios**, not absolute levels — see
   [ATPEffectReconciliation](../ATPEffectReconciliation/conclusions.md) for the
   targets (force ×1.36 rundown-corrected, `ktr` ×0.55 raw).
2. Fit **normalised shape + one free per-prep force scale** (§4).
3. **Exclude Baker amplitudes** entirely (§2); keep its `ktr`.
4. **Retire `ovrsht_*` and `ktr2_*`** from cost functions (§5).
5. **Detrend every recording by its own slope** before extracting features (§7).
6. Re-measure or de-weight 04/03, and get a per-day passive reference (§6).
