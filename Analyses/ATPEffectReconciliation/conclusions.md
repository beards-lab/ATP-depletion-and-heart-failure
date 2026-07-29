# ATPEffectReconciliation — what 2 mM ATP actually does, across all four preps

**Headline.** After rundown correction, four protocol days agree:

> **Low ATP (2 mM vs 8 mM) makes the muscle ~×1.36 STRONGER and ×0.55 SLOWER.**
> Force CV falls from 22 % to 11 %; `ktr` was already consistent at CV 8 %.

`RunATPReconciliation.m` → `results/atp_reconciliation.png`, `.mat`.
Correction from [RundownCorrection](../RundownCorrection/conclusions.md);
dataset provenance from [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md).

---

## The problem

Three days ran 8 mM → 2 mM; **04/10 ran 2 mM → 8 mM**. Raw, they disagree badly on
force (×1.18 to ×1.73) while agreeing tightly on `ktr` (×0.49 to ×0.59). That is
the signature of rundown, not biology: whichever condition is recorded *second* is
degraded, so the reversed-order prep carries the opposite bias.

## The correction

Damage accrues during the **earlier** run's activation, so the **later** run is the
degraded one and the fractional correction must be referenced to *its* force:

| day | order | loss (kPa) | F_later | frac of later |
|---|---|---|---|---|
| 03/27 M | 8→2 | 3.18 | 83.59 | 3.8 % |
| 04/03 F | 8→2 | 5.28 | 63.54 | 8.3 % |
| 04/10 M | **2→8** | 7.31 | 43.03 | **17.0 %** |

04/10's correction is much larger because its earlier run was the 2 mM one, which
decays ~2× faster within-run — so the reversed order is not merely a sign flip, it
is a *bigger* correction.

## Result

| feature | raw 03/27 / 04/03 / 04/10 | raw CV | **corrected** | corr CV |
|---|---|---|---|---|
| `A` | 1.184 / 1.207 / 1.725 | 22 % | 1.229 / 1.307 / 1.475 | **9 %** |
| `Am` | 1.178 / 1.177 / 1.747 | 24 % | 1.222 / 1.275 / 1.494 | 11 % |
| `steady` | 1.180 / 1.296 / 1.651 | 18 % | 1.225 / 1.403 / 1.412 | **8 %** |
| `peak1_y` | 1.254 / 1.092 / 1.959 | 32 % | 1.301 / 1.182 / 1.674 | 19 % |
| `peak2` | 1.179 / 1.178 / 1.855 | 28 % | 1.223 / 1.276 / 1.586 | 14 % |
| `vall_y` | 1.198 / 1.194 / 1.858 | 27 % | 1.243 / 1.294 / 1.588 | 14 % |
| **all force** | | **22 %** | **×1.36** | **11 %** |

`ktr`, deliberately **not** corrected (correcting it degrades consistency — see
RundownCorrection §5):

| 03/27 | 04/03 | 04/10 | Baker | mean | CV |
|---|---|---|---|---|---|
| 0.537 | 0.593 | 0.492 | 0.566 | **0.547** | **8 %** |

`ktr` is the most reliable number in the whole dataset: four preps, two ATP
orders, two eras, two different slack protocols.

![The ATP effect, reconciled](results/atp_reconciliation.png)

*Top left: force ratio per prep, raw vs corrected. Top middle: per feature (dashed =
raw, solid = corrected) — the three preps converge. Top right: `ktr`, already
consistent without any correction. Bottom left: the correction moves the preps
together, not apart.*

## The trustworthy picture

| observable | effect of 2 mM vs 8 mM | confidence |
|---|---|---|
| isometric / quasi-isometric **force** | **×1.36** | good — 6 features × 3 preps, CV 11 %, after a correction validated two independent ways |
| **`ktr`** (force redevelopment rate) | **×0.55** | very good — CV 8 %, uncorrected, 4 preps |
| force × `ktr` (turnover proxy) | ×0.74 | the slowed cycle outweighs the extra force |
| length dependence (`SLslack` slope) | unchanged in shape | see SlackDataAnalysis §5 |

Physiologically consistent: fewer ATP molecules → slower ADP release/detachment →
bridges dwell longer in force-bearing states → **more force per bridge-cycle but a
slower cycle**. It matches Beard 2022 in direction, with a larger force effect
(×1.36 vs the ~×1.18 previously quoted from uncorrected 03/27 data).

## Caveats

* **Baker is excluded from the force number.** Its recovery windows are 2–6×
  shorter than the 2026 protocol, so `A`/`Am` are truncated (Am/A = 0.70–0.83 vs
  0.94–1.05), and the truncation is worse at low ATP because `ktr` is slower there
  — which manufactures a spurious force *loss*. Its `ktr` is still usable, and it
  agrees.
* **04/10 remains the highest** (×1.54 vs ×1.24 / ×1.29) even after correction.
  Its 8 mM run is the weakest 8 mM recording in the whole set. Either that run was
  additionally compromised, or a first activation at 2 mM damages the fibre more
  than the correction accounts for. Both readings imply the ×1.36 mean is, if
  anything, an **over**-estimate; the two 8→2 preps alone give ×1.27.
* The correction assumes φ = 0.45 transfers between preps. It is calibrated on one
  bracket. **A second end-of-session repeat on another prep is the single most
  valuable next measurement.**

## For fitting

1. Target **force ×1.36** and **`ktr` ×0.55** as *ratios*, with a free per-prep
   force scale — never two absolute levels.
2. Use the **rundown-corrected** force ratio; use the **raw** `ktr` ratio.
3. Weight `ktr` above force: it is 8 % CV uncorrected versus 11 % CV after a
   model-dependent correction.
4. If a single-prep fit is wanted, use **03/27** — it is the only day with the
   bracket, so its correction is measured rather than assumed.
