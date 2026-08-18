# Transferring the low-ATP fit onto the realign4 8 mM baseline

**Question.** The 8 mM baseline moved on `Further-optim` (`realign4_opt`, the
load-realignment refit). The low-ATP fit was built on the *previous* baseline
(`rskR2_w025_opt`). Re-base it — new baseline, same ATP gates, no refit — and ask how
good the 2 mM fit is now.

**Verdict.** The transfer works and it is a large, real improvement: **2 mM feature cost
19.58 → 12.44** on one fixed yardstick, against a scatter floor of ~0.25. Nothing was
refitted at low ATP. The gain is **not** the gates doing something new — it is that
realign4 fixed a defect the two conditions *share*, and the 2 mM condition inherited it.
The specifically-ATP part of the residual is essentially untouched: the model's
`rsK` ratio is 0.588 against a measured 0.363, versus 0.600 before.

---

## 1. What was transferred

The low-ATP specifics are a cleanly separable block — everything else in
`lowATP_powergate2_opt.m` *is* the old 8 mM baseline. So the re-base is exactly:
`realign4_opt.m` verbatim, plus

- `UseAtpPowerGates = 1`, `MgATP_ref = 8`
- the 21 `pAtp_*` exponents
- the 2 mM condition itself (`MgATP = 2`, `Pi = 0.8`, 2 mM slack protocol)

written to **`params/lowATP_realign4_gates.m`**. The gates are exactly 1 at `MgATP_ref`
for any exponents, so **the 8 mM fit above them is untouched by construction**.

Exponents were taken from the optimiser state at **full precision**, not from the
committed `%.6g` snapshot — that rounding alone is worth ~0.60 of the 2 mM objective
([AtpGateHillForm](../AtpGateHillForm/conclusions.md) §3a).

Merging `Further-optim` was a prerequisite, not a convenience: `realign4_opt` sets
`UseLoadRealign`, `SRXFromR2HighStrain`, `UseMaxwellTensionOnly` and the
`FV_normpowerAvg` feature, none of which existed on `lowatp-2state`.

## 2. The numbers

**The two baselines were fitted against different objectives**, so the incumbent's
recorded 15.98 is *not* comparable to anything below. `realign4`'s list drops
`ovrsht_dy` (retired by [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md)),
reweights `rsK` 0.1 → 0.25, and rescales the passive terms. Everything here is
therefore scored **both ways**, and the like-for-like column is the one to read.

| configuration | OLD list (19 feats) | **NEW list (18 feats)** |
|---|---|---|
| OLD baseline + gates @ 2 mM | 15.985 | **19.575** |
| **NEW baseline + gates @ 2 mM** | 15.687 | **12.440** |
| NEW baseline, **gates OFF** @ 2 mM | 133.49 | **138.30** |
| OLD baseline @ 8 mM (floor) | 4.685 | 5.669 |
| NEW baseline @ 8 mM (floor) | 6.738 | **3.112** |

Three readings, in order of importance:

- **On the current objective the 2 mM fit improves 19.58 → 12.44 (−36 %)**, ~29× the
  measured peak-to-peak scatter of 0.25. The 8 mM floor improves 5.67 → 3.11, so the
  **excess over floor** goes 13.91 → 9.33 (−33 %). Roughly a third of the win is the
  shared baseline getting better and two thirds is better low-ATP behaviour — though
  that split is a framing, not a decomposition, since the features interact.
- **On the *old* objective the same change is worth only −0.30** — about one
  peak-to-peak, i.e. nothing. The entire improvement lives in `rsK`, which the old list
  weighted at 0.1 and which realign4 was built to fix. The two lists disagree about
  whether this transfer is a triumph or a wash; the new one is right, because `rsK` is
  where the model is actually wrong.
- **The transferred gates still carry the whole low-ATP effect on the new baseline**:
  turning them off costs 12.44 → 138.30. They did not become inert or harmful under
  the new parameters, which is the thing that could have gone wrong and did not.

> **Constant offset, read before quoting any absolute number.** The bound term
> `k2|…,ksrd|…,kmsrd|…` contributes **exactly 1.200 to every cell above**. A slack-only
> run does not produce `k2`/`ksrd`/`kmsrd`, so each sub-entry returns
> `MISSING_FEATURE_COST = 100` × its weight (0.01/0.001/0.001). It is not a bound
> violation and not a residual. It cancels in every comparison here, but the honest
> floors are 1.91 (new) and 4.47 (old), not 3.11 and 5.67.

## 3. Where the improvement is, feature by feature

2 mM, NEW list, old baseline → new baseline:

| feature | OLD 2 mM | **NEW 2 mM** | Δ |
|---|---|---|---|
| **`rsK`** | 10.595 | **2.756** | **−7.84** |
| `coolDownLS` | 2.020 | 1.616 | −0.40 |
| `steady` | 0.476 | 0.171 | −0.31 |
| `vall_y` | 0.978 | 0.798 | −0.18 |
| `peak2` | 0.330 | 0.240 | −0.09 |
| `vall2_dy` | 1.204 | 1.801 | +0.60 |
| `ktr` | 0.273 | 0.725 | +0.45 |
| `A` | 0.195 | 0.610 | +0.42 |
| `ktr_rmse` | 0.074 | 0.191 | +0.12 |
| (bound-term artifact) | 1.200 | 1.200 | 0 |
| **total** | **19.575** | **12.440** | **−7.14** |

`rsK` alone is **110 %** of the net improvement; six terms get modestly worse and pay
some of it back. `rsK` is still the largest single residual at 2 mM (2.76, 22 % of the
total), so it remains the thing to work on — it is just no longer catastrophic.

## 4. The physiology: is the *ATP effect* any better?

Mostly no, and this is the honest limit of the result. Ratios are 2 mM / 8 mM:

| observable | data | OLD model | **NEW model** | new model/data |
|---|---|---|---|---|
| `steady` (force) | 1.180 | 1.153 | **1.163** | 0.99 |
| `A` (amplitude) | 1.184 | 1.176 | **1.193** | 1.01 |
| `peak1_y` | 1.254 | 1.213 | **1.233** | 0.98 |
| `ktr` | 0.537 | 0.490 | **0.469** | 0.87 |
| `vall_y` | 1.198 | 1.085 | **1.072** | 0.89 |
| `peak2` | 1.179 | 1.064 | **1.049** | 0.89 |
| `rsK` | 0.363 | 0.600 | **0.588** | **1.62** |

- **Force amplitudes are excellent and were already**: force, `A` and `peak1_y` all land
  within 1–2 % of the measured ATP effect. The 2 mM force enhancement is a solved
  problem in this model.
- **`rsK`'s *ratio* barely moved** (0.600 → 0.588, still 1.6× too shallow) — but its
  *absolute* value improved a lot at both concentrations: 8 mM 102.6 → 67.1 (data 43.7,
  so ×2.35 → ×1.53) and 2 mM 61.6 → 39.5 (data 15.9, ×3.88 → ×2.48). That is exactly
  the shape of a defect **shared** with 8 mM, which is what
  [AtpGateHillForm](../AtpGateHillForm/conclusions.md) §5.4 predicted would be immune to
  anything done to the gates. Fixing it at 8 mM fixed it at 2 mM. The *residual* ATP
  asymmetry in `rsK` is a separate, still-open defect.
- **`ktr`, `vall_y` and `peak2` drifted slightly the wrong way** (m/d 0.87–0.89): the
  model now over-slows `ktr` at low ATP and under-produces the low-ATP enhancement of
  the valley and second peak. These are the terms that got worse in §3 and they are the
  natural targets for a refit.

## 5. What this does and does not license

- **Does:** adopt `params/lowATP_realign4_gates.m` as the low-ATP snapshot. It is
  strictly better than the incumbent on the current objective and the gates transfer
  intact.
- **Does not:** claim the low-ATP *mechanism* improved. No low-ATP parameter was
  refitted, and the ATP ratios are the same or marginally worse except through `rsK`'s
  shared component.
- **Next, and cheap:** the gates have never been refitted on this baseline. `rsK`,
  `vall2_dy`, `ktr` and `A` moved; a re-run of the power-gate optimiser from these
  exponents should recover most of the +1.6 that the six worsened terms cost — with
  accept thresholds at ~0.25, not `1e-4` (§3a of AtpGateHillForm).

## 6. Files

| file | what it does |
|---|---|
| `BuildCombinedSnapshot.m` | writes `params/lowATP_realign4_gates.m` = realign4 + the gate block, exponents at `%.17g` from the optimiser state |
| `RunTransferCheck.m` | the five slack simulations, scored under both feature lists; writes `results/realign4_transfer.mat` |

    matlab -batch "run('Analyses/LowATP_Realign4Transfer/BuildCombinedSnapshot.m')"
    matlab -batch "run('Analyses/LowATP_Realign4Transfer/RunTransferCheck.m')"
