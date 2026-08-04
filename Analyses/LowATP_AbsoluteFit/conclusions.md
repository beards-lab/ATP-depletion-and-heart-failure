# Low-ATP fit on an ABSOLUTE target (2026-08-04)

Supersedes the ratio-scored low-ATP fits (`LowATP_kmsrd`, `LowATP_2State`). Baseline is
`params/rskR2_w025_opt.m` (the rsK/rsR2-weighted high-ATP snapshot).
Driver: `RunOptimLowATPAbs.m`. Baseline diagnostic: `../../Workbench/CheckJointAtpBaseline.m`.

## Why absolute, and why JOINT

The ratio cost (2 mM / 8 mM) cancels the high-ATP baseline's own error, so it can score the
*change* the ATP levers produce but can never make the model reproduce the 2 mM data itself.
Scoring 2 mM absolutely alone is worse: most of the 2 mM residual is only reachable through
**shared** parameters, and moving those against a 2 mM-only cost would silently wreck 8 mM.

So the objective is the weighted sum of two ABSOLUTE evaluations:

    Et = w8 * E(8 mM, absolute) + w2 * E(2 mM, absolute)

Wired as `params0.UseAtpJointFit` in `Model/evaluateBakersExp.m` — a guarded recursive call
(inert when false), so `optimizeFeatures`, the per-feature decomposition and the iteration
gallery are reused unchanged. Per-feature labels come back prefixed `hi:` / `lo:`.

The two ATP gates are `MgATP_ref`-normalised, hence **exactly inert at 8 mM for any Km**:
`K_T2` and `K_srx` are free low-ATP levers that cost the 8 mM leg nothing. That asymmetry is
what makes an absolute low-ATP fit tractable.

## Where the low-ATP residual actually is

Joint baseline = **32.57** (hi 7.61 + lo 24.96), at `K_T2 = 1.97`, `K_srx = 13`:

| term | cost | | term | cost |
|---|---|---|---|---|
| **lo:rsK** | **9.72** | | hi:rsK | 3.39 |
| **lo:ovrsht_dy** | **5.81** | | hi:FV_fnorm | 0.52 |
| **lo:coolDownLS** | **3.80** | | hi:coolDownLS | 0.48 |
| lo:vall2_dy | 1.63 | | lo:ktr | 0.43 |
| lo:vall_t | 1.06 | | lo:A | 0.29 |
| lo:vall_y | 0.68 | | lo:peak1_y | 0.15 |
| | | | lo:steady | 0.11 |

**Low-ATP terms take 5 of the top 7 slots, and the low leg is 3.3× the high leg.**

Two clean conclusions:

1. **The amplitude/kinetic half of the low-ATP problem is SOLVED** by the two existing gates.
   Absolute steady force is essentially exact at both conditions — model 76.7 / 90.8 against
   data 77.3 / 91.3 — and `A`, `peak1_y`, `ktr` are all ≤ 0.43.
2. **~88 % of what is left is the post-restretch recovery cluster**
   (`rsK` + `ovrsht_dy` + `coolDownLS` + `vall2_dy` + `vall_t` = 22.0 of 24.96).

## rsK is the whole story, and it is a baseline defect amplified by low ATP

| | model | data | model/data |
|---|---|---|---|
| rsK @ 8 mM | 114.2 | 43.7 | **×2.61 too fast** |
| rsK @ 2 mM | 60.0 | 15.9 | **×3.77 too fast** |
| **rsK ratio (2/8)** | **0.53** | **0.36** | |

- The post-restretch recovery is too fast at **both** conditions — so most of the fix must come
  from **shared** parameters, which is exactly what the joint cost prices correctly.
- The ATP gates already deliver a genuine slowing (ratio 0.53) but under-shoot the data's 0.36.
- `rsK`'s ATP effect (**×0.36**) is the **largest ATP effect in the dataset** — larger than
  `ktr` (×0.55) or force (×1.36). It is the sharpest low-ATP observable available.

This is the same defect measured independently in
[`../LowATP_RestretchRecovery/conclusions.md`](../LowATP_RestretchRecovery/conclusions.md)
by exponential fits to the hold window (post-restretch τ 71.8 ms vs data-implied ≈ 63 ms;
ratio 0.39 there vs 0.36 from `rsK` — two independent routes agreeing).

## Lever structure the driver is built around

From [`../RestretchRecoveryFit/conclusions.md`](../RestretchRecoveryFit/conclusions.md) §6:

- **`k2` sets the rsK level**; **`eta_M` removes its amplitude dependence**
  (`k2`×0.6 + `eta_M`×0.3: rsK 137 → 48.7, slope +714 → +21, data ≈ 52 / −126).
- `K_T2` (`UseAtpDetachR2`) is an **[ATP] gate on `k2`** → a low-ATP-only rsK lever.
- `K_srx` gates `kmsrd` → the low-ATP-only ktr lever.

Hence `cfg.compulsory = {k2, kstiff2, ka, eta_M, K_T2, K_srx}` — every round can trade both
the shared rsK levers and the low-ATP-only gates. Pool = 78 params.

## Deliberately EXCLUDED by default

`UseR2Ceiling` and `UseD0State` are off (`ALLOW_CONTESTED = false`). They are falsified at high
ATP by the ktr overstretch peak — see [`../ApoPoolDetachment/conclusions.md`](../ApoPoolDetachment/conclusions.md).
Opting in adds `k2max, K_T2c, k_D0T, K_D0T, k2rip`; the 8 mM leg should then punish them, which
is the honest way to re-open the question.

## How to run

    matlab -batch "run('Analyses/LowATP_AbsoluteFit/RunOptimLowATPAbs.m')"

Validate the setup without spending the budget:

    DRYRUN_SETUP = true; run('Analyses/LowATP_AbsoluteFit/RunOptimLowATPAbs.m')

Knobs: `W2` / `W8` (leg weights), `TAGSUF`, `TIME_HRS` (default 18), `NWORKERS`,
`ALLOW_CONTESTED`. Resumes automatically from `params/lowatpabs_<TAGSUF>_state.mat`.
A joint evaluation is ~2× a single-condition one — budget accordingly.

## Caveats

1. **Rundown.** [`../ATPEffectReconciliation/conclusions.md`](../ATPEffectReconciliation/conclusions.md)
   recommends fitting *ratios* with a free per-prep force scale and warns against "two absolute
   levels", because 03/27's 2 mM run was recorded second and is ~3.8 % degraded. This driver
   fits absolute levels by request. 3.8 % on force is small next to a ×3.8 `rsK` error, so it
   does not change where the residual is — but any *force* conclusion drawn from the absolute
   fit carries that bias, and 03/27 is the right day to use (it is the only one with a measured
   bracket). A rundown-corrected variant would scale the 2 mM force features by ≈1.04.
2. Only 03/27 is fitted. 04/03 and 04/10 are now curated and present in `data/` — adding them as
   further legs is the natural extension (04/10 also has 2 mM ktr and FV2).
3. The snapshot's `fn` carried a duplicated `rsR2|0.4-1.0|0.5`, double-counting that guardrail;
   the driver de-duplicates. Worth fixing in `RunOptimRsK.m`, which de-dupes `rsK` but not `rsR2`.

## Infrastructure note

This worktree's `data/` was **stale** (gitignored, so the merge does not carry it): the 03/27
files predated `AddRestretchFeatsToData.m`, so `rsK`/`rsR2` were absent and silently scored as
the missing-feature sentinel (100 × 0.25 = **25.0 of a 29.2 "cost"**). Refreshed from the main
worktree. **Always re-sync `data/` into a worktree after pulling analysis work.**
