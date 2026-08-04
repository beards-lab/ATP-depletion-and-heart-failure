# Low-ATP fit on an ABSOLUTE target (2026-08-04)

Supersedes the ratio-scored low-ATP fits (`LowATP_kmsrd`, `LowATP_2State`). Baseline is
`params/rskR2_w025_opt.m` (the rsK/rsR2-weighted high-ATP snapshot).
Driver: `RunOptimLowATPAbs.m`. Baseline diagnostic: `../../Workbench/CheckJointAtpBaseline.m`.

---

# CAN THE OTHER MECHANISMS PUSH rsK THE RIGHT WAY? — No. (2026-08-04)

With refitting 8 mM now permitted, every remaining mechanism was scanned at **both** ATP levels.
`../ApoPoolDetachment` had scored these at 8 mM only; the ATP **ratio** was never asked.
Reproduce: `RunMechanismRatioScan.m`. Baseline = `rskR2_w025_opt` + gates (K_T2=1.8, K_srx=9).

Need: `rsK8` 114.2 → **43.7**, ratio 0.56 → **0.36**, holding `peak1_8` ≈ 89–94 and `steady8` ≈ 77.

| config | rsK8 | peak1_8 | steady8 | rsK2 | **ratio** | verdict |
|---|---|---|---|---|---|---|
| baseline (gates only) | 114.2 | 93.7 | 76.7 | 63.8 | 0.56 | — |
| `eta_M`×0.3 | **68.2** | **93.6** | **76.7** | 45.2 | 0.66 ↑ | level ✅ / ratio ✗ |
| `k2`×0.6 | 62.0 | 121.8 ✗ | 93.2 ✗ | 40.2 | 0.65 ↑ | level ✅ / peak+force ✗ |
| `k2`×0.6 + `eta_M`×0.3 | **46.2** | 121.7 ✗ | 93.1 ✗ | 32.2 | 0.70 ↑ | level ✅ / ratio ✗✗ |
| ceiling `k2max`=1500 | 95.9 | 101.6 | 78.1 | 67.1 | 0.70 ↑ | weak / ratio ✗ |
| ceiling `k2max`=800 | 106.8 | 109.4 | 79.3 | 67.1 | 0.63 ↑ | weak / ratio ✗ |
| ceiling `k2max`=1500, `K_T2c`=6 | 120.8 ↑ | 108.1 | 79.1 | 61.0 | **0.50 ↓** | ratio ✅ / level ✗ |
| **`K_T2`=6 (gate pushed hard)** | 114.2 | 93.7 | 76.7 | 46.0 | **0.40 ↓** | ratio ✅ / level ✗, force ✗ |
| `eta_M`×0.3 + `K_T2`=6 | 68.2 | 93.6 | 76.7 | 36.0 | 0.53 ↑ | **the weld, in one row** |
| `ceil400+D0+rip3000` | NaN | NaN | 26.5 ✗ | NaN | — | **model destroyed** |
| `ceil300+D0+rip10000` | NaN | NaN | 14.8 ✗ | NaN | — | **model destroyed** |
| *data* | *43.7* | *89.3* | *77.3* | *15.9* | ***0.36*** | |

## The result: a LEVEL ↔ RATIO weld

**Every lever that fixes the 8 mM `rsK` level makes the ATP ratio worse; the only levers that
improve the ratio do not fix the level.** 15 configurations, no exceptions.

The clearest single row is `eta_M`×0.3 + `K_T2`=6: `K_T2`=6 alone gives ratio **0.40** (nearly the
data's 0.36); adding the level fix pushes it straight back to **0.53**.

**Why:** `rsK` saturates in the detachment rate. At the baseline the model sits in the *steep*
region, where the multiplicative ATP gate produces a large fractional change in `rsK`. Fixing the
8 mM level (`eta_M`/`k2`) moves the operating point into the flat region, so the same gate buys a
smaller fractional change → ratio → 1. **The model's ATP effect on `rsK` is partly an artefact of
sitting in the steep region, and it is spent as soon as the level is corrected.**

## Secondary findings

- **`eta_M` is the clean level lever**: `rsK8` 114.2 → 68.2 with `peak1_8` (93.6) and `steady8`
  (76.7) **untouched**. The worry that refitting 8 mM must make things worse is unfounded for
  `eta_M` — it is nearly free. `k2` is the dirty one (`peak1_8` +28, `steady8` +16 uncompensated).
- **The D0 / tearing pool is dead on this baseline** — `steady` collapses to 15–27 and the cost
  explodes to 8×10⁴–1.4×10⁵. It was viable only on the older TNF + `UseMaxwellTensionOnly` +
  `tau_reg`=2 ms baseline. Not a candidate here at any setting tried.
- The ceiling remains weak in *both* directions even before the ktr-overstretch guardrail
  (`../ApoPoolDetachment`) is applied — which would only make it worse.

## What the data actually demands

`rsK` (×0.36) is suppressed by low ATP **more** than `ktr` (×0.55) and far more than force (×1.36).
Any lever acting on the shared cycle rate moves `rsK` and `ktr` together, so it cannot produce that
ordering. The missing mechanism must be **specific to the post-lengthening window** and
**ATP-gated** — i.e. it must slow re-attachment after a stretch, at low ATP, without touching the
isometric cycle.

That is exactly the shape of the one candidate never built: the **lengthening-driven registration
availability penalty** (`RegAvailShorteningOnly = 1` today, so lengthening does not modulate
availability, and the post-restretch window is precisely a post-lengthening window) — independently
recommended by [`../ApoPoolDetachment/conclusions.md`](../ApoPoolDetachment/conclusions.md) §3 on
entirely different grounds. Built `MgATP_ref`-normalised it would also be 8 mM-inert.

**Recommendation:** do not spend an 8 mM refit on `rsK` with the current mechanism set — the weld
above says the ratio will be paid for the level. Build the lengthening-availability lever first,
then refit.

---

# CURRENT ANSWER — 2 mM fit on the FROZEN 8 mM baseline (2026-08-04)

**Deliverable: `params/lowATP_only_opt.m` — `K_T2 = 1.80`, `K_srx = 9.00`.**
No shared parameter refitted; verified `max |shared param delta| = 0` against `rskR2_w025_opt`.
Driver: `RunFitLowATPOnly.m`.

## Observe — which mechanisms are even admissible

`dPUdT_CombinedTransitions.m` has exactly **four** `[MgATP]`-dependent paths. Only two are
`MgATP_ref`-normalised (`gATP = 1` at 8 mM for **any** Km), i.e. provably inert at 8 mM:

| path | line | normalised? | admissible with 8 mM frozen |
|---|---|---|---|
| `K_T2` (`UseAtpDetachR2`) — R2 detachment | 482 | **yes** | ✅ force amplitude |
| `K_srx` (`UseAtpKmsrd`) — SRX-ADP return | 539 | **yes** | ✅ ktr |
| `k2max`/`K_T2c` (`UseR2Ceiling`) | 472 | no | ❌ would move 8 mM (also falsified) |
| `k_D0T`/`K_D0T` (`UseD0State`) | 504 | no | ❌ would move 8 mM (also contested) |

So the justifiable lever set under this constraint is **exactly `{K_T2, K_srx}`** — a result read
off the code, not assumed.

## Orient — what those two levers can and cannot reach

The 2 mM `rsK` error **factorises**:

    x3.77 observed  =  x2.61 inherited 8 mM defect  x  x1.47 low-ATP specific

With 8 mM frozen only the ×1.47 is reachable. Confirmed by sweep: the best `rsK` anywhere in the
lever plane is 49.0 (vs data 15.9) and that corner wrecks everything else. **`rsK` is structurally
unreachable here** — it is an 8 mM baseline defect (model 114.2 vs data 43.7), not a low-ATP one.

## Decide — score three objectives from the same runs

`absFull` (all features), `absReach` (post-restretch recovery cluster removed: `rsK`, `rsR2`,
`coolDownLS`, `ovrsht_dy`, `vall2_dy`), and the legacy `ratio`. 5×5 lever sweep + local refine;
the 8 mM leg is computed **once** (the gates are inert there by construction).

| objective | optimum | K_T2 | K_srx |
|---|---|---|---|
| `absReach` | **5.033** (8 mM floor 2.723) | **1.80** | **9.0** |
| `ratio` (legacy relative) | 0.142 | 1.97 | 9.0 |
| `absFull` | 26.03 | 1.40 | 13.0 |

**`absReach` and the legacy `ratio` agree** (K_srx 9.0; K_T2 1.8 vs 1.97, adjacent grid points) —
two independent metrics selecting the same levers is the strongest validation available here.

**`absFull` disagrees, and that is the warning:** it buys an unreachable `rsK` by degrading the
force fit — `steady` 87.5 against data 91.3, versus 90.0–91.1 at the `absReach` optimum. Scoring a
defect the levers cannot reach actively corrupts the features they can. **Do not fit 2 mM against
`rsK` while 8 mM is frozen.**

## Fit — where it lands

| feature | model (2 mM) | data (2 mM) |
|---|---|---|
| `steady` | 92.1 92.9 93.4 93.9 77.7 | 95.8 95.7 94.8 94.2 75.7 |
| `ktr` | 30.2 27.9 26.0 24.0 23.0 | 32.3 28.9 26.0 23.0 21.9 |
| `peak1_y` | 122.2 122.9 120.9 119.0 93.5 | 121.3 121.3 114.6 109.6 92.8 |
| `rsK` | 64.2 62.5 64.1 70.6 57.4 | 14.8 15.2 15.5 15.0 18.9 ⟵ unreachable |

`ktr` tracks the per-segment slope essentially exactly; `steady` and `peak1_y` are close, both
slightly under/over on the middle segments. Everything except the post-restretch recovery is fit.

## Consequence

The post-restretch recovery **cannot** be fixed from the low-ATP side. It is an 8 mM defect that
low ATP amplifies, so it needs either a refit of 8 mM (the joint driver below) or a new mechanism —
and the only untested candidate that would not disturb either measured peak is the one
[`../ApoPoolDetachment/conclusions.md`](../ApoPoolDetachment/conclusions.md) recommends: a
**lengthening-driven registration-availability penalty** (`RegAvailShorteningOnly = 1` today, so
lengthening does not modulate availability — and the post-restretch window is exactly a
post-lengthening window). If built `MgATP_ref`-normalised it would also be 8 mM-inert and thus
admissible under this constraint.

---

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
