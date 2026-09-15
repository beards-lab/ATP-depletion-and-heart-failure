# MechanismAblation — nothing in the model is worth adding, and four things are worth deleting

**Run:** 2026-09-11 → 2026-09-13, 39.6 h, 2 instances × 5 workers. 30 of 34 arms
completed (2 skipped as dead switches, 4 failed — §5).
**Base:** `params/floored_seed.m` (= `rigaware_a_opt` under the one-prep-width deadband), cost **0.9081**.
**Scope as asked:** one mechanism at a time, 8 mM only.

---

## 0a. What can be deleted — the explicit list

Everything here is either measured dead (changes the cost by exactly 0) or
measured free (inside the 0.110 threshold). Grouped by what deleting it means.

### A. Dead code — delete the switch and its branches. No cost, no re-fit.

| what | where | evidence |
|---|---|---|
| `UseSerialStiffness` | `getParams.m:92` default only | flipping it changed cost by **exactly 0.0000**; read nowhere in `dPUdT_CombinedTransitions`. **Also delete its line from CLAUDE.md's "Key Behavioral Switches"** |
| `UseP31Shift` | legacy `dPUdTCa.m` only | no reader in the active ODE |
| `UseTORNegShift` | legacy `dPUdTCa.m` only | no reader in the active ODE |
| `UseAtpOnUNR` | legacy `dPUdTCa.m` only | no reader in the active ODE |
| `UseMutualPairingAttachment` | commented out in the active ODE | `Driver.m:77` even says "just breaks the sim" |
| `UseKstiff3` | `assertActiveParams.m:72` only | no reader in the active ODE |
| `UseTitinModel` | output-only in `evaluateModel` | not a mechanism |

### B. Inert as configured — delete the flag, or give the parameter a value.

| what | why inert | tested? |
|---|---|---|
| `UseDynamicPassive` | on, but `c_titin_visc = 0` | preflight: exactly 0.0000 |
| `UseTargetZoneSaturation` | superseded by `UseGlobalOccupancySaturation`, mutually exclusive in code | not run (would double-count) |

### C. Live mechanisms that cost nothing — delete mechanism and parameters.

| mechanism | parameters freed | cost |
|---|---|---|
| **`UseA2Reattaching`** | `kA2re`, `kA2hop`, `sA2hop` | **+0.0001** |
| **`UseMaxwellTensionOnly`** | `x_M_slack` | +0.037, and `rsK80` 0.06 **better** |

### D. Parameters that are dead inside a live mechanism — remove from the pool.

| param | why | when it would become live |
|---|---|---|
| `sigma1` | only in `RSR2PT = kmsr·exp(F_SR/sigma1)·P_SR`; `kmsr = 0` | if `kmsr > 0` — **tested, R1: −0.005, inside noise** |
| `sigma_srd2` | only in `RPD2SRD = ksrd·exp(−F_SR/sigma_srd2)·PD`; `ksrd = 0` | if `ksrd > 0` — **tested, R2: −0.014, inside noise** |
| `LXBpivot` | overwritten by `SL0` at `getParams.m:605` | **never tested — see the caveat below** |

**On "have the inert parameters been tested for when they become sensitive?"**
Partly, and the answer matters:

- `sigma1` and `sigma_srd2` were tested *through their gates*. R1 set `kmsr = 10`
  and R2 set `ksrd = 1`, which makes each parameter live, and the arm then had
  it compulsory in the draw. Both landed inside the noise. So these are not
  merely inert at the incumbent — **the transitions they shape have nothing to
  contribute even when switched on**.
- But each was tested at **one** revived value, not scanned. A gate revived at
  `kmsr = 10` explores one slice; the calibration line search moved it over
  ×0.1–×5, so the slice is a decade wide, which is decent but not a scan.
- `LXBpivot` was **not** tested. It cannot be, as things stand: `getParams`
  resets it, so the optimiser's modifier is discarded. It should simply leave
  the pool, and the reset should be made explicit in `getParams` rather than
  being an accident of the init block.

---

## 0. The short version

Thirty-five mechanisms were flipped, each given a scalar calibration plus a
budget-matched re-optimisation with its own relevant parameters made compulsory.

1. **No mechanism, added or removed, improves the incumbent beyond the noise
   floor.** The three best results (`limited_srd` −0.036, `srd_fixeddwell`
   −0.027, `a2shift` −0.021) are all inside a measured phantom floor of 0.055.
2. **Four things can be deleted for free**, the largest being
   **`UseA2Reattaching` — 3 parameters for +0.0001**.
3. **Two switches are dead code**: `UseSerialStiffness` and — as configured —
   `UseDynamicPassive`. `UseSerialStiffness` is listed in CLAUDE.md as a key
   behavioural switch and is read nowhere in the ODE.
4. **Three parameters in the optimiser pool are inert** and have been drawn and
   "optimised" for the whole campaign history: `sigma1`, `sigma_srd2`,
   `LXBpivot`.
5. The two mechanisms with the strongest prior case — fixed-dwell SRD drain
   (the rsK amplitude-slope defect) and the viscoelastic SE (the re-stretch
   overshoot) — **both failed to deliver**. C1 landed inside the noise; C2's
   apparent 0.8488 did not reproduce.

---

## 1. The threshold, and why it is 0.110 and not zero

Both surviving control replicates re-optimised the incumbent to **0.9081
exactly**: at this budget a round buys the unmodified model *nothing*. Their
spread is therefore 0, which would make any positive Δ "significant" — wrong.

The real floor is the objective's irreproducibility. `optimizeFeatures`
re-evaluates its best point before writing a snapshot and refuses if it has
drifted; that guard fired twice in this campaign:

| | search found | re-evaluated | phantom |
|---|---|---|---|
| control 3 | 0.8561 | 0.9081 | 0.052 |
| C2 `se_viscoelastic` | 0.8488 | 0.9042 | 0.055 |

Two independent measurements, on this machine, on this objective, of a gain
that was not real. **σ = 0.055, threshold = 2σ = 0.110.** This is consistent
with the σ ≈ 0.076 in `Docs objective_noise_floor`, and it is worth more because
it was measured here.

That the guard caught both is the campaign's own evidence that the design was
right to demand controls.

---

## 2. The table

`calib` = after the single-parameter line search; `reopt` = after the
budget-matched round. Ranked by cost against the control.

| id | arm | raw | calib | reopt | vs ctrl | ΔN | verdict |
|---|---|---:|---:|---:|---:|---:|---|
| C5 | limited_srd | 0.897 | 0.881 | **0.872** | −0.036 | +1 | inside noise |
| C1 | srd_fixeddwell | 269 | 0.913 | **0.882** | −0.027 | +1 | inside noise |
| C11 | a2shift | 6.15 | 1.088 | **0.888** | −0.021 | +3 | inside noise |
| R2 | ksrd_revive | 0.913 | 0.896 | 0.894 | −0.014 | +1 | inside noise |
| C14 | a2popping | 0.953 | 0.908 | 0.902 | −0.006 | +1 | inside noise |
| R1 | kmsr_revive | 1.10 | 0.922 | 0.903 | −0.005 | +1 | inside noise |
| C10 | a2recocking | 0.969 | 0.911 | 0.907 | −0.001 | +2 | inside noise |
| **B3** | **a2reattach_off** | 0.999 | 0.999 | **0.908** | **+0.0001** | **−3** | **FREE — delete** |
| R4 | dA2re_revive | 0.963 | 0.917 | 0.911 | +0.003 | +1 | inside noise |
| A3 | occupancy_langmuir | 1.39 | 0.986 | 0.931 | +0.023 | 0 | neutral |
| **B8** | **maxwell_bidirectional** | 0.991 | 0.991 | 0.945 | +0.037 | **−1** | **FREE — delete** |
| A5 | halfactive_SR | 2.19 | 1.036 | 0.945 | +0.037 | 0 | neutral |
| B10 | overlapfactor_off | 1.30 | 1.304 | 0.987 | +0.079 | 0 | neutral |
| C13 | r2ceiling | 3.96 | 1.015 | 1.000 | +0.092 | +2 | not worth it |
| C6 | direct_srx | 1.48 | 1.434 | 1.002 | +0.094 | 0 | neutral |
| B9 | lattice_off | 1.13 | 1.132 | 1.030 | +0.122 | **−6** | see §4 |
| C12 | catchbond | 7204 | 1.199 | 1.059 | +0.151 | +1 | reject |
| B6 | negkstiff_off | 1.16 | 1.162 | 1.155 | +0.247 | −2 | keep |
| C3 | a1kernel_on | 29.6 | 5.002 | 1.304 | +0.396 | +1 | reject |
| B11 | overlap_off | 16.1 | 16.07 | 1.363 | +0.455 | −3 | keep |
| B2 | loadrealign_off | 1.96 | 1.871 | 1.638 | +0.730 | −3 | keep |
| B7 | maxwell_off | 37.1 | 2.695 | 1.880 | +0.972 | −2 | keep |
| A1 | occupancy_off | 64.9 | 6.553 | 2.129 | +1.22 | −1 | keep |
| C4 | passive_for_SR | 466 | 7.151 | 4.077 | +3.17 | 0 | reject |
| B5 | srx_off | 37.9 | 18.18 | 11.01 | +10.1 | −8 | keep |
| C7 | stretch_activation | 389 | 21.03 | 19.92 | +19.0 | +1 | reject |
| B13 | r1d_gaussian | 2738 | 334.1 | 83.27 | +82.4 | −8 | keep (see §4) |
| A4 | regavail_bidirectional | 840 | 839.6 | 131.2 | +130 | 0 | keep directional |
| C8 | velgauss_for_regavail | 27850 | 835.8 | 747.1 | +746 | 0 | reject |
| B1 | regavail_off | 46007 | 1106 | 794.3 | +793 | −3 | **most load-bearing** |

All 30 completed arms are listed. Not in the table: **A2** and **B12** (dead
switches, §3); **B4** `srxadp_off`, **C9** `vernier_on`, **C2**
`se_viscoelastic` and control 3 (failed, §5).

---

## 3. Dead code and inert parameters

Found by the preflight and by mapping every switch to its reader.

**`UseSerialStiffness` is never read.** Flipping it changed the cost by exactly
0.0000. It appears only in a comment in `dPUdT_CombinedTransitions`, a
commented-out `rmfield` in `evaluateModel`, and its default in `getParams`. The
serial elastic element is unconditionally on. CLAUDE.md lists it under "Key
Behavioral Switches"; that line is wrong.

**`UseDynamicPassive` is on and inert**, because `c_titin_visc = 0`. Either
give it a value or drop the flag.

**Six more switches are dead in the current ODE** and survive only in the legacy
`dPUdTCa.m`: `UseP31Shift`, `UseTORNegShift`, `UseAtpOnUNR`,
`UseMutualPairingAttachment`, `UseKstiff3`, `UseTitinModel`.

**Three pool parameters have no effect on the trajectory**, and have been drawn
and optimised throughout this project's campaign history:

| param | why it is dead |
|---|---|
| `sigma1` (31.1) | only enters `RSR2PT = kmsr·exp(F_SR/sigma1)·P_SR`, and **`kmsr = 0`** |
| `sigma_srd2` (37.6) | only enters `RPD2SRD = ksrd·exp(−F_SR/sigma_srd2)·PD`, and **`ksrd = 0`** |
| `LXBpivot` | reset to `SL0` at `getParams.m:605`, which `RunBakersExp.m:133` triggers |

Consequence: the incumbent's thick-filament topology is a **one-way cascade** —
`PT →(ksr0=367)→ SRX ⇄(ksr2srd=128 / ksrd2sr=9.1)⇄ SRD →(kmsrd=31)→ PD` — with
no direct force-dependent SRX→DRX return at all, and no PD→SRD route.

Both revivals were run and both are inside the noise: **R1** (`kmsr`, which also
revives `sigma1`) −0.005, **R2** (`ksrd`, which also revives `sigma_srd2`)
−0.014. So the cascade is not costing anything and the three dead parameters are
dead in the useful sense too — restoring the routes they gate buys nothing. They
should still be removed from `cfg.pool`, where they currently consume draws.

---

## 4. What to actually do

**Delete now, free:**
- `UseA2Reattaching` → `kA2re`, `kA2hop`, `sA2hop` (**3 params, +0.0001**).
  Note its χ² movers: `vall2_dy` −0.13 *better*, `coolDownLS` +0.03 worse. This
  contradicts the earlier `a2reattaching` finding (3.22 → 3.12); under the
  pooled, floored objective the lever buys nothing.
- `UseMaxwellTensionOnly` → `x_M_slack` (1 param, +0.037; `rsK80` −0.06 better).

**The parsimony question — B9 `lattice_off`: 6 parameters for +0.122.** It is
2σ over the threshold by 0.012, i.e. barely "load-bearing", while carrying
`d10_ref`, `SL_ref_lattice`, `R_thick`, `R_thin`, `d_optimal`, `sigma_lattice`.
Six structural parameters for one tenth of a preparation-width is the weakest
per-parameter case in the model. Recommend a full-length run of this arm before
deciding — it is the one arm where the screen's 40-evaluation budget is
plausibly the binding constraint rather than the mechanism.

**Same question, larger: B13 `r1d_gaussian` (+82 for −8 params).** The free-form
R1D strain law *is* earning its 8 parameters, decisively. That is a real answer
to a parsimony question that has been open in this project.

**Load-bearing, confirmed, in order:** registration availability (B1, +793 —
and its directional form matters: making it bidirectional costs +130), the R1D
strain law (B13), the SRX pool (B5, +10 for 8 params), the passive Maxwell
element (B7), occupancy saturation (A1), load realignment (B2).

**Reject (adding them costs):** stretch activation (+19), velocity-Gaussian
attachment as a substitute for availability (+746), F_passive as the SRX driver
(+3.2), the A1 attachment kernel (+0.40), catch bonds (+0.15).

---

## 5. Caveats, honestly

**The `coolDownLS` problem.** It is an unbounded raw MSE on the cool-down trace
(`extractSlackAttributes.m:456`), not a χ² term, and it dominates the "better"
column of nearly every arm's movers — B11 −12.5, A4 −706, B13 −250. Much of
what the re-optimisation bought was that one guard term recovering, not data
terms improving. Any arm whose gain is mostly `coolDownLS` should be read as
"the search escaped a penalty", not "the mechanism was compensated".

**Four runs failed.** `B4 srxadp_off` and control 3 and `C2 se_viscoelastic` hit
the snapshot-divergence guard (§1); `C9 vernier_on` hit the 150 s wall-clock
watchdog. B4 and C2 are worth re-running — both were tracking ~1.4 and ~0.85
respectively when the guard refused them.

**The budget is a screen, not a result.** 10 calibration + 40 simplex
evaluations per arm. Arms that landed near the threshold (B9, B6, C6, B10)
deserve a full-length run before any of them is believed.

**One-at-a-time under-values overlapping mechanisms.** B3/C10/C11 are three
readings of A2 re-binding and all three landed inside the noise; that is
consistent with them substituting for each other, which this design cannot
distinguish. Same for B1/C8/C9 (velocity-dependent attachment).

---

## 7. Parameter coupling — which parameters are the same parameter

Asked directly: are some pool entries bound together as `X = A·B` (fit one, drop
the other) or as `X = k·(A+B+C)` (where a separate `k` *helps* the search)? Both
exist here, and the lattice mechanism is the extreme case.

### 7.1 Lattice spacing: six parameters, two identifiable. Exactly.

```
d10       = d10_ref · sqrt(SL_ref_lattice / SL)
d_surface = d10 − R_thick − R_thin
excess    = max(0, d_surface − d_optimal)
f_lattice = exp(−excess² / (2·sigma_lattice²))
```

Substitute and the structure falls out:

```
excess    = C₁/√SL − C₂        C₁ = d10_ref·√SL_ref_lattice
                               C₂ = R_thick + R_thin + d_optimal
f_lattice = exp( −(1/√SL − a)² / (2b²) )     a = C₂/C₁ ,  b = sigma_lattice/C₁
```

- **`d10_ref` and `SL_ref_lattice` are an `A·B` pair** — only their product
  `d10_ref·√SL_ref_lattice` can be identified. Fitting both is fitting a ridge.
- **`R_thick`, `R_thin` and `d_optimal` are an `A+B+C` triple** — only their sum
  enters. The optimiser cannot tell a fatter thick filament from a wider optimal
  gap.
- Then `C₁` cancels into the two ratios, so the whole six-parameter mechanism has
  **two** degrees of freedom, and both are physically readable:

| combination | incumbent | meaning |
|---|---|---|
| `SL_thr = 1/a²` | **2.090 µm** | the SL below which attachment starts being penalised |
| `b` | **0.0471** | roll-off width, in units of 1/√SL |

That is the whole mechanism: *attachment is attenuated below SL ≈ 2.09 µm.* At
`SL0 = 2.2` it is inert (`f_lattice = 1`) and it bites only on shortening —
`f_lattice ≈ 0.95` at SL 2.0, `≈ 0.52` at SL 1.8.

**This is why B9 looked weak at 6 parameters for +0.122.** Four of the six were
never doing independent work. Judged properly the mechanism costs **two**
parameters for +0.122, which is an entirely ordinary price — and the four
redundant ones are exactly the ones that are *observations*: `d10_ref` and
`SL_ref_lattice` are an X-ray measurement and the length it was taken at;
`R_thick` and `R_thin` are EM filament radii. They should be constants in the
source, not entries in a fitting pool. `d_optimal` and `sigma_lattice` are the
two that are genuinely unmeasured — and they are exactly the two the reduction
keeps. **Recommendation: fix the four, fit the two. Do not delete the mechanism.**

### 7.2 Overlap: three parameters, one shape plus one gain that `ka` already owns

In the operating regime (`L_thick/2 < SL/2`, and `SL−LSE < 2·L_thin − L_hbare =
2.233 µm`, i.e. everywhere these protocols go):

```
N_overlap = [ (L_thick − 2·L_thin) + (SL−LSE) ] / (L_thick − L_hbare)
```

The numerator identifies only the combination `L_thick − 2·L_thin`. The
denominator is a **pure gain**, and attachment is `RD1 = ka·PD·N_overlap·f_lattice`
— so `(L_thick − L_hbare)` is exactly degenerate with `ka`.
**Recommendation: drop `L_hbare` from the pool.** The degeneracy is exact below
SL−LSE = 2.233 µm and breaks only above it, where the `max(...)` switches branch.

### 7.3 `ka`, `A0`, `P_bound_max`: degenerate at steady state, separated only by dynamics

`ka_eff = ka·A_reg·(1−A_def)`, then `RD1 ·= f_sat(P_bound)`. At isometric steady
state `A_reg = A0`, `A_def = 0`, so isometric attachment ∝ `ka·A0·f_sat`. Three
gains multiplying to one number. They are distinguishable only through their
*dynamics* — `v_ref_reg`/`tau_reg` give `A_reg` its velocity dependence, and
`P_bound_max` makes `f_sat` state-dependent feedback rather than a constant.

The campaign measured this degeneracy without meaning to: **A1 deleted `f_sat`
entirely (64.87) and a line search on `ka` alone recovered it to 6.55** — 90 % of
the damage undone by the partner gain. Keep all three (the dynamics are real),
but expect them to trade, and do not read any one of them alone.

### 7.4 The case where redundancy HELPS: `xrate`

`updateRates.m` scales every turnover rate by a single `xrate`. That is an
*exactly* redundant coordinate — `xrate·ka` cannot be separated from `ka` — and
it is the good kind. The sloppiest direction of any rate model is "all rates
faster or slower together"; without a common-mode handle the optimiser has to
discover it by coordinating eight parameters through a simplex, which a 30-
evaluation round will not do. One parameter spans it in one step.

**`xrate` is not in the pool.** It should be, with the caveat written down: its
*value* is meaningless in isolation (it can always be traded back into the
individual rates), so it is a search coordinate, not a result. This is the
`X = k·(A+B+C)` case the question anticipated, and it is the one place where
adding a redundant parameter is the right move.

### 7.5 Net effect on the pool

| action | parameters | why |
|---|---|---|
| remove, inert | `sigma1`, `sigma_srd2`, `LXBpivot` | §0a.D |
| remove, deleted mechanism | `kA2re`, `kA2hop`, `sA2hop`, `x_M_slack` | §0a.C |
| remove, redundant | `d10_ref`, `SL_ref_lattice`, `R_thick`, `R_thin`, `L_hbare` | §7.1, §7.2 |
| **add** | `xrate` | §7.4 |

Twelve out for one in: a 70-entry pool becomes 59, and with `N_DRAW = 12` the
chance of a given parameter being drawn in a round rises by ~19 % for free.

---

## 8. `coolDownLS` is retired

It is an unbounded raw MSE (`extractSlackAttributes.m:456`), it is not a
preparation-width chi-2 like every other data term, and §5 showed it dominating
the movers of nearly every arm — B11 −12.5, B13 −250, A4 −706 against data terms
moving by 0.01. A term that can swing by 700 in a 0.9 objective is not measuring
the fit; it is a guard rail that the search spends its budget climbing.

It is replaced by **`rsK80`**, which measures the same thing — the post-re-stretch
force recovery — as a floored relative chi-2 with a measured between-preparation
width, and which the RestretchMetricAudit already established as the sound
estimator (`rsK80 = ln4/(t80−t20)` with the amplitude taken from steady state).
`rsK80` takes double weight in the lean objective to inherit the role.

---

## 9. PROPOSAL (not implemented): where to record that a knob is abandoned

There is currently nowhere to write "this was measured dead". The knowledge
lives in `conclusions.md` files, which nothing reads at runtime, so the same
dead parameter gets drawn by the next campaign. `sigma1` has been optimised for
months while having no effect on the trajectory — that is the failure this is
meant to stop.

**One registry, three thin hooks, one test.** Nothing here changes model
behaviour; it changes what the tooling will let you spend effort on.

**1. `Model/retiredKnobs.m`** — the single source of truth, a function returning
one entry per knob:

```matlab
R(end+1) = struct( ...
  'name',     'sigma1', ...
  'kind',     'inert', ...          % dead | inert | ineffective | redundant
  'gate',     'kmsr', ...           % what would make it live ([] if nothing can)
  'evidence', 'MechanismAblation R1: revived at kmsr=10, -0.005, inside noise', ...
  'measured', '2026-09-13', ...
  'action',   'exclude from optimiser pools');
```

Four kinds, because they need different handling: **dead** (no reader — delete
the source), **inert** (a live mechanism whose parameter is gated off),
**ineffective** (live and tested, buys nothing), **redundant** (identifiable only
in combination — §7).

**2. Hook it in three places, each one line of intent:**

- `getParams` — after defaults are filled, warn once if a retired switch is
  enabled or a retired parameter is non-zero: *"UseSerialStiffness was measured
  dead on 2026-09-13 (MechanismAblation B12); setting it has no effect."* A
  warning, never an error: the registry must not be able to block a deliberate
  re-test.
- **Pool construction** — `optimizeFeatures` filters `cfg.pool` through the
  registry and prints what it removed. This is the hook that actually saves
  compute, and it belongs in `optimizeFeatures` rather than in each campaign
  script so no future campaign can forget it.
- `assertActiveParams` — already a manual mirror of the ODE's switch logic;
  have it skip retired entries so physiology tiles stop asserting bounds on
  parameters that cannot move anything.

**3. `Workbench/Tests/TestRetiredKnobs.m`** — the part that matters most,
because a deprecation list that is never re-checked becomes folklore. For each
`kind`, a cheap assertion:

| kind | test |
|---|---|
| dead | flip the switch, assert Δcost == 0 exactly |
| inert | assert the gate parameter is still 0 |
| redundant | assert the algebraic identity still holds in the source (e.g. `L_hbare` appears only in the denominator) |
| ineffective | no runtime test — it is a measurement, carries its date, and goes stale honestly |

One evaluation per dead entry, ~8 evaluations total. Cheap enough for CI.

**4. Two documentation fixes that should ship with it:**
`CLAUDE.md` lists `UseSerialStiffness` under "Key Behavioral Switches" — it is
dead code and that line should go. `Docs/README.md` should point at the registry
so the knowledge map has an entry for "what we stopped fitting, and why".

**What I would not do:** put this in `parameterBounds.m`. That file is about
physiological priors on parameters that *are* being fitted; mixing "this value
is implausible" with "this parameter does nothing" would blur two different
claims. Keep them separate and let both be read.

---

## 6. What it cost, and the two traps that cost it

39.5 h, after two aborted attempts. Both aborts were the same misreading:

**`MaxRunTime` is a wall-clock watchdog on ONE SIMULATION, not one evaluation**
(`evaluateModel.m:102`), and one evaluation runs ~18 simulations (8 F-V points,
5 slack, 5 passive). At 300 s an arm had no effective bound — one ran 5.7 h.
`TIME_BUDGET_HRS` could not help because it is only checked *between* rounds and
the design used a single round.

**It also throws** (`evaluateModel.m:253`), and being wall-clock it trips on
machine contention rather than on a stiff model. Cutting it to 60 s while adding
a third MATLAB instance killed every arm including the unmodified control —
60 s evaluates the baseline fine on an idle box (0.9081 in 97 s). The fix was
150 s, two rounds so the wall-clock cap can fire, two instances rather than
three, and try/catch on every evaluation made outside `optimizeFeatures`.
