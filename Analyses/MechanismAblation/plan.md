# MechanismAblation — what does every switch actually buy?

**Status:** design only, nothing run. 2026-09-11.

The question: take the incumbent parametrisation as final, flip every implemented
mechanism (on→off, off→on), give each flip a *budget-matched* minor re-optimisation
so a mechanism that can be compensated parametrically is allowed to be compensated,
and report what each one is worth.

This file is the inventory: which switches are worth an arm, which are not, and for
each arm, **which parameters must be free** for the answer to mean anything.

---

## 0. The base point

| | |
|---|---|
| Snapshot | `params/floored_seed.m` (= `rigaware_a_opt` re-weighted with the one-prep-width deadband) |
| Objective | floored, pooled-8 mM target, `OptimizeOn='Feats'`, 25 `fn` entries |
| Cost here | **0.9081** (`campaignSpec('floored').seedExpect`) |
| Engine | `Analyses/FeatureWeighting/RunOptimSpreadW.m` via `optimizeFeatures` |
| Eval cost | ~38 s at 5 workers (one objective evaluation = FV + slack×5 + passive) |

Not `estiff125_opt`, not `spreadw_opt`, not `pooled8_opt` — all superseded; see
`Analyses/RecoveryOscillation/conclusions.md` §7 and §9.

---

## 1. Five traps that would silently invalidate the scan

These are properties of the existing engine, not hypotheticals. Each one turns an
arm into a null result that *looks* like a real answer.

**T1 — zero-valued parameters are dropped from the pool.**
`RunOptimSpreadW.m:117-129` filters `cfg.pool` with `params0.(f) ~= 0`. Every
currently-off mechanism seeds its parameters at 0 (`k_catch_bond=0`, `k_SA=0`,
`c_SE_visc=0`, `c_titin_visc=0`, `k2d=0`, `drmr=0`, `Srd_max=0`, `pop_s=-Inf`,
`SRDOutflowK=0`, `D0FromR2=0`). Flip the switch on and the mechanism is **on but
numerically inert, and undrawable** — the arm measures nothing and reports "no
effect". Every ON-side arm needs a hand-chosen non-zero seed value.

**T2 — `compulsory` is intersected with the pool.**
`optimizeFeatures.m:183`: `compulsory = compulsory(ismember(compulsory, pool))`.
A compulsory entry that T1 has already dropped disappears without a warning.

**T3 — an undrawn parameter is untested.**
`N_DRAW=10`, 6 of which are compulsory, drawn from a ~64-entry pool: a given
non-compulsory parameter is touched in roughly **1 round in 16**. In a one-round
arm, anything not in `cfg.compulsory` is frozen. Each arm's own mechanism
parameters *and* its named compensators must be compulsory for that arm.

**T4 — the deadband hides regressions.**
The floored objective pays nothing inside one preparation-width. A mechanism whose
only job is defending an already-solved term (F-V at χ²≈0.07, rsK80 at 0.42, ktr at
0.97) can be deleted for **zero scalar cost** while the fit visibly degrades. Every
arm must report the **raw χ² vector**, not only the floored scalar.

**T5 — the baseline is at a local minimum and the flips are not.**
Ranking a re-optimised flip against the *un-re-optimised* incumbent charges every
flip for the re-start. This is the `StairsPrediction` lesson (the incumbent was
where the search *stalled*, not a minimum). The control arm must be the incumbent
re-optimised under identical budget and identical draw protocol, repeated ≥3× with
different RNG seeds to give σ_reopt — the only honest significance threshold.

---

## 2. Three parameters the incumbent is already carrying dead

Found while mapping the switches. Independent of the scan, but they change what a
draw is worth:

| Param | In pool? | Why it is dead |
|---|---|---|
| `sigma1` (=31.1) | yes | only enters `RSR2PT = kmsr*exp(F_SR/sigma1)*P_SR`, and **`kmsr = 0`** |
| `sigma_srd2` (=37.6) | yes | only enters `RPD2SRD = ksrd*exp(-F_SR/sigma_srd2)*PD`, and **`ksrd = 0`** |
| `LXBpivot` | yes | `getParams.m:605` resets it to `SL0` in the init block, which `RunBakersExp.m:133` triggers (`updateInit=true`) |

So the incumbent's SRX topology is a one-way cascade —
`PT →(ksr0=367)→ SRX ⇄(ksr2srd=128 / ksrd2sr=9.1)⇄ SRD →(kmsrd=31)→ PD` —
with **no direct force-dependent SRX→DRX return at all**. That is a structural fact
about the current best fit, not a setting anyone chose, and it makes "revive `kmsr`"
a mechanism arm in its own right (§4, R1).

`UseGlobalOccupancySaturation` is on with `P_bound_max = 1.135`, against a
physiological prior of 0.12–0.15 — with `P_bound ≈ 0.15` that is a flat ×0.87 on
attachment, i.e. a mechanism that has drifted into being a constant absorbed by `ka`.
`UseDynamicPassive` is on with `c_titin_visc = 0` — exactly inert.

---

## 3. Excluded from the scan, and why

| Class | Switches | Why |
|---|---|---|
| Numerics / solver | `UseFastPPval`, `UseSpaceInterpolation`, `UseSpaceExtension`, `UseSpaceDiscretization`, `UseCompiledMex`, `UseJPattern`, `UseCalculatedN`, `UseSafeSRReset`, `UseStrictDetachmentAt` | not physiology; flipping changes accuracy/speed, not the model |
| Protocol / harness | `UseSlack`, `UseKtrProtocol`, `UseSLInput`, `UseCa`, `UsePooledTarget`, `UseForceOnsetShift` | select which experiment runs or how data is aligned |
| **Dead in the current ODE** | `UseP31Shift`, `UseTORNegShift`, `UseAtpOnUNR`, `UseMutualPairingAttachment`, `UseKstiff3`, `UseTitinModel` | no reader in `dPUdT_CombinedTransitions`; they survive only in the legacy `dPUdTCa.m` (or, for `UseKstiff3`, only in `assertActiveParams`). Flipping them is provably a no-op |
| **Provably inert at 8 mM** | `UseAtpDetachR2`, `UseAtpKmsrd` | both gates are normalised to 1 at `MgATP_ref = 8`; identically neutral in an 8-mM-only objective. Test them in the ATP worktree against 2 mM |
| Superseded | `UseTargetZoneSaturation` | deprecated per-bin occupancy; `UseGlobalOccupancySaturation` replaces it and they are mutually exclusive in the code |
| Own campaign, not an ablation | `NumberOfStates = 3`, `UseUniformTransitionFunc`, `UseWDetachment`, `UseNegativeForceRip` | alternative rate-law families / extra state: a full refit, not a minor re-optimisation. The 3-state case is the known open item (iron-law `XTOR = attached × k2eff`) |

A one-evaluation **pre-flight** (§6 stage 0) confirms the "dead" and "inert" claims
empirically before anything is spent on them.

---

## 4. The arms

`add` = parameters that become live and must be in `cfg.pool` **and**
`cfg.compulsory` for that arm. `drop` = parameters that become dead and must leave
the pool (otherwise the search wastes draws on ghosts). `compensators` = existing
parameters that must be free so the mechanism is allowed to be absorbed —
this is the part that makes it an ablation rather than a perturbation.

### Tier A — nominally on, plausibly inert. Cheapest possible wins.
Expect Δ ≈ 0; if confirmed, delete the mechanism and its parameters outright.

| # | Arm | Now | add | drop | compensators |
|---|---|---|---|---|---|
| A1 | `UseGlobalOccupancySaturation` → 0 | on, `P_bound_max=1.135` | — | `P_bound_max` | `ka`, `A0` |
| A2 | `UseDynamicPassive` → 0 | on, `c_titin_visc=0` | — | — | `eta_M`, `k_pas` |
| A3 | `OccupancyForm` `linear`→`langmuir` | linear | — | — | `P_bound_max`, `ka` |
| A4 | `RegAvailShorteningOnly` → 0 | 1 (directional) | — | — | `A0`, `v_ref_reg`, `tau_reg`, `ka` |
| A5 | `useHalfActiveForSR` → 1 | 0 | — | — | `sigma1`*, `sigma2`, `sigma_srd1`, `ksr0` |

\* dead unless R1 revives `kmsr`.

### Tier B — on and load-bearing. The real ablations.

| # | Arm | add | drop | compensators |
|---|---|---|---|---|
| B1 | `UseRegistrationAvailability` → 0 *(state dim −1)* | — | `A0`, `v_ref_reg`, `tau_reg` | `ka`, `k2`, `PieceWiseStrainDep2*` (F-V shape), `kA2re` |
| B2 | `UseLoadRealign` → 0 *(state dim −1)* | — | `k_cr`, `tau_cr`, `F_cr` | `kmsrd`, `tau_reg`, `kA2re`, `ksr2srd` |
| B3 | `UseA2Reattaching` → 0 | — | `kA2re`, `kA2hop`, `sA2hop` | `k1`, `k_1`, `kd`, `PieceWiseStrainDepR1D*` |
| B4 | `UseSuperRelaxedADP` → 0 | **`kmsr` (must be revived — see below)** | `kmsrd`, `ksrd`, `ksrd2sr`, `ksr2srd`, `sigma_srd1`, `sigma_srd2` | `ksr0`, `sigma1`, `sigma2` |
| B5 | `UseSuperRelaxed` → 0 | — | `ksr0`, `kmsr`, `sigma1`, `sigma2`, + all SRD params | `k1`, `k_1`, `ka`, `kd`, `k2` |
| B6 | `UseNegativeKstiff` → 0 | — | `kstiff1_n`, `kstiff2_n` | `kstiff1`, `kstiff2`, `estiff`, `dr` |
| B7 | `UseMaxwellDashpot` → 0 | `c_titin_visc` (seed ≠ 0) | `kSE_M`, `eta_M`, `x_M_slack` | `k_pas`, `gamma`, `Lsc0` |
| B8 | `UseMaxwellTensionOnly` → 0 (dashpot kept) | — | `x_M_slack` | `eta_M`, `kSE_M` |
| B9 | `UseLatticeSpacing` → 0 | — | `d10_ref`, `SL_ref_lattice`, `R_thick`, `R_thin`, `d_optimal`, `sigma_lattice` | `L_thick`, `L_thin`, `L_hbare`, `ka` |
| B10 | `UseOverlapFactor` → 0 | — | — | `ka`, `L_thick`, `L_thin` |
| B11 | `UseOverlap` → 0 | — | `L_thick`, `L_hbare`, `L_thin` | lattice params, `ka` |
| B12 | `UseSerialStiffness` → 0 | — | `kSE`, `ekSE` | `mu`, `mu_neg` |
| B13 | `UseStrainDep4R1D` → 1 (Gaussian R1D) | `kd_sigma` | 8× `PieceWiseStrainDepR1D*` | `kd`, `PieceWiseStrainDep*` |

**B4 is the one that breaks if flipped naively.** With `kmsr = 0` the only exit from
SRX is `ksr2srd → SRD → kmsrd → PD`. Turning off `UseSuperRelaxedADP` sets both
`RSRD2PD` and `RPD2SRD` to zero, so SRX becomes a black hole and the model empties
into it. The arm is only meaningful with `kmsr` seeded non-zero and compulsory.

B12 is expected to be catastrophic. It is kept as a **dynamic-range anchor**: an
ablation study with no large negative result has not demonstrated it can detect one.

### Tier C — off, candidate additions. Each needs a hand-picked seed value (T1).

| # | Arm | add (seed) | compensators | why |
|---|---|---|---|---|
| C1 | `SRDOutflowK > 0` *(not a `Use*` flag)* | `SRDOutflowK` ~ 0.05 | `kmsrd`, `sigma_srd1`, `k_cr` | **highest prior value.** Turns the SRD drain from first-order into fixed-dwell, so a bigger slug takes proportionally longer — the sign the post-restretch amplitude dependence needs, and the open structural defect (model's rsK amplitude slope has the wrong sign) |
| C2 | `UseViscoelasticSE` → 1 | `c_SE_visc` ~ 0.01 | `kSE`, `ekSE`, `mu` | SE dashpot resists rapid lengthening → directly targets the re-stretch overshoot |
| C3 | `UseA1AttachmentKernel` → 1 | — (`A1AttachmentWidth=0.006` already set) | `ka`, `dr`, `kstiff1` | distributed attachment; the only ON-side arm that needs no new value |
| C4 | `UsePassiveForSR` → 1 | — | `ksr0`, `sigma2`, `sigma_srd1`, `kmsrd` | changes *what force the thick filament senses* from F_total to F_passive. No new parameters |
| C5 | `UseLimitedSRDTransition` → 1 | `Srd_max` ~ 20, `Srd_n` | `kmsrd`, `ksr2srd` | currently PD→SRD flux is exactly zero (`ksrd=0`); this is the alternative way to give it one |
| C6 | `UseDirectSRXTransition` → 1 | — | `ksr0`, `kmsr`, `ksr2srd`, `kmsrd`, `sSRXrip` | routes **all** R2 detachment to SRX (now partial, `SRXFromR2HighStrain=0.213`) |
| C7 | `UseStretchActivation` → 1 | `k_SA` ~ 0.5 | `ka`, `kA2re` | ka rises with lengthening velocity |
| C8 | `UseVelGaussAttachment` → 1 | `v_att_*` (defaults exist) | `ka`, `A0` | velocity-dependent ka — a **substitute** for B1, not an addition; run as B1-off + C8-on |
| C9 | `UseVernierVelocity` → 1 | `alpha_vernier=0.3`, `v_ref_vernier=1` | `ka` | third form of the same idea; mutually exclusive with `UseTargetZoneSaturation` |
| C10 | `UseA2MechanicalRecocking` → 1 | `k2d` ~ 50, `drmr` ~ 0.01 | `kA2hop`, `k2` | ATP-free recocking of high-strain heads; overlaps `kA2hop` → substitution arm |
| C11 | `UseA2AttachmentShift` → 1 | `s_threshold_R/L` (`d_actin`, `slope` set) | `kA2hop`, `kd` | hop to the neighbouring actin site; overlaps C10 |
| C12 | `UseCatchBond` → 1 | `k_catch_bond` ~ 0.1 | `kd`, `k2`, `kA2re` | tried before: filled the valley but cost `peak1_y` (see `a2reattaching` note). Re-test under the floored objective |
| C13 | `UseR2Ceiling` → 1 | `k2max=1500`, `K_T2c=0.15` (set) | `k2`, `PieceWiseStrainDep2*` | ceiling on per-head detachment rate. **Not** ATP-normalised, so it is testable at 8 mM |
| C14 | `UseA2Popping` → 1 | `pop_s` ~ −0.05 | `kd` | hard removal at extreme negative strain; mostly a numerical guard — low priority |
| C15 | `UseD0State` → 1 *(state dim +1)* | `D0FromR2` or `k2rip` > 0 | `k_D0T`, `alphaRip`, `kd` | apo pool. Nearly inert at 8 mM (`k_D0T·8/8.15`); **belongs in the 2 mM worktree** |
| C16 | `UseTitinIdentifiedPassive` / `UseTitinInterpolation` → 1 | needs `TitinTable` | `k_pas`, `gamma`, `Lsc0` | alternative passive form — **blocked until the table is confirmed present** |

### Tier R — parameter revivals (structural, but not `Use*` flags)

| # | Arm | seed | compensators |
|---|---|---|---|
| R1 | `kmsr > 0` — restore the direct force-dependent SRX→DRX return | ~10 | `sigma1`, `ksr0`, `ksr2srd`, `kmsrd` |
| R2 | `ksrd > 0` — restore first-order PD→SRD | ~1 | `sigma_srd2`, `kmsrd` |
| R3 | `D0OutflowK > 0` | ~0.05 | `k_D0T` (only with C15) |
| R4 | `dA2re ≠ 0` — landing strain for re-attached heads | ~0.005 | `kA2re`, `sA2hop` |

R1 and R2 revive `sigma1` and `sigma_srd2`, which are in the pool being drawn and
optimised today while having no effect on the trajectory.

**Counts:** 5 (A) + 13 (B) + 16 (C) + 4 (R) = **38 arms**, plus ≥3 control replicates.

---

## 5. Scoring

Per arm, record:

1. **Floored scalar** — what the campaign optimises. Ranked against **σ_reopt** from
   the control replicates, never against 0.9081 (T5).
2. **Raw χ² vector** (all 25 `fn` entries, unfloored) — catches the regressions the
   deadband hides (T4). Report the largest 3 movers in each direction.
3. **Δ#parameters** — the parsimony ledger. A mechanism that buys less than σ_reopt
   while costing 4 parameters is a deletion even when Δcost is slightly negative.
4. **Physiology tiles separately** — `XTOR`, `SRX_ss`, `attached_ss`, `PT_ss`,
   `k2/ksrd/kmsrd` bounds are prior knowledge, not data. B5 (no SRX) is pre-judged by
   `SRX_ss[1]|0.01-0.40`; that is legitimate, but it must be visible as a separate
   column rather than buried in the scalar.
5. **Unscored predictions** — the ktr protocol and the staircase are not in the
   objective. A mechanism that improves them while neutral on cost is a *different*
   kind of win and should be logged (`RunKtr`/`RunStairs` on the arm's best point).

---

## 6. Execution

**Stage 0 — pre-flight, ~1 evaluation per arm (~40 min total, no optimisation).**
Flip each switch at frozen parameters and record Δcost. Three outcomes:
*no change* → dead switch or inert parameter, report and drop the arm;
*change* → proceed; *penalty/NaN* → the arm needs a seed value or a companion
revival (B4 is the known case). This is the cheapest possible filter and it
empirically confirms every "dead"/"inert" claim in §2 and §3.

**Stage 1 — budget-matched minor re-optimisation, surviving arms.**
Identical for every arm and every control replicate: 1 round, `N_DRAW = 10`,
`SIMPLEX_EVALS = 60`, `KICK_FRAC = 0.10`, `SURR_EVALS = 0`, `MAX_RUN_TIME = 300`.
`cfg.compulsory` = that arm's mechanism parameters + its named compensators,
truncated to 10. Pool = base pool ± the arm's `add`/`drop`.

≈60 evals × 38 s ≈ **38 min per arm**. 38 arms + 3 controls ≈ **26 h serial**,
≈ **6.5 h on four instances** (the engine saturates at 5–6 workers; extra cores buy
parallel arms, not faster evaluations — see `RunOptimRigAware.m` header).

**Stage 2 — promote.** Arms beating the control by > 2σ_reopt get a full-length
re-optimisation before being believed. One round is a screen, not a result: a
short-budget campaign can rank arms wrongly (`stairs_prediction_estiff`).

---

## 7. Open decisions

1. **Budget.** 38 arms × 1 round is the screen above. Halving `SIMPLEX_EVALS` doubles
   the arm count per hour but moves the screen closer to a fixed-parameter probe.
2. **Where the ATP-dependent arms run.** C15 and the two excluded ATP gates only mean
   something against 2 mM data, which lives in the `ATP-wt-lowatp` worktree. Either
   run a second scan there or leave them out and say so.
3. **Combinations.** Everything here is one-at-a-time. Mechanisms that overlap
   (B1/C8/C9 velocity-dependent attachment; B3/C10/C11 A2 re-binding; C1/B2 slow
   recovery) can compensate for *each other*, so a one-at-a-time scan will under-value
   each. The substitution arms above are the minimum; a full pairwise scan over the
   three families is ~20 more arms.
