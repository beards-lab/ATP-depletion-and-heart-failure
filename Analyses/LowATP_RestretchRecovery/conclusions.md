# Post-restretch recovery is ATP-limited — and the model has the asymmetry BACKWARDS (2026-07-27)

Base: `params/lowATP_TNF_best.m` (TNF baseline + `UseAtpDetachR2`/K_T2=1.97 + `UseAtpKmsrd`/K_srx=13).
Reproduce: `AnalyzeRecoveryData.m` (data only), `ModelVsDataRecovery.m` (model vs data).

## The protocol has two different recoveries
Per cycle: **slack** (force→0) → hold_short → **restretch ramp** → **hold_long (~0.28 s; 0.84 s on cycle 5)**.
So there are two force rises: the post-slack **redevelopment** (= ktr) and the post-restretch **recovery**.

## Finding 1 — in the DATA, post-restretch is SLOWER than post-slack, and the penalty is ATP-dependent

| | post-slack τ | post-restretch τ | restretch penalty |
|---|---|---|---|
| **2 mM** | 42.2 ms (ktr 23.7/s) | **71.8 ms** (13.9/s) | **1.70× slower** |
| **8 mM** | 23.6 ms (ktr 42.4/s) | 28.0 ms (35.8/s) | 1.18× slower |

Low ATP slows the post-restretch recovery (**2.57×**) *more* than it slows ktr (**1.79×**).
The user's read of the data is correct: after restretch force keeps creeping up slowly, more slowly
than redevelopment after a slack at the shorter length.

**Per-cycle structure is the tell:**
- post-slack τ **rises steeply** with deeper slack (2 mM: 32→38→45→54→61 ms)
- post-restretch τ is **flat** (2 mM: 69→67→72→78 ms; 8 mM: 28→30→28→26 ms)

The restretch always returns to the *same* length, so the post-restretch state is reset identically each
cycle → its own characteristic, shortening-history-independent timescale set by [ATP]. This is a
**distinct process**, not just "re-attachment of the heads the slack removed".

**It cannot be passive/viscoelastic or inter-sarcomere elasticity alone:** the timescale is 2.6×
ATP-dependent. Passive titin/SE creep does not care about [ATP]. It is an ATP-consuming process.

## Finding 2 — the MODEL is 4.4× too FAST and inverts the asymmetry

| | post-restretch τ data | model | valley drop data | model |
|---|---|---|---|---|
| 2 mM | 71.8 ms | **16.3 ms** (4.4× too fast) | 18.4 | **22.6 (too deep)** |
| 8 mM | 28.0 ms | **10.4 ms** (2.7× too fast) | 16.8 | 13.7 (too shallow) |

Restretch-penalty ratio: **data 1.70× / 1.18× SLOWER; model 0.31× / 0.35×** — the model recovers ~3×
*faster* after a restretch than after a slack. **The sign is wrong.** At 2 mM the model plunges into a
too-deep valley and snaps back within ~50 ms, while the data dips less and creeps for >200 ms.

## Why — the ATP gate is numerically irrelevant exactly where the restretch lives

Active branch (`UsePieceWiseStrainDep=1`), `dPUdT_CombinedTransitions.m:400`:
`R2 = k2·p2·ppval(PieceWiseStrainDep2, s+dr2-dr)` — then the ATP gate (`:451`) multiplies it **uniformly**
by `gATP_R2(2 mM) = 0.628`.

Measured strain profile of R2 (k2=101.2/s):

| strain | strain mult | R2 @8 mM | R2 @2 mM | τ @2 mM |
|---|---|---|---|---|
| ~0 (isometric) | 0.4 – 1.2 | 41–124 /s | 26–78 /s | **25–39 ms** |
| 0.02 (post-stretch) | 21 | 2120 /s | 1332 /s | 0.75 ms |
| 0.05 – 0.10 | 42 – 46 | 4246–4632 /s | 2667–2909 /s | **0.34 ms** |

At high positive strain the ATP gate changes the detachment time constant by **0.13 ms** — invisible
against a 72 ms recovery. So:

- **Near zero strain** (isometric: steady force, ktr) the strain multiplier is ~1 and the ATP gate has real
  leverage → this is exactly why K_T2/K_srx fit steady force and ktr so well.
- **At stretch strains** the strain multiplier is 20–46×, R2 ≈ 2000–4600 /s, so **every attached bridge is
  annihilated in <1 ms regardless of [ATP]** → valley too deep, nothing left "stuck", and the recovery is
  just fresh re-attachment at the fast intrinsic rate.

**The physiological error:** you cannot detach a strongly-bound head faster than ATP can bind it. After a
stretch the head sits in a post-stroke A·M·ADP state; release requires ADP release *then* ATP binding
(second-order, `k_on·[ATP]`). That imposes a **rate ceiling** that strain can approach but never exceed.
The model has no ceiling — `R2(s)` grows without bound with strain, and a multiplicative gate cannot
impose one.

## Proposal — replace the multiplicative gate with a series (ATP-binding) CEILING

    1/R2_eff = 1/R2_strain(s) + 1/R2_ATPmax,    R2_ATPmax = k2max · MgATP/(MgATP + K_T2)

- low strain → `R2_eff ≈ R2_strain` (isometric behaviour, steady force and ktr preserved)
- high strain → `R2_eff → R2_ATPmax` (a hard, [ATP]-set ceiling; heads stay attached and strained)

Predicted consequences, all matching the residuals: shallower valley (fixes `vall_y`, the largest
remaining 2 mM feature residual), a slowly-clearing stuck population → slow post-restretch creep, and a
restretch penalty that **grows** at low ATP.

Back-of-envelope from the data (2 mM/8 mM recovery-rate ratio 14/36 ≈ 0.39) gives `K_T2 ≈ 8–9 mM` for the
ceiling — i.e. **ATP is limiting even at 8 mM**.

## Honest caveat — this is NOT purely a low-ATP fix
The model is 2.7× too fast post-restretch at **8 mM** as well, and its 8 mM valley is too *shallow*. So the
missing rate ceiling is a **baseline structural deficiency**, amplified at low ATP. Implementing it properly
perturbs the high-ATP baseline and would require re-fitting TNF — unlike K_T2/K_srx, it cannot be made
inert at `MgATP_ref` while still doing its job. That is the real cost of this mechanism, and it should be a
deliberate decision.

Alternative not excluded: **inter-sarcomere redistribution** (a single half-sarcomere model structurally
cannot produce it). It is cross-bridge-mediated so it *can* be ATP-dependent — it may contribute to the
slow tail alongside the ceiling.

Figures: `results/recovery_data.png`, `results/recovery_model_vs_data.png`.
