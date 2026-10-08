# How much SRX is there during contraction — and can the model tell? (2026-09-27)

**Question.** The model's super-relaxed pool (`P_SR` + `P_SRD`) is load-gated
(entry `ksr0·e^(−F/σ₂)`, escape `kmsrd·e^(+F/σ_srd1)`). Through that gate it helps
shape the length dependence and the recovery transients. The mant-ATP chase defines the
SRX by a slow nucleotide-turnover phase. Its time constants (tens to hundreds of seconds)
are 3–4 decades slower than the model's rates, and "SRX is force-sensitive" is the
justification we give for the gap. Can the model estimate the SRX fraction *during
contraction*? How uncertain is it, what does a low vs a high SRX change, and which
experiment would decide?

**Short answer.**

1. **The time-scale gap is a speed gap, not a size gap, and load cannot close it.** At
   rest the model's pool is **0.27** of heads. The mant fit's relaxed SRX is **0.35**,
   and the two recruit/escape ratios agree within ×1.2 (6.4 vs 5.2). What differs is
   speed: every exchange rate is **×6 000–15 000** faster. Scaling the model's rest
   rates down by one factor, **µ ≈ 1.1×10⁻⁴**, reproduces the mant chase (slow-phase
   share **0.47 vs 0.48**). Load cannot supply that factor: moving the escape rate
   ×5 700 would take **≈200 kPa** at σ_srd1 = 22.8 kPa, about 3× the isometric stress.
   The force gate is therefore the *within-activation* modulation (×141 on the
   recruit/escape ratio from rest to isometric). The rest→active speed-up has to come
   from **activation itself** (Ca²⁺ / thin-filament / MyBP-C signalling), which the
   model does not represent: it only exists at saturating Ca.
2. **During contraction the model's "SRX" is a fast, load-gated switch, not the mant
   population.** Isometric SRX is **0.030** with a **10 ms** mean dwell. A virtual chase in
   the active state has no mode slower than 43 ms. At the model's own rest rates the
   chase would be a **single** exponential: a mant experiment on the model would report
   "no SRX" although 27 % of heads sit in it. The model's SRX is kinetically the
   **structural OFF state** of X-ray work (switches in ms, follows stress). It is not
   the biochemical SRX, and the two are only the same *pool* at rest, where both
   measures apply.
3. **No, the total pool is not identifiable from the current 8 mM data, at least between 50 %
   and 200 %** (§8). With the pool *pinned* at ×0.5 or ×2 and 9–12 other knobs refitted,
   the fibre returns to χ² **3.0–4.8** (17 features, data-SD units). The same refit applied
   to the incumbent pool gives **2.4–3.9**, and the incumbent itself is 6.5, so the pinned
   pools cost at most **Δχ² ≈ 1.4** against their control. The knobs that absorb it are
   cross-bridge ones (`ka`, `k2`, `kd`, `kah`, `k_cr`), with the SRX *kinetics* as an
   optional extra. What the data *do* require is that the pool **exchange fast**:
   slowing it ×10 at fixed size costs Δχ² ≈ 51 even after `ka` compensation (§3). A
   slowly exchanging pool on top is invisible (degenerate with head number). The earlier
   estimate of 0.015–0.08 from `ka`-only compensation (§3) was an artefact of allowing one
   knob and is superseded. On the project's own weighted objective, ×0.5 is equally
   hidden (1.8–2.2 vs control 1.2–2.1). ×2 is not, but almost only because of the
   objective's very tight F–V term (its variance ignores the 13.6 % reproducibility).
4. **Low vs high SRX at the same force:** more SRX slows `ktr` (44 → 26–38 s⁻¹), lowers
   ktr-protocol k80, depresses F–V power, and needs more `ka` (up to +55 %). Less SRX does
   the opposite. **It does not change the net ATPase** (17.24 s⁻¹ per head in all 16
   variants): at matched force the SRX takes heads only from the *detached* pool. The
   objective's `XTOR` suggested savings of up to 4×, but `XTOR` is the gross DRX hydrolysis
   flux and misses the hydrolysis that moves into the SRX path (§4).
5. **~~The experiment that decides it is~~ a release–hold–restretch series with short holds
   (5–80 ms) separates large, uncompensated SRX changes, but §9 shows it cannot pin the
   pool once the other knobs are free.** A fast SRX refills in slack within ~10 ms and drains within ~40 ms of
   force redevelopment. This produces a hold-duration signature in the redevelopment rate
   (k₆₃ suppressed ~20–30 % at holds ≤ 20 ms vs a frozen pool), separating every variant at
   **3–7 SD** in one preparation. The existing slack protocol's holds (50–150 ms, confounded
   with release amplitude) sit just past that window, which is why the current data pin the
   pool only weakly. An activation step from rest (caged Ca²⁺ / fast solution switch) reads
   the *rest* pool. Paired ktr and SL series are weak (≤ 2 SD). Time-resolved X-ray on the
   same release–hold protocol would read the OFF fraction directly (§5).
6. **What the SRX pool does to the features: it sets the redevelopment *kinetics*, not the
   force levels** (§6, uncompensated kicks of inflow and outflow ÷11…×11). Slack ktr,
   ktr-protocol k80/k63/F(20 ms) and F–V power move by **2–6 SD**. Steady force, the
   redeveloped force `A`, the re-stretch peak, the valley and its timing, and the
   post-restretch rate stay within **1 SD**. More inflow and less outflow are mirror
   images, and both act through the *return time* of heads parked while the fibre is
   unloaded, not through the isometric pool size: outflow ×11 barely changes isometric SRX
   (0.030 → 0.023) but speeds ktr 44 → 74 s⁻¹. Removing the detour (inflow ÷11 or outflow
   ×11) gives ktr ≈ 70 s⁻¹, ktr-protocol k63 ≈ 80 s⁻¹ and F–V power at 4 ML/s ≈ 0.70.
   **In the incumbent, the SRX detour is what brings those down to the data** (44, 43,
   0.41). The ktr protocol is the most sensitive read-out: its ~5 ms hold sits in the
   refill window.

7. **The candidate protocol does not rescue it, and a seed-noise test would have said it
   does** (§9). In a synthetic-data twin, a ±20 % activated pool is below one-preparation
   noise on *every* protocol even before compensation (χ² ≤ 1.7). Refitting 13 knobs
   hides it completely (χ² ≈ 0.01), whether or not the release–hold series (5–80 ms) is
   in the fit. The Fisher information puts the 95 % range at ×/1 000+ with all
   protocols combined. The obstacle is the **load sensitivity of the inflow (σ₂)**: it
   sets the pool under load but is invisible to protocols that probe unloaded refill.
   With σ₂ known the range shrinks to ×/3.5, and to ×/2.25 with both load sensitivities
   known (~20 preparations for ±20 %). The proposed "5 %-noised seed" test reports ×0.80
   ± 25 % for a truth of ×0.80, which looks like success. That precision comes from the
   seed, which sits next to the truth, not from the data.
8. **Mavacamten-type drugs: the model predicts a robust *relative* change in activated
   SRX but not the absolute level, and predicts that SRX stabilisation barely touches
   force** (§10). For fibres that are indistinguishable at control, "recruitment ×3"
   raises activated SRX **×2.45 ± 0.1** (×2.0–2.7 over pools ×0.5–×2), and "retention
   (outflow ÷3)" ×1.5. The absolute value under drug inherits the control uncertainty.
   The measured drug effects are the same for all these fibres (χ² 0.1–0.4 at ±20 %), so
   the drug does not identify the baseline pool either. Two things *are* identifiable
   and useful. **Mechanism:** both models raise the *resting* pool similarly (0.27 →
   0.46 / 0.51, what a mant assay under drug would report), but recruitment slows slack
   ktr to 36 s⁻¹ and retention to 22 s⁻¹. So mant plus ktr under the same dose separate
   them. **Force:** both reduce maximal isometric force by only **1–3 %**. If mavacamten
   is found to cut maximal Ca-activated force much more than that, it is not acting
   through this load-gated pool.

Base: `robust1_seed` (= `refit9_opt` parameters + robust feature reader, 8 mM cost
0.4240), 8 mM only, pooled target. The 2 mM leg was dropped on purpose: in this setup
ATP reaches the model only through the k2 gate, so SRX has no ATP dependence of its own
and a 2 mM leg would test how SRX changes the ATP *response*, not how much SRX there is.
Reproduce: [`SrxRateBridge.m`](SrxRateBridge.m) (2 min) → [`RunSrxFamily.m`](RunSrxFamily.m)
(16 variants, ~40 min on 5 workers, resumable) → [`SummarizeSrxFamily.m`](SummarizeSrxFamily.m) →
[`RunSrxEnergetics.m`](RunSrxEnergetics.m) → [`RunSrxProtocols.m`](RunSrxProtocols.m) (90 runs,
~1 h on 2 workers) → [`SummarizeSrxProtocols.m`](SummarizeSrxProtocols.m) →
[`RunSrxRateKicks.m`](RunSrxRateKicks.m) (14 kicks, ~30 min on 5 workers, resumable) →
[`SummarizeSrxRateKicks.m`](SummarizeSrxRateKicks.m). Helpers:
[`srxApplyVariant.m`](srxApplyVariant.m), [`srxIsoForce.m`](srxIsoForce.m), [`srxProtocol.m`](srxProtocol.m),
[`srxScaleRates.m`](srxScaleRates.m), [`srxIsoState.m`](srxIsoState.m) →
[`RunSrxCompensation.m`](RunSrxCompensation.m) (Jacobian + pinned refits, 28 battery runs, ~1 h,
cached/resumable) → [`PlotSrxCompensation.m`](PlotSrxCompensation.m). Helpers: [`srxApplyTheta.m`](srxApplyTheta.m),
[`srxResiduals.m`](srxResiduals.m), [`srxTraces.m`](srxTraces.m) →
[`RunSrxTwin.m`](RunSrxTwin.m) (protocol Jacobian + ±20 % twin profile + Monte Carlo, ~1.5 h) →
[`PlotSrxTwin.m`](PlotSrxTwin.m) → [`RunSrxDrug.m`](RunSrxDrug.m) (6 fibres × control/2 drugs, ~1 h) →
[`PlotSrxDrug.m`](PlotSrxDrug.m). Helper: [`srxProtoFeatures.m`](srxProtoFeatures.m).

---

## 1. The model's SRX and the mant chase, rate for rate ([`SrxRateBridge.m`](SrxRateBridge.m))

The detached pool of `dPUdT_CombinedTransitions` and the kinetic core of the mant-ATP model
(`../mantATP-SRX`, bioRxiv 10.1101/2025.08.06.667936; chase subsystem in
`identifiability-analysis/chase_identifiability.py`) have **the same topology**:

| sarcomere | mant | state |
|---|---|---|
| `P_SR` | S_T | SRX·ATP |
| `P_SRD` | S_D | SRX·ADP·Pi |
| `PT` (`PuATP`) | D_T | DRX·ATP |
| `PD` (`PuR`) | D_D | DRX·ADP·Pi, actin-ready |

`results/rate_bridge.txt`:

| step | model F = 0 | rest (9.4 kPa passive) | isometric (71 kPa) | mant | model(F=0)/mant |
|---|---|---|---|---|---|
| DRX→SRX recruitment (`ksr0` / kDT) | 367 | 260 | 27.3 | 0.0244 | 1.5×10⁴ |
| SRX→DRX escape (`kmsrd` / kSD) | 26.8 | 40.5 | 601 | 0.0047 | 5.7×10³ |
| DRX·D→SRX·D (`ksrd` / none) | 10 | 7.8 | 1.5 | — | — |
| hydrolysis in SRX (`ksr2srd` / kH) | 130 | 130 | 130 | 18 | 7.2 |
| hydrolysis in DRX (`kah` / kH) | 176 | 176 | 176 | 18 | 9.8 |

| | model rest | model isometric | mant relaxed |
|---|---|---|---|
| SRX fraction (SR+SRD) | **0.274** (SR 0.045, SRD 0.229) | **0.030** (SR 0.023, SRD 0.006) | **0.345** |
| recruit/escape ratio | 6.42 | 0.045 | 5.19 |
| mean SRX dwell | — | **10 ms** (flux 3.0 s⁻¹ per head) | 213 s (1/kSD) |

"Rest" means attachment switched off (`ka` = 0) at SL 2.2 µm, so only the passive force
loads the gate. Nothing in that number was fitted to rest: SRX was identified only at
saturating Ca. The agreement with the mant fraction is an independent check, and it
**passes**.

**Virtual chase.** The label is released from DRX·ADP·Pi at the mant `kleak`
(0.027 s⁻¹, basal product release; the model has no actin-free release step, so this is
borrowed). The initial condition is the aging steady state. Results:

| chase run on | slow rates (s⁻¹) | slow-phase share |
|---|---|---|
| mant rates | 0.0039, 0.027 | 0.48 |
| model rest rates | 0.019 (one phase) | — (single exponential) |
| model rest rates, exchange × µ = 1.1×10⁻⁴ | 0.0039, … | **0.47** |

With the model's fast exchange, heads interconvert thousands of times before they release
a nucleotide. The chase therefore sees one population turning over at the pool-averaged
rate. Two phases become resolvable (ratio ≥ 3, each ≥ 10 %) only below
**µ = 6×10⁻⁴**, i.e. escape ≤ 0.026 s⁻¹. **The existence of a two-phase chase at rest is,
in itself, proof that exchange at rest is slower than basal turnover.** The model is
consistent with that only if its exchange is ~10⁴× slower at rest than when active.

**Active chase.** The attached pool is lumped from the isometric fluxes (PD→A 43 s⁻¹,
A→PD 26.6 s⁻¹, A→PT 60.8 s⁻¹). The chase is dominated by one 23.5 s⁻¹ mode (the
cycle), with no mode slower than 43 ms. During contraction the SRX is not a separate
kinetic population.

![rate bridge](results/rate_bridge.png)

## 2. Reconciling the time scales: what the force argument can and cannot do

"The SRX is force sensitive" explains **why SRX empties on activation and refills in
slack**: the ratio swings ×141 between rest and isometric. It does **not** explain why the
exchange is 10⁴× faster in active muscle than in the chase, because the chase is run at
zero active force and the model's rates at F = 0 are *already* 10⁴× the mant rates.
Three readings are consistent with the numbers:

* **(a) Activation catalyses the exchange at nearly fixed equilibrium.** A uniform speed-up
  (µ ≈ 10⁴) takes the mant core to the model's rest state with no change in size. That is
  exactly what a Ca²⁺-/thin-filament-dependent release of the interacting-heads
  motif would do. The model then describes only the activated limit, and a twitch model
  needs `µ(Ca)`, not just `e^(F/σ)`. (CalciumTwitch §2.1 already found the diastolic SRX
  pool moves inside a 10 ms pulse, which requires the fast limit.)
* **(b) Two pools.** The mant SRX (slow, biochemical) and the model's SRX (fast,
  structural OFF) are different populations that coincide at rest. The model's pool
  would then map onto X-ray OFF-state switching (ms, stress-driven), and the biochemical
  SRX would be a slow sub-population that is effectively frozen over any mechanical
  protocol.
* **(c) The chase measures the wrong thing for mechanics.** Its slow phase is *nucleotide
  turnover* in SRX, bounded by escape. Escape could be fast in muscle and slow in
  the chase only if the chase conditions (relaxing solution, no Ca²⁺, often no
  thin-filament activation or blebbistatin) remove the signal that accelerates it.
  This is (a) seen from the assay side.

What the model can say: the *size* of the pool at rest is compatible with mant, the
*speed* in contraction must be ms-scale for the battery to be fitted (§3, family MA),
and nothing in the saturating-Ca data can distinguish (a) from (b). §3 adds why: a frozen
pool is equivalent to fewer heads, so (b)'s slow population cannot be seen mechanically at all.

## 3. Is the SRX fraction during contraction identifiable? ([`RunSrxFamily.m`](RunSrxFamily.m))

> **Read §8 for the answer on pool size.** This section lets only `ka` compensate. Its
> conclusions on *speed* stand, but its size range (0.015–0.08) does not survive a refit
> with the other cross-bridge knobs free.

Four one-parameter families on `robust1_seed`, each **ka-compensated** (secant on `ka`
so the isometric force at SL 2.2 is 71.03 kPa to ±0.05). SRXT_0/SRXD_0 start at the
variant's own isometric state:

| family | change | what it does |
|---|---|---|
| **L** | `ksr0`, `ksrd` × λ | more recruitment at every load (rest SRX rises too) |
| **S** | `σ₂`, `σ_srd1` × s | weaker load sensitivity: active SRX rises, rest SRX barely moves (0.27 → 0.32) |
| **M** | `ksr0`, `ksrd`, `kmsrd` × µ | entry/escape slowed; the pool *shrinks* (0.030 → 0.008) because SR→SRD hydrolysis still drains SR |
| **MA** | all five SRX rates × µ | pure speed change, pool size kept (0.030 → 0.034) |

The reference run through the same path scores 0.4259 (incumbent 0.4240; the difference
is the SRX initial condition). Distance from the data is measured **in data-SD units**,
independent of the objective's weights. The mechanics block is 8 features: `A`, `steady`,
`peak1_y`, `vall_y`, `ktrfit`, `rsK80`, `t0` against the pooled 8 mM target, plus
ktr-protocol `k80` against `pooledKtrTarget(8)`. The F–V block is normalised power at 0.5/1/2/4 ML/s
against the Baker sets, SD = max(across-set SD, 13.6 %). `results/family.txt`:

| variant | SRX iso | SRX_ss (obj.) | ka × | ktrfit | ktrP k80 | t0 ms | F–V pow @2 / @4 ML/s | Δχ² mech | Δχ² F–V |
|---|---|---|---|---|---|---|---|---|---|
| data | — | — | — | ≈44 | ≈43 | ≈11 | 0.595 / 0.413 | | |
| L × 0.25 | 0.008 | 0.020 | 0.985 | 59.4 | 57.9 | 9.3 | 0.760 / 0.629 | +9.2 | +18.6 |
| **ref** | **0.030** | **0.069** | 1 | 43.6 | 46.3 | 10.9 | 0.565 / 0.389 | 0 | 0 |
| S × 1.3 | 0.059 | 0.096 | 1.030 | 40.6 | 41.9 | 9.9 | 0.518 / 0.375 | +0.2 | +1.6 |
| S × 1.7 | 0.101 | 0.136 | 1.089 | 39.0 | 39.0 | 9.0 | 0.489 / 0.365 | +1.3 | +3.9 |
| S × 2.2 | 0.154 | 0.185 | 1.180 | 38.4 | 37.6 | 8.3 | 0.472 / 0.357 | +2.1 | +6.2 |
| S × 3 | 0.224 | 0.251 | 1.344 | 38.2 | 37.0 | 7.9 | 0.457 / 0.348 | +2.5 | +9.0 |
| S × 4 | 0.289 | 0.313 | 1.552 | 38.1 | 36.9 | 7.7 | 0.446 / 0.340 | +2.8 | +11.7 |
| L × 3 | 0.071 | 0.152 | 1.031 | 36.0 | 39.0 | 13.8 | 0.419 / 0.254 | +2.3 | +13.7 |
| L × 10 | 0.140 | 0.266 | 1.093 | 26.3 | 29.0 | 18.3 | 0.315 / 0.188 | +17.2 | +33.1 |
| L × 30 | 0.207 | 0.380 | 1.182 | 19.7 | 1.1 | 10.6 | 0.231 / 0.151 | +65.9 | +56.7 |
| MA × 0.1 | 0.034 | 0.038 | 1.032 | 70.3 | 47.3 | 6.7 | 0.769 / 0.710 | +18.7 | +32.1 |
| MA × 0.01 | 0.034 | 0.034 | 1.037 | 83.7 | 70.1 | 9.7 | 0.866 / 0.798 | +43.4 | +57.9 |
| MA × 0.001 | 0.036 | 0.034 | 1.037 | 80.3 | 71.6 | 9.4 | 0.857 / 0.789 | +41.6 | +55.0 |
| M × 0.001 | 0.008 | 0.008 | 0.991 | 80.1 | 71.4 | 9.5 | 0.863 / 0.795 | +40.7 | +57.0 |

("SRX iso" = isometric steady state at SL 2.2. "SRX_ss (obj.)" = the objective's reading in
the window before the *first* slack, which is not yet steady. Cycles 2–5 read 0.034 for
ref, matching the isometric value.)

![family](results/family_cost.png)

**Findings.**

* **The pool is pinned from below by ktr and F–V together.** Less SRX (L × 0.25) makes ktr
  +36 % and stops the high-velocity power from falling: Δχ² 28.
* **From above, how tightly it is pinned depends on the direction.** Along L (recruitment),
  ktr and t0 degrade quickly: Δχ² 16 already at iso SRX 0.07. Along S (weaker load
  sensitivity), **the mechanics hardly care**: Δχ² ≤ 2.8 from 0.03 to 0.29. The F–V
  power is what resists (Δχ² 1.6 → 11.7). Combined Δχ² = 4 is crossed between S × 1.3 and
  S × 1.7, i.e. iso SRX **≈ 0.07–0.08**.
* **The exchange must be fast.** Keeping the pool size but slowing every SRX rate ×10 (MA)
  costs Δχ² 51. Frozen (×10⁻³) costs 97: ktr jumps to ~80 s⁻¹ (data ~44) and F–V
  power stays high at 4 ML/s (0.79 vs 0.41). A pool that can refill in slack and in fast
  shortening, and drain as force redevelops, is *load-bearing* for ktr and for the F–V
  shape. M (a shrunk pool) and MA (a kept pool) give the same answer once frozen, which
  shows the damage comes from speed, not size.
* ~~Estimate: isometric SRX ≈ 0.03, 2σ range ≈ 0.015–0.08~~. **Superseded by §8:** this
  range comes from letting only `ka` compensate. With the other cross-bridge knobs free,
  ×0.5 and ×2 are both within Δχ² ≈ 1.4 of the control, so the *size* is not pinned. The
  *speed* conclusion (ms-scale exchange) stands.

**How much to trust it.**

1. The compensation is `ka` only. A full refit at a fixed SRX level (a real profile
   likelihood, several hours per point) could absorb part of each penalty, so the range
   above is a *lower* bound on the true uncertainty.
2. Two thirds of the constraint rides on the **Baker F–V curve**, borrowed from another
   lab. It is transferable at 8 mM (DataPooling §6: 0.21 sd), but with 13.6 %
   reproducibility, and its power is what `FV_normpowerAvg` (weight 1804) already made
   the most expensive term in every family. Drop it, and along S the pool is
   unidentifiable from our own slack/ktr data.
3. The constraint acts **through the model's mechanism**: SRX is where heads go when load
   falls. Another mechanism that depresses F–V power and slows ktr (e.g. registration
   availability `A_reg`, which acts on the same windows) could take over part of that
   role. The estimate is conditional on this structure.
4. **A slowly exchanging pool cannot be seen at all.** During every protocol it only
   removes heads (MA/M × 10⁻³ SRX stays flat through slack, §5 panel B), which is what a
   smaller `ka` or head number does. A mant-type slow SRX of, say, 10–20 % *in addition*
   to the fast pool is therefore degenerate with the force scale. The model estimates the
   fast, exchanging pool, not the total. (Inference from the frozen variants; an explicit
   frozen-plus-fast variant was not run.)

## 4. What low vs high SRX does, at the same force ([`RunSrxEnergetics.m`](RunSrxEnergetics.m))

| | less SRX (L × 0.25) | incumbent | more SRX via S (× 4) | more SRX via L (× 10) |
|---|---|---|---|---|
| iso SRX | 0.008 | 0.030 | 0.289 | 0.140 |
| `ka` to keep F_iso | × 0.985 | 1 | × 1.55 | × 1.09 |
| attached fraction | 0.284 | 0.284 | 0.284 | 0.284 |
| slack ktr / ktr-protocol k80 | 59 / 58 | 44 / 46 | 38 / 37 | 26 / 29 |
| t0 (ms) | 9.3 | 10.9 | 7.7 | 18.3 |
| F–V power at 2 ML/s (data 0.595) | 0.76 | 0.57 | 0.45 | 0.32 |
| objective's `XTOR` (gross DRX hydrolysis) | 18.4 | 15.0 | 7.7 | 4.3 |
| **net ATPase R2T** (s⁻¹ per head) | **17.24** | **17.24** | **17.24** | **17.24** |
| hydrolysis inside SRX (SR→SRD) | 0.8 | 3.0 | 10.7 | 13.3 |

* **More SRX slows everything kinetic** (ktr, k80, F–V power), and the pool has to be
  offset by faster attachment to keep force. Less SRX does the opposite.
* **SRX size does not change the energy cost of isometric force in this model.** The
  steady ATPase is the strong-detachment flux, and at matched force the attached pool and
  its turnover are identical to the digit: SRX competes only with *detached* heads. The
  model's SRX has only one route to saving ATP: lowering force (fewer attached heads). A
  real ATP saving at equal force would need the SRX to suppress a *basal* (actin-free)
  turnover, which the model does not have (§1: no `kleak`).
* **Side finding on the objective.** `XTOR` = `out.RTD` = `kah·PT` counts only DRX
  hydrolysis, gross of reversal. The incumbent's true net turnover is **17.2 s⁻¹ per head,
  above the 3–15 bound** the `XTOR` term is meant to enforce. The term also *rewards*
  moving hydrolysis into the SRX path (L × 0.25 pays +1.1 for `XTOR` 18.4 with identical
  net ATPase). Scoring `R2T` (or `RTD − RDT + RSR2SRD − RSRD2SR`) instead would fix both.
  Not changed here.

## 5. Which experiment would decide it? ([`RunSrxProtocols.m`](RunSrxProtocols.m))

Five ka-compensated variants: ref (iso SRX 0.030), L × 10 (0.140), S × 3 (0.224), M × 10⁻³
(0.008, frozen), MA × 10⁻³ (0.035, frozen). Four candidate protocols at 8 mM were run directly
through `evaluateModel`. Each difference from ref is divided by an assumed
within-preparation CV (rates 8 %, force ratios 5 %, rise times 10 %); ≥ 2 counts as detectable in
one preparation. `results/protocols.txt`:

| protocol | read-out | L × 10 | S × 3 | frozen (MA × 10⁻³) | verdict |
|---|---|---|---|---|---|
| **B** release 0.15 ML, hold T = 5–320 ms, restretch in 2 ms | k₆₃ after restretch vs T (absolute / shape) | 3.4 / 4.2 | **7.2** / 3.3 | 4.3 / 3.3 | **best: separates all** |
| **A** activation step from the relaxed steady state | t50 / t90 of force rise | **6.3** / 2.5 | 1.5 / 1.5 | (not at rest equilibrium) | reads the *rest* pool |
| **C** two ktr manoeuvres, interval 30–500 ms | k₂/k₁ | 0.4 | 2.0 | 0.5 | weak |
| **D** isometric SL 1.9–2.3 µm | F(SL)/F(2.2) | 1.8 | 0.8 | 0.2 | weak |

k₆₃ (s⁻¹) in protocol B:

| hold T (ms) | 5 | 10 | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|---|---|
| ref | 56.2 | 62.5 | 71.0 | 55.9 | 47.4 | 48.9 | 49.0 |
| L × 10 | 43.0 | 47.4 | 59.2 | 68.5 | 58.7 | 51.7 | 53.1 |
| S × 3 | 37.0 | 38.2 | 40.1 | 41.0 | 43.8 | 46.7 | 46.9 |
| MA × 10⁻³ (frozen, same size) | 72.8 | 88.3 | 74.4 | 49.9 | 51.1 | 53.1 | 53.4 |
| SRX at restretch, ref | 0.37 | 0.41 | 0.31 | 0.07 | 0.05 | 0.05 | 0.05 |

![protocols](results/protocols.png)

* **Why B works.** In slack, the fast pool fills from 0.03 to ~0.5 within ~10 ms
  (entry `ksr0` ≈ 370 s⁻¹ at F = 0). Once force redevelops at the short length, it drains
  again by ~40 ms. A restretch inside that window meets a fibre with its heads parked, and
  redevelopment is slower. **ref vs MA × 10⁻³ is the clean contrast: same size, only the
  speed differs**. It gives k₆₃ 20–30 % lower at holds ≤ 20 ms, converging by 40 ms. The
  position and depth of the dip read the refill/drain kinetics. Its height relative to the
  long-hold plateau reads the pool size (L × 10 shifts the peak to 40 ms; S × 3 flattens
  the curve). The frozen curves are not flat either (a 10 ms peak from the other
  transient states), so the SRX reading rests on the *contrast* and on the rest of the
  model being right. Fit the whole k(T) curve, not one point.
* **Why the current data barely see it.** The slack protocol's holds are 50–150 ms, and each
  is paired with a different release amplitude. By then the pool has drained (0.05 at
  restretch for T ≥ 80 ms). The protocol samples the plateau, not the dip.
* **Recommended protocol:** from isometric SL 2.2, release 0.15 ML (≥ 1 ms), hold
  **T = 5, 10, 20, 40, 80 ms** at one fixed amplitude, restretch in ≤ 2 ms, record 400 ms.
  Randomise T, bracket with a reference T = 80 ms, and score k(T)/k(80 ms) within a
  preparation, so rundown and cross-section scale cancel. It is the existing ktr manoeuvre
  with a variable hold, so the rig can already do it.
* **A (activation from rest)** is the protocol that links the model to the mant number. Its
  rise time is set by how many heads start parked at rest (L × 10: rest SRX 0.70, t50
  16 vs 8.6 ms). It needs a fast activation (caged Ca²⁺ photolysis or a sub-10 ms
  solution switch). With a mant chase on the *same* preparation it would test whether the
  rest pool the mant assay sees (~0.35) is the pool that has to leave on activation,
  i.e. readings (a) vs (b) of §2. The model cannot simulate the chase itself (no
  actin-free turnover), and the frozen variants cannot be brought to a rest equilibrium
  in 5 s, so they are left out of A.
* **Beyond what the model can simulate:** time-resolved X-ray (M3/M6 spacing, I₁₁/I₁₀) on
  protocol B would read the OFF fraction directly through the dip. Pharmacology that moves
  the SRX in a known direction (mavacamten up, dATP down), run through protocol B, would
  test whether the dip moves the way the model says. A relaxed/basal ATPase would test
  the missing `kleak`. An isometric ATPase would **not** help: at matched force it is
  SRX-blind in this model (§4).

## 6. What the SRX pool does to the features: inflow and outflow kicks ([`RunSrxRateKicks.m`](RunSrxRateKicks.m))

§3 asked what the data *allow*, with force held fixed. This section asks what the pool
*does*. SRX is read as one pool. **SRX+** scales the inflow: `ksr0` (PT → SRX·ATP) and
`ksrd` (PD → SRX·ADP·Pi) together. **SRX−** scales the outflow: `kmsrd` (SRX·ADP·Pi → PD;
`kmsr`, SRX·ATP → PT, is 0 in the incumbent). Each is run at ÷11, ÷3, ×3, ×11 (×11 = +1000 %),
**uncompensated**, so force is free to move. Each run gets the full 8 mM battery, a traced
slack–restretch cycle 1 and a traced ktr protocol. Shifts are in **data-SD units**
(shift ÷ pooled between-preparation SD). `results/kicks.txt`:

| feature | data | incumbent | SRX+ ÷11 | SRX+ ÷3 | SRX+ ×3 | SRX+ ×11 | SRX− ÷3 | SRX− ×3 | SRX− ×11 |
|---|---|---|---|---|---|---|---|---|---|
| isometric SRX | — | 0.030 | 0.003 | 0.011 | 0.073 | 0.155 | 0.045 | 0.025 | 0.023 |
| isometric force (kPa) | — | 71.0 | 71.7 | 71.5 | 69.9 | 67.7 | 70.2 | 71.3 | 71.4 |
| steady force | 65.5 | 68.0 | +0.1 | +0.0 | −0.1 | −0.4 | −0.1 | +0.0 | +0.0 |
| redeveloped force A | 49.8 | 49.2 | +0.1 | +0.1 | −0.2 | −0.9 | −0.3 | +0.1 | +0.1 |
| re-stretch peak | 73.8 | 74.4 | +0.2 | +0.1 | −0.2 | −0.9 | −0.3 | +0.1 | +0.1 |
| valley / valley time | | | 0.0 / −0.1 | 0.0 / 0.0 | −0.1 / 0.0 | −0.3 / −0.4 | −0.1 / −0.2 | 0.0 / −0.1 | 0.0 / −0.1 |
| post-restretch rate rsK80 | 47.5 | 47.8 | −0.1 | −0.1 | +0.3 | +1.4 | +0.1 | 0.0 | 0.0 |
| **slack ktr (fit)** s⁻¹ | 44.0 | 43.6 | **+3.3** (69.8) | +1.5 | −1.0 | **−2.4** (24.4) | **−2.8** (21.6) | +2.6 | **+3.9** (74.2) |
| slack t₀ | 11.0 ms | 10.9 ms | −0.6 | −0.6 | +1.0 | **+2.5** (18 ms) | (−3.5)* | +0.1 | −0.3 |
| **ktr protocol k80** s⁻¹ | 43.6 | 46.3 | +2.9 | +1.4 | −1.2 | **−3.1** (26.7) | **−4.4** (18.7) | +2.4 | +3.2 |
| **ktr protocol k63** s⁻¹ | 43.2 | 51.3 | **+5.0** (80) | +2.7 | −2.0 | **−4.6** (25) | **−4.1** (27.6) | +3.5 | **+5.0** (80) |
| **ktr protocol F(20 ms)/F** | 0.56 | 0.60 | +2.8 | +1.7 | −2.1 | **−5.6** (0.34) | −2.8 | +2.2 | +2.9 |
| **F–V power, 2 ML/s** | 0.60 | 0.56 | +3.2 | +2.0 | −1.8 | −3.1 | −2.6 | +2.2 | +3.0 |
| **F–V power, 4 ML/s** | 0.41 | 0.39 | **+6.0** (0.72) | +3.5 | −2.3 | −3.5 | −1.5 | +2.8 | **+5.5** (0.70) |

\* The t₀ reading at SRX− ÷3 (0.98 ms) is an extraction artefact, not a real delay. SRX− ÷11
is left out of the table: heads get trapped (SRX 0.6–0.9 through the slack hold), force
redevelops only to ~26 kPa (`A` −68 %), and the slack features stop being extractable
(ktr "345 s⁻¹", no re-stretch peak). It is outside any physiological range; it is shown
in the figures as the limit case.

![dose response](results/kicks_dose.png)

![heat-map](results/kicks_heatmap.png)

![traces](results/kicks_traces.png)

**What is affected, and what is not.**

* **Affected: every redevelopment rate and the F–V power, by 2–6 SD at ×/÷11.** These
  include slack ktr (fit and 80 %), the ktr protocol (k80, k63, F at 20 ms), and F–V power
  at 2 and 4 ML/s. Delay t₀ responds to inflow only (+2.5 SD at ×11: parked heads delay
  force onset).
* **Not affected (< 1 SD even at ±1000 %): the force levels and the re-stretch
  transient.** These are steady force, `A`, re-stretch peak, valley, valley time, and
  post-restretch rate rsK80. The exception is SRX+ ×11: rsK80 +1.4, `A`/peak −0.9. This
  agrees with RestretchSRXRecoil: the SRX brake disengages under load, so the
  post-restretch window, which runs at high force, cannot feel the pool.
* **Inflow up and outflow down are mirror images.** Both slow every rate and depress F–V
  power; the reverse speeds them. The pool acts as a **detour whose delay is paid when the
  fibre is unloaded**: in slack or fast shortening, heads enter within milliseconds, then
  wait for the escape `kmsrd·e^(F/σ_srd1)` (≈27 s⁻¹ at F ≈ 0, i.e. τ ≈ 37 ms) before they
  can re-attach. The traces show it: in the incumbent, SRX jumps 0.03 → 0.48 within ~10 ms
  of the slack, then drains over ~50 ms as force comes back. With SRX+ ×11 it reaches 0.83,
  and force redevelops visibly later in both the slack and the ktr protocol.
* **The features read the detour, not the isometric pool.** SRX− ×11 leaves isometric SRX
  almost unchanged (0.030 → 0.023), yet ktr goes 44 → 74 s⁻¹ and k63 51 → 80. The other
  way, SRX+ ÷11 empties the isometric pool (0.003) and gives the same kinetics. Both are
  the "no detour" limit. The side run with inflow *and* outflow ×11 raises the isometric
  pool to 0.11 but makes the kinetics *faster* (ktr +3.8 SD). A big but fast-exchanging pool
  behaves like a small one. This is why §3 found the pool size weakly determined along
  some directions: the data see the *return time at low load*.
* **How much of the fitted kinetics is the SRX?** In the "no detour" limit, ktr ≈ 70 s⁻¹,
  ktr-protocol k63 ≈ 80 s⁻¹ and F–V power at 4 ML/s ≈ 0.70. The incumbent (and the data)
  are at 44, 43–51 and 0.39–0.41. **The SRX detour accounts for roughly a third to a half
  of the slowing of redevelopment, and most of the high-velocity power decline, in the
  current fit.** It is a load-bearing mechanism, not decoration.
* **The ktr protocol is the most SRX-sensitive measurement we already have** (k63 ±5 SD,
  F(20 ms) −5.6 SD). Its hold (~5 ms at 20 % release) sits right in the refill window.
  The slack cycles (holds of ~300 ms) see the detour only through the first redevelopment,
  and F–V through repeated unloading during shortening.
* **The response is steep toward slower outflow.** ÷3 already halves ktr (21.6 s⁻¹); ÷11
  traps the heads. Toward faster outflow it saturates (×3 → ×11 changes little). The
  incumbent sits near the knee.
* Side runs, single inflow rates: `ksrd` ×11 alone (PD → SRX·ADP·Pi) carries most of the
  kinetic effect (ktr −2.5, k63 −4.5 SD) with the isometric pool barely changed (0.047).
  `ksr0` ×11 alone mainly raises the isometric pool (0.136) with milder kinetic effects.
  The ready (ADP·Pi) heads parking directly is the fast lane into the detour.

## 7. What this changes / next steps

* **Say "SRX" carefully.** The model's pool is a fast, load-gated OFF state, ≈3 % during
  isometric contraction and ≈27 % at rest. It matches the mant *size* at rest, but it
  is not the mant *kinetic* population. The `SRX_ss|0.10–0.40` literature bound (now
  `0.01–0.40`) was comparing a fast pool with biochemical numbers. Leaving it relaxed is
  right.
* **Any twitch / Ca model needs an activation-dependent exchange speed** `µ(Ca)` spanning
  ~10⁴. Force alone gives ×141 and would need ~200 kPa to do the rest.
* **Fix `XTOR`** to score net turnover (§4). This touches the incumbent's objective, so it is
  a separate change and was not made here.
* **Low ATP was deliberately left out.** If the SRX question is taken to 2 mM, it needs an
  SRX that feels ATP (`UseAtpKmsrd`, or ATP-dependent re-priming). The model's rest pool
  is mostly SRX·ADP·Pi (0.23 of 0.27), which is where nucleotide type would act.
* If protocol B is run, fit λ, s and µ to the k(T) curve jointly with the battery, as a
  proper profile (refit at fixed SRX) instead of the one-knob compensation used here.

## 8. Identifiability, done properly: pin the pool, refit everything else ([`RunSrxCompensation.m`](RunSrxCompensation.m))

**Question.** If the total SRX pool during contraction were really half (×0.5) or double (×2)
the incumbent's, could other knobs compensate so that the 8 mM data cannot tell?

**Design.** The isometric pool (SR+SRD at SL 2.2) is **pinned** by the inflow (`ksr0` &
`ksrd` together): 0.030 → 0.015 / 0.060. It is held there by a constraint residual
(SD 5 % of target). Two sets of knobs are then refitted, all as log multipliers:

* **set A, cross-bridge only:** `ka kd k1 k_1 k2 kah kstiff2 A0 k_cr` (+ inflow to hold the pool);
* **set B, A + the SRX kinetics:** outflow `kmsrd`, its load sensitivity `σ_srd1`, and the
  inflow's `σ₂`.

The target is the 17 features of §3/§6 (slack levels, rates and timing, ktr protocol,
F–V power) in data-SD units: **χ² = Σz²**. The Jacobian (13 knobs, step ×1.15) is taken
once at the incumbent. Then two damped Gauss–Newton steps are run (ridge 0.5, step ≤ ×4),
each verified by a full battery run, and the better point is kept. **Control:** the
incumbent pool (×1) goes through the identical refit. The honest comparison is
pinned-and-refitted vs control-refitted, not vs the incumbent. The incumbent is not at the
minimum of this χ², so part of any refit's gain is available at ×1 too.

![chi2](results/comp_chi2.png)

**Left: the verdict.** After refitting, every pinned pool fits the data about as well as the
refitted control:

| | uncompensated | set A | set B |
|---|---|---|---|
| pool ×0.5 | 26.3 | **3.3** | **3.8** |
| **control ×1** | 6.5 (= incumbent) | **3.9** | **2.4** |
| pool ×2 | 15.0 | **4.8** | **3.0** |

Worst gap to its control: **+1.4** (×0.5, set B). With 17 residuals and one pinned
parameter, a difference of ~4 would be needed to call the size identified. **Neither
halving nor doubling the pool can be detected.** Set A alone, with no SRX-kinetic knob,
already hides ×0.5 (3.3 < 3.9) and nearly hides ×2 (4.8 vs 3.9).

**Right: the project objective.** This is the weighted, floored sum the optimiser uses; the
refit did not minimise it. Even the control rises (0.43 → 1.2–2.1) when the refit targets
χ², so again compare within a column. At ×0.5 the pinned fits (1.8–2.2) match the control.
At ×2, set B reaches 4.3 and set A 57.6. The excess is almost entirely `FV_normpowerAvg`
(3.1 and 54.4), whose inverse-variance weight uses only the spread *between* the Baker
datasets. In reproducibility units (13.6 %) those same points are within 1 SD (left figure,
F–V rows). The rest is small floored terms (`rsTfitE`, `vall_y`, `peak1_y`) and the
artefactual gross-`XTOR` bound (§4). So the objective resists doubling mainly through its
own F–V weighting choice.

![residuals](results/comp_residuals.png)

**Feature by feature.** Pinning alone (orange squares) mostly throws the *rates*: ktr, ktr
protocol k63/F(20 ms) and F–V power, as §6 predicts. After refitting (blue/green), every
feature returns inside or at the edge of the one-preparation band, in all three panels,
and no single feature stays flagged. The largest leftovers (valley level, ktr-protocol
k63, valley time, rsK80, ≈ 1 SD) are the same ones the control keeps. They are
incumbent-model misses, not SRX signatures.

![knobs](results/comp_knobs.png)

**What does the compensating.** Two things are mixed here.

* **A common move present in all six refits, including the control:** `ka` ×1.4–1.5,
  `k_1` ÷2–3, `A0` ×0.66–0.85. It is an improvement the incumbent could take anyway, so it
  says nothing about SRX.
* **The SRX-specific part is what differs from the control's diamonds:**
  * *Pool ×2:* `kah` up (×1.3–1.45) and `k2` up (set A ×1.32) restore the fast
    redevelopment. In set B, `σ₂` ×1.49 instead makes the inflow *less* load-sensitive, so
    the extra pool sits in isometric contraction without adding much refill in slack. That
    is the near-degenerate direction the Jacobian flagged: σ₂ moves the pool the most
    (d log SRX/d log σ₂ = +2.0) with the smallest feature change of any knob (|dz| = 1.6).
  * *Pool ×0.5:* `k_1` ÷3 and `k_cr` ÷3 (set A), or `kmsrd` ×0.69 in set B, slow the
    return a little to replace the missing detour.

**The pool can be traded against the cross-bridge rates, because what the data see is
the redevelopment kinetics (§6), not the pool itself.**

![traces](results/comp_traces.png)

**The demonstration.** On the left, force in the slack–restretch cycle and the ktr protocol:
the compensated ×0.5 and ×2 fibres (solid) lie on the incumbent, and the uncompensated
ones (dotted) are visibly off. On the right, the SRX pool in the same runs differs by a factor of 4
throughout, at rest, at the slack peak (0.38 vs 0.54) and during redevelopment.
**Identical force, different SRX.** Nothing in these protocols measures the pool directly.

**Limits.** This is a linearised refit (one Jacobian, two steps), not a converged profile.
A full optimiser could only *lower* the pinned χ², which strengthens the verdict. Two
targets were tested (×0.5, ×2); beyond that range the refit was not tried. The pool
constraint drifted to ×2.4 in set A at ×2 (0.072), which is no weaker a test. Everything
here is 8 mM. The protocols that *could* break the degeneracy are §5's short-hold
release series (it reads refill and drain timing, which the compensation had to move) and
any direct structural read-out (X-ray OFF fraction).

## 9. Would the proposed identification protocol pin a ±20 % pool? A twin test ([`RunSrxTwin.m`](RunSrxTwin.m))

**The proposed test and what it can and cannot show.** The proposal: reduce or raise the
pool ±20 %, noise a parameter subset ~5 %, refit the control, then fit the low and high legs
from the same noised seed. That is a synthetic-data (twin) experiment, which is the right
frame. It has three problems as specified:

1. **Without measurement noise, it is too optimistic.** A noiseless twin recovers the truth
   whenever the pool is identifiable at all, however weakly. The scatter it shows is the
   optimiser's, not the data's.
2. **A noised seed tests the optimiser's basin, not identifiability.** In the linearised
   problem the least-squares answer does not depend on the start. A seed matters only if
   the optimiser stays near it, and then it pulls the estimate toward the truth, which it
   is 5 % away from. The real uncertainty of the other knobs is far larger: in §8 the
   compensating knobs moved 30–200 %.
3. **Refitting the incumbent against real data mixes in model misfit** (§8: the incumbent
   is not at the χ² minimum). In a twin the truth sits exactly at χ² = 0, so every residual
   is identifiability.

**What was run instead.** Truth = `robust1_seed` (pool 0.030). Synthetic data = the model's
own features. Noise = one preparation: the pooled between-preparation SD for the battery;
8 % for protocol B's k₆₃(T) and 5 % for its F₂₀(T) (release 0.15 ML, hold T = 5, 10, 20,
40, 80 ms, restretch); 10 % for protocol A (activation step, t50/t90). The same 13 knobs as
§8 were used. Three analyses:

![fisher](results/twin_fim.png)

**(a) Fisher information: the best precision the data could give, all 13 knobs free.**
Protocol B shrinks the pool uncertainty ~2.5× over the battery, and A + B ~3×. The
95 % range is still **×/1 000+** in every case, and even with every knob known within ×2
it is ×/10. The right panel shows why: the pool estimate is almost perfectly correlated
with **σ₂** (the load sensitivity of the inflow), with `A0`, `inflow` and `kmsrd` next. σ₂
changes the pool the most per unit change (d log SRX/d log σ₂ = +2.0) and the protocols the
least (|dz| 1.6 battery, 2.1 B, 0.4 A). It acts only under load, and every protocol here
reads the *unloaded* refill. Taking the load sensitivities as known (the two bottom rows)
gives ×/3.5 (σ₂ known) and ×/2.25 (σ₂ and σ_srd1 known). That still needs ~20–50
preparations for ±20 % (middle panel).

![profile](results/twin_profile.png)

**(b) Twin profile at ±20 %, the nonlinear check.** With the pool pinned at ×0.8 or ×1.2 and
**nothing refitted**, the fibre is already only χ² 1.7 / 1.0 (battery), 0.5 / 0.4 (B) and
0.2 / 0.1 (A) from the truth, all below one-preparation noise. A 20 % pool change is simply
too small to see in one preparation. Refitting the 13 knobs to the battery takes the battery
to ≈ 0.01. **Protocol B, which was not used in that fit, stays at 0.1–0.2: it does not
expose the hidden pool.** Fitting battery + B together hides it in both (≈ 0.01).

![protocol B](results/twin_protocolB.png)

The release–hold curves of every pinned and refitted fibre lie inside the truth's
one-preparation band at every hold time.

![mc](results/twin_mc.png)

**(c) The proposed noised-seed test, linearised Monte Carlo** (4 000 draws, one preparation
of measurement noise). Top row, data only: the estimate is unbounded (95 % range
×/4 000 to beyond ×/10⁹). Most draws hit the plot limits, and adding B does not change the
picture. Bottom row, the same fits when the optimiser stays near a seed that is the truth
plus 5 % noise: **the estimate comes back unbiased at ×0.80 / ×1.20 with ×/1.5–1.7**. That
is exactly the "good estimate" the proposed test would have reported, and all of it comes
from the seed. **Do not use a truth-centred noised seed as an identifiability test.** Use
measurement-noise Monte Carlo with free starts, a profile, or the Fisher information.

**Verdict.** In this model, neither the current battery nor the proposed release–hold series
(nor the activation step) can resolve the activated SRX pool to ±20 %, or even to a factor
of 2. A protocol that could would have to load the pool: vary the force at which heads
park, not the unloaded interval. Candidates: isotonic load steps to intermediate forces;
release–hold with partial rather than full unloading; an SL series at fixed activation
where passive load changes σ₂'s effect. Or measure the pool directly (X-ray OFF-state
reflections during contraction). Not tested here.

## 10. Does a mavacamten-type drug make the activated SRX quantifiable? ([`RunSrxDrug.m`](RunSrxDrug.m))

**Question.** Your target is the activated pool, and whether a drug that shifts heads into
SRX (mavacamten) has a quantifiable effect on the model-predicted activated SRX. A drug is
measured *within* a preparation (control vs drug), so it could succeed where absolute
measurements fail. The test uses fibres that fit the control data equally well but carry
different activated pools: the twin refits at ×0.8 and ×1.2 plus the truth, and the §8
real-data refits at ×0.5 and ×2 plus the ×1 control refit. The same drug is applied to
each. Two drug models:

* **D1 recruitment:** SRX inflow ×3 (`ksr0` & `ksrd`), i.e. mavacamten as a stabiliser of
  heads entering SRX from DRX;
* **D2 retention:** SRX outflow ÷3 (`kmsrd`), heads held in SRX, also under load.

![drug srx](results/drug_srx.png)

**The predicted activated SRX.** The *fold change* is robust: D1 **×2.45** (×2.36–2.54
for the ±20 % fibres, ×2.0–2.7 over ×0.5–×2); D2 **×1.5** (×1.40–1.51). The *absolute*
value under drug is not: it inherits the control uncertainty (D1: 0.041–0.125). So "the drug
increases activated SRX by ~2.5×" is a model prediction you can quote. "Activated SRX
under drug is 7 %" is not.

![drug effects](results/drug_effects.png)

![drug discrimination](results/drug_discrim.png)

**Does the drug response identify the baseline pool?** For the twin fibres (identical at
control, χ² 0.14–0.24), the drug effects on every measurable feature coincide: χ² of the
difference **0.1–0.4** in within-preparation units. **The drug does not break the
degeneracy.** The ×0.5/×2 real-data fibres differ more (χ² 1.4–7.2), but they already differ
at control on features they were not fitted to (protocol B and force: χ² 4.5 and 12.6). Most
of that gap is not drug-specific.

**What the drug *does* make quantifiable:**

* **The mechanism.** Both drug models raise the **resting** pool about equally (0.27 → 0.46
  for D1, 0.51 for D2). That is the shift a mant-ATP chase under the drug would report, and
  the two cannot be told apart that way. Under activation they differ strongly: D1 slows
  slack ktr 44 → 36 s⁻¹ and ktr-protocol k₆₃ −1.5 SD, while D2 slows slack ktr to 22 s⁻¹, with
  k₆₃ −3 SD and the release–hold k₆₃ at 10 and 40 ms −2 to −4 SD. **A mant chase (resting
  shift, i.e. drug strength) plus ktr (mechanism) on the same dose pins the drug's action,
  and with it the predicted activated-SRX fold change.**
* **The force prediction is sharp and falsifiable.** Both mechanisms reduce maximal
  isometric force by only **1–3 %**, because under full activation the load-gated escape
  empties the pool however it was filled. If mavacamten turns out to cut maximal
  Ca-activated force substantially in these preparations (a literature check is needed),
  that loss cannot come through this SRX pool. It needs either an SRX exit that is not
  force-accelerated (a Ca/activation-dependent pool, cf. §2(a)) or a direct action on
  cycling heads (e.g. product release). Measuring force *and* ktr under the drug tells the
  two apart.

**Caveats.** Two drug mechanisms, one dose, 8 mM, saturating Ca. The model has no Ca
dependence, so partial activation (where mavacamten acts most strongly) is out of reach.
The mechanistic conclusion is conditional on the load-gated SRX structure of
`dPUdT_CombinedTransitions`.
