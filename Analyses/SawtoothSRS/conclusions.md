# Porting the skeletal short-range-stiffness sawtooth to cardiac

**Can the rat-soleus SRS protocol be run on this cardiac model, what does it predict, and
does the answer survive the translation?**

Repo: `ATP-depletion-and-heart-failure` · Model: `dPUdT_CombinedTransitions` ·
Params: `params/estiff125_opt.m` (current best 8 mM) · Date: 2026-09-08

Source protocol: van der Zee TJ, Simha SN, Milburn GN, Campbell KS, Ting LH, De Groote F.
*Cross-bridge model for predicting muscle short-range stiffness during movement.*
PLOS Comput Biol 2026, [10.1371/journal.pcbi.1014748](https://doi.org/10.1371/journal.pcbi.1014748).

---

## 0. TL;DR

The protocol ports mechanically but **not dynamically**. The published numbers (4 % L₀ in
100 ms) are almost the same *normalised velocity* as the skeletal experiment (0.40 vs
0.45 L₀/s), so the verbatim port reproduces their measurement window faithfully — but
cardiac is a **5.8× faster muscle by V_max and >20× faster by k_tr**, so the same input
lands at **5.3 % V_max where they sat at 35 %**. Matching the relative dynamics instead
requires a **15.2 ms ramp at 2.625 ML/s**, and that changes essentially every outcome:
peak force 1.17 → 1.41 F₀, hysteresis 0.29 → 1.87 kPa·ML, elastic return 90 % → 43 %.

Two things worth a lab's attention came out of it, neither of them about SRS:

1. **A symmetric length input produces an asymmetric force response, and at slow speeds
   the model delivers net work *out* over the cycle** (hysteresis −0.252 at 0.1 ML/s).
2. **Force overshoots the isometric level by up to 22 % *after* the muscle is back at L₀**
   (1.223 F₀ at 0.1 ML/s; 1.177 F₀ at 32 ms in the verbatim case).

3. **On a conditioning/test pair the model gets the paper's sign wrong**: the second ramp's
   onset stiffness is **enhanced up to 2×** at 5 ms recovery, crossing to depression only
   at ~40 ms (minimum 0.782 at 70 ms) and gone by 350 ms (§6).

The first two are **entirely** attributable to directional registration availability
(`RegAvailShorteningOnly`) — ablating it drives the overshoot to exactly 1.000 at every
speed. That mechanism was introduced to fix the *steady-state* F–V shoulder
([FV_Shoulder](../FV_Shoulder/conclusions.md)), so the sawtooth is a cheap, nearly
independent test of it. **This is the recommended experiment.**

---

## 1. The translation

| | rat soleus (source) | cardiac (this model) |
|---|---|---|
| reference length | L₀ ↔ SL **2.60 µm** | L₀ = 1.0 ML ↔ SL **2.0 µm** (`ML` = 2 µm) |
| amplitude | 3.83 % L₀ | **4 % L₀** = 0.04 ML = **0.08 µm SL = 40 nm/HS** (SL 2.0 → 2.08) |
| ramp rate | 0.45 L₀/s | see below |
| V_max | ~1.3 L₀/s (stated) | **7.50 ML/s** (measured, §2) |
| k_tr | **< 2 s⁻¹** (stated) | **40.3 s⁻¹** (measured, §2) |
| SRS definition | regression slope of F/F₀ on length over the **first 10 ms** | same recipe, plus a strain-matched variant |

SL 2.0 µm is the reference because every other runner in this repo hardcodes it as its
isometric reference (`runStairsExperiment`, `runKtrExperiment`, the F–V isometric point).

**The two cases.** Amplitude is identical; only the ramp rate differs.

| case | v (ML/s) | ramp | v/V_max | rationale |
|---|---|---|---|---|
| **A — unmatched / verbatim** | ±0.400 | 100.0 ms | **0.053** | the published numbers applied as-is |
| **B — kinetically matched** | ±2.625 | 15.24 ms | **0.350** | the paper's *own* normalisation (0.45 L₀/s ≡ 35 % V_max) |
| *(turnover-matched)* | *8.06* | *5.0 ms* | *1.07* | *same number of XB turnovers per ramp — **not reachable**, it exceeds V_max* |

The third row is the sharpest statement of the mismatch: to make the ramp span the same
fraction of a cross-bridge cycle as it does in soleus, cardiac would have to be stretched
**faster than its own unloaded shortening velocity**. There is no ramp speed at which this
protocol is simultaneously amplitude-matched, velocity-matched and turnover-matched.

## 2. Model constants measured for the port (`RunCharacterize.m`)

- **F₀ = 65.73 kPa** at SL 2.0 µm (65.42 on the F–V grid — a 0.5 % grid sensitivity;
  independently corroborated by `simulateForceLengthEstim`'s 65.5).
- **V_max = 7.50 ML/s** (15.0 µm/s SL), where *total* force reaches zero.
- **k_tr = 40.3 s⁻¹** (τ = 24.8 ms), via `ktrProbe`.
- The **F–V shoulder is present**: F(−0.25 ML/s) = 67.83 kPa **exceeds** isometric 65.42.
  This is the fingerprint of the availability mechanism and it is what §5 is about.

> ⚠️ **V_max is the soft number here.** Active force *asymptotes* to ~0.2 kPa rather than
> crossing zero, so 7.50 ML/s is where the serial element goes slack, and it sits **outside
> the fitted F–V range** (`FV_velocities` stops at 4 ML/s). A 5 %-F₀ criterion gives
> **2.69 ML/s** instead, which would put case B at 0.94 ML/s / 43 ms — a **2.8× spread on
> case B's velocity**. This is precisely why §4 sweeps velocity rather than reporting two
> points: the sweep brackets the ambiguity (the 0.8 ML/s row is the alternative case B).

## 3. What the two cases predict

| | **A verbatim** (0.400 ML/s, 100 ms) | **B kin-matched** (2.625 ML/s, 15.2 ms) |
|---|---|---|
| SRS, 10 ms window | **32.62** F₀/L₀ | **8.60** ← meaningless, see below |
| SRS, strain-matched | **31.10** F₀/L₀ | **37.15** F₀/L₀ |
| peak F/F₀ | 1.171 at 19.5 ms (**7.8 nm/HS**) | 1.407 at 5.98 ms (**15.7 nm/HS**) |
| F/F₀ at apex | 1.143 | 1.081 |
| min F/F₀ during release | 0.980 | **0.245** |
| **post-release overshoot** | **1.177 F₀** at 32.3 ms | 1.040 F₀ at 58.7 ms |
| hysteresis (kPa·ML) | 0.288 (**90 %** elastic return) | 1.866 (**43 %**) |
| attached: iso → apex → back at L₀ | 0.250 → 0.225 → **0.391** | 0.250 → 0.129 → 0.210 |
| passive share of apex force | 1.7 % | 7.2 % |
| RFE after 1 s | −0.000 | +0.000 |

**Case A is a nearly-elastic, oscillatory response**; case B is a large yield followed by a
near-slack release. In A the force *oscillates* during the 100 ms ramp — peak 1.171 at
19.5 ms, trough 1.06 at ~45 ms, recovering to 1.14 by 75 ms — a ~55 ms period against
τ_ktr = 25 ms, i.e. strain-driven detachment followed by reattachment at a refreshed
register. That oscillation is invisible in case B because the ramp is over before it starts.

**Passive force is negligible** (1.7 % of apex force in case A). At this amplitude the
protocol reads cross-bridges, not titin.

## 4. The 10 ms SRS window does not survive the port (two ways)

**(a) It stops being "short range."** At 0.45 L₀/s the skeletal 10 ms window covers
5.9 nm/HS. At 0.400 ML/s the cardiac verbatim window covers **4.0 nm/HS** — faithful. At
2.625 ML/s it covers **26.2 nm/HS**, which is **1.7× past that case's own yield point**
(15.7 nm/HS). The number it returns (8.60) is a chord drawn across the yield and down the
far side, not a stiffness — see the orange dashed line in `results/traces.png`.

The two definitions cross at v = 0.45 ML/s **by construction** (that is where 10 ms and
0.0045 L₀ are the same excursion), which is 6 % V_max. So:

> **Use the strain-matched window whenever the ramp speed changes.** The literal 10 ms
> recipe is only transferable at ~0.4–0.45 ML/s, which the verbatim port happens to hit.

The velocity sweep shows how badly the two diverge (`results/sweep.png`):

| v (ML/s) | % V_max | SRS 10 ms | SRS strain-matched | peak F/F₀ | yield (nm/HS) | hysteresis | post-release max |
|---|---|---|---|---|---|---|---|
| 0.10 | 1.3 | 53.39 | 21.29 | 1.159 | 39.6 *(= apex, no yield)* | **−0.252** | **1.223** |
| 0.20 | 2.7 | 41.01 | 27.22 | 1.162 | 39.9 *(= apex, no yield)* | **−0.081** | 1.216 |
| **0.40** | **5.3** | **32.62** | **31.10** | 1.171 | **7.8** | +0.288 | 1.177 |
| 0.80 | 10.7 | 27.60 | 33.85 | 1.233 | 10.0 | +0.766 | 1.137 |
| 1.50 | 20.0 | 21.06 | 35.74 | 1.313 | 12.9 | +1.344 | 1.086 |
| **2.62** | **35.0** | 8.60 | **37.15** | 1.407 | 15.7 | +1.866 | 1.040 |
| 4.00 | 53.3 | 3.68 | 37.92 | 1.494 | 18.3 | +2.049 | 1.018 |
| 6.00 | 80.0 | 9.73 | 38.52 | 1.591 | 21.5 | +1.967 | 1.005 |

The strain-matched SRS rises monotonically 21.3 → 38.5 and saturates; the 10 ms SRS is
non-monotonic and is reporting window geometry, not muscle.

**(b) 10 ms is not "instantaneous" for cardiac.** With k_tr = 40.3 s⁻¹ the window spans
**0.40 k_tr time constants**; for soleus (< 2 s⁻¹) the same window spans **< 0.02**. Cardiac
SRS is *inherently* turnover-contaminated — a ≥20× difference in how much cycling happens
inside the measurement. Any cross-species SRS comparison built on a fixed 10 ms window is
comparing different mixtures of elasticity and kinetics.

**The verbatim port sits exactly on the yield threshold.** Yield (a force peak before the
apex) first appears between 0.2 and 0.4 ML/s. At 0.4 ML/s the yield strain is **7.8 nm/HS,
the minimum of the entire sweep**, climbing to 21.5 nm/HS at 6 ML/s. So the verbatim cardiac
protocol lands at the single most sensitive point of the whole velocity axis — a poor choice
of operating point if you want a robust measurement, and a good one if you want a
discriminating one.

## 5. The finding worth an experiment: symmetric input, asymmetric output

Two signatures, both maximal at **slow** ramp speeds — i.e. both accessible with the
verbatim protocol and ordinary equipment:

1. **Negative hysteresis.** At 0.1 and 0.2 ML/s the work loop runs the wrong way round:
   hysteresis −0.252 and −0.081 kPa·ML, elastic return **109 %** and **103 %**. The muscle
   delivers **net work over a symmetric length cycle**. This is not a thermodynamic problem
   (it is an ATP-fed motor), but it is a strong, easily-falsified claim.
2. **Post-release force overshoot.** Force does not relax back to F₀ after the muscle
   returns to L₀ — it **overshoots**, monotonically with slowness: 1.223 F₀ at 0.1 ML/s down
   to 1.005 at 6 ML/s, and **1.177 F₀ peaking 32 ms after the return** in the verbatim case.
   The mechanical cause is visible in the state variables: the release leaves the attached
   fraction at **0.391 vs 0.250 isometric**, and those extra bridges are strained up the
   moment the shortening stops.

**Attribution is demonstrated, not inferred** (`RunMechanismTest.m`). Re-running with
`UseRegistrationAvailability = 0` and nothing else changed:

| v (ML/s) | hysteresis ON → OFF | attached at L₀ ON → OFF | post-release max ON → OFF |
|---|---|---|---|
| 0.100 | **−0.252 → +0.856** | 0.357 → 0.475 | **1.223 → 1.000** |
| 0.200 | **−0.081 → +1.215** | 0.390 → 0.474 | **1.216 → 1.000** |
| 0.400 | +0.288 → +1.544 | 0.391 → 0.454 | **1.177 → 0.999** |
| 2.625 | +1.866 → +2.462 | 0.210 → 0.417 | 1.040 → 0.999 |

With the mechanism off the overshoot is **exactly 1.000** at every speed, the hysteresis is
positive everywhere, and the attached fraction at the end of the release is
indistinguishable from its isometric value (0.454–0.475 vs 0.455) — no memory at all. The
mechanism also *reverses the sign of the force response during the shortening half*: ON,
force **rises** to 1.24 during the release; OFF, it **falls** to 0.885
(`results/mechanism.png`).

> ⚠️ **The ablation is an attribution test, not a rival model.** Switching the flag off
> without refitting raises isometric force **65.73 → 119.33 kPa** and isometric attached
> fraction **0.250 → 0.455**, because depressing isometric attachment is exactly what the
> mechanism is for. Do not read the OFF column as "the model without availability."

Why this is a good experiment: `RegAvailShorteningOnly` was identified from **steady-state**
F–V data. The sawtooth turns it into a **large transient** with an unmistakable sign — a
force overshoot after a return to the starting length, and a work loop that runs backwards
at slow speed. Neither has been used to constrain it.

## 6. First vs second ramp-up: the effect has the *wrong sign* at short recovery

`RunConditioningTest.m` runs the paper's full structure — conditioning sawtooth, isometric
recovery hold, identical test sawtooth — and regresses the onset slope of **both** ramps,
sweeping the recovery interval. Both events see exactly the same length input, so any
difference is history alone. Table: `results/condtest_table.csv`; figure:
`results/condtest.png`.

**Verbatim case** (0.400 ML/s, 100 ms ramp). SRS on the first ramp is 31.10 F₀/L₀ throughout:

| t_rec (ms) | SRS₂ (F₀/L₀) | **SRS₂/SRS₁** | F/F₀ at 2nd onset | attached at 2nd onset | ratio per unit onset force |
|---|---|---|---|---|---|
| 5 | 61.72 | **1.985** | 1.049 | 0.3873 | 1.891 |
| 10 | 54.45 | 1.751 | 1.101 | 0.3803 | 1.591 |
| 20 | 44.17 | 1.420 | 1.156 | 0.3619 | 1.229 |
| 30 | 36.14 | 1.162 | 1.177 | 0.3395 | 0.988 |
| **40** | 30.71 | **0.988** ← crossover | 1.171 | 0.3158 | 0.844 |
| 50 | 27.31 | 0.878 | 1.147 | 0.2932 | 0.766 |
| **70** | 24.31 | **0.782** ← max depression | 1.092 | 0.2627 | 0.716 |
| 85 | 26.84 | 0.863 | 1.051 | 0.2511 | 0.821 |
| 100 | 29.52 | 0.949 | 1.023 | 0.2469 | 0.928 |
| 150 | 30.76 | 0.989 | 1.006 | 0.2490 | 0.983 |
| 200 | 30.98 | 0.996 | 1.001 | 0.2494 | 0.995 |
| ≥ 350 | 31.10 | **1.000** | 1.000 | 0.2498 | 1.000 |

**Kinetically-matched case** (2.625 ML/s, 15.2 ms ramp), SRS₁ = 37.15 F₀/L₀:
ratio **1.378 → 1.294 → 1.208 → 1.079 → 1.005** at t_rec = 5, 10, 20, 50, 100 ms, and
**never dips below 0.999**. It is enhancement-only, with no depression phase at all.

Three things follow.

**(a) The model gets the sign wrong at short recovery.** The paper reports *depression* of
the test-stretch SRS. This model predicts **up to 2× enhancement** (1.985 at t_rec = 5 ms).
Dividing out the elevated onset force does not explain it away (1.891) — it is a genuine
stiffness change. That is a clean, falsifiable disagreement with the source phenomenon.

**(b) There is a depression window, but it is narrow and slow-case-only.** The verbatim case
crosses unity at **~40 ms** and bottoms at **0.782 (−22 %) at 70 ms**, recovering to 1.000 by
350 ms. The kinetically-matched case never depresses. So whether you reproduce the paper's
*direction* depends entirely on a protocol choice the paper does not vary.

**(c) The two effects have different mechanisms.** The enhancement tracks the leftover
attached fraction — 0.387 at 5 ms against 0.250 isometric, i.e. the §5 availability transient
leaving extra bridges bound. The depression does **not**: at its 70 ms minimum the attached
fraction is still **above** isometric (0.2627 vs 0.2498) while the slope is 22 % **below**.
More bridges, less stiffness ⇒ it is the *strain distribution*, not the bridge count — the
surviving bridges sit at strains that contribute little to the onset slope.

**Against the §6 prediction made before this was run:** the "≥ 0.5 s shows nothing" half is
confirmed exactly (1.000). The "any depression lives below ~350 ms" half was right about the
window but missed that the effect **reverses sign** below 40 ms — the single-sawtooth
recovery time constant did not reveal that, and it would not have been guessable from it.

**Reality of the dip** (`RunCondTestVerify.m`): the refined 12-point grid through 20–100 ms
is smooth and monotone into the minimum (no stray point); the ratio moves **−0.38 %** under
64× denser solver output at the dip, and **+0.10 %** at the 5 ms enhancement point.

> The 10 ms window misbehaves here exactly as §4 predicts: in the kinetically-matched case it
> reports ratios up to **4.667**, because the depressed onset force at short recovery moves
> the yield inside the window. Use the strain-matched column.

## 7. Verification (`RunVerify.m`) — why the numbers above are trustworthy

- **Repeatability:** bit-identical across repeat runs (identical time vectors, |ΔF₀| = 0).
- **Output-sampling convergence:** forcing 64× denser solver output (by subdividing the ramp
  into sub-segments of identical velocity) moves `SRS_t10` by **0.10–0.17 %**, `SRS_eps` by
  **0.01–0.19 %**, peak force by **0.015 %**. The SRS windows initially held only 21–23 raw
  solver points, so this needed checking; the slopes are data, not interpolation.
- **Warm start:** the 2 s isometric hold is settled to **0.0000 %** over its final 0.5 s
  (attached and SRX pools flat to 5 decimals). Every ratio above has a valid F₀.
- **One defect found and fixed:** the ramp force peak is *broad*, so a raw argmax over
  solver points was unstable — perturbing the commanded velocity by **one ulp**
  (`FRAC/0.100` vs literal `0.4`) moved the reported peak time 0.5 ms while leaving peak
  amplitude identical to 4 digits. Peak time is now taken from a uniform resample, and is
  determined to **±0.4 ms**; peak *amplitude* to < 0.02 %. Do not quote peak timings tighter
  than that.

## 8. Limitations

- **8 mM only.** The 2 mM ATP gate code lives only in the `ATP-wt-lowatp` worktree; the main
  tree would silently run 8 mM. The obvious extension — is the §5 signature ATP-sensitive? —
  needs that worktree.
- **Prediction, not validation.** No cardiac sawtooth data exists. The protocol files written
  to `data/` carry the *model's own* output as their `datatable`, which makes them a
  round-trip self-test of the plumbing, not a target.
- **Single sawtooth**, per the brief — no conditioning/test pair, so §6 is inference from the
  single-event recovery rather than a measurement.
- **V_max is definition-sensitive and extrapolated** beyond the fitted F–V range (§2). Case B
  is uncertain by ~2.8× in velocity; the sweep is the mitigation.
- SRS is reported in **F₀/L₀**. The paper does not state its units explicitly, so absolute
  cross-comparison of SRS magnitudes with theirs is not attempted here — only the internal
  velocity dependence and the window-validity argument.

## 9. Files

| file | what it does |
|---|---|
| `buildSawtoothProtocol.m` | the velocity table: warm-start → stretch → release → recovery. Ramp duration is *derived* from amplitude and speed |
| `extractSawtoothAttributes.m` | SRS (both windows), yield, work loop, post-release overshoot, state snapshots |
| `RunCharacterize.m` | measures V_max and k_tr, derives the two cases → `results/characterize.mat` |
| `RunSawtoothCardiac.m` | **main**: both cases + velocity sweep + figures → `results/sawtooth.mat`, `traces.png`, `workloops.png`, `sweep.png` |
| `RunVerify.m` | repeatability, sampling convergence, warm-start convergence → `results/verify.mat` |
| `RunMechanismTest.m` | availability ON/OFF attribution → `results/mechanism.mat`, `mechanism.png` |
| `buildConditioningTestProtocol.m` | two identical sawtooths separated by a recovery hold — the paper's full conditioning/test structure |
| `RunConditioningTest.m` | **§6**: onset slope of the 1st and 2nd ramp-up vs recovery interval → `results/condtest.mat`, `condtest_table.csv`, `condtest.png` |
| `RunCondTestVerify.m` | is the 40–150 ms depression real? refined grid + convergence → `results/condtest_verify.mat` |

Generated protocol files (gitignored, in `data/`):
`protocol_sawtooth_cardiac_verbatim.mat`, `protocol_sawtooth_cardiac_kinmatched.mat` —
each holds `velocitytable`, `datatable`, `protocolMeta`, normalised exactly as
`runStairsExperiment` expects.

**The round trip is verified lossless.** Setting
`params0.stairs_velocitytableonfile` to either file and calling `runStairsExperiment`
reproduces the trace exactly: `E_stairs = 0.000e+00`, `max|ΔF/F₀| = 0.000e+00`. That is
expected (the runner hardcodes the same `SL0 = 2.0`, `Slim_l = 1.4`, `Slim_r = 2.3` grid this
analysis uses) and it proves the build → save → load → run plumbing adds nothing, so any
residual against *real* data later will be physics.

## 10. Next

1. **Run the slow sawtooth in the lab** (0.4 ML/s, 4 % L₀, SL 2.0 µm). One trace tests two
   independent predictions: post-release overshoot ≈ 1.18 F₀ at ~30 ms, and a work loop that
   is nearly closed (90 % elastic return) rather than dissipative.
2. **Add the 0.1–0.2 ML/s point** — that is where the model claims *net work output* over a
   symmetric cycle, the most falsifiable statement in this analysis.
3. ~~Conditioning/test pair~~ — **done (§6)**, and it is the sharpest disagreement with the
   source paper in this analysis: the model predicts **enhancement** where they report
   depression, crossing over only at ~40 ms. The experiment to run is a **recovery-interval
   sweep from 5 to 350 ms** at 0.4 ML/s; the sign of SRS₂/SRS₁ below 40 ms decides it, and
   the model commits to ≈2× at 5 ms.
4. **Repeat at 2 mM** in the `ATP-wt-lowatp` worktree. Availability is a sliding-speed
   mechanism and V_max is ATP-sensitive, so §5 should scale with ATP — an independent
   handle on the ATP-depletion question.
