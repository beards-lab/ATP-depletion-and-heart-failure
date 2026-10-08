# 1 Hz paced twitch: `twoState` thin filament driven by an analytic Ca transient

**Date:** 2026-09-24; **re-run at SL 2.0 µm on 2026-09-29** (§0). **Branch:** `CaActivation`
(worktree `ATP-wt-ca`). **Script:** [`RunPaced1Hz.m`](RunPaced1Hz.m). **Data:** `results/paced1Hz.mat`
(SL 2.0). The SL 2.2 run is kept as `results/paced1Hz_SL22.mat` and `results/paced1Hz_SL22_*.png`.

Numbers are measured unless marked *literature*. Drug mechanisms in §5 are *literature*,
recalled rather than re-checked for this report. Treat them as pointers to verify, not as
citations.

## 0. Why the diastolic force was so large, and what changed

The first version of this report ran at the snapshot's own SL, 2.2 µm. Its diastolic force was
**8.95 kPa, 98 % of it passive**. That value was wrong for a twitch, for two reasons.

1. **The passive curve is the high-Ca curve.** `ExampleHighATP` computes passive force as
   `k_pas·(SL − LSE − Lsc0)^γ`, with `k_pas` 12104, `Lsc0` 1.80 µm and `γ` 7.82. The curve was
   identified against the passive trace recorded at pCa 4.51 with cross-bridges blocked
   (PNB/mavacamten). In that trace titin is in its **Ca-stiffened** state, about 10 kPa at the
   operating SL. The low-Ca passive on the same preparations is about **3 kPa**
   ([PassiveForceID](../PassiveForceID/conclusions.md)). Subtracting the high-Ca trace is correct
   for the saturating-Ca fits. Diastole, however, is at low Ca. **At SL 2.2 the model was carrying
   roughly 3× the diastolic passive force the preparation actually has.**
2. **The curve is extremely steep near 2.2 µm.** With `γ` ≈ 7.8, the curve gives the following
   (series element neglected):

   | SL (µm) | 1.9 | 2.0 | 2.1 | 2.2 | 2.3 |
   |---|---|---|---|---|---|
   | passive (kPa) | 0.000 | **0.040** | 0.96 | **9.20** | 52.8 |

   The 2.2 µm operating point sits where 0.1 µm changes passive force tenfold.

**Fix: the twitch now runs at SL 2.0 µm** (`SL_TWITCH`, default 2.0, at the top of the script;
any other value writes files with a `_SLxx` suffix). `getParams` rebuilds `PU0` and the strain
pivot at that SL. No parameter was changed.

At SL 2.0 the passive force is **0.040 kPa**. This is probably an *under*estimate: low-Ca passive
at 2.0 µm in skinned mouse myocardium is of order 0.5–1 kPa (*estimate*). But it is now small next to the
developed force, and that was the problem to solve. A Ca-dependent titin (low-Ca curve in
diastole, high-Ca curve at saturating Ca) would be the proper fix. It would need the low-Ca
passive curve identified over SL, which has not been done.

**What else changes between 2.2 and 2.0 µm** (both measured; same parameters otherwise):

| | SL 2.2 | SL 2.0 | why |
|---|---|---|---|
| diastolic force (kPa) | 8.95 (passive 8.76) | **0.19** (passive 0.04) | passive curve, above |
| saturating developed force (kPa) | 70.5 | **64.3** | filament overlap 0.977 → **0.851** |
| twitch developed force (kPa) | 3.47 | **2.67** (−23 %) | overlap (−13 %) plus fewer primed heads (below) |
| diastolic SRX fraction | 0.147 | **0.204** | the SRX exit is force-driven, and the force used is **total** force, passive included |
| primed `PD` at rest | 0.833 | 0.781 | the counterpart of the SRX shift |
| TTP / RT50 / RT90 (ms) | 70.0 / 88.0 / 250 | 69.5 / 90.0 / 255 | shape essentially unchanged |

**The model has a passive-force route to length-dependent activation.** The force-sensitive SRX
exit (`kmsrd·exp(F_SR/σ)`) uses `F_SR = F_total` (`UsePassiveForSR = 0`,
`useHalfActiveForSR = 0`). Passive force therefore pulls heads out of SRX. At SL 2.2 the (too
large) passive tension recruited 5.7 % of heads compared with SL 2.0. That is a Frank–Starling-like
thick-filament mechanism. Its strength currently rests on the high-Ca passive curve, so it is
overstated in diastole.

## Setup

| item | value |
|---|---|
| parameter set | `params/ExampleHighATP.m`: the frozen best, 2-state, 8 mM MgATP; **SL0 overridden to 2.0 µm** |
| boundary | isometric with fixed ends; the series elastic element (SE) is in the loop, so sarcomeres shorten internally |
| Ca model | `CaModel = 'twoState'`: Tewari/Campbell/Marzban cooperative N state with bound-head feedback; Marzban 2020 rat constants (`tf_k_on` 101.2 µM⁻¹s⁻¹, `tf_k_off` 723.9 s⁻¹, `tf_K_coop` 9.68) |
| Ca input | `CaInput = 'analytic'`, `rice2008` (rat, 24 °C: 0.09 → 1.45 µM, rise τ 20 ms, decay τ 110 ms, 5 ms onset), `CaPeriod = 1` s |
| protocol | 2 s at diastolic Ca from an all-detached start, then 6 beats; the last beat is analysed |
| integrator | `ode15s`, RelTol 1e-3, AbsTol 1e-4, MaxStep 1 ms, restarted at every beat onset and at the +5 ms onset kink; output every 0.5 ms |
| cost | 24 s wall for 8 s simulated (2.9 s per beat); 275 s for the whole script including 10 perturbation runs |

Reproduce, from the worktree root:

```matlab
addpath(genpath(pwd)); run('Analyses/CalciumTwitch/RunPaced1Hz.m')          % SL 2.0
SL_TWITCH = 2.2; run('Analyses/CalciumTwitch/RunPaced1Hz.m')                % -> *_SL22.*
```

## Summary (SL 2.0 µm)

| measure | value |
|---|---|
| diastolic force | **0.19 kPa**: 0.04 kPa passive, 0.15 kPa cross-bridge |
| developed twitch force | **2.67 kPa**, i.e. **4.1 %** of saturating-Ca developed force (64.48 − 0.19 = 64.3 kPa) |
| time to peak (TTP) / half-relaxation from peak (RT50) / 90 % relaxation (RT90) | **69.5 / 90.0 / 255 ms** |
| steady force at the transient's peak Ca (1.45 µM, pCa 5.84) | 3.16 kPa (developed 2.98 kPa). The twitch reaches **90 %** of it |
| attachment gate: diastole → peak | 0.0012 → **0.026** |
| attached heads: diastole → peak | 0.0006 → **0.012** |
| SRX (`P_SR + P_SRD`): diastole → minimum | 0.204 → 0.194 (at 92 ms) |
| primed heads `PD` (M·ADP·Pi) at rest | **0.78** |
| ATP use (net hydrolysis) | **0.103 ATP per head per beat**; 0.080 of it in the first 325 ms; resting rate 0.026 per head per s |
| beat-to-beat change | developed force identical from beat 2 (Δ ≤ 1e-4 kPa) |

## 1. Last beat: Ca, force and internal states

![last beat](results/paced1Hz_lastbeat.png)

- **A.** The prescribed Ca transient (input).
- **B.** Total force, which equals the force in the series element and is what a transducer
  measures, split into its cross-bridge and passive parts.
- **C.** Thin filament: the permissive fraction of sites `1−N`, and the attachment gate
  `Pu/(1−bound)` that multiplies `ka_eff`.
- **D.** Myosin head pools.
- **E.** Fluxes, in heads per s: attachment `PD→p1`, detachment `p2→PT`, SRX entry
  `PT→SR` and SRX exit `SRD→PD`.
- **F.** The slow attachment gates `A_reg` and `A_def` (two more gates on `ka_eff`).

**Order of the peaks** (ms after the stimulus): Ca **46.5** → gate **49.5** → permissive
sites **56.0** → attached heads and force **69.5**.

![whole run](results/paced1Hz_overview.png)

The whole run. Force (A), gate (B), and SRX and attached heads (C). The SRX transient in
the first 0.1 s comes from the all-detached start and has settled long before beat 1.

![force-Ca loop](results/paced1Hz_loop.png)

The force–Ca loop of the last beat. It runs counter-clockwise: Ca rises first with little
force (lower branch), force peaks after Ca has begun to fall, and during relaxation force is
higher at a given Ca than it was during the rise.

## 2. Observations

1. **The twitch is small: 4.1 % of saturating developed force.** The cause is Ca
   sensitivity, not kinetics. The steady force at the transient's peak Ca is only 4.6 % of
   saturating developed force, and the twitch reaches 90 % of even that. The Marzban
   constants on this core put the rat 24 °C transient's peak (pCa 5.84) on the foot of the
   force–pCa curve (pCa₅₀ 5.33 at SL 2.2, [conclusions §2.2](conclusions.md)). Intact trabeculae
   develop twitch stresses far above 5 % of the maximal Ca-activated stress (*literature*).
   The thin-filament constants are therefore not calibrated for this core. That was
   expected, and force–pCa data would fix it (conclusions §6.2).

2. **Myosin is not what limits the twitch.** At rest 78 % of heads are primed in `PD` and
   ready to attach, and only 1.2 % attach at the peak. Every lever on head availability is
   therefore weak at this activation level:
   - tripling SRX entry (`ksr0 ×3`) costs only −5 % developed force;
   - tripling SRX exit (`kmsrd ×3`) gains only +13 % (§4).

   Per-head kinetics are strong levers by comparison: `ka ×0.5` gives −52 % and `k2 ×0.5`
   gives +52 %.

3. **Bound heads hold about a third of the permissive sites at the peak** (0.012 of 0.036).
   This is the cooperative cross-bridge feedback at work. It keeps sites on after Ca has started to
   fall (the gate peaks at 49.5 ms, the permissive fraction at 56 ms), and it produces the
   counter-clockwise force–Ca loop. It is also why making heads stay attached longer
   (`k2 ×0.5`) raises the gate peak and increases the twitch disproportionately.

4. **In this regime `k_off` sets amplitude, not relaxation.** Doubling `tf_k_off`:
   - cuts developed force by 57 %;
   - slightly **lengthens** RT50 (90 → 97 ms).

   The thin filament (k_off 724 s⁻¹) is much faster than the Ca decay (τ 110 ms), so N sits
   close to quasi-steady and relaxation follows Ca. Relaxation responds to the Ca side
   instead: decay τ ×0.6 shortens RT50 by 27 % and RT90 by 34 %. This corrects the
   expectation in conclusions §2.3 that `k_off` would mainly govern relaxation. For
   lusitropy, this model puts the lever on Ca removal (SERCA), not on troponin off-rate.

5. **The SRX pool is a fast-exchanging OFF pool, not a mant-ATP-type SRX.** At rest 20 % of
   heads are in SRX with a mean dwell of order 30 ms (32 ms at SL 2.2). Mant-ATP SRX lifetimes are of order 100 s
   (*literature*). During a twitch the force-sensitive exit (`kmsrd·exp(F/σ)`) draws SRX down by
   only 1 percentage point. The diastolic SRX fraction was never identified (it is extrapolated
   from saturating-Ca fits, conclusions §6.5), and §0 shows it also depends on the passive
   curve. Any drug mechanism that works by stabilizing slow SRX (§5) is therefore
   misrepresented in its time scale.

6. **Resting ATP turnover is 0.026 per head per s,** close to the basal DRX turnover in
   mant-ATP chase experiments (about 0.03 s⁻¹, *literature*). The SRX loop
   (PT→SR→SRD→PD→PT) consumes no net ATP here: the SR→SRD hydrolysis is undone by
   PD→PT reversal. About 25 % of the ATP used per beat is this basal turnover. Economy
   (developed force per ATP used) is 26 kPa per (ATP/head/beat).

7. **The series element is now a negligible buffer.** The SE lengthens by 0.6 nm at peak force.
   At SL 2.0 the passive curve is flat, so passive force drops by only 0.001 kPa. Measured and
   cross-bridge developed force agree (2.665 vs 2.666 kPa). At SL 2.2 the same internal shortening
   dropped passive force by 0.20 kPa, 5 % of the twitch. That was an artefact of the high-Ca
   passive curve (§0).

8. **There is no beat-to-beat memory.** Beat 2 already equals beat 6, because no myofilament
   state is slower than about 0.2 s. With a prescribed transient, the model has no
   myofilament route to force–frequency or post-rest potentiation: any such effect must
   come from the Ca input.

9. **The attachment rate is gated three times over.** `A_reg` sits at its isometric floor
   A0 = 0.211 throughout. It rises only under shortening (`RegAvailShorteningOnly`). So
   isometric attachment runs at 0.211 × gate × ka. Both `A_reg` and the Ca gate were fitted
   only at saturating Ca, as conclusions §6.5 warned. **`A_def` dips slightly below zero**
   while bound mass grows, because its smooth rectifier lets small negative values through.
   The effect is tiny and harmless here, but not physical.

## 3. Internal states at a glance (SL 2.0)

| quantity | diastole | peak or minimum in beat | at time (ms) |
|---|---|---|---|
| Ca (µM) | 0.090 | 1.45 | 46.5 |
| attachment gate | 0.0012 | 0.0257 | 49.5 |
| permissive sites `1−N` | 0.0018 | 0.0362 | 56.0 |
| attached heads `p1+p2` | 0.0006 | 0.0115 | 69.5 |
| cross-bridge force (kPa) | 0.145 | 2.81 | 69.5 |
| passive force (kPa) | 0.040 | 0.039 (min) | 69.5 |
| primed `PD` | 0.781 | 0.776 (min) | ~70 |
| SRX | 0.204 | 0.194 (min) | 92.5 |

## 4. Illustrative single-parameter perturbations (SL 2.0)

These are **not drug models.** Each run changes one parameter by an arbitrary factor, to
show the direction and relative strength of each lever. The protocol is 1.5 s diastolic
equilibration and 4 beats, analysing the last. All runs are periodic (beat-to-beat drift
≤ 2e-5 kPa).

![perturbations](results/paced1Hz_perturb.png)

| change | stands for (§5) | F_dia (kPa) | F_dev (kPa) | ΔF_dev | TTP (ms) | RT50 (ms) | RT90 (ms) | gate peak | SRX dia → min | ATP per beat | ΔATP | F_dev / ATP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | — | 0.19 | 2.67 | — | 69.5 | 90.0 | 255 | 0.026 | 0.204 → 0.194 | 0.103 | — | 25.9 |
| `tf_k_on ×3` | TnC Ca sensitizer | 0.49 | 19.98 | **+650 %** | 75.0 | 59.4 | 176 | 0.216 | 0.203 → **0.127** | 0.519 | +404 % | 38.5 |
| `tf_k_off ×2` | TnI phosphorylation (PKA) | 0.11 | 1.14 | −57 % | 68.5 | 96.9 | 269 | 0.011 | 0.205 → 0.200 | 0.047 | −54 % | 24.1 |
| `k1 ×0.5` | slower Pi release / powerstroke | 0.17 | 2.35 | −12 % | 68.5 | 89.8 | 254 | 0.026 | 0.204 → 0.195 | 0.081 | **−21 %** | 28.9 |
| `ka ×0.5` | fewer weak→bound attachments | 0.11 | 1.27 | −52 % | 68.5 | 92.0 | 259 | 0.024 | 0.205 → 0.200 | 0.050 | −51 % | 25.4 |
| `ksr0 ×3` | SRX stabilization | 0.18 | 2.54 | −5 % | 69.0 | 90.5 | 255 | 0.026 | **0.246** → 0.239 | 0.098 | −5 % | 25.9 |
| `kmsrd ×3` | SRX recruitment | 0.20 | 3.02 | +13 % | 69.5 | 89.4 | 254 | 0.026 | **0.106** → 0.100 | 0.116 | +13 % | 26.0 |
| `k2 ×0.5` | slower detachment / ADP release | 0.28 | 4.05 | **+52 %** | **83.0** | 100.6 | 267 | 0.027 | 0.203 → 0.182 | 0.089 | −14 % | **45.7** |
| Ca peak ×1.5 | larger transient | 0.19 | 4.68 | +76 % | 70.0 | 83.2 | 239 | 0.045 | 0.204 → 0.187 | 0.154 | +50 % | 30.4 |
| Ca decay τ ×0.6 | faster Ca removal (SERCA) | 0.18 | 2.55 | −4 % | **60.0** | **65.9** | **168** | 0.025 | 0.204 → 0.195 | 0.079 | −23 % | 32.2 |

The ranking and signs are the same as at SL 2.2 (`results/paced1Hz_SL22_perturb.png`). Every
relative change is within 4 percentage points of its SL 2.2 value, except `tf_k_on ×3`
(+650 % against +727 %). What stands out:

- **Force responds supra-linearly to both the thin-filament and the Ca side,** because the
  twitch sits on the steep foot of the force–pCa curve:
  - `k_on ×3` gives 7.5× the developed force;
  - Ca peak ×1.5 gives 1.8×.

  With `k_on ×3`, SRX empties to 0.127, because the force-sensitive exit recruits heads
  once force is high. Relaxation also becomes *faster* (RT50 59 ms): with a steeper
  operating point, force falls off sharply as Ca decays. This is the opposite of the
  slowed relaxation usually attributed to sensitizers, and it follows from observation 4.
- **Slower detachment gives more force for less ATP.** `k2 ×0.5` adds 52 % force, uses
  14 % less ATP, and prolongs the time to peak by 13.5 ms. Longer duty time feeds the
  bound-head feedback. Economy improves by 76 %.
- **Slowing Pi release saves more ATP than force** (`k1 ×0.5`: −12 % force, −21 % ATP).
  Halving attachment instead (`ka ×0.5`) scales force and ATP down together.
- **SRX levers are weak** at this activation level (observation 2).

## 5. Which model parts common drugs would act on

The mechanism column is *literature*, recalled rather than re-checked for this report; the
model columns are this model's parameters. "Not representable" marks what the current
model cannot express. The Beard 2022 and Baker sub-maximal pCa sets were recorded with
cross-bridge inhibitors (conclusions §3.4), so blebbistatin and mavacamten are also data
conditions here, not only therapies. Perturbation numbers are the SL 2.0 values of §4.

| drug | class | primary mechanism (*literature*) | model layer | parameters it maps to | what the perturbations suggest | not representable now |
|---|---|---|---|---|---|---|
| **Mavacamten** | cardiac myosin inhibitor (HCM) | slows Pi release, so less actin-activated ATPase; stabilizes the SRX / interacting-heads (IHM) OFF state; lowers maximal force and Ca sensitivity in skinned muscle | myosin cycle + thick filament | `k1`↓, `ksr0`↑ and/or `kmsrd`↓, possibly `ka`↓ | the `k1` part saves ATP out of proportion to force; the SRX part barely moves a low-activation twitch (−5 % for `ksr0 ×3`) | slow-SRX time scale (model dwell ~30 ms against ~100 s); drug binding occupancy and kinetics |
| **Aficamten** | cardiac myosin inhibitor (HCM), a different allosteric site | slows Pi release and stabilizes a weak actin-binding (pre-powerstroke) state; the SRX contribution is reported smaller and is debated | myosin cycle | `k1`↓, `ka`↓ (weak→bound), perhaps `kd`↑ | force and ATP fall together (`ka ×0.5`: −52 % / −51 %); via `k1` alone, economy improves | a separate weak-bound state (the model's `p1` stands in for it) |
| **Omecamtiv mecarbil** | cardiac myosin activator | speeds entry into strongly bound states but slows the working stroke and ADP release, so heads stay attached longer; raises Ca sensitivity; lowers ktr; prolongs systolic ejection | myosin cycle, and the thin filament indirectly through bound-head feedback | `k2`↓, `dr`↓ (stroke size), `ka`/`k1`↑, `kmsrd`↑ | `k2 ×0.5` reproduces the phenotype qualitatively: +52 % force, TTP +13.5 ms, better economy, and a higher gate from bound-head feedback. The model makes OM's Ca sensitization emerge from `Joff ∝ Pu`, with no change in TnC constants | the reduced step size needs the `dr` × `kstiff` tradeoff; not exercised here |
| **Danicamtiv** | cardiac myosin activator | recruits heads from the OFF state; raises Ca sensitivity with less prolongation than OM | thick filament + myosin cycle | `kmsrd`↑, `ka`↑ | weak (`kmsrd ×3`: +13 %) because 78 % of heads are already primed | slow-SRX time scale (as for mavacamten) |
| **Levosimendan** (also pimobendan) | Ca sensitizer + PDE3 inhibitor | stabilizes Ca-bound cTnC; PDE3 inhibition raises cAMP (Ca input side) at higher doses | thin filament (+ Ca input) | `tf_k_on`↑ / `tf_k_off`↓; Ca peak↑, τ↓ | by far the strongest lever (`k_on ×3`: 7.5×), with the side effect of slightly higher diastolic force (+0.30 kPa) | extra TnC Ca buffering reshaping the transient (the input is prescribed) |
| **β-agonists** (isoproterenol, dobutamine), PDE3 inhibitors (milrinone) | PKA activation | TnI Ser23/24 phosphorylation (faster Ca dissociation); MyBP-C phosphorylation (heads released from OFF); phospholamban relief (faster SERCA); larger L-type Ca current; titin phosphorylation (less passive stiffness) | all layers | `tf_k_off`↑, `kmsrd`↑, Ca τ↓ and peak↑, `k_pas`↓ | the pieces pull opposite ways: TnI −57 % force, larger Ca +76 %, faster SERCA gives relaxation (RT50 −27 %). Only the Ca side speeds relaxation (observation 4) | rate and heart-rate coupling; only one fixed frequency |
| **Digoxin** (also istaroxime) | Na/K-ATPase inhibitor (istaroxime also activates SERCA) | more SR Ca load, so a larger transient and slightly higher diastolic Ca | Ca input | `Camax`↑, `Cadia`↑ (istaroxime also τ↓) | Ca peak ×1.5 gives +76 % force and +50 % ATP | — (fully input-side) |
| **L-type Ca channel blockers** (verapamil, diltiazem, nifedipine) | Ca-entry block | smaller transient | Ca input | `Camax`↓ | the mirror of Ca peak ×1.5: supra-linear loss of force | — |
| **Blebbistatin / para-nitroblebbistatin** | myosin ATPase inhibitor (research tool; used in `data/`) | traps myosin as M·ADP·Pi with low actin affinity | myosin cycle | `ka`↓ (strongly) | `ka ×0.5`: −52 % force and ATP | — |
| *(context)* **low MgATP / ATP depletion** | metabolic | slower ATP-dependent detachment; the project's main question | myosin cycle (+ SERCA in reality) | `k2`↓ (and any ATP gate on SRX exit, `UseAtpKmsrd`) | like `k2 ×0.5`: more twitch force, longer TTP | ATP-dependent SERCA (Tran 2009): at low ATP a real Ca transient would decay more slowly |

## 6. Caveats

- **Species and temperature are mixed:** a mouse skinned-trabecula core (rig temperature
  unknown), rat 24 °C thin-filament constants, and a rat 24 °C transient. Absolute twitch
  amplitude and timing are not physiological predictions. Directions and relative lever
  strengths are the useful output.
- **Passive force is still the high-Ca curve** (§0). At SL 2.0 it is too small to matter, but
  any twitch run above about 2.1 µm needs a low-Ca passive curve first.
- **The Ca input is prescribed:** there is no TnC buffering feedback, and drugs can change
  the transient only through the parameters shown.
- **The perturbation factors are arbitrary.** They show direction and relative strength, not
  dose.
- **SRX kinetics are fast** in this parameter set (observation 5), so any statement about
  SRX-targeting drugs is limited by that.

## Files

| file | content |
|---|---|
| [`RunPaced1Hz.m`](RunPaced1Hz.m) | simulation, references, perturbations, figures; `SL_TWITCH` selects the SL (default 2.0) |
| `results/paced1Hz.mat` | SL 2.0: time courses `B`, per-beat metrics, perturbation results `PV` |
| `results/paced1Hz_lastbeat.png`, `_overview.png`, `_loop.png`, `_perturb.png` | SL 2.0 figures above |
| `results/paced1Hz_SL22.mat`, `results/paced1Hz_SL22_*.png` | the original SL 2.2 run (2026-09-24) |
