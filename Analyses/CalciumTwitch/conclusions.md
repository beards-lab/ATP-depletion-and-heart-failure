# Adding calcium handling to simulate twitches: complexity, driving signals, models

**Status:** scoping plus a measured probe, 2026-09-23. **No model code was changed.** The probe
wraps the production RHS (`dPUdT_CaGate.m` → `dPUdT_CombinedTransitions`). With the gate held
at 1 it reproduces the core bit-exactly (max |Δf| = 0, max |Δoutputs| = 0). Numbers are
measured on this machine (R2023a, `params/fastgate_best.m`, the set TriSegCoupling used,
8 mM MgATP, isometric at SL 2.2 µm with fixed ends, so the series elastic element is in the
loop) unless marked *literature* or *estimate*.

Reproduce: `RUN = [0 1 2 3 5 6 4]; run('Analyses/CalciumTwitch/RunTwitchProbe.m')` (~4 min).
Figure: [`results/twitch_probe.png`](results/twitch_probe.png).

---

## Verdict

**The engineering is small. The identification is not.**

| | |
|---|---|
| **The cross-bridge core is fast enough for a twitch.** | With activation clamped on from rest, force reaches 50 % of steady state in **7.9 ms**. A **10 ms** activation pulse develops **56 %** of steady force, and relaxation after switch-off has **RT50 10.7 ms**. Force tracks activation. Twitch *shape* will come from the Ca/thin-filament layer, not from the calibrated rates. §2.1 |
| **The insertion point already exists, and it is not the `NP` slot.** | Uncommenting the legacy `NP` branch would gate the wrong transition (§1.2). The correct form is a factorized scalar gate on `ka_eff`, the same structure `A_reg` already uses. That form is ~30 lines, and at pCa 4.51 it keeps **98.5 %** of the calibrated force, so the existing campaign survives it. §1.3, §2.2 |
| **Literature thin-filament constants do not transfer.** | The Marzban 2020 rat constants on this core give **pCa₅₀ 5.33, n_H 3.33**. A physiological 1.5–2 µM transient then opens only **3–4 %** of the gate, and the twitch comes out at **5–6 % of F_ss**. §2.2 |
| **Twitch amplitude is one steep knob that no current data constrains.** | Multiplying `k_on` by **3** moves the twitch from 5 % to **37 %** of F_ss. By **10** it reaches **90 %**, but the twitch goes flat-topped, TTP and RT50 double, and diastolic force rises 1.8 kPa. All of this project's data is at pCa 4.51, and the sub-maximal pCa sets in `data/` are cross-bridge-inhibited (PNB/mavacamten). §2.3, §3.4 |

**Cost is a non-issue for isometric twitches:** 2.0–2.9 s wall clock per Ca-driven beat, against
TriSegCoupling's 119 s for a moving sarcomere. At fixed ends the strain window never slides.

**Recommended order** (§6):

1. Promote the probe's gate into the RHS behind `UseCa` (a day or two).
2. Measure steady **force–pCa at 8 and 2 mM ATP** on the same rig. This is the cheapest data
   that pins the steepest knob, and the gate makes a falsifiable prediction about it (§6).
3. Choose the Ca driving signal by the scientific question, not by availability (§3).

---

## 1. What already exists in this repo

### 1.1 Remnants and infrastructure

| item | where | state |
|---|---|---|
| `NP` state slot | `PU(Ns*ss+2)`, [dPUdT_CombinedTransitions.m:65-66](../../Model/dPUdT_CombinedTransitions.m#L65) | read as a hard-coded `NP = 0`; `dNP = 0` at [:840-841](../../Model/dPUdT_CombinedTransitions.m#L840); the Ca branch is commented out |
| legacy Ca kinetics | [dPUdTCa.m:322-326](../../Model/dPUdTCa.m#L322) | Tewari/Marzban `Jon`/`Joff` on `NP`, with `Ca_i = params.Ca` held constant |
| `UseCa`, `Ca = 1000` | [getParams.m:81-84](../../Model/getParams.m#L81) | `UseCa = false` in every snapshot; `params.Ca` is not read by the active RHS |
| bounds for `K_coop`, `k_on`, `k_off` | [parameterBounds.m:438-452](../../Model/parameterBounds.m#L438) | cited (Davis 2007, Tikunova & Davis 2004, Dobesh 2002); gated behind `UseCa` in [assertActiveParams.m:147-151](../../Model/assertActiveParams.m#L147) |
| measured driving signals | [preprocessDrivingSignal.m](../../DataCuration/preprocessDrivingSignal.m) → `sig` struct → [sigLookup.m](../../Auxiliary/sigLookup.m) | one slot, `params.driveSig`, **already used for measured passive force** ([getParams.m:667-671](../../Model/getParams.m#L667)) |

### 1.2 Why the legacy `NP` branch cannot simply be uncommented

In `dPUdTCa.m` there is one head pool: `P = (1 − NP) − bound`, and attachment draws on `P`. The
current RHS splits unattached heads into five pools: `PT` (M·ATP), `PD` (M·ADP·Pi), `P_SR`,
`P_SRD` (SRX) and `D0`. **Attachment draws only on `PD`**
(`RD1 = ka_eff*PD*N_overlap`, [:379](../../Model/dPUdT_CombinedTransitions.m#L379)), but `NP`
is subtracted only from `PT` ([:222](../../Model/dPUdT_CombinedTransitions.m#L222)). Turning
`NP` back on would therefore:

- gate the **hydrolysis bookkeeping** (`PT`) and not attachment;
- put thin-filament state (`NP`) and thick-filament state (`P_SR`, `P_SRD`) into **one simplex**.
  A head in SRX facing a switched-off actin site then has no consistent place to be.

### 1.3 The correct insertion point is already in the code

`A_reg` ([:368-373](../../Model/dPUdT_CombinedTransitions.m#L368)) is a scalar state that
multiplies `ka_eff`. Ca activation needs the same structure. This is the factorized
MyoSim/Campbell form, in which the thin filament and the myosin states are separate
populations:

```
Jon  = k_on·Ca·N·(1 + K_coop·(1−N))          N = non-permissive fraction of sites
Joff = k_off·Pu·(1 + K_coop·N),   Pu = max(0, 1 − N − bound)
gate = Pu / (1 − bound)          ->  ka_eff *= gate
```

Sites held by bound heads cannot switch off, which is the cooperative cross-bridge feedback.
At `N = 0` the gate is exactly 1. [`dPUdT_CaGate.m`](dPUdT_CaGate.m) implements this as a
wrapper. Promoting it into the RHS means reusing the `NP` slot, dropping `NP` from the `PT`
formula, and multiplying `ka_eff` by the gate. With `UseCa = 0` nothing changes.

### 1.4 A latent defect for any measured Ca signal

`sigLookup` returns a **nearest-sample** value and then `return`s
([sigLookup.m:16-18](../../Auxiliary/sigLookup.m#L16), last changed in `a811fc0`, 2026-04-27). The Catmull-Rom
Hermite interpolation below that line is dead code. A measured Ca(t) passed through it becomes
a 1 kHz staircase, and every step is a derivative discontinuity that `ode15s` must resolve.
That is harmless for the passive-force signal it serves today, because that signal enters the
force balance only. It is not harmless for an input that drives a stiff kinetic state. A Ca
signal needs its own slot and the Hermite path restored.

---

## 2. Measured probe ([`RunTwitchProbe.m`](RunTwitchProbe.m))

### 2.1 Activation clamp: the bandwidth of the cross-bridge core (E1)

The gate is prescribed directly (0 → 1 → 0), with no Ca model, starting from diastolic
equilibrium (gate 0 for 2 s, drift over the last 0.5 s 5×10⁻¹²). This bounds what *any* Ca
model can get out of the current core.

| protocol | developed force | timing |
|---|---|---|
| rest | F = 9.15 kPa (passive, borne by the SE) | SRX pool (`P_SR+P_SRD`) = **0.139** |
| gate step on, 1 s | F_ss − F₀ = **69.7 kPa** | t10 **1.0 ms**, t50 **7.9 ms**, t90 **48.1 ms**; SRX 0.139 → 0.022; bound fraction 0.282 |
| gate step off from steady state | → rest | RT50 **10.7 ms**, RT90 **18.6 ms** |
| 10 ms pulse | **56 %** of F_ss | TTP 10.5 ms, RT50 16.1 ms, SRX min 0.064 |
| 20 ms pulse | **72 %** | TTP 20.0 ms, RT50 16.1 ms |
| 50 ms pulse | **91 %** | TTP 50.0 ms, RT50 13.5 ms |
| 100 ms pulse | **98 %** | TTP 100 ms, RT50 11.5 ms |

- **Force peaks exactly when the gate closes.** The core is fast enough to follow any
  physiological activation time course. This includes a mouse twitch at 37 °C (*literature*:
  the Marzban 6 Hz rat transient peaks at 22.5 ms and decays with RT50 20.8 ms, §3.1).
- The rise from rest (t50 7.9 ms) is much faster than `ktr` (40 s⁻¹, τ ≈ 25 ms,
  [SawtoothSRS](../SawtoothSRS/conclusions.md)). This is physically right: at rest the heads are
  primed in `PD`, whereas a `ktr` release–restretch first redistributes them. **Twitch rise and
  `ktr` are different observables, and a twitch fit will not simply inherit `ktr`.**
- **The SRX pool takes part in the twitch.** A 10 ms activation already pulls SRX from 0.139 to
  0.064, through the force-sensitive exit `kmsrd·exp(F/σ)` (31 s⁻¹ at F = 0). This is the
  mechanism Mijailovich 2021 found necessary for rat twitches (*literature*, §4). The
  diastolic SRX fraction (0.139) is an **extrapolation**: SRX was identified only at saturating
  Ca, where force is never low. Compare it with mant-ATP SRX data (`../mantATP-SRX`).

### 2.2 Backward compatibility and force–pCa with literature constants (E2)

Thin filament: Marzban 2020 rat (`K_coop 9.6846`, `k_on 101.185 µM⁻¹s⁻¹`, `k_off 723.852 s⁻¹`,
from `Rat-Cardiac-Energetic/dXdT_cardiovascular_mechanics.m:75-77`; the file itself notes
`k_off` was "manually tuned"). Each pCa was run for 1.5 s from rest. Drift over the last
100 ms was ≤ 10⁻⁶ kPa.

| pCa | 7.0 | 6.5 | 6.2 | 6.0 | 5.8 | 5.6 | 5.4 | 5.0 | 4.51 |
|---|---|---|---|---|---|---|---|---|---|
| F (kPa) | 9.34 | 9.77 | 10.46 | 11.41 | 13.31 | 17.76 | 33.86 | 74.43 | 77.83 |
| gate | 0.001 | 0.004 | 0.009 | 0.016 | 0.031 | 0.066 | 0.214 | 0.882 | 0.975 |

- **At pCa 4.51 (the activating solution of Beard 2022, 30.9 µM), developed force is 98.5 %
  of the gate = 1 model.** Adding the gate costs the existing calibration 1.5 % of force. That
  can be absorbed by normalisation or left as is.
- Hill fit: **pCa₅₀ 5.33, n_H 3.33**. There is a shallow foot at pCa 6.2–5.8 that a single
  Hill curve misses. Constants identified with a different cross-bridge model need not give the same
  force–pCa on this core, because the cooperative feedback runs through *this* model's bound
  fraction. **Treat thin-filament parameters as non-portable between cross-bridge models.**

### 2.3 Ca-driven twitches (E3, E5)

Isometric, diastolic pre-equilibration for 1 s, then 3–6 beats. The last beat is reported.
Beats are periodic by beat 3 (beat-2→3 change in developed force < 3×10⁻⁴ kPa).

| driver | thin filament | F_dia | F_dev | % of F_ss | TTP | RT50 | gate peak | wall / beat |
|---|---|---|---|---|---|---|---|---|
| Rice 2008, 1 Hz | Marzban | 9.32 | 3.21 | **5 %** | 66.5 ms | 87.3 ms | 0.026 | 2.9 s |
| Marzban 2020, 6 Hz | Marzban | 9.47 | 3.93 | **6 %** | 35.0 ms | 26.6 ms | 0.043 | 2.0 s |
| Rice 2008, 1 Hz | `k_on` ×3 | 9.67 | 25.71 | **37 %** | 74.0 ms | 52.2 ms | 0.250 | |
| Rice 2008, 1 Hz | `k_on` ×10 | **11.15** | 63.02 | **90 %** | 107.5 ms | 115.4 ms | 0.935 | |

- With literature constants, **a physiological transient barely activates this core**. The
  transient peaks at pCa ≈ 5.7–5.8, on the foot of the curve in §2.2.
- **Amplitude is super-linear in one ratio.** `k_on/k_off` ×3 (Δlog = 0.48) gives **8.0×** more
  developed force. ×10 saturates. The ×10 twitch is **flat-topped**, relaxes late (RT50
  115 ms), and leaves **1.8 kPa of diastolic force** at 1 Hz. So scaling sensitivity alone
  cannot deliver both amplitude and shape. Amplitude is set by `k_on/k_off`; relaxation should be governed mostly
  by `k_off` and steepness by `K_coop` (expected, not scanned here). Each needs its own data (§3.4).

### 2.4 Numerical check (E6)

The 10 ms pulse and the Marzban 6 Hz twitch were re-run at `RelTol 1e-5`, `AbsTol 1e-6`,
`MaxStep 0.2 ms` on a 0.1 ms grid. Developed force changed by **0.1 %** and 0.03 %, and
TTP/RT50 by **≤ 0.2 ms**. The probe tolerances (`RelTol 1e-3`, `MaxStep 1 ms`) are adequate.
The production default (`RelTol 1e-2`, `evaluateModel.m:161`) was not tested for twitches.

---

## 3. Driving signals

A twitch has **two** inputs: the Ca²⁺ signal and the mechanical boundary condition.

### 3.1 Prescribed analytic Ca(t) (feed-forward; no state, no feedback)

All four are implemented in [`caTransient.m`](caTransient.m). Metrics come from E0.

| signal | origin | species, T | onset → peak (µM) | TTP | RT50 | where it already lives |
|---|---|---|---|---|---|---|
| `rice2008` bi-exponential, τ 20/110 ms | Rice 2008 CellML defaults | rat, 24 °C | 0.09 → 1.45 | 46.7 ms | 98.0 ms | CellML repository |
| `marzban2020` log-normal in phase, 2–10 Hz table | Marzban 2020 | rat, in-vivo rates | 0.16 → 2.13 (6 Hz) | 22.5 ms | 20.8 ms | `Rat-Cardiac-Energetic`, `Multiscale-Hx-P21-Model`, `multiscale-ventricular-interactions-model` |
| `tewari_tanh` tanh − exp + Kang restitution | Piacentino 2003 shape, Kang 2016 restitution | human, 37 °C | 0.04 → 0.59 (1 Hz) | 38.5 ms | 253.6 ms | **your** `cross-bridge-model` (`XB.mo` block `Ca2`; `InVitroTwitch/dXdT.m`; write-up `20220228 Calcium transients.docx`) |
| `fibersim2c` two-compartment SR pulse | FiberSim demo defaults | generic | 0.03 → 0.77 | = stimulus (10 ms) | 34.7 ms (τ = 1/k_serca) | `FiberSim` characterize `twitch` |

**Pros:** trivial, smooth, cheap, reproducible. **Cons:** species and temperature are fixed by
the source. None match mouse at the rig's temperature: rat 24 °C is the nearest. There is no
myofilament buffering feedback. And the transient **cannot respond to ATP depletion**, even
though SERCA is ATP-driven.

### 3.2 Measured Ca(t) (fura-2 / indo-1 / fluo-4)

This is exact for its preparation. It requires an intact preparation, indicator calibration
(Kd in-cell, buffering by the dye itself), and matched temperature. The Baker rig is skinned,
so it cannot produce one. The natural external source is intact mouse RV trabeculae at 37 °C
with fura-2 from the Janssen lab (Stull 2002, 2–12 Hz). **Janssen and Campbell are both
co-authors of Beard 2022.** Plumbing: a new `CaSig` slot plus the `sigLookup` fix (§1.4).

### 3.3 Model-generated Ca(t) (EP + Ca handling)

This is the only class that can carry **feedback**: troponin buffering, and SERCA's ATP
dependence. Cost is 40–70+ extra stiff states.

- **Mouse:** Bondarenko 2004 (coupled to Rice 2008 by Mullins & Bondarenko 2013, 25 °C,
  Fortran; no public code stated); Li 2010 (used by Land 2012); Morotti 2014.
- **Rat:** Gattoni 2016 (openCARP `Gattoni16Land17_DynSarc`); Pandit 2001.
- **Human:** ToR-ORd-Land; **T-World** (local clone `../TWorld`: MATLAB + CellML + Modelica,
  Land contraction with two-way TnC buffering, `model_TWorld.m:350,390`).
- **For this project specifically:** Tran 2009's thermodynamic **SERCA** model (CellML) makes
  Ca removal depend on ATP, ADP and Pi. Without something like it, an "ATP depletion → twitch"
  prediction freezes the Ca half of the answer.

### 3.4 Skinned-preparation twitch surrogates (what the Baker rig can supply)

| experiment | what it identifies | notes |
|---|---|---|
| **steady force–pCa**, at 8 and 2 mM ATP | `k_on/k_off` (pCa₅₀), `K_coop` (n_H), and the ATP shift | standard; hours per prep; **not in the dataset** (below) |
| pCa step by solution exchange | only steady states | diffusion-limited over seconds in a ~100–200 µm trabecula |
| caged Ca (NP-EGTA) / caged chelator (diazo-2) photolysis | activation and relaxation rates after a ms Ca jump | needs a UV flash; works on the same skinned trabeculae |
| myofibril rapid solution switching (Poggesi/Tesi/Stehle) | kACT; relaxation = slow linear phase (≈ detachment) + fast exponential | different preparation, same chemistry, ATP set exactly; rigs exist that impose "artificial cardiac cycles" on myofibrils |

**What `data/` actually holds at sub-maximal Ca:** the `Avg*pCa*.csv` sets and
`PassiveCaSrc2/*/…pCa_5.5_PNB…` are **cross-bridge-inhibited** (para-nitroblebbistatin, later
plus mavacamten). They measure Ca-dependent titin, not active force. The 2023-05-18 and
2023-09-19 sessions have runs named without an inhibitor (pCa 4.4/5.5/5.75/6/6.25). Their
logs, however, sit at only **2–6 kPa** at short length even at pCa 4.4, which suggests
inhibition throughout. **Ask the Baker lab.** Some 2024 sessions were also run at **10 °C**
(`…_10C.txt`), so the rig can hold temperature and a Q10 could be measured on the same preps.

### 3.5 The mechanical boundary condition

| mode | status here |
|---|---|
| fixed ends (SE in the loop) | works today; what the probe used |
| sarcomere-length clamp | the rig records SL optically (raw ASI600A exports carry `SL in (um)`); the model would need SL as an input (the old `UseSLInput` path is commented out) |
| afterloaded / work loop | needs a **force-controlled** boundary; the RHS only takes velocity segments (`Vums`) |
| in vivo | TriSeg, see [TriSegCoupling](../TriSegCoupling/conclusions.md) |

---

## 4. Available models and implementations

Only the **thin-filament layer** is needed: the cross-bridge core stays this project's own. So
the question for each model is what its activation block would cost to insert as the gate.

| model | activation layer | species, T | extra states | implementation | fit to this project |
|---|---|---|---|---|---|
| **Tewari 2016 / Campbell 2018 / Marzban 2020** two-state cooperative | N with `k_on`, `k_off`, `K_coop` + bound-head feedback | rat/human | **1** | **local:** `dPUdTCa.m`, `Rat-Cardiac-Energetic`, `Multiscale-Hx-P21-Model`, `cross-bridge-model/XB.mo`, MATMyoSim `update_*` (`J_off ∝ N_on − N_bound`) | **the probe's form. Start here.** Campbell 2018's force-dependent OFF→ON is the same idea as this model's `exp(F/σ)` SRX exits |
| **MyoSim / FiberSim** (Campbell) | two-state thin + three-state thick filament | any | 1–2 (MATMyoSim); spatially explicit (FiberSim) | **local:** `../MATMyoSim` (MATLAB), `../FiberSim` (C++/Python, twitch demos) | reference implementation and cross-check; FiberSim is too costly to embed |
| **Rice 2008** | TnC low/high affinity + N/P Markov, SL- and Q10-dependent | rat, 24 °C default | ~5 regulatory | CellML | only published block with **built-in Q10**; heavier than needed |
| **Tran 2010** | Rice 2008 + MgATP/MgADP/Pi/H⁺ thermodynamics | rat | as Rice | CellML (+ Tran 2009 SERCA, CellML) | the only published pair acting on **both** sides of the twitch under metabolite change |
| **Land 2012** | TRPN + tropomyosin, length/velocity feedback | **mouse, 37 °C** | 2–3 | CellML exposure `e/27d` | **right species**; physiological temperature |
| **Land 2017** | TRPN + TmBlocked; XS/XW with distortion | human, 37 °C | 2 (+4 XB) | **local** in `../TWorld`; ToR-ORd-Land; openCARP `Stress_Land17` | human; the XB part is redundant here |
| **Regazzoni RDQ18/RDQ20** | mean-field regulatory units with nearest-neighbour cooperativity | human | 20 (RDQ20-MF) | GitHub, MATLAB/C++/Python, MIT | most rigorous cooperativity; overkill for step 1 |
| **MUSICO** (Mijailovich) | spatially explicit; Ca-dependent parked state | mouse/rat/human, 25 °C | n/a | availability not checked | Prodanovic 2022 fits **mouse** twitches (Chung et al. data, 25 °C) |

**Beard-lineage precedents for the whole pipeline:** Tewari 2016 *JMCC* (rat Ca + tension time
courses, metabolite effects on twitch) and Lopez 2020 *Function* (rat TAC failure,
Pi-impaired cross-bridge cycle). Your own `cross-bridge-model/InVitroTwitch/` has already done
this once, for human: Mashali 2021 twitch data (37 °C, 0.5–3 Hz, TTP/RT50 by SL and
frequency), the tanh transient, a moment-closure cross-bridge model, and a `fminsearch` fit
(`fit_human_XB.m`, `XBInVitro_objective.m`). **That objective and its data table are directly
reusable.**

---

## 5. Complexity estimate

*Estimates, except where they cite §2.*

| phase | work | effort | risk |
|---|---|---|---|
| **0 — plumbing** | gate into the RHS behind `UseCa` (reuse `NP` slot, drop it from `PT`, multiply `ka_eff`); `caTransient` → `Model/`; a time-based Ca driver in `evaluateModel` (knots at beat onsets, as the probe does); `CaSig` slot + `sigLookup` Hermite fix; regression `UseCa = 0` byte-identical; twitch features (F_dev, TTP, RT50, RT90, F_dia) | **2–4 days** | low: the probe already runs; cost 2–3 s/beat |
| **1 — steady Ca sensitivity** | force–pCa at 8 and 2 mM, same preps → `k_on/k_off`, `K_coop`; add as a feature block | **1–2 weeks + one data session** | medium: n_H vs the model's two-component foot (§2.2) |
| **2 — twitch kinetics** | choose the driver (§3); identify `k_off` and the diastolic SRX fraction; **temperature** (rig T unknown, manuscript §2.1.1 TODO) and **species** (mouse core vs rat/human transients) | **1–3 months** | **high**: none of these are in the current data |
| **3 — ATP-dependent Ca handling** | Tran 2009 SERCA, or at least an ATP-dependent transient decay, for "ATP depletion → twitch" | **2–6 weeks** | medium; literature models exist |
| **4 — boundary modes** | SL clamp as an input; force-controlled boundary for afterloaded twitches | **1 week each** | low |

- **A Ca-driven twitch that runs: ~1 week** (phase 0).
- **A twitch whose amplitude means something: ~1 month**, gated on force–pCa data (phase 1).
- **A twitch whose time course means something: 2–4 months**, gated on the driving-signal and
  temperature decisions (phase 2). This is the same critical path that
  [TriSegCoupling](../TriSegCoupling/conclusions.md) §7 found, now with numbers.

---

## 6. Recommendations

1. **Promote `dPUdT_CaGate`'s factorized gate into the RHS.** Do not resurrect the legacy
   `NP` branch (§1.2).
2. **Request force–pCa at 8 and 2 mM ATP.** The gate makes a testable prediction here. Bound
   heads hold sites on (`Joff ∝ Pu`), so the larger bound fraction at low ATP should **shift
   pCa₅₀ leftward**, independent of any change in `k_on`/`k_off`. If the data shows no shift, the
   bound-head feedback is wrong for this core. *Prediction not yet run: it needs the 2 mM
   parameter set, which lives only in the `ATP-wt-lowatp` worktree.*
3. **Decide the question before the driving signal.**
   - *Mechanism, same preps:* a skinned surrogate (caged Ca, or myofibril Ca jumps) at 8 vs
     2 mM, with no transient to argue about.
   - *Physiology:* a literature mouse transient at matched temperature against intact mouse
     twitches, e.g. Janssen-lab data.
   - *Heart failure and energetics:* an ATP-dependent Ca removal (Tran 2009 SERCA) is
     required. A prescribed transient answers only half the question.
4. **Pin the rig temperature.** Every twitch comparison in §3 depends on it, and the rig can
   run 10 °C, so a within-prep Q10 is within reach.
5. **Watch two interactions:**
   - The isometric attachment already runs at `A_reg = A0 = 0.211` of `ka`. Under shortening,
     `RegAvailShorteningOnly` then raises it. That gives two multiplicative gates on one
     transition, one of them fitted only at saturating Ca.
   - The diastolic SRX fraction (0.139) has never been constrained.

## Files

| file | role |
|---|---|
| [`dPUdT_CaGate.m`](dPUdT_CaGate.m) | probe RHS wrapper: one thin-filament state, gate on `ka`; identity at gate = 1 |
| [`caTransient.m`](caTransient.m) | Rice 2008, Marzban 2020, Tewari/Piacentino, FiberSim two-compartment, constant pCa |
| [`RunTwitchProbe.m`](RunTwitchProbe.m) | E0 signals, E1 clamp, E2 force–pCa, E3 twitches, E5 `k_on` scan, E6 numerics, E4 figure |
| `results/probe_E*.mat`, [`results/twitch_probe.png`](results/twitch_probe.png) | outputs |

## Sources

- Beard et al. 2022, *Biophys J* 121:3213 — [PMC9463691](https://pmc.ncbi.nlm.nih.gov/articles/PMC9463691/) (mouse C57BL/6J RV trabeculae, pCa 4.51)
- Rice, Wang, Bers, de Tombe 2008, *Biophys J* — [CellML](https://models.cellml.org/exposure/a31519a27f4c2fe6158e04fd40eeda98/view) (transient τ 20/110 ms, 0.09→1.45 µM, 24 °C)
- Tran, Smith, Loiselle, Crampin 2010, *Biophys J* 98:267 — [PubMed 20338848](https://pubmed.ncbi.nlm.nih.gov/20338848/), [CellML](https://models.cellml.org/e/4b6); SERCA 2009 — [Biophys J](https://www.sciencedirect.com/science/article/pii/S0006349509002124), [CellML](https://models.cellml.org/exposure/cfa7684fb084ad748bf3061569d99334/view)
- Land et al. 2012, *J Physiol* 590:4553 — [doi](https://physoc.onlinelibrary.wiley.com/doi/10.1113/jphysiol.2012.231928), [CellML](https://models.cellml.org/e/27d/view)
- Campbell, Janssen, Campbell 2018, *Biophys J* 115:543 — [PubMed 30054031](https://pubmed.ncbi.nlm.nih.gov/30054031/)
- Mijailovich et al. 2021, *J Gen Physiol*, doi:10.1085/jgp.202012604 (rat twitch; nonlinear SE + parked state required)
- Prodanovic et al. 2022, *IJMS* 23:1135 — [PMC8835009](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8835009/fullTextXML) (mouse/rat/human twitches at 25 °C, MUSICO)
- Mullins & Bondarenko 2013, *PLoS One* — [doi](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0063141) (mouse Bondarenko + Rice, 25 °C, Fortran)
- Stull, Leppo, Marbán, Janssen 2002, *JMCC*, doi:10.1006/jmcc.2002.2065 (mouse RV trabeculae, 37 °C, fura-2, 2–12 Hz)
- Tewari et al. 2016, *JMCC*, doi:10.1016/j.yjmcc.2016.04.003; Lopez et al. 2020, *Function*, doi:10.1093/function/zqaa018
- Regazzoni cardiac-activation — [GitHub](https://github.com/FrancescoRegazzoni/cardiac-activation) (RDQ18/RDQ20, MATLAB/C++/Python, MIT)
- openCARP LIMPET (`Stress_Land17`, `Gattoni16Land17_DynSarc`, `Tomek20Land17_DynSarc`) — [changelog](https://opencarp.org/download/changelog)
- Myofibril relaxation kinetics (slow linear + fast exponential; kREL,slow ≈ detachment) — [PMC13073943](https://pmc.ncbi.nlm.nih.gov/articles/PMC13073943/); myofibril cardiac-cycle instrumentation — [PubMed 27047763](https://pubmed.ncbi.nlm.nih.gov/27047763/)
