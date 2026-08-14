# Can A2 detachment recoil into SRX instead of DRX? (2026-08-12)

**Question.** Two variants, both about where a detaching strong bridge *goes*:

- **A —** switch the A2 detachment to **SRXT** instead of DRXT **completely**, so the
  super-relaxed pool becomes an in-series step of every cycle rather than a side
  pool fed only from `PT`.
- **B —** leave the normal path alone and route only the **tear detachment at high
  strain** into SRX, so the SRX jump exists **only during the restretch**.

Which features get better, which worse; can other parameters rescue it?

**Short answer.** Both mechanisms engage cleanly, and **both move the post-restretch
recovery the wrong way**: `rsK` goes from ×2.35 of data to ×2.24–2.82 in all 12
variants, and **no variant beats the baseline objective** (best L2 7.59 vs 6.51).
The reason is a single measured fact — *the SRX exit is force-**accelerated***
(`kmsrd·exp(F/σ_srd1)` = 238 s⁻¹ at the post-restretch force against `kah` = 102 s⁻¹
on the normal route), so **SRX is a fast lane exactly where the restretch happens and
a slow lane exactly where the other two recoveries happen.** Flattening that gate is
the one parameter that flips the sign, and it is load-bearing for isometric force:
`σ_srd1 → ∞` costs F_iso 76.2 → 59.7 kPa and L2 6.5 → 422.

Variant B is nevertheless the first SRX manipulation in this repo that leaves `ktr`,
`steady` and FV **untouched to three digits** — it solves the *selectivity* problem
that killed every earlier SRX test (§4). What it lacks is a destination whose exit is
slow at high force. §5 and §6 say what that has to look like.

**Part 2 (§8–§9) builds the two mechanisms §6 specified and tests them.** The
**load-engaged availability deficit works**: `rsK` ×2.35 → **×1.12** with `ktr` ×1.04
and isometric force unchanged to the digit — the first mechanism in this repo with
full authority over that window. Its bill is the restretch transient shape, and it is
partly payable. The **fixed dwell time is falsified**: on the shared SRX pool it
brakes post-slack instead (`ktr` → 0.63), and on a private rupture-fed pool it is
worth **0.01** in `rsK`, because the recovery never waits for the torn heads at all.

Base: `params/rskR2_w025_opt.m`, protocol 03/27 8 mM (same base as
[`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) §6b, so the numbers
are directly comparable).
Reproduce in order: [`Baseline.m`](Baseline.m) → [`RunSRXRecoil.m`](RunSRXRecoil.m) →
[`RunSRXGate.m`](RunSRXGate.m) → [`RunDestControl.m`](RunDestControl.m) →
[`RunParadoxProbe.m`](RunParadoxProbe.m) → [`PlotParadox.m`](PlotParadox.m).
~20–70 s per evaluation; run one at a time (`MaxRunTime` guards every eval, so two
concurrent MATLAB jobs make healthy points score the timeout penalty).

---

## 1. What was wired

All default-OFF and **exactly inert at the defaults** (ΔL2 = 0, asserted in
`Baseline.m`). `Model/dPUdT_CombinedTransitions.m` + 6 defaults in `Model/getParams.m`:

| param | meaning |
|---|---|
| `SRXFromR2` | fraction of the WHOLE detachment flux (`R2` in 2-state, `R3` in 3-state) landing in SRX instead of `PT`. **Variant A**; `= 1` makes SRX obligatory |
| `SRXFromR2HighStrain` | same, but weighted by a ramp over `s ∈ [sSRXrip, sSRXrip+dSRXrip]` — **variant B**, the rupture branch only |
| `sSRXrip`, `dSRXrip` | that ramp |
| `SRXRecoilToD` | destination `P_SRD` (SRX·ADP) instead of `P_SR` (SRX·ATP). A *mechanically ruptured* bridge never released its ADP, so `P_SRD` is the nucleotide-consistent landing state and `PT` is the one destination it cannot have |
| `SRXRecoilToD0` | destination `D0` instead — the **control** of §5, identical except that its exit is force-*independent* |

Two design points matter for interpreting the result:

- This is a **destination switch, not a rate cap.** `R2(s)` is untouched, so the
  restretch peak and the ktr overstretch peak that pin it at both ends of the strain
  range are unaffected *by construction*. That is the objection that killed
  [`ApoPoolDetachment`](../ApoPoolDetachment/conclusions.md), and it does not apply here.
- The legacy `UseDirectSRXTransition` flag is now `SRXFromR2 = 1` (identical in
  2-state; in 3-state it previously double-counted, adding `R2` to `P_SR` while `dp3`
  also received it — now it routes `R3`, the flux that actually leaves the attached pool).

`rates` gained a 17th element (`RdetSRX`, appended last) surfaced as `out.RdetSRX`, so
the diverted mass is directly measurable.

---

## 2. Results — every variant, both families

L2 is what the optimiser minimises; `slope` is d`rsK`/d(amplitude) on the real
5-cycle protocol. **Lower `rsK ×` is better; baseline is already ×2.35 too fast.**

| variant | L2 | `rsK ×` | `ktr ×` | steady | slope |
|---|---|---|---|---|---|
| **baseline** | **6.51** | **2.35** | 1.07 | 76.2 | −186 |
| A1 all→SRXT f=1.00 | 9.32 | 2.49 | **0.94** | 75.0 | −155 |
| A2 all→SRXT f=0.50 | 8.86 | 2.77 | **1.00** | 75.6 | +89 |
| A3 all→SRXT f=0.25 | 8.58 | 2.79 | 1.03 | 75.9 | +44 |
| A4 all→SRXD f=1.00 | 23.16 | 2.24 | 1.00 | 82.4 | +295 |
| A5 all→SRXD f=0.50 | 9.30 | 2.39 | 1.03 | 79.2 | −49 |
| B1 rip→SRXT s₀=.015 | 7.84 | 2.66 | 1.06 | 76.2 | +39 |
| B2 rip→SRXD s₀=.015 | 8.19 | 2.69 | 1.07 | 76.2 | −83 |
| B3 rip→SRXD s₀=.010 | 9.20 | 2.80 | 1.06 | 76.2 | −36 |
| **B4 rip→SRXD s₀=.020** | **7.59** | 2.59 | 1.06 | 76.2 | −115 |
| B5 rip→SRXD s₀=.005 | 10.94 | 2.82 | 1.06 | 76.2 | +683 |
| B6 rip→SRXD s₀=.010 f=.5 | 8.26 | 2.73 | 1.06 | 76.2 | −95 |
| B7 rip→SRXT s₀=.005 | 8.24 | 2.69 | 1.06 | 76.2 | +90 |
| *DATA* | — | *1.00* | *1.00* | *~80* | *−94* |

### What gets better

| improves | where | by |
|---|---|---|
| **`ktr`** | A-family | cost 0.321 → **0.061** (A2), ratio 1.07 → **1.00**. Making SRX obligatory slows the *zero-force* windows, and `ktr` was 7 % too fast — this is a real, correctly-signed win |
| **`vall_y`** (restretch valley) | B-family with SRXD | 0.147 → **0.013**–0.120, monotone in how much mass is diverted |
| `vall2_dy`, `coolDownLS` | B-family, mild s₀ | 0.081 → 0.056; 0.449 → 0.291 |
| `peak1_y` | A2/A3 | 0.165 → 0.116 |

### What gets worse

| degrades | where | by |
|---|---|---|
| **`rsK`** — the target | **everywhere** | 2.317 → 2.03–4.48 |
| `FV_fnorm` | A-family | +0.23 → **+1.14** (A1). Routing the whole cycle through SRX changes the duty ratio, and FV pays |
| `ovrsht_dy` | B-family, aggressive s₀ | +0.06 → +1.42 (B5) |
| `coolDownLS` | A4 (all→SRXD, f=1) | +0.45 → **+13.4**, with steady 82.4 kPa — catastrophic |
| `XTOR`, `steady` | A1/A4 | ATP turnover and isometric force both move once SRX is in series |

**No variant is a net win.** The best is B4 at L2 7.59 against a baseline 6.51.

![the paradox in three panels](results/paradox.png)

---

## 3. The paradox: fewer cycling heads, yet a *faster* recovery

The intuition — *heads parked in SRX means less cycling, so the recovery must slow* —
fails for three reasons, all now measured on the real sequential protocol
([`RunParadoxProbe.m`](RunParadoxProbe.m)). State at the start of the post-restretch
window, cycle 1:

| | F₀/F_iso | bound (p1+p2) | **PD** | PT | SR | SRD | diverted | k₆₃ |
|---|---|---|---|---|---|---|---|---|
| baseline | 0.839 | 0.160 | **0.642** | 0.103 | 0.067 | 0.027 | 0 | **71.1** |
| rip→SRXD .015 | 0.852 | 0.167 | **0.654** | 0.091 | 0.058 | 0.030 | 0.076 | **79.6** |
| rip→D0 .010 k=15 | 0.857 | 0.164 | **0.550** | 0.086 | 0.050 | 0.020 | 0.164 | **65.5** |

### (a) For the SRXD destination, the detour is a **shortcut** at the force where it happens

*(This subsection is about the `P_SRD` route only. It is true there, and it is **not**
the general explanation — see (a′), which measures the slower `P_SR` route and finds
the real mechanism.)*

Transit time out of `p2`, at the *measured* window-start force F₀ = 66.4 kPa:

| route | at F₀ = 66 kPa | at F = 0 |
|---|---|---|
| normal `p2→PT→PD` (`1/kah`) | **9.84 ms** | 9.84 ms |
| `p2→SRD→PD` (`1/kmsrd·e^{F/σ_srd1}`) | **4.20 ms** ← 2.3× **faster** | 50.55 ms |
| `p2→SR→SRD→PD` | 13.25 ms | 59.61 ms |

`kmsrd·exp(F/σ_srd1)` = **238 s⁻¹** at F₀ versus **19.8 s⁻¹** at zero force. The
routed heads do not queue — they take an express lane straight into the
attachment-ready reservoir, **bypassing the `kah` hydrolysis step entirely**. The
bookkeeping is visible in the table: `PT` falls (0.103 → 0.091) while **`PD` *rises*
(0.642 → 0.654)**. Sending detachment "into the super-relaxed state" *refills the
ready pool faster* than the normal route does.

That is also why the sign flips between windows. At F = 0 the same exit is 50.6 ms,
5× slower than `kah` — which is exactly why the A-variants *slowed* `ktr`
(×1.07 → 0.94) while *speeding* post-restretch. **The mechanism has the right sign in
the wrong window.**

### (a′) CORRECTION — the p2 → SRXT route *is* the slow one, and it still fails

The first version of §3a explained the wrong-way result with the SRXD transit time.
That is right for SRXD and **wrong as the general explanation**, because the route
originally asked for is `p2 → SRXT` (`P_SR`), which is genuinely the *slower* one:

| route to PD | at F₀ = 66 kPa |
|---|---|
| normal `p2 → PT → PD` | 9.84 ms |
| `p2 → SRD → PD` (SRXD) | 4.20 ms — faster |
| **`p2 → SR → SRD → PD` (SRXT)** | **13.30 ms — 1.35× SLOWER** |

So the transit-time argument predicts SRXT should slow the recovery. **It does not:**
k₆₃ 71.1 → **75.4** (still faster). The pool census through restretch cycle 1 says why
(t relative to the ramp start; `detach` = PD+PT+SR+SRD):

| t (ms) | | PD | PT | SR | SRD | detach | bound |
|---|---|---|---|---|---|---|---|
| 20 | baseline | **0.662** | 0.126 | 0.065 | 0.022 | 0.875 | 0.125 |
| 20 | **rip→SRXT** | **0.658** | 0.105 | **0.084** | 0.029 | 0.875 | 0.125 |
| 20 | rip→SRXD | **0.678** | 0.108 | 0.057 | 0.028 | 0.872 | 0.128 |
| 20 | rip→D0 k=15 | **0.592** | 0.095 | 0.049 | 0.016 | **0.751** | 0.130 |

The SRXT route **does** park heads — `SR` rises 0.065 → 0.084, and it takes them out
of `PT` (0.126 → 0.105). The delay is real and measurable. But **`PD` is unchanged to
three decimals** (0.662 → 0.658), and `PD` is the pool attachment actually draws
from. The detour drains the *pass-through* (`PT`, which nothing waits on) and leaves
the *buffer* (`PD`) full. Total detached mass is identical, 0.875 in both.

Contrast `D0`, the only variant that slowed the recovery: it takes 0.124 out of the
**counted detached pool entirely** and `PD` falls 0.662 → 0.592.

> **The real law is not transit time, it is `PD`.** SRXT, SRXD and PT are all *inside*
> the same loop, so routing between them only shuffles mass between upstream pools;
> `PD` is refilled faster than it is drained and never notices. `D0` is *outside* the
> loop, so it actually withholds mass.

And the law is quantitative. Two independent measurements agree:

| | Δ`PD` | Δk₆₃ | slope |
|---|---|---|---|
| this analysis, `D0` k=15 | −0.070 | −5.6 | **+80 s⁻¹ per unit PD** |
| [`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) Part 4, swap **H2** | −0.176 | −12.7 | **+72 s⁻¹ per unit PD** |

**That closes the whole return-path class arithmetically.** To bring k₆₃ from 71 to the
measured 46 s⁻¹ by depleting `PD` alone needs Δ`PD` ≈ −0.33 — `PD` from 0.64 to 0.31,
i.e. **half the entire myosin population held refractory after every stretch**. Which
is exactly the `max D0` = 0.59–0.74 that the old apo-pool fit needed, and exactly why
that fit was rejected on energetics. The number was already in Part 4's H2 row; nobody
had read it as a scaling law.

### (b) The mass diverted is small next to the reservoir it comes from

Even the aggressive force-flat case diverts 0.164 of head mass and depletes `PD` only
0.642 → 0.550 — still **3.4× the bound population it has to rebuild** (0.16). The
attachment flux `RD1 = ka_eff·PD·N_ov·f_lat` barely notices. A 14 % `PD` depletion
buys an 8 % slowdown (71.1 → 65.5 s⁻¹). *Availability is not what limits this window*,
so removing availability cannot slow it.

### (c) What the diversion *does* buy, it buys at the other end of the same flux

Within the single-knob rip→SRXD family (`sSRXrip` the only thing changing),
**r = −0.993** between the two features:

| `sSRXrip` | .005 | .010 | .015 | .020 | none |
|---|---|---|---|---|---|
| `rsK` cost | 4.483 | 4.130 | 3.710 | 3.311 | **2.317** |
| `vall_y` cost | **0.013** | 0.039 | 0.093 | 0.120 | 0.147 |

Divert more mass and the mid-ramp valley `vall_y` deepens toward the data (good), but
the diverted heads return via the express lane, so by the time the ramp ends the
population has already partly rebuilt: the **post-ramp minimum rises** (0.839 → 0.852
→ 0.857 F_iso), the excursion left to recover is smaller, and its apparent first-order
rate is larger. One flux, two features, opposite signs.

*(Across the mixed 20-variant set the same scatter is uncorrelated, r = 0.04 — this is
a property of that one axis, not a universal. The `rip→D0` case proves the point from
the other side: it has the **highest** post-ramp minimum, 0.857, and yet the
**slowest** rate, 65.5 — because its pool holds the mass for 67 ms instead of 4.)*

---

## 4. Why the model makes post-restretch its fastest window — extending the earlier answers

### What was already established

- [`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) §Part 4:
  `kmsr = 0` makes SRX a **one-way force-gated trap**; entry/exit ratio **15.6** at
  F ≈ 0 versus **0.31** at F = 0.86 F_iso, τ_slow 51.9 vs 9.7 ms. **37 % of heads
  transit the trap during post-slack, 0.1 % during post-restretch.** The window is
  fast because "it never leaves the activated state". State swaps: the surviving bound
  pool (H1, −16.9) and the primed reservoir (H2, −12.7) carry it; the **SRX pool size
  at t₀ is irrelevant** (H3, +4.4). Verdict: *"No parameter fix is in sight."*
- [`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) §6b: **half** the
  excess is the Maxwell dashpot (×2.35 → ×1.53 with it off). And the design rule that
  motivated this analysis: *"slowing an existing rate does not clamp the recovery,
  because the model always has a faster parallel route. Only an obligatory in-series
  step can."*

My probe reproduces Part 4's gate numbers independently (entry 329 s⁻¹ at F=0 vs
66 s⁻¹ at F₀; exit 19.8 vs 238 — in/out 16.6 vs 0.28).

### Blind spot 1 — "an obligatory in-series step" was an untested premise

I made the step obligatory (for the rupture branch, §2 variant B; for the whole cycle,
variant A) and **it still does not clamp.** The premise omits *where* in the loop the
step sits. An in-series step slows a cycle only if it is **at the bottleneck**, and the
bottleneck here is not the return path — `PD` is a pre-filled buffer 4× the bound mass
(§3b). The corrected rule:

> Only an obligatory in-series step **upstream of attachment and downstream of a pool
> that is not already full** can clamp this window. Every state on the return path
> (`PT` 0.10, `PD` 0.64) is a reservoir, so the return path has no authority over the
> rate no matter how slow it is made.

This also retro-explains why the earlier `D0` fit "worked": it needed `max D0` =
0.59–0.74, i.e. it had to strip **59–74 % of all heads** before the depletion was deep
enough to matter. That was not a fitting detail — it is the *minimum* the reservoir
argument permits.

### Blind spot 2 — the SRX brake **auto-disengages under load**

Part 4 treated the force gates as levers to *tune* (σ₂, σ_srd1 sweeps) and the pool as
a state to *seed* (H3). Neither framing states the structural consequence: entry
∝ `exp(−F/σ₂)` and exit ∝ `exp(+F/σ_srd1)` **point the same way**, so the in/out ratio
swings **60×** (16.6 → 0.28) across the physiological force range. Therefore

> the model's force-redevelopment rate is a monotonically **increasing function of the
> force it is redeveloping at**, and the coupling is positive feedback: more force →
> gate opens → more available heads → faster rise → more force.

Post-restretch redevelops entirely above 0.84 F_iso; post-slack and ktr start at zero.
**The model is structurally obliged to make post-restretch its fastest window.** That
upgrades Part 4's empirical "no parameter fix is in sight" into a reason, and it
predicts what happened here: *any* pool reached through those gates inherits the
auto-disengagement, which is why both recoil hypotheses backfired and why the σ_srd1
rescue (§5) cannot be free — the same gate sets the isometric duty ratio.

### Blind spot 3 — `vall_y` and `rsK` are one lever, not two features

§3c. Both prior analyses score them independently (weights 10 and 0.25). Any mechanism
that improves the restretch valley by changing the `p2` detachment flux will **worsen**
`rsK`. This is a hard constraint on mechanism design and it removes a whole class of
candidates from consideration before they are coded.

---

## 5. Can other parameters rescue it? — the one lever, and its bill

`σ_srd1 → ∞` removes the force acceleration, making the SRX exit a genuine trap at
high force. It is the only parameter in the SRX system that flips the sign. Its cost
([`RunSRXGate.m`](RunSRXGate.m)); C0* are the **gate-alone controls**:

| variant | L2 | `rsK ×` | `ktr ×` | steady | SRX_ss |
|---|---|---|---|---|---|
| baseline | **6.51** | 2.35 | 1.07 | **76.2** | 0.089 |
| C0a gate only σ=100 | 43.6 | 2.13 | 0.95 | 68.5 | 0.214 |
| C0b gate only σ=1e6 | 422.0 | 2.39 | 1.05 | **59.7** | 0.346 |
| C1 rip.020→SRXD σ=100 | 44.7 | 2.05 | 0.93 | 68.5 | 0.215 |
| C2 rip.020→SRXD σ=1e6 | 439.9 | 2.37 | 1.05 | 59.7 | 0.346 |
| C4 rip.010→SRXD σ=1e6 | 498.9 | **1.96** | 1.05 | 59.7 | 0.346 |
| C5 + `kmsrd`=5 | 21147 | 7.50 | 0.84 | 30.7 | 0.735 |
| C6 all→SRXD σ=1e6 | 1880 | 3.21 | 1.19 | 51.4 | 0.521 |

**Rejected.** Isometric force collapses 76.2 → 59.7 kPa, the SRX pool balloons
0.089 → 0.346, and `coolDownLS` alone goes 0.45 → 375. Decisively, **the routing adds
almost nothing on top of the gate alone** (C0b 422.0 vs C2 439.9): what little `rsK`
moves is the gate, not the mechanism. This is the *third* independent falsification of
tuning SRX kinetics for this window (after Part 3A's `xsrx` scaling and Part 4's σ
sweeps).

### The control that isolates the cause

[`RunDestControl.m`](RunDestControl.m) holds variant B's entry route fixed and changes
**only the destination's exit law** — `D0` is identical to `P_SRD` as a refractory pool
except that `k_D0T` is force-*independent*:

| destination | exit at F₀ | L2 | `rsK ×` | `ktr ×` | steady |
|---|---|---|---|---|---|
| `PT` (baseline) | `kah` 102 /s, flat | **6.51** | 2.35 | 1.07 | 76.2 |
| `P_SRD` s₀=.015 | 238 /s, accelerated | 8.19 | 2.69 | 1.07 | 76.2 |
| `D0` s₀=.015 k=500 | 491 /s, flat | 8.00 | 2.71 | 1.07 | 76.2 |
| `D0` s₀=.015 k=40 | 39 /s, flat | 7.56 | 2.46 | 1.07 | 76.2 |
| `D0` s₀=.010 k=40 | 39 /s, flat | 9.04 | 2.33 | 1.07 | 76.2 |
| **`D0` s₀=.010 k=15** | **15 /s, flat** | 10.24 | **2.26** | 1.06 | **76.1** |

**Confirmed: the exit force gate is the whole story.** As soon as the exit is
force-flat and slower than `kah`, the sign flips — `rsK` finally moves *down*
(2.35 → 2.26), `rsK` per cycle 122/90/93/102/106 → 92/100/109/98/96, and **isometric
force, `ktr` and FV are untouched** because entry is still gated to high strain. The
magnitude is small (§3b caps it), and the bill lands almost entirely in `ovrsht_dy`
(0.41 → 4.55).

**What this buys the SRX framing.** The energetic objection that killed the `D0`
apo-pool — a stretch-induced ATPase burst, against Abbott & Aubert's cheap lengthening
contraction — **does not apply to an SRX·ADP recoil**: a head torn off before ADP
release folds back nucleotide-bound and needs no second hydrolysis. So the physiology
of variant B is sound; only its rate constant is wrong, and the rate constant it needs
is one it does not own.

---

## 6. How to compensate — postulates, ranked

The requirement is now precise: **something whose braking action grows with force (or
at least is force-independent), engaged during a high-force redevelopment, acting
upstream of attachment rather than on the return path.**

1. **A load-engaged availability deficit (recommended).** Drive an attachment gate from
   *bound mass disrupted while the filament stays loaded* — `deficit ∝ |d(bound)/dt| · F`,
   relaxing with τ. Post-restretch: `p2` 0.19 → 0.05 while F holds at 0.84 F_iso → large
   deficit. Post-slack and ktr: `p2` collapses but F → 0 → **no deficit**. This is the
   only sign structure found that separates the windows, it is *upstream* of the
   reservoir so it escapes blind spot 1, and it does not touch `R2(s)` so both pinned
   force peaks survive. It is also the **compliant-realignment** candidate that
   [`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) §"What remains
   structural" and [`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md)
   §C converged on independently (Daniel 1998; Tanner, Daniel & Regnier 2007) — this
   analysis adds the mechanistic *reason* it is the right class, and the arithmetic
   target: `d ln rsK` ≈ −0.85 with the amplitude slope kept negative. Cost: one scalar
   state, 2 params.
2. **Keep variant B as a supporting term, with its own exit.** Give the routed heads a
   private, force-flat, slow return (~15–40 s⁻¹) instead of `P_SRD`'s shared
   force-accelerated one — i.e. split `P_SRD` into a PD-fed branch (current gate,
   load-bearing for F_iso) and a rupture-fed branch. Worth ~10 % of the needed
   `d ln rsK` at L2 +1.1 (B4) to +3.7 (`D0` k=15), with the damage localised in
   `ovrsht_dy`. Not a fix on its own; the reservoir argument caps it.
3. **Pay the mechanical half on the mechanical side.** §6b already has it: `eta_M`×0.3
   gives ×2.35 → ×1.58 for +9.8 feature cost, leaving ~×1.5 for (1) to close.
4. **Do not** revisit the SRX force sensitivities or a global scaling of SRX rates.
   Three independent falsifications now (Part 3A, Part 4, §5 here).

**A discriminating test that costs one run.** If (1) is right, a synthetic restretch
that is immediately released to zero force should recover at the *post-slack* rate; if
the auto-disengaging gate (blind spot 2) is the dominant term, it should stay fast.
That separates the two candidate explanations of the same 2.35× on existing machinery.

---

---

# Part 2 — the two mechanisms §6 called for (2026-08-12)

[`RunLoadRealign.m`](RunLoadRealign.m), [`RunRealignCompensate.m`](RunRealignCompensate.m),
[`PlotRealign.m`](PlotRealign.m). Both wired default-OFF; the baseline regression
(ΔL2 = 0) still passes.

## 8. M1 — the load-engaged availability deficit

### What it is

Bound cross-bridges set the **register** of the neighbouring actin sites through
filament compliance (Daniel 1998; Tanner, Daniel & Regnier 2007) — a myosin head can
only attach where a monomer happens to face it, and the bound population itself
shifts where those monomers sit. Strip a large bound population **while the filament
stays loaded** and the local register is left disturbed; re-equilibrating it takes
time, during which fewer sites are attachable. One scalar state:

```
dA_def/dt = k_cr · max(0, −d(bound)/dt) · F/(F + F_cr)  −  A_def/tau_cr
ka_eff   *= (1 − A_def)
```

Three design choices carry the whole result:

- **Driven by the NET change in bound mass**, so it is identically zero at any steady
  state. The mechanism cannot rescale `ka` or move the operating point; it acts only
  on transients. (Confirmed empirically: `steady` = 76.2 kPa in *every* M1 variant,
  to the digit.)
- **Weighted by force** — this is the selectivity, and it is the entire point:

  | window | bound mass stripped | force during | ⇒ deficit |
  |---|---|---|---|
  | post-restretch | `p2` 0.19 → 0.05 | **holds at 0.84 F_iso** | **large** |
  | post-slack | `p2` collapses | **falls to 0** | ~none |
  | ktr | `p2` collapses | **falls to 0** | ~none |

- **Direction-agnostic** — it reads disrupted mass, not the sign of velocity. That is
  what [`RestretchVsKtrRecovery`](../RestretchVsKtrRecovery/conclusions.md) §C
  demanded when it rejected a lengthening penalty on `A_reg` as "a description of a
  fit, not of a mechanism".

It also never touches `R2(s)`, so both force peaks that pin it survive untouched.

### It works — and it is the first mechanism here with full authority over `rsK`

| variant | L2 | `rsK ×` | `ktr ×` | steady |
|---|---|---|---|---|
| baseline | **6.51** | 2.35 | 1.07 | 76.2 |
| `k_cr`=1, τ=20 ms | 7.67 | 2.43 | 1.06 | 76.2 |
| `k_cr`=3, τ=20 ms | 12.36 | 2.01 | 1.06 | 76.2 |
| `k_cr`=6, τ=20 ms | 47.11 | 1.41 | 1.06 | 76.2 |
| `k_cr`=12, τ=20 ms | 114.3 | **1.23** | 1.07 | 76.2 |
| `k_cr`=6, τ=40 ms | 223.2 | **0.76** ← overshoots | 1.06 | 76.2 |
| `k_cr`=6, τ=20, `F_cr`=60 | 22.58 | 1.66 | 1.06 | 76.2 |
| **`k_cr`=6, τ=20, LOAD-BLIND** | 102.2 | 1.23 | **0.83** | 76.2 |

`rsK` spans **×2.43 → ×0.76**: this is not a weak lever that needs help, it is a lever
strong enough to overshoot. And `ktr` sits at 1.06–1.07 and `steady` at 76.2 across
the whole range — the mechanism is *silent* on the two observables every previous
candidate broke. `FV_fnorm` even improves slightly at high gain (0.518 → 0.438).

**The load-blind control is the proof.** Set `F_cr` → 0 so the deficit no longer cares
whether the filament was loaded: `rsK` improves exactly as much (×1.23) but now
**`ktr` breaks, 1.07 → 0.83** (cost 0.321 → 1.476). The force weighting is doing the
separating, not the gain.

![the compliant-realignment result](results/realign.png)

### The bill: restretch transient SHAPE, and it is partly payable

The damage is entirely in five features, all of them the shape of the restretch
transient — the deficit suppresses re-attachment during and just after the ramp, so
the valley goes too deep and the recovery starts too low:

| `k_cr`=6, τ=20 | Δ L2 |
|---|---|
| `vall_y` | +12.6 |
| `coolDownLS` | +12.4 |
| `vall2_dy` | +10.4 |
| `ovrsht_dy` | +4.3 |
| `peak2` | +3.1 |
| **`rsK`** | **−2.08** |

Those features have their own levers. `kA2re` (PT → p2 re-attachment during
lengthening, zero in the baseline) and `eta_M`×0.3 were tested as compensators:

| variant | L2 | shape L2 | `rsK ×` | `ktr ×` | steady |
|---|---|---|---|---|---|
| baseline | 6.51 | 1.12 | 2.35 | 1.07 | 76.2 |
| M1 k=3 | 12.36 | 7.96 | 2.01 | 1.06 | 76.2 |
| M1 k=3 + `kA2re`=30 | **9.84** | **3.59** | 2.15 | 1.05 | 76.2 |
| M1 k=3 + `eta_M`×0.3 | 22.54 | 18.46 | 1.33 | 1.06 | 76.2 |
| M1 k=3 + both | 16.54 | 11.27 | 1.34 | 1.04 | 76.2 |
| **M1 k=6 + both** | 61.61 | 56.57 | **1.12** | 1.04 | 76.2 |
| `kA2re`=30 alone (control) | 9.68 | 1.57 | 2.58 | 1.05 | 76.2 |

`kA2re` pays down **more than half** the shape bill (7.96 → 3.59) for a small `rsK`
give-back, and it does nothing useful on its own (control: `rsK` ×2.58, *worse*) —
it is specifically a compensator for this mechanism's damage. `eta_M` stacks on the
`rsK` axis as [`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) §6
predicted.

**The best combination lands on the data:**

| | cyc 1 | 2 | 3 | 4 | 5 | ratio |
|---|---|---|---|---|---|---|
| baseline | 122.2 | 90.1 | 92.8 | 102.0 | 106.0 | ×2.35 |
| **M1 k=6 + `kA2re`30 + `eta_M`×0.3** | **48.7** | **48.0** | **48.1** | **49.2** | **50.6** | **×1.12** |
| data 03/27 8 mM | 46.3 | 37.8 | 41.6 | 42.9 | 50.0 | ×1.00 |

`rsK` cost 2.317 → **0.025**, with `ktr` ×1.04 and `steady` 76.2 kPa. That is the
target [`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md) set
(`d ln rsK` = −0.85) reached *without* `k2`, i.e. without its force and `ktr`
couplings — which is what made every earlier route expensive.

### Stated honestly

- **L2 is 61.6 against a baseline 6.51.** The rate is right and the trace through that
  window is not: `coolDownLS` +25.2 and `vall2_dy` +14.1 say the force dips too far
  after the ramp even though it then recovers at the right rate. `rsK` and
  `coolDownLS` measure the same window with different metrics and they now
  **disagree** — the model gets the rate right by getting the shape wrong.
- **Nothing was jointly refitted.** These are 2–3 hand-set knobs on top of a
  76-parameter optimum; the damaged features have many other levers (`kstiff1/2`,
  `kSE`, `kA2hop`, `dr2`, the piecewise strain functions). Unlike
  [`ApoPoolDetachment`](../ApoPoolDetachment/conclusions.md), **there is no
  falsification blocking a refit here**: no parameter bound is violated, no measured
  force peak is broken, and `ktr`/`steady`/FV are untouched. This is the first
  post-restretch mechanism in this repo that earns one.
- The model becomes **too flat across cycles** (spread 2.6 s⁻¹ against a data 12.2).
  The deficit saturates. `Adef_max` and the `F_cr` weighting are the knobs that would
  restore cycle-to-cycle variation.
- `k_cr` = 6–12 is a phenomenological gain with no independent measurement behind it.
  It should be bounded by what compliant realignment can plausibly deliver before the
  mechanism is claimed as physiology rather than as a working hypothesis.

## 9. M2 — the fixed dwell time (constant outflow)

### What it is

A first-order rate gives an **exponential** residence time: a head that arrived 1 µs
ago is exactly as likely to leave as one that has waited 50 ms. A *dwell* is the
opposite. Implemented as the user's suggestion — saturate the outflow in the pool
occupancy so it becomes a **constant flux** once the pool is stocked:

```
flux = k_legacy·P · K/(P + K)   →  k_legacy·P   as P → 0   (legacy exactly)
                               →  k_legacy·K    as P >> K  (CONSTANT FLUX)
```

so `K` is the crossover occupancy and the only new parameter. The pool then drains
**linearly**, a slug of size *m* takes *m/J* to clear — a dwell, not a rate — and a
**bigger slug takes proportionally longer**, which is the sign the amplitude
dependence wants. Two placements: `SRDOutflowK` on the shared `P_SRD` pool (the
literal "heads must stay in SRX for a fixed time") and `D0OutflowK` on `D0`, which is
fed **only** by the rupture recoil.

### Result: falsified, both placements, for two different reasons

| variant | L2 | `rsK ×` | `ktr ×` | steady |
|---|---|---|---|---|
| baseline | **6.51** | 2.35 | 1.07 | 76.2 |
| M2a shared SRXD dwell K=0.30 | 13.81 | 2.47 | **0.74** | 76.1 |
| M2a shared SRXD dwell K=0.10 | 40.83 | 2.35 | **0.63** | 75.9 |
| M2a shared SRXD dwell K=0.03 | 77761 | NaN | NaN | **12.8** |
| M2b recoil + D0 dwell K=0.30 | 7.67 | 2.47 | 1.07 | 76.2 |
| M2b recoil + D0 dwell K=0.10 | 7.99 | 2.48 | 1.07 | 76.2 |
| M2b recoil + D0 dwell K=0.03 | 8.14 | 2.47 | 1.07 | 76.2 |
| M2b recoil + D0 dwell K=0.01 | 7.85 | 2.48 | 1.06 | 76.0 |

**M2a — on the shared pool, a dwell is a post-slack brake.** `ktr` 1.07 → 0.74 → 0.63,
and at K = 0.03 the model collapses (F_iso 12.8 kPa). `rsK` never moves. The reason is
§4's occupancy asymmetry running the other way: **the SRX pool fills far more during
post-slack (+0.373 excursion) than during post-restretch (+0.001)**, so a saturating
outflow bites hardest exactly where the model is already right. Same structural trap
as every other shared-SRX lever.

**M2b — on the private rupture-fed pool, a dwell changes nothing at all.** Compare
like-for-like against §5's first-order control at the same settings
(`D0` s₀=.015, `k_D0T`=40): first-order gives `rsK` ×2.46, the dwell gives ×2.47.
**Converting the exit from a rate to a fixed dwell is worth 0.01 in `rsK`.**

That is a clean, independent confirmation of §4's blind spot 1 and it answers the
question the dwell was asked to answer: *the residence-time distribution of the
detached heads is irrelevant, because the recovery never waits for them.* `PD` = 0.64
is a reservoir 4× the bound population being rebuilt, so how long the torn heads are
held — exponentially, linearly, or for a hard fixed time — cannot matter. Making the
delay "properly fixed" rather than exponential does not rescue a return-path mechanism;
**nothing on the return path can work.**

That is also the sharpest available statement of why M1 succeeds where M2 fails: M1
acts on **attachment** (upstream, where the rate actually is set), M2 on **return**
(downstream of a full buffer).

## 10. Critical re-reading of every mechanism postulated for the recovery-window problem

Six mechanisms have been proposed across three analyses. Re-read with what §3a′ and
§8–9 now establish, three of the six were tested with a lever too blunt to answer the
question they were asked, and **two accepted recommendations were never actually
applied to the current parameter set.**

| # | mechanism | how it was tested | verdict then | verdict now |
|---|---|---|---|---|
| 1 | **Maxwell force floor** unique to ktr (`UseMaxwellTensionOnly` unset ⇒ the dashpot only charges, and the ktr protocol's 1.05→1.00 release cannot discharge it; the trace starts at 0.18 F_iso, apparent rate 192 s⁻¹) | set the flag | *"Strictly an improvement — take it before anything else"* (ktr 192→85.5, `E_slack` 183→173) | **sound, and NEVER APPLIED.** `rskR2_w025_opt.m` still has `UseMaxwellTensionOnly = 0`. §11 |
| 2 | **`tau_reg` sits inside the ktr timescale** (27 ms vs a 6.5 ms manoeuvre, so `A_reg` remembers which protocol it was) | free `tau_reg`, expect 1–5 ms | fix; but it trades against slack shape, minimum at ~10 ms | **sound, and NEVER APPLIED** — still 27.1 ms, because the current objective sets `RunKtr = 0`, so nothing counterweights it. §11 |
| 3 | **Slower SRX** as a memory that distinguishes the windows | scale **all four** SRX rates by one factor `xsrx` | falsified: *"One SRX cannot be both the fast brake that sets `ktr` and the slow memory that sets post-restretch"* | **conclusion right, test too blunt.** A uniform scaling cannot separate entry from exit or one window from another. §2's B-family shows **selective** SRX entry IS achievable (`ktr`/`steady`/FV untouched to 3 digits) — the statement survives for a sharper reason: the *destination*, not the *rate* |
| 4 | **R2 ceiling + tearing into an apo pool** | `k2max`/`k2rip` scan | falsified by the **ktr overstretch peak** (data 0.72–0.83 F_iso; uncapped model 0.746, every capped variant 1.2–2.2) | **falsification stands and is the strongest in the file.** Its surviving half — that `R2(s)` conflates chemical cycle completion with mechanical rupture — is what §1 acted on, as a destination switch with no cap, so both pinned peaks survive |
| 5 | **Directional `A_reg`** (lengthening suppresses availability) | reviewed, not built | rejected — *"a description of a fit, not of a mechanism"*; the site-sweeping argument is direction-blind | **rejection correct.** §8's M1 is the version that review pointed at: driven by disrupted mass under load, not by the sign of velocity |
| 6 | **State swaps** H1/H2/H3 (Part 4) | seed the post-restretch window with donor pools | *"No parameter fix is in sight"*; H3 ⇒ SRX pool size is irrelevant | **the most under-read result in the file.** H2 is not a null — it is a *calibration*: `PD` 0.646→0.470 bought −12.7 s⁻¹. §3a′ reproduces the slope independently and turns it into the arithmetic that closes the entire return-path class |

### What this buys — three usable things

1. **Two free, unapplied fixes** (#1, #2). Both were called improvements and neither is
   in the current snapshot, because the objective stopped running the ktr protocol and
   the terms that would have held them in place went silent. §11 tests them.
2. **The `PD` scaling law** (#6 + §3a′): **≈ 72–80 s⁻¹ per unit `PD`**, measured twice
   by independent methods. Any future mechanism that works by withholding heads can be
   costed on the back of an envelope before it is coded — and every one of them fails
   that envelope, which is why M1 (which gates attachment instead) is the one that works.
3. **"Blunt lever" is the recurring methodological failure** (#3, and #5 before review).
   Three of six were tested by scaling something global and concluding the mechanism
   class was dead. In two of the three the class was fine and the *lever* was wrong.
   The lesson for #6's successor: test selectivity explicitly before declaring a class
   falsified.

### What is still unexploited

- **The ktr overstretch peak as an explicit feature.** ApoPool recommended promoting it
  out of the `E_ktr` trace MSE; it is still buried, and it is the only measurement that
  constrains `R2(s)` at large strain.
- **The ktr protocol is not in the current objective at all** (`RunKtr = 0`), which is
  the direct cause of #1 and #2 drifting back. Any claim about protocol-independence
  made on this snapshot is unconstrained by the protocol it is about.

## 10b. `UseMaxwellTensionOnly` — still a free win, but for a different reason

[`RunMaxwellTensionOnly.m`](RunMaxwellTensionOnly.m). §10 item 1 flagged that this
never got applied. Tested at `params/rskR2_w025_opt.m`, slack protocol, 19 scorable
features:

| | slack L2 | `rsK` (s⁻¹) | `ktr` | steady | `vall_y` |
|---|---|---|---|---|---|
| **flag OFF (current)** | **5.869** | 111.0 ×2.54 | 52.63 | 76.11 | 67.57 |
| **flag ON, `x_M_slack`=0.01** | **5.672** | **101.8 ×2.33** | 52.83 | 76.11 | 67.57 |
| flag ON, `x_M_slack`=0.003 | 6.689 | — | — | — | — |

| feature | OFF | ON | Δ |
|---|---|---|---|
| **rsK** | 3.111 | **2.297** | **−0.814** |
| **coolDownLS** | 0.452 | **0.268** | **−0.184** |
| `ovrsht_dy` | 0.573 | 1.363 | +0.790 |
| ktr | 0.337 | 0.354 | +0.016 |
| everything else | | | < 0.006 |
| **TOTAL** | 5.869 | **5.672** | **−0.197** |
| **TOTAL without `ovrsht_dy`** | 5.296 | **4.309** | **−0.987 (−19 %)** |

**Confirmed as an improvement**, and a bigger one under the objective with
`ovrsht_dy` dropped (§11), since that feature is the *only* thing it costs.
`x_M_slack` = 0.003 is too aggressive; the 0.01 default is right.

### But the original rationale does not apply here — the mechanism is different

The flag was introduced because a stretch-then-release (the ktr protocol's
1.05 → 1.00 ML) leaves an undischarged force floor. **That cannot be what is
happening on the slack protocol.** Standing Maxwell stress, flag OFF:

| event | F_total | F_Maxwell | share |
|---|---|---|---|
| cyc1 restretch ramp END | 77.11 | **8.59** | 11.1 % |
| cyc1 +40 ms | 79.26 | 2.07 | 2.6 % |
| **cyc2 NEXT RELEASE** (where the flag was supposed to act) | 78.30 | **0.000** | **0.0 %** |
| cyc4 restretch ramp END | 85.37 | 14.78 | 17.3 % |
| cyc4 +40 ms | 81.44 | 3.46 | 4.3 % |

The hold between restretch and the next release is 279 ms = **10.4 × `eta_M`**, so by
the time the fibre shortens there is nothing left to dump. The flag is inert at every
release in this protocol.

**Where it actually acts is inside the post-restretch recovery**, through the serial
element. `dSLc = vel − dLSEdt`: the imposed velocity is zero during the hold, but
active force is *rising*, so `LSE` lengthens and the contractile part therefore
**shortens** — `dSLc < 0` with no external shortening at all. The flag reads that as
titin going slack and discharges the stored stress:

| t after ramp (ms) | F_Maxwell OFF | ON | Δ |
|---|---|---|---|
| 0 | 8.592 | 8.592 | 0.000 |
| 10 | 6.126 | 5.754 | −0.372 |
| **20** | **4.438** | **3.074** | **−1.364 (−31 %)** |
| 40 | 2.070 | 1.248 | −0.821 |

This is precisely the term [`RestretchRecoveryFit`](../RestretchRecoveryFit/conclusions.md)
§6b identified as half the `rsK` excess — *"the post-restretch net recovery is the
small difference of two large opposing terms, a +22 kPa active rise against a −8 kPa
passive decay"*. The flag makes that passive decay track the contraction instead of a
fixed 27 ms clock, which is why `rsK` and `coolDownLS` are the two features that move.

**So the recommendation to set it stands, but the reason recorded in
`RestretchVsKtrRecovery` is not the reason it helps here.** That matters, because it
means the benefit does *not* depend on the ktr protocol being in the objective.

**Caveat:** only the slack protocol was scored (the run had to share the box with a
live optimiser). The FV and passive terms were not evaluated. The code comment claims
passive cycles are untouched because they *"restretch monotonically and hold"* — that
is plausible (during a passive hold, force decays, `LSE` shortens, so `dSLc > 0` and
the release term stays off) but it is **unverified at this snapshot**. Check on the
full battery before adopting.

**The currently running refit does not have this flag set** — it inherits
`UseMaxwellTensionOnly = 0` through `realign1_opt`. Turning it on is worth roughly
−0.99 on the reduced objective by itself, so the next run should start from it.

## 11. The reoptim — does the mechanism survive a joint refit?

[`RunOptimRealign.m`](RunOptimRealign.m) → `params/realign1_opt.m`; assessed with
[`ReportSnapshot.m`](ReportSnapshot.m). Seed = §8's `k_cr`=6 + `kA2re`=30 +
`eta_M`×0.3 point, plus `SRXFromR2HighStrain`=0.2 so the §2 routing could be explored.
29-parameter pool (mechanism + series-viscoelastic + restretch + SRX system + core),
`k_cr` **in** the pool so the optimiser could switch the mechanism off if it did not
pay. **Objective unchanged** from the baseline snapshot, so L2 is comparable.

Budget 1.3 h → **3 rounds**. Cost trajectory **63.22 → 11.91 → 11.45 → 10.96**, still
descending when the clock ran out. (The baseline snapshot itself came from an 18 h
run. This is not a converged optimum and must not be read as one.)

### The optimiser kept the mechanism

`UseLoadRealign` stayed on, at **`k_cr` 6 → 1.56** — it dialled the gain down 4× and
paid for the rest with `ka`↑ (216→242), `kstiff2`↓ (24722→23316), `eta_M` relaxed back
up (0.0081→0.0128), `kSE_M`↑, `mu_neg`↓, `sigma_srd1`↑ (26.7→31.8), `ksr0`↑, `kmsrd`↓.
It did **not** walk `k_cr` to zero, which was available to it in every round.

### Where it landed

Slack-protocol features only (19 of 26 fn entries; the FV and passive terms are not
run by this report and are excluded rather than scored as missing):

| feature | baseline | M1 k=6 hand | **reoptim** |
|---|---|---|---|
| **rsK** | 3.111 | **0.031** | **0.639** |
| **ovrsht_dy** | 0.573 | 4.688 | **4.639** |
| coolDownLS | 0.452 | 25.320 | **0.271** |
| ktr | 0.337 | 0.357 | **0.233** |
| peak1_dSL | 0.313 | 0.419 | **0.299** |
| t0_crossing | 0.213 | 0.258 | 0.271 |
| peak1_y | 0.169 | 0.408 | 0.386 |
| vall_y | 0.146 | 9.019 | **0.025** |
| A | 0.106 | 0.807 | 0.470 |
| vall2_dy | 0.078 | 14.116 | 0.194 |
| doublePeak | 0.040 | 0.000 | 0.540 |
| peak2 | 0.037 | 2.983 | 0.041 |
| **TOTAL (slack)** | **5.869** | 58.819 | **8.396** |
| **TOTAL minus `ovrsht_dy`** | **5.296** | 54.131 | **3.757** |
| *full battery (optimiser's own number)* | *6.509* | *63.218* | *10.964* |

| observable | DATA | baseline | M1 k=6 hand | **reoptim** |
|---|---|---|---|---|
| `rsK` (s⁻¹) | 43.74 | 111.04 ×2.54 | 49.41 ×1.13 | **74.55 ×1.70** |
| `ktr` (s⁻¹) | 49.21 | 52.63 | 51.83 | **50.77** ← best |
| steady (kPa) | 77.35 | 76.11 | 76.11 | 75.60 |
| `vall_y` (kPa) | 71.11 | 67.57 | 41.57 | **71.89** ← near-exact |
| `peak2` (kPa) | 78.42 | 77.35 | 51.63 | 80.71 |
| `ovrsht_dy` (kPa) | 1.33 | 1.41 | 0.05 | 0.05 |

![timecourse](results/report_timecourse.png)

The traces are the clearest statement: the hand-set M1 (red) obviously overshoots — it
digs the valley to 42 kPa where the data is at 68. The reoptim (yellow) sits on the
data through the ramp, the double peak and the valley, and follows the recovery where
the baseline (blue) runs above it. `coolDownLS` — the least-squares error over exactly
that window — confirms it: **0.452 → 0.271**.

### Are we better? — the honest answer is "conditionally, and not yet on the objective"

- **On the objective as written: no.** 10.96 against 6.51.
- **The regression is one feature.** `ovrsht_dy` alone is **+4.07 of the +2.53 net**
  (161 % of it): improvements elsewhere total −2.89, degradations other than
  `ovrsht_dy` total +1.35. Strip that single term and the reoptim **beats** the
  baseline, **5.296 → 3.757 (−29 %)**, while taking `rsK` ×2.54 → ×1.70 and `ktr`
  ×1.07 → ×1.03.
- **And that feature is one this repo already voted to retire.**
  [`SlackDataAnalysis`](../SlackDataAnalysis/conclusions.md) §5 lists `ovrsht_dy` under
  **"do not fit"** — ratio CV **115 %** with the sign flipping across preps
  (+1.35 / +0.34 / +0.37 / −0.01) — and its recommendation 4 is *"Retire `ovrsht_*`
  and `ktr2_*` from cost functions."* It still carries weight 1 in `params0.fn`.
  So the single term blocking this result is a term with no reproducible sign, and
  the mechanism is being penalised for failing to reproduce a 1.3 kPa (1.7 % of F_iso)
  feature that three preps do not agree on.

That is a real caveat in both directions and it should not be waved away: **retiring a
feature after seeing that it blocks your mechanism is exactly how a fit gets talked
into itself.** The defensible statement is the conditional one, plus the note that the
recommendation to retire it predates this analysis by weeks and was made on
data-reproducibility grounds alone.

### What this run did NOT test

`SRXFromR2HighStrain` came out at exactly its seed value 0.2 — it was never drawn in
3 rounds of 8-parameter subsets from a 29-parameter pool. **The p2→SRXT-under-stretch
route still has not had a joint test**, only the one-at-a-time ones in §2. That is the
first thing a longer run should cover.

### Recommended next steps

1. **Re-run with `ovrsht_dy` removed from `fn`** per SlackDataAnalysis's standing
   recommendation, and with a real budget (8–18 h, not 1.3). The current trajectory was
   still descending at the cut-off.
2. **Apply the two never-applied fixes first** (§10): `UseMaxwellTensionOnly = 1` and
   free `tau_reg` with the ktr protocol actually in the objective. Both interact with
   the same series-viscoelastic levers this refit is using, so doing them afterwards
   would waste the refit.
3. Bound `k_cr` in `parameterBounds.m` before the mechanism is described as physiology
   rather than as a working hypothesis.

## 12. Files

| script | what it does |
|---|---|
| [`RunLoadRealign.m`](RunLoadRealign.m) | §8/§9 — M1 gain/τ/`F_cr` scans incl. the load-blind control, and both M2 dwell placements. Chunked: set `SEL` before running |
| [`RunRealignCompensate.m`](RunRealignCompensate.m) | §8 — can `kA2re`/`eta_M` pay the shape bill |
| [`PlotRealign.m`](PlotRealign.m) | `results/realign.png` |

## 7. Files (Part 1)

| script | what it does |
|---|---|
| [`Baseline.m`](Baseline.m) | reference point + the inertness regression (ΔL2 = 0) |
| [`srxRecoilCase.m`](srxRecoilCase.m) | one variant → L1/L2 feature costs, `rsK` per cycle, amplitude slope |
| [`RunSRXRecoil.m`](RunSRXRecoil.m) | §2 — the 12 A/B variants + the per-feature better/worse breakdown |
| [`RunSRXGate.m`](RunSRXGate.m) | §5 — the `σ_srd1` rescue attempt, with gate-alone controls |
| [`RunDestControl.m`](RunDestControl.m) | §5 — the force-flat-exit control (`D0` destination), `k_D0T` scan |
| [`RunParadoxProbe.m`](RunParadoxProbe.m) | §3 — diverted mass, `PD` depletion, transit-time arithmetic on the real protocol |
| [`PlotParadox.m`](PlotParadox.m) | `results/paradox.png` |

## Caveats

- The MEX path (`UseCompiledMex`, `dPUdT_core_mex.c`, `packMexVals`) does **not**
  implement the routing. It is off in this paramset; anything enabling it must
  regenerate the MEX or the routing will be silently ignored.
- `assertActiveParams.m` is a manual mirror of the switch structure and was not
  extended, so the new params are not physiology-scored. They have no entries in
  `parameterBounds.m` either — deliberate, since none of them is a rate.
- The amplitude `slope` column is measured on the real 5-cycle protocol, where release
  depth and re-stretch amplitude are decoupled; the +714 figure quoted in
  `RestretchRecoveryFit` §6 is from a *synthetic* sweep where they move together. On
  the real protocol the baseline slope is already −186 against a data −94, so **the
  slope is not the defect on this base — the level is** (`rsK` ×2.35).
