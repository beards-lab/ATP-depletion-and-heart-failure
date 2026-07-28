# ktr vs post-slack vs post-restretch — do the three force rises share one rate? (2026-07-28)

Data only, 8 mM ATP, three protocol days (03/27, 04/03, 04/10).
Reproduce: `CompareRecoveryRates.m` → `results/recovery_rates_traces.png`, `results/recovery_rates_summary.png`.

## Answer: yes, the rates coincide. No, you may not compare the normalised traces.

`k_63` = model-free, 1/(t63 − t_on) against the **measured** plateau. Cycles 1–4 only
(cycle 5 restretches to 1.00 ML instead of 1.10 ML and is excluded).

| k (s⁻¹) | 03/27 | 04/03 | 04/10 | mean | vs post-slack |
|---|---|---|---|---|---|
| post-slack     | 46.4 | 40.4 | 37.3 | **41.4** | — |
| post-restretch | 49.2 | 43.7 | 48.3 | **47.1** | 1.14× |
| ktr protocol   | 49.4 | 42.1 | 38.1 | **43.2** | 1.04× |

The between-window spread (41–47) is smaller than the between-day spread (37–49). **The ktr
protocol's redevelopment is the same process as the post-slack redevelopment** — ktr/post-slack
= 1.04, and the per-day ktr value tracks that day's post-slack value (03/27 highest for both,
04/10 lowest for both). That is the expected result: the ktr manoeuvre ends on a *release*
(1.05 → 1.00 ML), so it starts from zero force at the reference length, exactly like the slack hold.

## Correction to the earlier claim

I previously reported post-restretch as **1.18× slower** than post-slack (τ 28.0 vs 23.6 ms at 8 mM).
That came from a free-offset 3-parameter exponential fitted from the *end of the restretch ramp*.
It is wrong, for the reason below. With the baseline anchored at the true minimum, post-restretch
is if anything slightly **faster**. The low-ATP `LowATP_RestretchRecovery` conclusion used the same
estimator and needs re-running before its 2 mM numbers are trusted.

## Which minimum — this is what broke the earlier fit

The restretch produces **two** valleys, not one:

```
peak1  (~1.20 F_iso, during the ramp)
  → vall_y   first dip, still during the ramp        ← this is extractSlackAttributes' `vall_y`
  → peak2    at the ramp end, ~1.00 F_iso
  → vall2    ~0.81 F_iso, ~7 ms AFTER the ramp ends  ← the recovery actually starts here
  → monotone rise to steady
```

Starting the fit at the ramp end means the first ~7 ms are still *falling*. A free-offset fit
absorbs that descent into the baseline and reports a slower rate. `vall_y` (0.96 F_iso) is the
wrong minimum in the other direction — it is not the start of the recovery and it understates the
amplitude by 4×. Use **vall2 = min over the 80 ms after the ramp ends**.

## Why the normalised overlay is not a like-for-like test

| | post-slack | ktr | post-restretch |
|---|---|---|---|
| starts from | F = 0, no bridges | F = 0, no bridges | 0.80 F_iso, full bridge population just stretched |
| ends at | 0.74 F_iso (short length) | 1.00 F_iso | 1.00 F_iso |
| **amplitude A/F_iso** | **0.74** | **1.00** | **0.22** |
| length during the rise | 0.94–1.02 ML | 1.00 ML | 1.10 ML |
| single-exp r² | 0.97–0.99 | 0.97–0.99 | **0.46–0.86** |

1. **Amplitude.** Post-restretch is a 0.22 F_iso excursion; normalising rescales it 4.6× to sit on
   top of a full redevelopment. Noise is ±0.03 F_iso, so its SNR is ~7 against ~30 for the others —
   that is why the red traces in the normalised panel are so much noisier.
2. **Endpoints differ.** Post-slack plateaus at a *short-length* force (0.79 → 0.63 F_iso from
   cycle 1 to 5); ktr and post-restretch return to full F_iso. Normalising erases that.
3. **Length differs.** Overlap and lattice spacing are not the same in the three windows, so even
   the rate constant is not a pure comparison of the same f_app + g_app.
4. **Post-restretch is not a single exponential** (r² 0.46–0.86 vs 0.97–0.99). Its "rate" is a
   summary of ringing + a slow tail, not one rate-limiting step. Do not over-read it.

Where normalisation *is* defensible: comparing the **rate constants**. For a relaxation, the
eigenvalue is amplitude-independent, so a small perturbation and a full redevelopment should share
it — provided each window really is single-exponential. That holds for post-slack and ktr; it does
not hold for post-restretch.

## Estimator caveat for the ktr protocol

The ktr window is only 58.5 ms and the trace has reached 0.93 F_iso by its end, i.e. it is
**truncated before the plateau**. Consequences:

- Do **not** let A float in `A(1−exp(−k·t))` on a ktr recording. Refitting the post-slack windows
  over a matched 58.5 ms span with A free returns only **86%** of the true rate (post-restretch 92%);
  on 04/03 cycle 5 it collapses to 19 s⁻¹ against a true 30 s⁻¹.
- Instead **fix A at F_iso** (the muscle returns to the same length, so the plateau is known), or use
  the crossing estimator `k_63` with P = F_iso. Both are span-insensitive here because t63 ≈ 23 ms
  sits well inside the window.
- `extractSlackAttributes`' `ktr` fit anchors F0 = 0 but leaves A free. That is safe on the ~275 ms
  slack windows; it would not be safe on a 58.5 ms ktr window.

---

# Part 2 — why the model gets this wrong, and what to do (2026-07-28)

`ModelVsDataRecovery.m`, base `params/ThursdayNightFever.m`. The test is the **SL-matched pair**:
slack cycle 2 holds at exactly SL 2.0 µm, which is the length the ktr protocol redevelops at. Same
length, same starting force (zero), same endpoint — so the two rates must be equal, and in the data
they are.

| | post-slack #2 k63 | ktr k63 | ratio | E_FV | E_ktr | E_stairs | E_slack |
|---|---|---|---|---|---|---|---|
| **DATA 03/27** | 50.0 | 49.4 | **0.99** | | | | |
| V1 TNF as-is | 60.7 | **192.0** | **3.16** | 0.016 | 0.984 | 0.010 | 183.37 |
| V2 `+UseMaxwellTensionOnly` | 61.1 | 85.5 | 1.40 | 0.016 | **0.225** | 0.010 | **173.38** |
| V6 V2 `+tau_reg=2 ms` | **45.4** | **53.0** | **1.17** | 0.020 | **0.054** | 0.010 | 180.06 |
| *V4 A_reg pinned (diagnostic)* | 52.7 | 59.7 | 1.13 | **0.168** | 0.050 | 0.010 | **236.53** |
| *V5 no SRX entry (diagnostic)* | 97.3 | 131.3 | 1.35 | **0.174** | 0.451 | 0.012 | **353.59** |

## Defect 1 — a titin-dashpot force floor that only exists in the ktr protocol

`UseMaxwellTensionOnly` is **never set** in ThursdayNightFever.m, so it takes the `getParams`
default `false`. The Maxwell element then only ever charges, and the ktr protocol's final
1.05 → 1.00 ML release cannot discharge it. At the first sample of the redevelopment window:

```
Force = 12.85   |   active FXB = -2.47   |   passive = 16.30   of which F_Maxwell = 16.26
```

The trace starts at 0.18 F_iso where the data starts at 0.000, and the apparent rate is 192 s⁻¹.
The ODE comment at [dPUdT_CombinedTransitions.m:284](../../Model/dPUdT_CombinedTransitions.m) describes
this exact protocol as the reason the flag exists.

**Fix: set the flag.** ktr 192 → 85.5 s⁻¹, `E_ktr` 0.98 → 0.225, and `E_slack` 183.4 → **173.4**
(also better), FV and stairs unchanged. Strictly an improvement — take it before anything else.

## Defect 2 — `tau_reg` sits inside the ktr timescale

`A_reg` is charged by shortening (`RegAvailShorteningOnly = 1`) and relaxes with `tau_reg = 27.4 ms`.
The ktr manoeuvre shortens for 4.8 + 1.7 = **6.5 ms**; the slack release shortens for **2.2 ms**. With
a 27 ms memory the model still remembers the difference when redevelopment starts:

| at the start of the rise | PD | SR+SRD | p2_0 | **A_reg** |
|---|---|---|---|---|
| post-slack #2 | 0.517 | 0.153 | 0.147 | **0.284** |
| ktr | 0.508 | 0.358 | 0.053 | **0.341** |

20% more availability → 20% more `ka_eff` → a faster ktr, and the memory decays over exactly the
same ~27 ms the redevelopment takes.

**Two alternative explanations were tested and rejected:**
- *SRX charging during the low-force dwell* (my first hypothesis). **V3** stretched the ktr hold at
  0.80 ML from 5 ms to 275 ms so it matched the slack dwell. PD went 0.508 → 0.680 and SR+SRD
  0.358 → 0.187 — the pool state changed a lot — but the rate barely moved (85.5 → 92.3 s⁻¹, in the
  *wrong* direction). **The model's ktr rate does not depend on how long it sat at zero force.**
- *SRX as the source of the asymmetry*. **V5** (`ksr0 = 0`) speeds up *both* windows (61 → 97 and
  85 → 131) and leaves the ratio at 1.35. SRX is a global brake, not a differentiator.

**Fix: make `A_reg` fast, not slow.** At `tau_reg` ≈ 1–2 ms the gate tracks velocity
quasi-instantaneously — the FV shoulder survives, because that is a *steady-state-during-a-ramp*
property, not a memory — but nothing is remembered once the isometric hold starts, so the
redevelopment rate becomes protocol-independent by construction.

At `tau_reg = 2 ms` (with defect 1 fixed):

| post-slack k63 per cycle (SL 2.04 → 1.88 µm) | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| DATA | 48.8 | 50.0 | 43.5 | 43.5 | 43.5 |
| V2 (tau_reg 27 ms) | 62.6 | 61.1 | 58.0 | 54.0 | 55.4 |
| **V6 (tau_reg 2 ms)** | **47.1** | **45.4** | **44.4** | **44.0** | **42.0** |

— the absolute rate *and* the SL dependence now match, and the amplitude A per cycle is untouched
(0.877/0.849/0.811/0.769/0.692, same as V2). `tau_reg = 1 ms` goes further: ratio **1.05**,
`E_ktr` **0.034**.

**The trade-off, stated honestly:** `E_slack` goes 173.4 → 180.1 (2 ms) → 184.5 (1 ms), and its
minimum is at `tau_reg` ≈ 10 ms (169.1). So **`tau_reg` is precisely the parameter that trades slack
trace shape against ktr protocol-independence** — and because defect 1 kept `E_ktr` pinned at 0.98
(dominated by the dashpot floor, insensitive to kinetics), the fit had no counterweight and settled
at 27 ms. Free `tau_reg` in the next round with the ktr trace properly weighted.

**Do not** pin `A_reg` (`v_ref_reg → ∞`). It reaches ratio 1.13 but costs `E_FV` 0.016 → 0.168 and
`E_slack` 173 → 237: the velocity boost is load-bearing for force-velocity.

## Defect 3 — post-restretch is still ~2× too fast, and neither fix touches it

| post-restretch k63 per cycle | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| DATA | 44.4 | 50.0 | 43.5 | 58.8 | 37.7 |
| V2 | 72.6 | 78.6 | 79.6 | 91.1 | 78.2 |
| V6 (tau_reg 2 ms) | 80.8 | 89.6 | 100.0 | 100.2 | 73.5 |

Independent of defects 1 and 2, and slightly *worsened* by the tau_reg fix. This is the
unbounded strain-dependent detachment: `R2 = k2·p2·ppval(PieceWiseStrainDep2, ...)` reaches
2100–4250 s⁻¹ at the strains the restretch creates, so every stretched bridge is gone in <0.5 ms and
the recovery is pure fast re-attachment. Needs a series ATP-binding ceiling on R2 — see
`ATP-wt-lowatp/Analyses/LowATP_RestretchRecovery`.

## Recommended order

1. `params0.UseMaxwellTensionOnly = 1` in ThursdayNightFever.m. Free; improves ktr **and** slack.
2. Add `tau_reg` to `params0.mods` and refit with the ktr trace weighted. Expect it to land in the
   1–5 ms range. Do not touch `A0` (moving it to 0.17 blew `E_slack` to 1111).
3. Separately, the R2 detachment ceiling for defect 3.

---

---

# Part 3 — candidate mechanisms for defect 3

The **SRX-kinetics** test is below (it belongs here — it is a parameter test on the existing model).
The **ATP-ceiling / apo-pool** mechanism grew into its own analysis, including a critical review of
its physiological justification and a refit driver:
→ [`../ApoPoolDetachment/conclusions.md`](../ApoPoolDetachment/conclusions.md).

## A. Slower SRX — mechanism engages, but it cannot pay for itself (2026-07-28)

Baseline throughout: TNF + `UseMaxwellTensionOnly` + `tau_reg = 2 ms` (defects 1 and 2 fixed).
Target: post-restretch 49.2 s⁻¹; baseline gives 92.7 s⁻¹.

## A. Slower SRX — mechanism engages, but it cannot pay for itself

Hypothesis: if the SRX loop were slow enough to be a *memory* rather than a fast equilibrium, it
could distinguish the post-restretch window (which follows a 275 ms hold at reduced force, so SRX
loads) from the post-slack window (2 ms after full activation). Test: scale **all** SRX rates
(`ksr0, ksr2srd, ksrd2sr, kmsrd`) by `xsrx`.

| xsrx | post-slack #2 | post-restretch | ktr | E_FV | E_ktr | E_slack |
|---|---|---|---|---|---|---|
| 1 | **45.4** | 92.7 | 53.1 | 0.020 | 0.054 | **180** |
| 0.3 | 48.1 | 82.5 | 48.9 | 0.037 | 0.091 | 854 |
| 0.1 | 72.4 | 72.1 | 75.6 | 0.049 | 0.116 | 1613 |
| *(data)* | *50.0* | *49.2* | *49.4* | | | |

The SRX pool excursion during the short hold does grow as intended (+0.030 at xsrx = 0.3, +0.066 at
0.1 vs +0.00 at 1), and post-restretch does slow — but only to ~72 s⁻¹, and `E_slack` rises 5–9×.
Steady state is nearly intact (F_iso 76.5 → 70.9, attached 0.270 → 0.247) and the post-slack
*plateaus* actually improve, so the blow-up is in the transients, not the operating point.

**Verdict: falsified as a fix at fixed other parameters.** A slow SRX stops braking the fast window
(post-slack 45 → 72, heading for the no-SRX limit of ~97) at the same rate as it starts braking the
slow one. One SRX cannot be both the fast brake that sets ktr and the slow memory that sets
post-restretch. Caveat: no compensating refit was run, so this rules out the *lever*, not every
possible slow-SRX model.

## C. Directional (one-sided) registration availability — REVIEWED, NOT RECOMMENDED

I proposed extending `A_reg` so that *lengthening* suppresses availability, as a way to slow the
post-restretch recovery without touching detachment. On review the justification does not hold, and
the honest finding is that the **existing** `RegAvailShorteningOnly = 1` is also a fitting device.

**What `A_reg` encodes.** `A_inf = A0 + (1-A0)·v/(v+v_ref)`, gating `ka_eff = ka·A_reg`. The static
part is well founded: actin's 5.5 nm monomer / ~36 nm helical repeat is mismatched against myosin's
14.3 nm crown / 43 nm repeat, so only a subset of heads face a favourably oriented monomer — the
*target zone* (Squire; Daniel, Trimble & Chase 1998, compliant realignment of binding sites). `A0 < 1`
is a fair statement of that.

**The velocity-dependent part is much weaker, and may be double-counting.** In a Huxley-type
strain-distributed model, attachment already occurs continuously at s≈0 while sliding advects the
distribution — the flux of sites past a head is *implicit in the advection*. Multiplying `ka` by a
velocity-rising factor on top is not derived from site availability; it is a phenomenological
correction that happens to fix the FV shoulder.

**And one-sidedness has no mechanism behind it.** The site-sweeping argument that motivates the boost
is direction-blind: sites sweep past just as fast when lengthening. To make availability *fall*
during lengthening requires a *different* physical process than the one that makes it rise during
shortening — one state with one time constant would be standing in for two unrelated things. The
code comment is candid about what it is actually for: *"restretch/lengthening does not boost
attachment -> preserves post-restretch undershoot while keeping the FV shoulder"*. That is a
description of a fit, not of a mechanism.

Arguments considered and rejected:
- *Attachment during lengthening is unproductive* (the head is strained the wrong way immediately) —
  true, but the model **already** penalises it: a head attaching at s≈0 during lengthening is
  advected straight into the high-strain region where R2 is ~700–3000 s⁻¹. Adding a `ka` penalty
  double-counts.
- *Cooperative thin-filament activation* — predicts availability **higher** during lengthening
  (force rises), i.e. the opposite sign.

**The one argument that survives is not about direction at all.** Compliant realignment (Daniel 1998;
Tanner, Daniel & Regnier 2007) says bound bridges shift the register of neighbouring sites through
filament compliance. A perturbation that strips a large bound population leaves the local register
disturbed, and re-equilibration takes time. That predicts a deficit driven by **how much bound mass
was disrupted while the filament stayed loaded** — direction-agnostic, but naturally larger after a
restretch (p2 0.19 → 0.05 while force stays at 0.80 F_iso) than after a slack (p2 collapses, but so
does force, to zero). This is the version worth building, if any.

**Recommendation:** do not add a lengthening penalty to `A_reg`. If an attachment-side mechanism is
to be tested at all, drive it from disrupted bound mass under load, not from the sign of velocity —
and label the current `RegAvailShorteningOnly` as the empirical device it is when interpreting fits.

## B. ATP-limited detachment + tearing into an apo pool

Reaches the measured post-restretch rate (44.8 s⁻¹ against a data 49.2) — the only mechanism tested
that does. But only at settings that fall outside the repo's own α-MHC bounds, and three of its four
physiological justifications do not survive review. Full results, the critical review, and the refit
driver: → [`../ApoPoolDetachment/conclusions.md`](../ApoPoolDetachment/conclusions.md).

The wiring itself is recorded below because it touches `Model/`, not the analysis folder.

## Breaking-change assessment — the wiring is cheap, the tuning is not

**Wiring (done, on this branch):**

| file | change |
|---|---|
| `Model/dPUdT_CombinedTransitions.m` | read `D0`; subtract it from the `PT` complement; `RD0T`; `dD0` balance; append to `f`; the `UseR2Ceiling` block | 
| `Model/getParams.m` | 6 defaults; `expected_PU0` term; `PU0` append | 
| `Model/evaluateModel.m` | `out.D0` diagnostic (2 lines) |

Nothing else needed, because **every** PU consumer indexes head-relative (`Ns*ss+1..8`) — nothing
reads from the tail — and `UseSpaceExtension` only slides mass *within* the p1/p2/p3 blocks, never
resizing the vector.

**Verified** ([Workbench/Tests/TestD0State.m](../../Workbench/Tests/TestD0State.m)):
- PU0 length 128 → 129 when enabled.
- With both flags off, the full battery is **exactly** the pre-change result:
  `E = [0.016 0.984 0.010 183.37]`.
- `UseD0State = true` with `k2rip = 0`: max |ΔForce| = 0.064% of full scale (adaptive-step noise —
  the extra state changes ode15s's error norm; the RHS is identical).
- Probability conserved to 1e-6 with D0 carrying flux.
- Known gap: `Auxiliary/diagnostics/ktrProbe.m` hand-builds `PU0 = zeros(1, 2*ss+7)` and so already
  does not support `A_reg` either. It needs updating before use with any tail-state feature.

**Tuning: a genuine breaking change.** The ceiling cannot be made inert at `MgATP_ref` while still
doing its job, and it shifts the isometric operating point. Getting from "mechanism reaches 44.8 s⁻¹"
to "model still fits everything else" requires a joint refit of at least `k2, k2max, k2rip, alphaRip,
dr3, k_D0T, kstiff2` — that belongs on its own branch, not on top of the current fit.

## Secondary structure worth keeping

- **post-slack k falls with slack depth** on two of three days (04/03: 48→30 s⁻¹ over cycles 1–5;
  04/10: 42→31). 03/27 is flat (49→44). The dead time t0 grows monotonically everywhere
  (03/27 2.3→12.7 ms; 04/03 4.1→17.9) — that is the unloaded-shortening take-up, not the rate.
- **post-restretch k is flat across cycles** (the restretch always returns to the same length).
- **cycle 5** (restretch to 1.00 ML, 0.84 s hold) gives 71–74 s⁻¹ on 04/03 and 04/10 with r² ≈ 0.47 —
  a different, much smaller and noisier event. Keep it out of any mean.
