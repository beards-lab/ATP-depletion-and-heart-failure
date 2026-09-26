# Coupling the strain-discretized cross-bridge model into a TriSeg wall

**Status:** feasibility assessment, 2026-09-23. No coupling code written. All numbers below
are measured on this machine (R2023a, `params/fastgate_best.m`) unless marked *estimate*.

---

## Verdict

**Feasible, and computational cost is not the blocker.** The measured cost is bad
(**119 s of wall clock for one 0.6 s cardiac-like length excursion, one sarcomere**), but
~95 % of it is recoverable with two levers that already exist in skeleton form in this
repo. The blocker is **scientific**: the model has no calcium activation and is calibrated
entirely on maximally-activated skinned trabeculae.

Three findings drive everything else:

| | |
|---|---|
| **The mechanical interface already matches TriSeg.** | Not "could be adapted" — the series-elastic + viscous contractile-element structure is *already the same equation*. See §1. |
| **The premise "big position change -> steep dxdt -> small steps" is wrong for this model.** | The model integrates in a **Lagrangian (co-moving) frame**. There is no advection term and no CFL limit. The 11.9x cost penalty for a moving sarcomere is real but comes from three *nameable and fixable* sources, none of them advective stiffness. See §2. |
| **The expensive path may not be the right path.** | `../../../Multiscale-Hx-P21-Model` already runs TriSeg + a 9-moment cross-bridge model, fast. Porting *this project's ATP-dependence* into that is weeks, not months. Whether the 60-bin distribution buys anything at the organ level is untested — and cheaply testable. See §6. |

---

## 1. The interface already matches (the good surprise)

TriSeg (`Multiscale-Hx-P21-Model/model/model.m`) drives its sarcomere as:

```matlab
sigma_LV = -Kse*(Lsc_LV - Ls_LV);                       % series elastic element
dLsc_LV  = (sigma_LV - sigma_pas_LV - sigma_XB_LV)/eta; % viscous contractile element
```

`Ls` comes from wall geometry (algebraic, from `xm`/`ym`); `Lsc` is an ODE state.

This model (`Model/dPUdT_CombinedTransitions.m:272-307`) does:

```matlab
Force   = params.kSE*LSE;              % series elastic element      <-> Kse*(Ls - Lsc)
velHS   = (Force - F_total)/params.mu; % viscous contractile element <-> eta
dLSEdt  = vel - velHS;
```

The correspondence is exact: `LSE` <-> `(Ls - Lsc)`, `params.kSE` <-> `Kse`, `params.mu` <-> `eta`.
Current values `kSE = 4775`, `mu = 0.0192`; force is already in **kPa**, the same unit TriSeg's
`sigma` uses.

Better still, the cross-bridge block is **already driven by the contractile-element velocity,
not the imposed one**. Line 192 sets the strain grid from `(SL - LSE)`, whose derivative is
`velHS`. So the internal coupling is already the physically correct one.

**What actually has to change is one thing:** `vel = params.Vums` (line 36) is a *parameter*,
constant per velocity segment, and `SL` is an ODE state with `dSL = vel`. Under TriSeg, `SL`
becomes an *algebraic output of the wall geometry* and must be passed in. That is a small,
well-defined edit — not a rewrite.

---

## 2. Why a moving sarcomere costs more — and it is not what you think

**Measured** (`params/fastgate_best.m`, `ss = 60`, `dS = 0.002 um`, 129 states, 2-state model):

| run | wall clock | per simulated second |
|---|---|---|
| isometric, 1.0 s | 16.7 s | 16.7 s/s |
| SL 2.2 -> 1.75 over 0.2 s, back over 0.4 s | **119.3 s** | 199 s/s (**11.9x**) |

Single RHS evaluation: **2343 us**. So the beat-like run made **~51 000 RHS calls in 0.6 s of
simulated time.**

### There is no advection term

The strain grid is a *lab-frame array whose coordinate labels move*:

```matlab
s = params.s - (-(SL - LSE) + params.LXBpivot)/2;    % line 192
```

An attached bridge stays in the bin it attached into; its *strain* changes because the grid
label slides under it at `velHS/2`. Every attached bridge sees the same `dstrain/dt`, so this
is exact. There is no `-v*dp/ds` anywhere in the RHS (confirmed: no `diff`, `gradient`,
`circshift`, or neighbour-differencing in the 989-line file).

**Consequence: rapid length change does not, by itself, force small time steps.** No CFL
condition. This is a genuinely good design decision and it is what makes TriSeg coupling
tractable at all.

### Where the 11.9x actually comes from

**(a) Window-slide restarts — 9 of them in 0.6 s.** The grid spans only
`60 * 0.002 = 0.12 um` of half-sarcomere strain. A 0.45 um SL excursion is 0.225 um of
half-sarcomere strain — nearly **2x the window width**. So `movingWindow` fires, terminates
`ode15s`, shifts the arrays by `nds = floor(ss/3) = 20` bins, and restarts
(`Model/evaluateModel.m:210-250`). Measured, in order: 5 slides right during ejection, 4 left
during filling. Each restart drops BDF to order 1, takes a tiny first step, and refactors the
Jacobian.

There is a correctness cost too: the slide **zeroes** the vacated bins
(`PU0(ss-nds+1:ss) = 0`). Any attached mass carried off the edge is *deleted*. Under TriSeg
that is a force discontinuity the ventricular pressure will see.

**(b) The attachment point sweeps across bins, with a kink at every crossing.**
`UseA1AttachmentKernel = 0` in the current best set, so attachment uses 2-bin linear
interpolation (`attachmentPoint`, lines 818-820). That is C0 but not C1 — a derivative kink
every time the attachment point crosses a bin boundary. Over 0.225 um at `dS = 0.002` that is
**~110 kink events per beat**, each one a step rejection for `ode15s`'s error estimator.

**(c) The Jacobian is re-formed by full finite differences — 129 RHS calls each.**

---

## 3. The Jacobian is the single biggest lever (measured)

Numerical Jacobian at a populated, shortening state:

```
density: 5.1 % nonzero (852 / 16641)
DENSE ROWS (>20 nz):  121 P_SR (115)   124 LSE (121)   125 PD (35)   129 D0 (85)
DENSE COLS (>20 nz):  123 SL (125)     124 LSE (125)
```

The structure is: **a sparse, near-banded 120x120 cross-bridge block, plus exactly four dense
rows and two dense columns.**

- The **two dense columns** (`SL`, `LSE`) are the grid-shift coupling — harmless for column
  grouping, they just each need their own group.
- The **four dense rows** (`P_SR`, `LSE`, `PD`, `D0`) are the ones that kill grouping: their
  derivatives are *global sums* over all bins (force, total attached mass). They are also
  **analytically trivial** — `Force = kstiff*dS*sum(s.*p)` differentiates by hand in a few lines.

Supplying just those four rows analytically and grouping the rest:

```
RHS calls per Jacobian:  full FD = 129  ->  grouped = 8   (16x)
```

**16x for perhaps a day of work**, on what is almost certainly the dominant cost term.

---

## 4. The MEX core exists but has bitrotted

`Model/dPUdT_core_mex.c` + `.mexw64` are present and built. It does not currently run:

```
MEX FAILED: packMexVals: expected 187 elements but got 192. Check pchip segment counts.
```

The flat parameter array is config-specific (`packMexVals.m` says so in its header) and the
active switch set has drifted since it was generated. This is regeneration work, not new work.
For a vectorized-MATLAB RHS over 60-element arrays — where per-operation interpreter overhead
dominates — expect **20-50x**, *to be re-measured, not assumed*.

---

## 5. Cost estimate for the coupled model

Composing the measured levers (*estimate*, single sarcomere, one beat):

| configuration | per beat |
|---|---|
| as-is, today | **119 s** |
| + analytic-row Jacobian (16x measured) | ~7 s |
| + working MEX (20-50x, unverified) | **~0.2-1 s** |

Scaling to TriSeg — 3 walls, ~3x129 + ~40 TriSeg/circulation states ~= **430 states**:

- RHS cost ~3x (three independent cross-bridge blocks).
- Jacobian: the three blocks are structurally independent, so groups are *shared* across them.
  ~8 (cross-bridge) + ~10-15 (TriSeg + circulation) ~= **20-25 groups vs 430 columns — a 17-21x
  saving, larger than in the single-sarcomere case.**

| | per beat | to steady state (15-30 beats) |
|---|---|---|
| **unoptimized** | 10-20 min *(est.)* | **5-10 h** — single demo run only |
| **both levers** | **5-25 s** *(est.)* | **2-12 min** |

At 2-12 min per simulation, forward runs are comfortable and a fitting campaign is a
pool-week, not impossible. **Computational cost is a solved engineering problem here.**

### On widening the grid instead of sliding

Covering a full 0.5 um SL excursion plus strain tails needs ~200-250 bins (`ss = 60 -> 250`),
i.e. ~510 states per wall, ~1550 total. RHS cost goes up ~4x, but the **number of Jacobian
groups is set by bandwidth, not by `ss`** — so it stays ~8-10. In exchange you delete all 9
restarts per beat and the mass-discarding discontinuity. Likely **net neutral or better, and
far more robust**. This is the recommended direction.

---

## 6. The strategic question: do you need the distribution at all?

`../../../Multiscale-Hx-P21-Model` already couples TriSeg to a cross-bridge model — using
**moment closure** (`P1_0, P1_1, P1_2, P2_0, ...`, 9 moments + `N` + `U_NR` per wall, ~46 states
total). It already carries MgATP / MgADP / Pi dependence:

```matlab
kd_LV  = kd*(Pi_LV/K_Pi)/(1.0 + Pi_LV/K_Pi);
km2_LV = km2*(MgADP_LV/K_D)/(1.0 + MgADP_LV/K_D + MgATP_LV/K_T);
k3_LV  = k3*(MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);
```

Moment closure is cheap but requires the strain-dependence of the rates to be **low-order
polynomial in strain**. This project's distinctive mechanisms — `XB_Ripped`, catch bonds,
`UseA1AttachmentKernel`, `UseRegistrationAvailability`, occupancy saturation — are precisely
the ones that *cannot* be moment-closed.

**So the real question is whether those mechanisms change anything a ventricle can feel.**
That is cheap to answer and should be answered before committing months:

> **Decisive test (days, not weeks).** Extract SL(t) for one wall from a TriSeg beat. Run the
> full 60-bin model open-loop on that waveform (already demonstrated — the probe in §2 is
> exactly this). Run a moment-truncated version of the same parameter set. Compare force.
> If they agree within the objective's noise floor, use moments and the job is weeks.
> If they diverge, you have *earned* the expensive path and can say why.

---

## 7. What is actually missing (the real work)

Ordered by how much they threaten the result, not by effort.

**1. No calcium activation.** `UseCa = 0`; the primary RHS hardcodes `NP = 0` (line 65) and the
Ca branch is commented out. The Ca path survives only in `Model/dPUdTCa.m`, an older variant.
Every calibration in this project is at **fixed, maximal Ca in skinned trabeculae**. A ventricle
needs a twitch. Two consequences, the second worse than the first:

- Thin-filament activation must be added (the `N` state with `k_on`/`K_coop` is already
  parameterized in `Model/parameterBounds.m:277` and implemented in the P21 model — portable).
- **Re-calibration is unavoidable.** The model currently attributes the *entire* force time
  course to cross-bridge kinetics because activation was constant. Once activation has its own
  time course, `ktr`, the re-stretch rate `rsK80`, and the SRX time constants will all shift —
  and the skinned data cannot identify the activation parameters, because it has no Ca
  transient. You need literature values, new force-pCa data, or an explicit decision to fix
  them from Campbell/Tewari and say so.

**2. Temperature and species.** Skinned trabeculae near room temperature vs. an in-vivo
ventricle. Rate constants differ several-fold. Currently unaddressed.

**3. Passive force double-counting.** This model has titin-based passive identified on the same
preparation (`UsePassive`, `UseTitinIdentifiedPassive`, `UseMaxwellDashpot = 1`). TriSeg has
its own `sigma_pas = k_passive*(Ls - Lsc0)^gamma`. Exactly one of them must survive.

**4. Consistent DAE initialization.** TriSeg's `triseg.m` root-finds `xm`/`ym` for consistent
initial conditions. That has to be extended to include a steady-state cross-bridge
distribution, or the first step will be a transient the ventricle does not deserve.

**5. Mass-conserving window slide.** If the sliding window is kept rather than widened, the
discard at the edge must become a redistribution.

---

## 8. Timescale

*Estimates. Phases 0, 1 and 3 are engineering on a well-understood codebase; phase 2 is
science and is the critical path.*

| phase | work | effort |
|---|---|---|
| **0 — numerics** | regenerate MEX + verify parity; analytic Jacobian rows + grouping; smooth attachment kernel (`UseA1AttachmentKernel = 1`, re-fit the 1-2 params it moves); widen grid or make the slide conservative; benchmark | **1-2 weeks** |
| **1 — open-loop + the §6 test** | drive the model with a TriSeg-derived SL(t); moment-truncation comparison; **go/no-go on the whole distribution-in-TriSeg idea** | **~1 week** |
| **2 — calcium activation** | port thin-filament activation; decide the identifiability question; re-calibrate against the existing four protocols; temperature/species | **1-3 months** (critical path) |
| **3 — closed-loop DAE** | `SL` from wall geometry; reconcile `kSE`<->`Kse`, `mu`<->`eta`, passive; consistent init; steady-state beat loop | **2-4 weeks** |
| **4 — calibration & validation** | organ-level fitting, PV loops, ATP-depletion predictions | **months** |

- **A coupled model that runs and produces a PV loop: ~2 months.**
- **A coupled model whose output means something: ~6-12 months**, and that estimate is
  dominated entirely by phase 2, not by any of the numerics.

---

## Recommended order

1. **Phase 1 first, before phase 0.** The §6 moment-truncation test is days of work and can
   make phases 0, 2 and 3 unnecessary. Do not optimize a path you may not take.
2. If the distribution is justified: phase 0 (cheap, and it makes everything else pleasant).
3. Then treat calcium activation as **its own project with its own data question** — not as a
   step inside the coupling work. It is the only part of this that is genuinely hard.

## Reproducing the numbers

Probes used are in the session scratchpad, not committed (they are ~20 lines each):
timing (`ss`, RHS us, isometric vs. beat-like), MEX parity attempt, and the Jacobian
structure/grouping count. All ran against
`getParams(loadParams('fastgate_best'), [], true, false)`.
