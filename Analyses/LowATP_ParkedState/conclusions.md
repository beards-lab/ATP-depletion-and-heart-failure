# Forced-detachment "parked" state (PK) — implemented, tested, REFUTED (2026-07-14)

## What was built (branch `lowatp-2state`, worktree `ATP-wt-lowatp`)
A new scalar ODE state `PK`, flag `UseForcedDetachmentPark` (OFF by default → baseline byte-identical):
- **Entry** `p2 → PK`: high-strain tearing rip, strain-gated only (both ATP levels):
  `R_park = kpark·p2·exp(alphaRip·(s+dr3))`. Subtracted from `dp2` (`-R_park`), collapsed to the scalar
  pool `dPK = sum(R_park)·dS − R_unpark`.
- **Exit** `PK → PT` (DRX-ATP): ATP-gated `R_unpark = kunpark·MgATP/(MgATP+K_T_park)·PK` (fast at high ATP
  → PK≈0; slow at low ATP → PK accumulates). Returns to `PT` by conservation.
- Wiring: `getParams.m` (flag/rate defaults, `expected_PU0`, `PU0` tail after `A_reg`), `dPUdT_Combined
  Transitions.m` (read at `Ns*ss+8+UseRegistrationAvailability`, `PT`+overlap conservation, flux block,
  `dp2 -R_park`, `f`-append), `evaluateModel.m` (`out.PK`). Runs clean; `numel(PU0)` grows by +1 when on.

## Result — hypothesis REFUTED (Step 2, ratio 2mM/8mM vs data)
Data ktr ratio slopes DOWN across segments: **0.64 0.59 0.53 0.47 0.46**.

| config | cost | ktr ratio (per segment) | peak1R | steadyR |
|---|---|---|---|---|
| no parking | 0.473 | 0.76 0.77 0.77 0.80 0.77 | 1.25 | 1.18 |
| kpark=40  | 0.661 | 0.83 0.82 0.79 0.81 0.86 | 1.22 | 1.14 |
| kpark=120 | 1.105 | 0.86 0.85 0.83 0.82 0.85 | 1.17 | 1.09 |
| kpark=300 | NaN   | (unstable) | NaN | 1.04 |

- **ktr ratio goes UP (0.77→0.86) — the WRONG direction** (data wants ↓ to ~0.54).
- **No per-segment slope** emerges (stays flat).
- Cost worsens (0.47→1.1); `peak1R` and `steadyR` both fall away from data.

## Why (the mechanistic lesson)
Parking **removes force-bearing heads = an AMPLITUDE effect**; it does NOT slow the **cycling RATE** of the
remaining heads. **ktr is a RATE.** So parking lowers force (steady/peak1), not ktr — and in the ratio the
rip disrupts the 8 mM redevelopment *more* than the 2 mM one, so the ktr ratio *rises*. The ktr residual
needs a **rate/timescale lever** (ADP-trap slows R2; SRX `kmsrd` slows the SRX-return timescale), **not a
head-count lever**. This is exactly why `kmsrd` (a return-RATE) breaks the force↔ktr weld and parking
(a head-count) cannot.

## 8 mM absolute is also perturbed (Step 1)
Activating the rip (baseline had `k2rip=0`) removes `p2` during the high-strain restretch at BOTH ATP:
dSteady ~9.5, dKtr ~10, dPeak1 ~15 at kpark=50 (grows with kpark). `kunpark` 500 vs 3000 barely matters —
the **p2-loss during the transient**, not the PK pool size, is what perturbs. So the mechanism would also
require an 8 mM re-fit — not worth it given it fails the ratio direction test.

## Status
Implementation committed with the flag **OFF by default** → baseline/`optfull` untouched. Kept for the
record and any future variant, but refuted as a ktr-slope lever. Reproduce: `RunParkPreTest.m`.
