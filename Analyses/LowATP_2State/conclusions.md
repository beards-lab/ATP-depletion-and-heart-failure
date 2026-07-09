# Best 2-state low-ATP parametrization (on the Further-optim baseline)

**Date:** 2026-07-09 · worktree `ATP-wt-lowatp`, branch `lowatp-2state` (off `Further-optim` b399eb9).
**Deliverable:** `params/params_2state_lowATP.m` — the best 2-state low-ATP fit we can get now.

## What it is
- **Baseline:** `params/optfull_opt.m` (your newest global optim, 2-state, high-ATP feature cost **2.82**).
- **Mechanism (2-state, physiological):** low-ATP condition = **ADP-trap** (`UseAdpTrap`: R2 P2-detachment × g2(MgADP,K_D)) + **Pi-force** (`UsePiForce`: kstiff2 × 1/(1+Pi/K_Pi)). No rigor/3rd state.
- **Tuned constants:** `K_D = 0.35`, `K_Pi = 6`; concentrations `MgATP=2, MgADP=1.6, Pi=0.8`.
- **ONE file, both conditions:** as-is → low-ATP; set `MgATP=8, MgADP=0, Pi=0` → high-ATP.
  The mechanism is **identically inert at 8 mM** (g2=1, fPi=1) → **high-ATP unaffected, verified max|Δ|=0**.
- Driver: `Workbench/DriverLowATP_2State.m`.

## Fit: relative-ratio cost 0.47 (default constants 1.29)

| feature | model 2/8 | data 2/8 | verdict |
|---|---|---|---|
| steady | 1.159 | 1.179 | ✅ |
| A | 1.187 | 1.182 | ✅ |
| Am | 1.201 | 1.176 | ✅ |
| peak2 | 1.135 | 1.180 | ✅ |
| vall_y | 1.176 | 1.200 | ✅ |
| ktr (mean) | 0.67 | 0.54 | ~ high, slope flat |
| **peak1_y** | **1.375** | 1.252 | ✗ over |
| **t0 (mean)** | **1.90** | 1.26 | ✗ onset over-delayed |

**Force amplitude is solved** (all ~+18% from physiological [ADP]/[Pi]). The residual is the **restretch
transient**: onset (t0) over-delayed and peak1 over-lifted — the 2-state force↔onset weld (R2 gates both).
`K_D×K_Pi` sweep confirms the weld: any K_D giving correct force forces t0 ≥ 1.9× (data 1.26).

## Plan to tune further (2-state, physiologically justified)

**1. ATP-dependent SRX return `kmsrd` — the answer to "will SRX-ATP help?"** *(most promising for the residual)*
   - Current config: `kmsr=0`, `UseAtpOnUNR=0` (ATP→SRX mobilization inert); the live SRX lever is
     `kmsrd=20.87` (SRX-ADP return), which moves **ktr while holding force** (max-Ca SRX pool is
     mechanosensitively suppressed, so it barely touches force). That is exactly the *kinetic-only*
     lever the ADP-trap lacks.
   - Test: make `kmsrd` (and/or `ksrd` parking) **[ADP]/[ATP]-dependent** (e.g. `kmsrd_eff = kmsrd·h([ADP])`),
     physiologically grounded (SRX/DRX equilibrium is nucleotide-set). Target: pull t0 & the ktr slope
     toward data **without disturbing the solved force**. This is the one lever that can break the weld
     on the kinetic side inside a 2-state model.
   - Caveat: the *naïve* SRX-ATP (more parking at low ATP) lowers force and slows onset further — wrong
     way. Use the return-rate/mechanosensitivity, not the pool size.

**2. Strain-dependent ADP-trap.** Make g2 act only where P2 bears force (operating-strain window), so it
   traps force without gating the near-zero-strain heads that set onset — decouples force from t0 in 2-state.

**3. peak1 damping** (shared with the high-ATP residual): kstiff2 strain-shape / c_SE_visc so accumulated
   P2 contributes steady force but less transient first-peak.

**4. Classic-vs-slack ktr side quest** (SRX timescale separation) — an extra constraint if useful.

## Reproduce
`cd(root); addpath(genpath('.')); DriverLowATP_2State`  (needs `data/protocol_03_27_2026_{2,8}mM_slack.mat`).
Fit/verify scripts logged in the session scratchpad; the ADP-trap+Pi wiring is committed to `Model/`.
