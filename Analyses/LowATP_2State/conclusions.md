# Best 2-state low-ATP parametrization (on the Further-optim baseline)

**Date:** 2026-07-09 · worktree `ATP-wt-lowatp`, branch `lowatp-2state` (off `Further-optim` b399eb9).
**Deliverable:** `params/params_2state_lowATP.m` — the best 2-state low-ATP fit, PHYSICAL mechanisms only.

## What it is
- **Baseline:** `params/optfull_opt.m` (newest global optim, 2-state, high-ATP feature cost **2.82**).
- **Mechanism: ADP-trap ONLY** (`UseAdpTrap`): elevated [ADP] blocks ADP-release/detachment, so the
  P2→detachment rate is scaled `R2 × g2(MgADP,K_D)`. Physically well-established (Cooke & Pate; Dantzig).
- **Tuned constant:** `K_D = 0.7`; concentrations `MgATP=2, MgADP=1.6, Pi=0.8`.
- **ONE file, both conditions:** as-is → low-ATP; set `MgATP=8, MgADP=0` → high-ATP.
  ADP-trap is `g2=1` (inert) at MgADP=0 → **high-ATP identical to optfull, verified max|Δ|=0**.
- Driver: `Workbench/DriverLowATP_2State.m`.

## Fit: relative-ratio cost 0.47

| feature | model 2/8 | data 2/8 | verdict |
|---|---|---|---|
| steady | 1.177 | 1.179 | ✅ |
| A | 1.201 | 1.182 | ✅ |
| Am | 1.209 | 1.176 | ✅ |
| **peak1_y** | 1.255 | 1.252 | ✅ (now exact) |
| peak2 | 1.096 | 1.180 | ~ slightly under |
| vall_y | 1.115 | 1.200 | ~ slightly under |
| t0 (mean) | 1.50 | 1.26 | ~ closer |
| **ktr (mean)** | **0.774** | 0.536 | ✗ the residual |

Force amplitude is captured; the dominant residual is now **ktr** (kinetics too fast at low ATP).

## Why NOT a Pi force-trim (correcting the earlier version)
The earlier file scaled `kstiff2 × 1/(1+Pi/K_Pi)` — **unphysical**: stiffness is the mechanical spring
constant of an attached head; [Pi] cannot soften an attached bridge. Pi's real action is **mass action on
the power stroke** (`A·M·ADP·Pi → A·M·ADP + Pi` releases Pi), so it belongs on the **rate `R12`**, not the
stiffness. I wired the honest version (`UsePiReversal`: `R12 × 1/(1+Pi/K_Pi)`, shifts p2→p1, force falls via
POPULATION). **Result: it does not improve the fit** — the K_D×K_Pi grid drives K_Pi→∞ (Pi off) and prefers
a weaker ADP-trap alone (cost 0.473 either way). At [Pi]=0.8 mM the stroke-inhibition is small. So the
best PHYSICAL config is ADP-trap alone; `UsePiReversal` stays wired and available (e.g. for higher [Pi]).
The `UsePiForce` kstiff2 scaling is removed/inert.

## The 2-state weld (why ktr can't also be fit here)
ADP-trap strength (K_D) sets BOTH force and ktr, oppositely: strong trap (K_D≈0.2) → ktr≈0.65 (good) but
force overshoots ~1.4; weak trap (K_D=0.7) → force≈1.18 (good) but ktr≈0.77. The old kstiff2 hack faked a
way around this (strong trap for ktr + fake softening for force). Physically you can't — R2 gates both.
**This is the structural floor (~0.47) for a single-lever 2-state model.**

## Plan to tune further (physical)
**1. ATP-dependent SRX return `kmsrd` — the dedicated kinetic knob** *(fixes the ktr residual)*.
   `kmsr=0`, `UseAtpOnUNR=0` here, so ATP→SRX mobilization is inert; the live lever is `kmsrd` (SRX-ADP
   return). At max-Ca the SRX pool is mechanosensitively suppressed, so `kmsrd` moves **ktr while holding
   force**. Make `kmsrd` [ADP]/[ATP]-dependent (SRX/DRX equilibrium is nucleotide-set) → pull ktr 0.77→0.54
   without disturbing the solved force. This is the one physical lever that breaks the weld in 2-state.
**2. Strain-dependent ADP-trap** — let g2 act only in the force-bearing strain window, decoupling force
   from the near-zero-strain heads that set ktr/onset.
**3. Pin K_D** — 0.7 vs literature 0.194 (Yamashita). The fit wants a weak trap to keep force realistic;
   the "correct" K_D over-traps. That gap is itself the weld — relieved by lever 1, not by K_D.

## Reproduce
`cd(root); addpath(genpath('.')); DriverLowATP_2State`  (needs `data/protocol_03_27_2026_{2,8}mM_slack.mat`).
