# Part C — Coupled ATP mechanism on the CLEAN baseline (opt2state_v2, cost 3.22)

**Date:** 2026-07-08 · worktree `ATP-wt-lowatp` (LowATP + merge of better baseline 236a950).

## Setup
- Baseline: `params/opt2state_v2_opt.m` (2-state regavail, recorded feature cost **3.15**, recomputed
  on merged code **3.22** — merge verified non-destructive; ktr now [56 53 52 50 51] vs data [51 49 49 49 48]).
- Overlay: honest coupled ATP mechanism — `UseAtpDetach + UseAdpTrap + UsePiForce` (NOT PiReversal).
  Constants (getParams defaults): `K_T_detach=4.0`, `K_Pi=4.007`, `K_D=0.194`, `K_T1=0.5`.
- **k3 compensated ×(8+K_T_detach)/8 = ×1.5** so the 8 mM effective k3 = baseline (verified: 8 mM
  steady reproduces the clean baseline exactly).
- 2 mM condition: `MgATP=2, MgADP=1.6, Pi=0.8`. Scored with `atpRatioCost` (relative ratio 2/8 vs data).

## Result: cost 1.34, cleanly localized

| cluster | features | w·cost | verdict |
|---|---|---|---|
| **Force amplitude** | steady .001, peak2 .007, vall_y .011, Am .029, A .048 | **0.096** | ✅ SOLVED — all ~+18%, match data |
| **peak1_y** | model ratio 1.57 vs data 1.25 | **0.507** | ✗ rigor over-amplifies first restretch peak |
| **ktr** | model ~flat 0.6 vs data slope 0.64→0.46 | 0.092 | ~ mean right, slope wrong (2-state wall) |
| **t0 (onset)** | model ratio 2.5–3 vs data 1.0–1.4 | **0.645** | ✗ mechanism over-delays onset, wrong slope sign |

**86% of the cost is the restretch transient (t0 0.645 + peak1 0.507 = 1.15).** The nucleotide
mechanism captures the isometric/steady force gain physiologically; what it gets wrong is the
**restretch ONSET timing and first-peak height** at low ATP.

## Key findings
1. **Force amplitude is a solved problem.** Coupled rigor(AtpDetach)+ADP-trap, tempered by honest
   Pi-force, reproduces the flat +18% force gain across all 5 segments from physiological
   concentrations. steady ratio 1.175 vs 1.179 (cost 0.001).
2. **Pi now HELPS (reversal of the old finding).** On the clean baseline rigor+ADP over-shoot force
   (A/Am ratio 1.22–1.24 > target 1.18); `UsePiForce` (fPi=0.83 at Pi=0.8) pulls it back to ~1.18.
   On the OLD baseline Pi pulled an already-too-low force further down — hence "wrong direction."
   The clean baseline flips this: **honest Pi is beneficial here.**
3. **t0 over-delay is the headline residual.** Model onset ratio 2.5–3× vs data 1.0–1.4×, and the
   SLOPE SIGN is wrong (model decreases with slack, data increases). The rigor+ADP-trap slows the
   P2→P3→PD→reattach onset far more than the data shows. This is over-aggressive kinetic slowing.
4. **peak1 over-amplification.** Accumulated rigor is too stiff on the first restretch peak
   (ratio 1.57 vs 1.25). Same restretch-transient physics as the high-ATP peak1 residual, amplified.
5. **ktr slack-slope** is unreproducible in 2-state mean-field (0.092, structural — 3rd-state / desync).

## Proposed further steps (physiologically bounded)
- **A. Bounded fit of the ATP CONSTANTS only** (`K_T_detach`, `K_D`, `K_Pi`, plus k3 comp) against
  `atpRatioCost`, holding the 8 mM baseline fixed. Low-dim, physiological. Target: pull t0 & peak1
  down without touching force amplitudes. Expected cost 1.34 → ~0.5–0.7.
- **B. Rigor stiffness at low ATP** — `kstiff3` may over-bear peak1; test a modest reduction /
  strain-dependence so accumulated rigor contributes steady force but less transient peak.
- **C. t0 slope-sign** — the mechanism delays onset most at small slack; data most at large slack.
  Points to a slack/strain-dependent onset term (attachment availability during re-lengthening),
  possibly the same lever as the high-ATP t0. Structural piece → 3rd low-force state (breaks ktr
  slope + t0 slope together).
- **D. ktr side quest (classic vs slack)** still pending — SRX timescale separation test.

## Part A (2026-07-08) — constant fit + the decisive force↔onset weld

**Correction to Part C attribution:** `opt2state_v2_opt.m` is **NumberOfStates=2** — there is NO rigor
(P3) state, so `UseAtpDetach` (scales R3=k3·p3) is **completely inert** (a K_T_detach grid gave identical
cost at 4/6/8/12). The low-ATP mechanism here is therefore **`UseAdpTrap` (slows P2 detachment R2 by
g2=(MgATP/K_T1)/(MgADP/K_D+MgATP/K_T1)) + `UsePiForce`**. The peak1/t0 residuals are the ADP-trap, NOT
rigor. `UseAdpTrap` at g2=0.33 ≈ a **physiological k2×0.33** → it inherits the k2-frontier bimodality.

**K_D × K_Pi grid (`results/partA_force_t0_weld.png`, `RunAdpTrapSweep.m`):** best compromise cost
**0.69** (K_D=0.35, K_Pi=6), down from 1.34. But the sweep exposes the wall (K_Pi=4 slice):

| K_D | steady ratio (force) | t0 ratio (onset) | peak1 ratio |
|---|---|---|---|
| 0.10 | 1.28 | 4.2 | 1.9 |
| **0.194** | **1.18 ✓** | 2.75 | 1.57 |
| 0.35 | 1.09 | 2.0 | 1.35 |
| 0.60 | 1.02 | 1.55 | 1.19 |
| **1.0** | 0.96 | **1.29 ✓** | 1.09 |

*(targets: force 1.18, onset ~1.26)*. **Force and onset are welded to one knob and anti-correlate.**
The trap strength that gives +18% force (K_D≈0.19) forces onset to 2.75× (data 1.26); the K_D that fixes
onset (≈1.0) gives ZERO force gain. They want K_D **5× apart** — no constant combo satisfies both.

**Structural conclusion (the headline).** A 2-state cross-bridge cycle uses the SAME detachment step (R2)
to (a) hold force and (b) gate redevelopment onset — so raising low-ATP force via ADP-trapping *necessarily*
over-slows onset. Data demands force UP with onset only mildly changed → the cycle cannot do both. The
resolution is a **3rd low-force, slow-detaching state** (rigor/A.M-like) that bears force WITHOUT gating
the P1↔P2 turnover, decoupling force from onset. This is [[mech-tradeoff]]'s "3rd slow-detaching state",
now proven necessary from the ATP data itself (not just the 8 mM iron law). Next: build/fit a 3-state
baseline (`params/params_3state_DRAFT.m`) and re-run the coupled ATP with rigor engaged (UseAtpDetach live).

## Reproduce
`cd(root); addpath(genpath('.'))`, run `Analyses/LowATP_Mechanisms/RunCoupledCleanBaseline.m`
(needs `data/protocol_03_27_2026_{2,8}mM_slack.mat`). Outputs `results/lowatp_partC.mat`.
Part A sweep: `RunAdpTrapSweep.m` → `results/partA_force_t0_weld.png`.
