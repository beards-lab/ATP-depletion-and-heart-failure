# Lab diary — physiologically-justifiable low-ATP mechanisms

Re-do of the mechanism search with ONLY ATP-justifiable levers (user constraints: `ka` is not
ATP-justifiable; `kstiff3 < kstiff2` is unphysical — a rigor state must be stiff). The prior
3-state win (`../LowATP_3rdState`, cost 0.39) used a "wobbly" low-stiffness p3 + ad-hoc `ka`
reduction — rejected on physiology. Scoring unchanged: relative ratio vs per-segment target
(`atpRatioCost.m` + `../LowATP_k2Frontier/results/atp_target.mat`).

## Decisions & progress (2026-06-29)

### Part 1 — infra (Sonnet agent)
- `Auxiliary/loadParams.m` (sandbox loader — kills the `params0` variable-name footgun: every
  `params/*.m` hardcodes `params0`, so a driver naming its struct anything else silently kept
  defaults — the `base=` NaN bug). `Auxiliary/refreshPool.m` (restart pool so workers pick up
  edited code). Verified: `foo=loadParams(...)` loads overrides under a non-`params0` name.
- Also fixed earlier: dup-timestamp `interp1` crash in `Model/extractSlackAttributes.m`.

### Part 3 — wired nucleotide factors into the ACTIVE pchip path (were inert / dead branches)
`Model/dPUdT_CombinedTransitions.m`, behind flags, all **baseline-preserving** at (MgADP=0,Pi=0):
- `UseAdpTrap`: `R2 *= g2` (g2=ATP-ready fraction; ↑ADP traps force-bearing P2).
- `UsePiReversal`: `R12 *= f2`, `R21 *= (1+Pi/K_Pi)` — Pi is a power-stroke PRODUCT → inhibits
  forward stroke + drives reverse → **lowers force** (Cooke&Pate). [First tried `f2` on R2
  (detachment) → force went UP (wrong); moved it to R12 (forward stroke) → force down ✓.]
- `UseAtpDetach`: `R3 *= MgATP/(MgATP+K_T_detach)` — ATP-limited rigor P3→detach; low ATP →
  rigor piles up. P3 = rigor: `kstiff3 = kstiff2` (physiological min of ≥), `drp3 = dr`.
- getParams defaults: `UseAdpTrap/UsePiReversal/UseAtpDetach=false`, `K_T_detach=4`.
- **Validation (single levers, baseline-preserving when off):** flags-off = exact baseline ✓;
  ADP-trap MgADP=2 → steady +18%, ktr −27% ✓✓; Pi=4 → steady −8% ✓; rigor MgATP 8→2 →
  steady +11%, peak2/vall up (restretch) ✓.

### Part 4 — mechanism scoring (concentration-driven; `RunMechanisms.m`)
| mechanism | lever | cost | steady | pk2/vall | ktr | t0 |
|---|---|---|---|---|---|---|
| ADP-trap | MgADP=2 | 0.80 | 1.18✓ | 1.05/0.98 ✗ | 0.73 | 1.49 |
| Pi only | Pi=4 | 8.9 | 0.92 (DOWN) | up | up | — (tempering term, not a driver) |
| Rigor only | MgATP→2 | 2.54 | 1.04 (flat) | 1.15/1.36 ✓ | 0.93 | 1.59 |
| **COUPLED** | MgADP=1.5,Pi=3,MgATP=2 | **0.555** | 1.16 | **1.17/1.24 ✓✓** | 0.64 | 2.10 |

**RESULT:** the **coupled metabolic state (↓ATP + ↑ADP + ↑Pi together)** reproduces the FULL
force-amplitude signature — steady/A/Am/peak2/vall all +16–24%, INCLUDING the restretch
transients that broke the 2-state wall — from **physiological concentrations alone** (no `ka`,
no unphysical stiffness). Division of labour: **ADP-trap** → isometric force (traps stiff P2);
**rigor** → restretch transients (stiff P3 accumulates, bears force when strained); **Pi** is
ESSENTIAL — it tempers the ADP+rigor force overshoot (1.41→1.16). Numerically 0.555 > the
unphysical 0.39, but it is the *justifiable* answer.

**RESIDUAL (next frontier) = kinetic TIMING:** `t0` (onset 2.1× vs data ~1.3), `peak1`
(instantaneous stiffness over), `ktr` (level 0.64 vs 0.54, FLAT vs data sloping 0.64→0.46). The
stiff rigor makes onset too slow + restretch too stiff, and no scalar concentration produces the
ktr slack-SLOPE. Candidates: catch-bond / cooperative restretch detachment (mech E), lattice-
spacing SL-dependence (mech F), or refining the rigor detachment kinetics / a bounded re-fit of
the 8 mM 3-state baseline (peak2 baseline ~106 vs data 78).

**Implied physiology:** at 2 mM bulk ATP the fit implies local [ADP]≈1.5 mM, [Pi]≈3 mM —
a plausible energetic-stress state (elevated ADP/Pi from active ATPase under ATP limitation).

### 2026-06-29 (cont) — raw traces + residual-fix attempts
- **Raw force-time traces** (`PlotRawOutputs.m`; `results/raw_trace_{8mM,2mM_coupled}.png`): the
  steady plateaus match data well (the amplitudes the mechanism explains); the restretch PEAKS
  overshoot (sharp elastic spikes, model ~120-160 vs data ~80-95). This overshoot is **pre-existing**
  (the 2-state best already overshot peak1 ~10%), amplified by the stiff rigor — a baseline-fit
  issue separate from the ATP question (the relative-ratio amplitude result is unaffected).
- **Re-fit attempt (transient rigor, k3=2500-4000):** WORSE (cost 0.81-1.15 vs 0.555); baseline
  peak2 stays ~100 (the overshoot is p2, not rigor). k3=1200 coupled (0.555) remains best.
- **Maxwell dashpot (`UseMaxwellDashpot=1`, default kSE_M=100/eta_M=0.1):** negligible
  (0.568→0.555) — default damping too weak to tame the peak1 spike.
- **Conclusion:** the t0/peak1/ktr-slope residual is STRUCTURAL — not fixed by rigor tuning,
  transient re-baselining, or the default dashpot. Real fixes (next phase): **tuned** viscoelastic
  damping of the restretch (kSE_M/eta_M or UseViscoelasticSE), catch-bond / cooperative restretch
  detachment, and SL-dependence for the ktr slope; plus a dedicated bounded re-fit of the
  pre-existing restretch-peak overshoot. The coupled metabolic mechanism is the settled answer to
  *which nucleotide effects drive the ATP difference* (force/amplitude); the *restretch kinetics*
  are a distinct, pre-existing model-fit frontier.
