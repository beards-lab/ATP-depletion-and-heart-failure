# Low-ATP mechanisms: physiologically-justifiable drivers

Full decision log + validation: [`labdiary.md`](labdiary.md). Reproduce: `run
Analyses/LowATP_Mechanisms/RunMechanisms` (figure `results/mechanisms_compare.png`).

## Question
Which **ATP-justifiable** mechanism reproduces the 2 mM-vs-8 mM signature (force +18% uniform,
ktr ×0.54 sloping, restretch transients up)? Constraints (user): no `ka` lever (attachment isn't
ATP-controlled); a rigor state must be stiff (`kstiff3 ≥ kstiff2`, not the earlier "wobbly" p3).

## Method
Low ATP = simultaneous **↓ATP + ↑ADP + ↑Pi** (coupled). Three nucleotide factors wired into the
ACTIVE pchip path of `dPUdT_CombinedTransitions.m` (previously inert), behind flags, baseline-
preserving at MgADP=0/Pi=0: `UseAdpTrap` (R2·=g2), `UsePiReversal` (R12·=f2, R21·=(1+Pi/K_Pi)),
`UseAtpDetach` (R3·=MgATP/(MgATP+K_T_detach); P3 = rigor, kstiff3=kstiff2, drp3=dr). Driven by
**concentrations**; relative-ratio scoring vs the per-segment target.

## Result — the coupled metabolic state reproduces the force signature
| mechanism | lever | cost | steady | pk2 / vall | ktr | t0 |
|---|---|---|---|---|---|---|
| ADP-trap | MgADP=2 | 0.80 | 1.18 ✓ | 1.05 / 0.98 ✗ | 0.73 | 1.49 |
| Pi only | Pi=4 | 8.9 | 0.92 ↓ | up | up | (tempering term, not a driver) |
| Rigor only | MgATP→2 | 2.54 | 1.04 flat | 1.15 / 1.36 ✓ | 0.93 | 1.59 |
| **COUPLED** | MgADP=1.5, Pi=3, MgATP=2 | **0.555** | 1.16 | **1.17 / 1.24 ✓✓** | 0.64 | 2.10 |

**The coupled state (↓ATP+↑ADP+↑Pi) reproduces the full force-amplitude signature — steady, A,
Am, peak2, vall all +16–24%, including the restretch transients that broke the 2-state wall —
from physiological concentrations alone, no `ka`, no unphysical stiffness.** Division of labour:
ADP-trap → isometric force (traps stiff P2); rigor → restretch transients (stiff P3 bears force
when strained); **Pi is essential** — it tempers the ADP+rigor force overshoot (1.41→1.16). It is
numerically worse than the rejected unphysical fit (0.39) but is the *justifiable* answer.

## Residual = kinetic TIMING (the next frontier)
`t0` (onset 2.1× vs data ~1.3), `peak1` (instantaneous stiffness over), `ktr` (level 0.64 vs
0.54, FLAT vs data sloping 0.64→0.46). The stiff rigor makes onset too slow and the restretch
first-peak too stiff; no scalar concentration produces the ktr **slack-slope**. Candidates:
catch-bond / cooperative restretch detachment (mech E), lattice-spacing SL-dependence (mech F,
targets the slope), or a bounded re-fit of the 8 mM 3-state baseline (peak2 baseline ~106 vs 78).

## Implied physiology
At 2 mM bulk ATP the fit implies local **[ADP]≈1.5 mM, [Pi]≈3 mM** — a plausible energetic-stress
state (elevated ADP/Pi from active ATPase under ATP limitation; the heart-failure metabolic
phenotype).
