# Lab diary — Low-ATP mechanism fit + 3rd cross-bridge state

Goal: reproduce the measured 2 mM-vs-8 mM signature (relative-ratio scoring, per-segment
target from `../LowATP_k2Frontier/results/atp_target.mat`). Two phases:
- **P1 — quantify the 2-state wall:** bounded multi-lever fit (k2, ka, k_1, kstiff…) as ATP
  multipliers; find the best 2-state cost and pin the `peak2` residual as a hard number.
- **P2 — break the wall:** add a 3rd cross-bridge state (low-force, slow-detaching post-stroke)
  and re-fit; target = let `peak2` rise with `steady` while `ktr` slows.

Scoring recap (from `../LowATP_k2Frontier/`): force-amplitude features flat ×1.18 (scalar
targets); ktr (0.64→0.46) & t0 (1.01→1.41) slope with slack (per-segment vector targets);
restretchSlopeStart/ktr2_overshoot noisy (de-weight). `model(2mM)/model(8mM)` vs data ratio.

Prior finding (lever map): k2 owns ktr+force but overshoots force; the best 2-lever is k2+ka
(steady/A/ktr on target) but **peak2 collapses to 0.96** (data +18%). peak2 is the binding
residual → motivates the 3rd state.

---

## Decisions log

### 2026-06-26 — kickoff
- New folder `Analyses/LowATP_3rdState/` holds both phases + this diary.
- Decision: relative-ratio cost in **log-space**, weights from `atp_target.weights`
  (steady/A high, ktr-vector high, t0 medium, noisy ~0). Vector features scored per-segment.
- Open question: does the existing `NumberOfStates=3` path give a usable low-force slow-detaching
  state, or do I need to re-wire p3 force/rates? Investigating before committing.
- **Finding (infra):** the 3-state path ALREADY exists. Force = `kstiff3·p3_1 + kstiff2·p2_1 +
  kstiff1·p1_1`; cycle p1→p2→p3→PT with `R2`=p2→p3 (rate k2), `R3=k3·p3`=p3→detach (k3=44.3,
  k3m=0), `kstiff3=14000` (p3 is lowest-force). Current best = `NumberOfStates=2`. So enabling
  3-state + ATP-limited final detachment (k3↓) should pile heads into a low-force slow-detaching
  p3: modest isometric force (low kstiff3) but attached during restretch (keeps peak2/vall_y up)
  while the cycle slows (ktr↓). Decision: reuse the existing p3 infra first (try-both: custom p3
  only if it fails).
- **Key simplification:** relative-ratio scoring means the 3-state 8 mM baseline need NOT be a
  perfect refit — only sane. Score the 2 mM/8 mM ratios. So Phase 2 = (a) sane 3-state baseline,
  (b) ATP lever via k3 (±k2), (c) check it breaks the vall_y/peak2/ktr wall.

### 2026-06-26 — Phase 1 result: the 2-state wall
Reusable cost `atpRatioCost.m` (log-space, per-segment vectors, weights from `atp_target`).
Grid over (k2,ka) ATP-multipliers, relative-ratio cost vs per-segment target.
**2-state floor cost = 0.704 @ k2×0.62, ka×1.00.** Saved `phase1_2state.mat`. Residuals (w*cost):
- **vall_y 0.287** — model 1.01 vs data 1.20 (restretch valley flat; data +20%). DOMINANT.
- **ktr 0.216** — model flat ~0.71 vs data sloping 0.64→0.46 (under-slowed AND no slack-slope).
- A 0.078, peak2 0.031 (model +9% vs data +18%), steady/Am ~0.03 each (well fit), t0 0.019.
Complementarity seen: ka↓ trims steady but kills peak2/vall_y; Pi (k_1↑) raises peak2/vall_y to
~1.15-1.18 but kills steady (→1.05). 2-state can fit {steady,A} OR {peak2,vall_y}, not both.
**Conclusion:** the binding residuals are restretch-transient forces (vall_y, peak2) +
ktr(level+slope). These motivate a slow-detaching attached state (p3) that bears force during
restretch without inflating isometric steady. → Phase 2.

### 2026-06-26 — Phase 2: 3-state breaks the restretch-transient wall
Reused existing p3 infra (`NumberOfStates=3`). Untuned baseline broke (ktr 6, A 161) — p2→p3→PT
is a slow 2-step bottleneck. Re-baselined with **k3=600** (fast p3→detach so p3 stays transient
at 8 mM): steady 77, A 62, vall 74, ktr 63 (ktr a bit high but relative scoring tolerates).
- **k3↓ ALONE overshoots:** raises vall/peak2 strongly but barely raises steady (p3 is near-zero-
  strain at isometric → a *restretch*-force state). So p3 ≠ isometric force.
- **2-lever within 3-state (k2↓ steady, k3↓ restretch+ktr) WORKS:** grid min **cost 0.448 @
  k2×0.65, k3×0.50** — beats the 2-state wall 0.704 by 36%. steady 1.19✓, A 1.21✓, **peak2 1.15✓
  (was the wall), vall 1.31** (now slightly over, was 1.01 under). Restretch-transient residual
  SOLVED.
- **New dominant residual = ktr** (w*c 0.209): model flat ~0.71, data slopes 0.64→0.46 — under-
  slowed + no slack-slope. Also t0 too slow (1.7 vs 1.0-1.4). Decision: add a 3rd lever (k1/ka)
  to slow ktr; the ktr *slope* (slack-dependence) may need a strain-dependent ATP effect (try
  later). Saved nothing yet — exploring before locking a script.
