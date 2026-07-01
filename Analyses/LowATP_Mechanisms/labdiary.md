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

### 2026-06-30 — Part A/B: restretch trade-off diagnosed (kSE / kstiff2 / A2-shift)
8 mM baseline vs data: `restretchSlopeStart` 1256(2-state)/1285(rigor) vs **data 1588** (LOW → kSE
under-set, confirmed — it's frozen during the `ka=0` passive fit); `peak1` 97/114 vs 89; `peak2`
92/106 vs 78; `ktr` 66/62 vs 49.
- **kSE sweep** (fig `results/restretch_tradeoff.png`): kSE×1.5 (≈4800) matches the data slope
  (1611) BUT raises `peak1` (114→117) AND speeds `ktr` (66→79). kSE↑ fixes the slope, worsens
  peak1 + ktr.
- **kstiff2↓** reduces peak1 (the velHS lever) but drops `steady` (it is the force scale): ×0.70 →
  peak1 87 ✓ but steady 60 (data 77).
- **c_SE_visc** (`UseViscoelasticSE`) damps peak1 (88→74) + improves vall/ktr BUT **crashes
  restretchSlope** (1504→268).
- **`UseA2AttachmentShift` AMPLIFIES, not damps** — peak1 99→110, vall 86→102, peak2→119. The
  implementation re-attaches strained heads at a lower-strain site (they stay bound) → traps force
  (matches RestretchMechanisms "valley fills in"). The "hop off to relieve" intuition doesn't match
  the code (it re-binds). RULED OUT as a damper.
- **Verdict:** the restretch peak is a genuine 4-way trade-off (slope↑kSE, peak1↓kstiff2,
  steady↑kstiff2, ktr) — no single lever fixes it. Needs the full `RestretchMechanisms`-style
  bounded re-optimization (free ka/kd/k1/k2 to compensate) = scoped Part B-full. The pre-existing
  **ktr-too-fast (66 vs 49)** is the deeper coupled residual underneath.

### 2026-06-30 (cont) — restretch deep-dive: valley / overshoot / low-ATP kSE
Zoom diagnostic — model restretch = a TALL spike (peak1 115@8mM / 163@2mM vs data ~76/95) + a
RINGING 2nd bump that settles ABOVE steady; data = small spike + smooth recovery.

![restretch zoom](results/zoom_restretch.png)

**Thread 2 — the post-restretch peak (`ovrsht_dy`):** the visible artifact is the 2nd bump =
`peak2` (106) **≫ `steady`** (80); the data has peak2≈steady (78≈77). It's a damped-oscillation
overshoot from synchronized DRX reattachment after the high-strain-detachment valley. Feature-wise
the model's `ovrsht_dy`/`vall2_dy` ≈ 0 (over-damped settle) while the data has real ± dynamics — but
the dominant eye-catching error is peak2≫steady. Fix = bring the spike down + damp/spread the
reattachment.

**Thread 1 — k2-strain-dep + A2-shift (user is right):** my earlier "A2-shift amplifies" was tested
WITHOUT first resolving the spike. The correct sequence: strong **high-strain k2-detachment** (the
R2 / `PieceWiseStrainDep2` curve) removes the spike → deep valley; **A2-shift then fills the valley**
(strained heads hop to a lower-strain site, STAY attached → sustained plateau) = trades the sharp
peak for duration. NEXT EXPERIMENT: steepen `PieceWiseStrainDep2` at high +s (or a velocity-dependent
R2 slip bond) to drop peak1 to ~data, then enable `UseA2AttachmentShift` (slope sweep,
`s_threshold_R=0.0046` locked) to fill the valley; target peak1≈90, shallow valley, plateau≈steady.

**Thread 3 — does low ATP change serial elasticity?** Data `restretchSlope` ratio 1.26 vs model
1.10. Experiment = scale kSE at 2 mM: **kSE×1.2 matches restretchSlope (1.24≈1.26)** — but most
likely this is *more attached rigor XBs* stiffening the series path (attached XBs sit in series),
not a true kSE change. CRUCIALLY **t0 is DECOUPLED**: t0 ratio stays ~2.1 across all kSE (data 1.31).
So "kSE influences t0" does NOT hold in the model — t0 is a rigor-onset/kinetic delay, not a
compliance effect. The model under-stiffens at 2 mM (1.10 vs 1.26) → rigor pool/stiffness slightly
under-represented; closing restretchSlope is easy, but t0 needs a kinetic fix.

![low-ATP kSE vs t0](results/lowatp_kSE_t0.png)

### 2026-06-30 (cont) — Thread 1 result: R2-strain-dep can't isolate the restretch peak
Steepening `PieceWiseStrainDep2` at the restretch zone (knot3, s≈0.026) drops peak1 (114→83) BUT
**crashes steady (80→56→42) and explodes ktr (62→157)** — raising R2 anywhere speeds the whole
cycle, and the restretch-peak strain **overlaps the isometric force-bearing strain**. A2-shift can't
rescue a crashed baseline. ⇒ R2-strain-dependence does NOT isolate peak1 from steady/ktr. With
kSE/kstiff2/c_SE_visc/A2-shift/R2-strain-dep ALL trading off, **no single lever fixes the restretch
peak → it is structural (mean-field synchrony)**. Full synthesis written to `SYNTHESIS.md`.

### 2026-07-01 — Pi without the ktr artifact (user: "Pi pulls it the wrong direction")
Confirmed: the model's Pi wiring (R12↓ forward + R21↑ reverse stroke) tempers steady but **SPEEDS
ktr** (0.55→0.63) — a spurious kinetic side effect (the reverse stroke fast-equilibrates P1↔P2).
Tries:
- **Pi-FREE (ADP+rigor rebalanced):** steady overshoots 1.24–1.35 (ADP-trap iron law, nothing
  tempers); ktr nailed 0.54 at MgADP=1.4; cost 1.03. Key: **rigor does NOT slow ktr** (0.95) — the
  ktr slowdown comes ENTIRELY from ADP-trap, so it can't be decoupled from steady without a temper.
- **FIX 1 — ktr-neutral Pi** (new flag `UsePiReverseStroke=false`, now default): Pi via forward-
  stroke inhibition only (R12↓). ktr stays honest (0.56 not 0.63) but tempers steady weakly (1.23);
  cost 1.07.
- **FIX 2 — Pi-as-force** (new flag `UsePiForce`): Pi weakens P2 binding force
  (kstiff2_eff = kstiff2/(1+Pi/K_Pi)), ktr-neutral, spares P3 rigor. At **MgADP=1.6, Pi=0.8: steady
  1.20 ✓, ktr 0.56 ✓ (BOTH honest), cost 0.98.** Residual now = vall/peak1/t0 (the structural
  restretch / rigor-onset walls) — no longer MASKED as the old Pi did. (Pi≥3 over-tempers: steady
  crashes <1.0, since fPi=0.57 at Pi=3.)
**Conclusion:** user was right — the old Pi faked ktr. The honest Pi (force reduction) fits
steady+ktr cleanly; the higher cost vs 0.555 is because it exposes rather than masks the restretch/
timing residuals. New flags: `UsePiReverseStroke` (default off = ktr-neutral), `UsePiForce` (P2
force reduction). Best HONEST coupled = ADP-trap + rigor + `UsePiForce` (MgADP 1.6, Pi 0.8).
