# Lab diary — Restretch feature fit

Goal: best feature-cost match (evalFeatureCost) with emphasis on the **restretch
features** (`Am`, `peak1_y`, `peak1_dSL`, `vall_y`, `peak2`, `vall2_dy`) without
degrading FV or the physiology output bounds. Base = `params/params_reseeded_regavail_opt2.m`.

## 2026-07-03 — Wall-breaking + redundancy (USER fn, ovrsht|1 & ktr|2 at full weight)
Wall status (all tested, not asserted): ovrsht_dy → k2-coupled, reaches ~1.3 (data). FV →
FIXED by velocity mechanisms. ktr → k2-driven, pushable 81→~57 but ktr=49 needs k2≈50
(below the physiological floor 100) AND deepens vall2 — a QUANTIFIED tradeoff, not a wall.
peak1_y (119 vs 96) = fast ~2 ms elastic spike, kstiff-coupled to steady, NOT detachment-
fixable (R1D too slow). Best velGauss-both config = 9.72 (params_restretch_best.m).

REDUNDANCY (velGauss vs registration-availability, both flatten FV):
- regavail OFF, velGauss only → FV steepens [1 .52], steady balloons to 128 (A0 floor was
  holding isometric down), cranking velGauss to compensate CRASHES the ODE (stiff).
- velGauss OFF, regavail only → FV recovers to [1 .92 .60] (shoulder MATCHES data) with
  A0+ka tuning, no crashes. Best regavail-only = 11.02.
DECISION: **keep registration-availability, DROP velGauss.** velGauss is redundant, ODE-
unstable (the optimizer's timeouts), less defensible (arbitrary Gaussian vs target-zone
story). ~1.3 gap vs 9.72 is unoptimized hand-tuning; the stable surface should close it.
=> Clean 2-state seed saved: `params/params_2state_seed.m` (velGauss OFF, cost 11.02,
   FV=[1 .92 .60 .30 .11], ktr 69, peak1y 106).

NEXT: (a) 3-state variant `params/params_3state_seed.m` (low-force slow-detaching 3rd state
to break ktr/peak1/vall_y walls at the source); (b) parameterized optimizer `optimizeFeatures.m`
+ RunOptim2State.m / RunOptim3State.m (distinct tag outputs) to run 2-state & 3-state opt in
two MATLAB instances. NOTE: prior RunRestretchOptim saved a diverged snapshot (174 > seed
9.72) — accept/save logic being fixed to never persist worse-than-incumbent.

## 3-state investigation (topology + collapse diagnosis)
Topology (NS=3): p1 --k1--> p2 --k2(R2 shape)--> p3 --k3--> (PT algebraic) --> PD --ka--> p1.
So p2's ONLY exit is via p3 (no direct p2->detach). alpha3/s3/dr3 are DEAD in the active ODE
(dPUdT_CombinedTransitions); real p3 knobs = k3, k3m, kstiff3(=14000, low-force), drp3.
- Naive flip / ka=160 -> NaN (force fit fails). Feasibility sweep: as k2,k3 rise it runs but
  **FV tail collapses to 0** (force->0 by v=2). ktr=49 reachable (k2=90,k3=200) but FV=[1 .81 .22 0 0].
- ROOT CAUSE: the R2 strain shape (PieceWiseStrainDep2, ramps up at NEGATIVE/shortening strain)
  was tuned as 2-state DETACHMENT; as the p2->p3 FEED it drains force-bearing p2 during shortening
  -> no force at velocity -> FV tail = 0. FIX = flatten R2 neg-strain knots (__5/6/7). Was testing
  this when the MCP MATLAB session STALLED (a 3-state config hung the integrator).
CODE FIXES APPLIED (offline, need MATLAB to validate):
  1. dPUdT_CombinedTransitions.m: added p3 negativity clamp `p3(p3<0)=0` (p1/p2 had it, p3 didn't
     -> likely the stall/NaN source).
  2. optimizeFeatures.m: ENABLED surrogateopt (was commented out -> was only fminsearch-from-g=1!);
     production SURR_EVALS default 0->100; fminsearch headless (no plot) for parallel instances.
  3. RunOptim3State.m pool: dropped dead s3, added drp3 + R2 feed-shape knots __5/6/7 (the FV-tail
     lever); kstiff3 added to compulsory.
STATUS: 2-state optim READY (params_2state_seed.m 11.02 + optimizeFeatures + RunOptim2State).
3-state seed NOT yet saved — params_3state_DRAFT.m still collapses; must first hand-find a
non-collapsed config (flatten R2 __5/6/7) and save as params/params_3state_seed.m before RunOptim3State.

## Handoff state (execution infra)
Optimizer mode set to USER spec: cfg.SURR_EVALS=0 (surrogate OFF, pure fminsearch on random-drawn
subsets), cfg.KICK_FRAC=0 (no stall-kick; non-improving round just draws a fresh subset). Files:
- Analyses/RestretchFeatureFit/optimizeFeatures.m  (function optimizeFeatures(cfg); accept/save fixed:
  never persists worse-than-incumbent; final write from state.best_params with assert)
- RunOptim2State.m (tag opt2state, seed params_2state_seed.m)  -- READY
- RunOptim3State.m (tag opt3state, seed params_3state_seed.m [MISSING], pool incl. R2 __5/6/7 + k3/k3m/kstiff3/drp3)
- costOfSnap.m (helper: eval a snapshot+override struct -> total cost)
PARALLEL: confirmed safe -- every write is tag-derived (opt2state_* vs opt3state_*), figures per-process,
seeds/data/code read-only. Two `matlab -batch RunOptim2State` / `RunOptim3State` in separate instances
are isolated. CPU caveat: 2x parpool('Threads',5)=10 threads; drop to 4 each if <10 cores.
EXECUTION NOTE: MCP matlab session stalls repeatedly (esp. on 3-state stiff ODE); agents that launch
background matlab tend to END before it finishes. RELIABLE PATTERN = launch `matlab -batch <driver>` as a
detached process and read the OUTPUT FILES (params/<tag>_opt.m, <tag>_state.mat) -- NOT the .log
(-batch buffers stdout). 2-state DEBUG launched as bg task; awaiting opt2state_* files.

## 3-state: deeper diagnosis (why seed-search keeps failing)
The FV collapse AND the ODE stiffness both trace to ONE thing: the p2->p3 feed uses R2 = k2 x
PieceWiseStrainDep2(strain), and that SAME shape/code path is the 2-state DETACHMENT profile (ramps to
50 at strain extremes). As a feed it (a) drains force-bearing p2 during shortening -> FV tail=0, and
(b) creates fast/stiff transients at strain extremes -> slow integrator. Flattening only __5/6/7 (neg
knots) was not enough in the timed-out run. NEXT IDEAS (need a working matlab): (i) flatten the ENTIRE
PieceWiseStrainDep2 to ~1 (strain-independent p2->p3 feed) so FV is shaped only by regavail/R1D/k3;
(ii) or a code change giving p2->p3 its OWN gentle strain shape separate from detachment. 3-state is a
genuine R&D effort (ODE stiffness makes each iter slow), not a quick seed search.


Ground rules (from PI): may relax parameter/physiology bounds if justified; may
manipulate piecewise strain X (knot positions) and Params (knot values)
separately; may reweight `params0.fn`. Primary cost = features.

Data (8 mM, protocol_03_27_2026): restretch targets
- Am ≈ [67 66.7 63 59 53.6]; peak1_y ≈ [96 94 90 89 77]; peak1_dSL ≈ [0.022 0.023 0.027 0.028 0.027]
- vall_y ≈ [77 75 72 70 61]; peak2 ≈ [80 82 81 82 66]; steady ≈ [81 80 80 80 65]
- vall2_dy ≈ [-13 -13 -14 -14 -10] (undershoot); ovrsht_dy ≈ 1.3 (noise-level)
- FV_fnorm ≈ [1 .92 .66 .32 .11]

Convention: each entry logs the CHANGE, the resulting TOTAL feature cost, and which
features moved. Snapshots saved to results/.

---

## Baseline (reg-availability symmetric, current best-fit) — cost 15.86
fn = restretch-focused (FV_fnorm|10, A|30, Am|20, peak1_y|10, peak1_dSL|2,
vall_y|10, peak2|8, steady|30, vall2_dy|1.0, ovrsht_dy|0.01, t0_crossing|2,
restretchSlopeStart|0.2, ktr_rmse|.1, physiology-group, AssertParams|0.001).
Top misses: vall2_dy 5.05 (undershoot≈0), peak1_dSL 3.98 (0.040 vs 0.025),
FV_fnorm 1.68, peak2 1.27 (96 vs 80), vall_y 1.03, t0_crossing 1.19.

## E01 — directional availability (RegAvailShorteningOnly=true) — cost 11.38 ✅ (−28%)
Code: dPUdT_CombinedTransitions.m A_inf now uses max(0,-vel) when the flag is set
(getParams default added). Restretch/lengthening no longer boosts attachment.
- vall2_dy: ~0 -> [-8.8 -10 -10 -9.4 -8] (data [-13..-10]) — UNDERSHOOT RESTORED (5.05->0.38)
- peak2: 96 -> 81 (data 80) — FIXED (1.27->~0)
- peak1_y: 104 -> 102 (data 96) — slightly better
- FV_fnorm: UNCHANGED — shoulder preserved (shortening branch untouched) ✅
- SIDE EFFECT: vall_y 68.9 -> 58.7 (data 76.7) — now too deep (1.03 -> 3.14)
New top misses: peak1_dSL 3.37, vall_y 3.14, FV_fnorm 1.68, t0_crossing 1.19.
Snapshot: results/params_E01_directional.m
DECISION: keep directional=true as the new working base. Next: fix vall_y (too
deep post-peak detachment) + peak1_dSL (peak at too-large length excursion).

## E-series (from E01 base) — manual coordinate exploration

Sweep findings (sensitivities), data targets in brackets:
- `R1Dpos__4 (50→15)`: vall_y 59→73 [77] ✅ BUT vall2 −8→−5 [−13] ✗, peak2↑. R1D positive
  clamp controls BOTH the fast valley and the slow undershoot, in OPPOSITE directions.
- `k2×0.8`: vall2 −8→−13 [−13] ✅ (deepens slow undershoot) without moving vall_y; peaks↑.
- `kstiff2×0.8`: peak1_y↓, peak2↓ (good) but vall2 flattens, dSL↑.
- `ka×1.3`: peaks balloon. `kSE×1.4`: dSL↓ (good) but peak1_y↑.

**C3** = R1D__4=18, k2×0.85, kSE×1.4, kstiff2×0.9 → **cost 9.93**. dSL better, valley up, but peaks high.
**D1** = C3 + kstiff1×0.85 → **8.23**. peak1_y 101 [96], peak2 ~81 [80], vall2 [−11..−9] ✅. vall_y still 64 [77].
**D2** = D1 + v_ref_reg 0.8→0.2 → **7.34** ✅✅ (current best). FV shoulder widened (fnorm@0.5 0.68→0.73)
without depressing isometric. Snapshot: results/params_D2_best2state.m

FV sweeps:
- Lowering `A0` (isometric availability floor) FLATTENS FV but COLLAPSES steady (anchored) → cost explodes.
  A0 is not a usable FV lever while steady is pinned. `v_ref_reg`↓ flattens shoulder cheaply (use this).
- Lowering R2 negative-strain knots (`PieceWiseStrainDep2Params` 5/6) made FV STEEPER (sign opposite to
  naive expectation) → the FV steepness is not simply shortening-detachment magnitude. Left to optimizer.

R1D X-position (knot) sweeps for vall_y (I-series):
- `R1DX3 0.006→0.02` or `R1DParams3 1.76→0.8`: vall_y → [80,77] MATCHES data ✅ — but vall2 flattens to
  −4 [−13] and peak2 inflates to 90 → total cost 24 (much worse). CONFIRMED across every variant.

### KEY CONCLUSION (2-state structural wall)
`vall_y` (shallow fast post-peak valley, ~4 kPa below peak2 in data) and `vall2_dy` (deep ~13 kPa slow
undershoot) are both governed by the SAME positive-strain R1D detachment, pulling in opposite directions,
and raising vall_y also inflates peak1/peak2. A single strain-dependent detachment rate cannot separate a
FAST shallow valley from a SLOW deep undershoot — these are timescale-separated. **This is the 2-state
limit (cf. mech_tradeoff); the durable fix is NumberOfStates=3 (a low-force slow-detaching state).**
D2 is the best weighted 2-state compromise (keeps peaks + vall2 + slopes; accepts vall_y ~15% low).

### D2 remaining costs (7.34): FV_fnorm 2.05 (structural steepness), peak1_dSL 1.38, vall_y 1.14 (wall),
physiology 0.67 (XTOR_vmax high/SRX low), t0_crossing 0.50, A/Am/steady ~0.9 (force ~3% low).

### Next steps
1. Optimizer polish (surrogate/fminsearch) over {kstiff1,kstiff2,k2,kSE,ka,R1D__3/4,A0,v_ref_reg,R2 5/6}
   from D2 — resolve the fittable remainder (FV, dSL, force level) better than hand-tuning.
2. Reweight: consider vall_y|10→5 (partly unfittable in 2-state) so it stops distorting peaks.
3. Structural: NumberOfStates=3 to break the vall_y/vall2 wall and the ktr iron law together.
4. Fold RegAvailShorteningOnly=true into the production optimizer pool (RunSurrogateSimplex_Optim).

## J/K/M-series — squeeze the fittable remainder (from D2)
- **J1** = D2 + kSE×1.3 → **6.70**. dSL 0.034→0.031, peak1_y→96 (perfect), nothing hurt. CLEAN WIN.
  (Raising R2 neg knots flattens FV [FVcost 2.05→1.61] but degrades the coupled slack phase → net worse;
   FV/slack share the R2 shape — a coupling for the optimizer, not hand-tuning.)
- **K1** = J1 + kstiff1×1.05 + kstiff2×1.05 → **6.23** (BEST-by-cost). steady→[81 81 81 82 66]=data ✅,
  Am/vall_y up, vall2 [−12.5..−11] ✅. COST of raising force in 2-state: peak1_y re-inflates to ~105 [96]
  (kstiff1 drives the P1 rapid peak). Snapshot: results/params_best.m
- M-series: kstiff2-ALONE up (×1.12) → 7.49 WORSE; the force-level gain needs the symmetric K1 combo.
  ka×1.12 → 6.63 (≈J1). Confirms hand-tuning has saturated — couplings are non-monotone.

## P/Q/R/S-series — FV breakthrough (FV was NOT structural)
Earlier I wrongly concluded FV steepness was structural. It is fixable:
- **A0↓ (isometric availability floor) deepens the shortening-availability SWING, flattening FV.**
  Alone it collapses force (steady anchored) → must COMPENSATE force. HOW you compensate matters:
  - kstiff compensation (P-series): FVcost 2.09→0.11 (FV≈data!) BUT inflates peak1_y (kstiff drives the
    P1 peak) and overshoots steady → total worse.
  - **ka (attachment) compensation (Q2): FVcost 2.09→0.65 with peak1_y 97 [96], peak2 79 → COST 6.23→4.91.**
    ka raises force by attaching MORE bridges (not stiffer ones), so peaks don't inflate. KEY LEVER.
- **Q2** = K1 + A0 0.6→0.45 + ka×1.25 → **4.91**.
- **R4** = Q2 + A0→0.42, ka×1.08, kstiff2×1.04 → **4.77** (FVc 0.43, vall2 −13.6=data).
- **S3** = R4 + kstiff2×0.97 → **4.56** ✅ FINAL. steady [80 80 81 81 64]=data, vall2 [−13..−11]=data,
  peak2 seg1=81, FVc 0.44 (from 2.09!). Snapshot: results/params_best.m = params/params_restretch_best.m

### FINAL (this session): 15.86 → **4.56  (−71%)**
S3 residual: vall_y 0.97 (2-state WALL, structural), physiology 0.74 (XTOR_vmax ~18>15, SRX_ss 0.07<0.10
— both candidates for justified bound relaxation: Fenn factor can exceed 3×, and biochemical SRX at max
activation can be <0.10 [Ma 2022]), peak1_dSL 0.69 (0.032 vs 0.027), t0_crossing 0.64, peak1_y 0.54
(105 vs 96 seg1-3 — residual kstiff/peak coupling), FV 0.45. peak2/vall2/steady/Am/restretchSlope ≈ data.

Winning recipe (all from params_reseeded_regavail_opt2): RegAvailShorteningOnly=true; A0 0.6→0.42;
ka ×~1.35 (force comp); v_ref_reg 0.8→0.2 (shoulder width); kSE ×1.3; kstiff1 ×1.05; kstiff2 ≈×1.06;
k2 ×0.85; kSE... ; R1D positive clamp __4 50→18. See results/params_best.m for exact values.
