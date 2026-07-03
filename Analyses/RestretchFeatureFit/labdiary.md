# Lab diary — Restretch feature fit

Goal: best feature-cost match (evalFeatureCost) with emphasis on the **restretch
features** (`Am`, `peak1_y`, `peak1_dSL`, `vall_y`, `peak2`, `vall2_dy`) without
degrading FV or the physiology output bounds. Base = `params/params_reseeded_regavail_opt2.m`.

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
