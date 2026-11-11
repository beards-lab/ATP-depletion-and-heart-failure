# The slack `ktr` ladder: why the model's profile is flat

2026-09-22. Baseline `params/refit4_opt.m` (objective 0.2909).

## 1. The five slack segments are a length ladder, not five repeats

| seg | SL held during the ktr window | release velocity |
|-----|------------------------------|------------------|
| 1 | 2.041 um | -36 um/s |
| 2 | 2.001 | -45 |
| 3 | 1.960 | -53 |
| 4 | 1.920 | -62 |
| 5 | 1.879 | -70 |

`extractSlackAttributes` fits `ktr` to the 297 ms isometric hold **at the shortened
length**, so the feature is `ktr(SL)` sampled at five lengths. The release duration is
2.2 ms in every segment, so velocity and depth move together.

## 2. The data falls; the model does not

| source | seg1 | seg2 | seg3 | seg4 | seg5 | slope/seg | span |
|---|---|---|---|---|---|---|---|
| 03/27 | 50.66 | 48.91 | 49.14 | 49.15 | 48.20 | -0.47 | -4.9 % |
| 04/03 | 57.47 | 56.11 | 53.07 | 49.43 | 44.06 | -3.35 | -23.3 % |
| 04/10 | 46.70 | 43.49 | 41.41 | 40.39 | 39.61 | -1.73 | -15.2 % |
| **pooled** | 51.61 | 49.50 | 47.87 | 46.32 | 43.96 | **-1.85** | **-14.8 %** |
| refit4_opt | 47.23 | 47.19 | 47.07 | 48.15 | 48.44 | **+0.34** | +2.6 % |

All three preparations fall. The sign is not in question; only the magnitude is
(the three disagree by 5x on how steeply).

## 3. The model already has a 40 % SL gradient on the attachment flux

`dPUdT_CombinedTransitions.m:379` -- `RD1 = ka_eff*PD*N_overlap*f_lattice`:

| factor | seg1 | seg5 | seg5/seg1 |
|---|---|---|---|
| N_overlap | 0.822 | 0.673 | 0.818 |
| f_lattice | 0.972 | 0.707 | 0.727 |
| **product** | 0.799 | 0.475 | **0.595** |

The prefactor falls 40 % across the ladder and `ktr` does not move at all. The
ingredients for a length dependence are present and something cancels them.

## 4. What does move the slope (need d(slope) <= -2.19)

Structural variants, one change at a time:

| change | d(slope) | note |
|---|---|---|
| `P_bound_max` x0.25 | **-3.00** | biggest mechanism lever |
| `kmsrd` x0.5 | -1.31 | SRD -> PD, the surviving force-gated SRX exit |
| `tau_reg` = 20 ms | -1.15 | **conflicts** with the ktr-protocol fast-gate fix |
| `UsePassiveForSR` = 1 | -0.71 | canonical thick-filament mechanosensing, currently OFF |
| `UseOverlapFactor` = 0 | -0.20 | |
| `sigma1` x0.5 or x2 | **0.000** | dead -- see below |
| `useHalfActiveForSR` = 1 | +0.64 | wrong way |

Lean pool at +15 %, most negative:
`PieceWiseStrainDep2X_logOffset` -3.77, `ekSE` -1.71, `dr2` -1.43, `L_thin` -0.74,
`PieceWiseStrainDepR1DX_logOffset` -0.64. **1 of 76 reaches the data slope alone**, and it
takes ktr(1) to 68 (data 51.6), so no single parameter is the answer.

## 5. Two things worth naming

**`P_bound_max` = 1.208.** The physiological binding-site ceiling is ~0.15
(`Analyses/BindingSiteOccupancy`). At 1.208 the gate `f_sat = 1 - P_bound/P_bound_max`
sits at ~0.9 and never bites, so the SL-dependent site count only sets the *target* of
force redevelopment, never its *rate*. Restoring a real ceiling is what converts the
40 % prefactor gradient into a rate gradient -- it is the mechanism, not a fudge.

**`sigma1` is dead code at this operating point.** `RSR2PT = kmsr*exp(F_SR/sigma1)*P_SR`
and `kmsr = 0`, so the direct force-gated SRX exit is switched off entirely; the only
route out of SRX is `ksr2srd` (force-independent) then `kmsrd*exp(F_SR/sigma_srd1)`.
That is why `sigma1` scores exactly zero and `kmsrd` scores -1.31.

## 6. Why the optimiser was free to do this

`ktr` costs **nothing**: pooled s 0.1807 gives fl = 0.0544 while the raw cost is 0.021 --
inside the deadband by 2.6x. Every one of the parameters above could drift wherever the
rest of the objective wanted. Same shape of blind spot as `tau_reg` and the ktr protocol
(`Analyses/KtrProtocolVsSlackKtr`).

## 7. The force block's deadbands are built from the wrong variance

`RunSteadySpread.m`. Each feature's between-preparation spread, split into the part that
is each heart's own force level and the part that is shape:

| feature | s used | s level | s shape |
|---|---|---|---|
| steady | 0.172 | 0.233 | **0.021** |
| A | 0.179 | 0.264 | **0.021** |
| peak1_y | 0.175 | 0.237 | **0.016** |
| vall_y | 0.180 | 0.244 | **0.026** |
| peak2 | 0.175 | 0.211 | **0.024** |
| ktr | 0.181 | 0.104 | **0.041** |

The three hearts' plateau forces are 80.6 / 68.9 / 50.6 kPa; divided by their own means
they agree to 1-2 %. `fl` goes as `s^2`, so these deadbands are ~100x too wide, and
refit4_opt sits inside every one of them while running steady +5.0 %, vall_y -9.7 %,
peak2 +11 % at segment 4, peak1_y -6.5 % at segment 5.

**Fix applied** (`BuildRefit5Seed.m`): for the level-dominated force features
`s := 2 * s_shape` (the factor 2 concedes that a shape spread from n=3 is itself only
known to ~2x), and a new level-free ratio feature `ktr_shape = ktr/ktr(1)` scored on
elements 2:5, alongside the existing `A_over_steady`. The wide `steady[1]` boundary term
is left alone -- that one is the level guard and 52.5-87.7 kPa is the honest heart range.

refit4_opt scores **0.2909 -> 0.7988** on the corrected objective. Written to
`params/refit5_seed.m`. No campaign started.

## 8. Side finding: the staircase overshoots every step

Zooming the held-out staircase to the manoeuvre (it was being plotted with 5.5 s of
isometric warm-up) shows the model peaking 1.15-1.28 on each 1 % step against 1.05-1.12
in the data, 1.51 vs 1.23 on the final stretch, and settling at 1.33 vs 1.19. This is
the **same re-stretch overshoot** as the spikes in the slack trace and the 1.9x
overstretch peak in the ktr protocol -- one defect, now visible on three protocols, two
of them held out. The staircase `E = 0.039` badly understates it.

## 9. CORRECTION + the P_bound_max scan: the ceiling is welded to F-V

Section 5 said `P_bound_max = 1.208` is "8x above the ~0.15 physiological ceiling" and that
the gate "never bites". **Both statements were wrong in an important way, and measuring it
changes the conclusion.**

`P_bound_max` is a **ceiling**, not an operating point. Measured on the slack protocol at
`refit5_seed`, the model's own bound fraction is

    P_bound = p1_0 + p2_0 :   median 0.105,  max 0.327

so the model already sits **inside** the literature's ~0.12-0.15 max-Ca range at the
plateau. The gate is not inert either -- `f_sat = 1 - P_bound/P_bound_max` = **0.91 at the
plateau, 0.73 at the transient peak**, i.e. it attenuates attachment by 9-27 %. What is
true is that its *gradient* is shallow: `d f_sat/d P_bound = -1/P_bound_max = -0.83`, eight
times flatter than it would be at a ceiling of 0.15. That shallowness, not inertness, is
why the 40 % SL gradient never becomes a rate gradient.

Setting the ceiling **at** the literature bound fraction is therefore not "restoring
physiology" -- it puts the operating point on top of the ceiling and shuts attachment off.
The ceiling has to stay above `max(P_bound) = 0.327`.

Scan on the refit5 objective (seed = 0.7988), everything else fixed:

| P_bound_max | total | of which F-V | of which rsK80 | ktr slope |
|---|---|---|---|---|
| 0.12 | 599.5 | 505.6 | 1.33 | -9.15 |
| 0.15 | 476.9 | 418.6 | 1.33 | -8.72 |
| 0.20 | 318.4 | 287.2 | 1.33 | -6.24 |
| 0.30 | 138.0 | 121.6 | 1.33 | -2.73 |
| 0.50 | 27.0 | 9.69 | 13.30 | **-0.68** |
| 0.80 | 1.34 | 0.00 | 0.37 | +0.09 |
| **1.208** | **0.799** | 0.00 | 0.01 | +0.32 |
| data | | | | **-1.85** |

**The ceiling is welded to F-V.** Every value that buys a meaningful slope costs 10-500
units of F-V, and F-V is currently a solved, zero-cost term (it was already known to be a
1749x lever, `Analyses/RefitStrategy`). `P_bound_max = 0.5` is the only interesting point
-- slope -0.68, a third of the way -- and it costs **27 units** against a `ktr_shape` term
worth **0.032**. No optimiser will ever make that trade, so seeding the campaign there
would simply walk it back to 1.2.

**Decision: do not pre-pull it.** `P_bound_max` stays in the lean pool and stays free, but
the ktr slope has to come from the cheap levers instead --
`PieceWiseStrainDep2X_logOffset` (-3.77 at +15 %), `ekSE` (-1.71), `dr2` (-1.43) -- which
is what the refit5 campaign is now pointed at.

## 10. refit5 result: the force block moved, the ktr ladder did not

6.3 h, 16 rounds, 9 snapshots. **0.7988 -> 0.6885** recorded; reloading
`params/refit5_opt.m` from disk gives **0.6950** (+0.95 %, the known serial-vs-parallel
chunking noise; the repeat delta inside one process is exactly 0).

| feature | seed | opt | delta |
|---|---|---|---|
| vall_y | 0.1939 | 0.0847 | **-0.1092** |
| steady | 0.0609 | 0.0341 | **-0.0268** |
| peak1_y | 0.0907 | 0.0727 | -0.0180 |
| ktr_shape | 0.0324 | 0.0255 | -0.0069 |
| rsK80 | 0.0072 | 0.0033 | -0.0039 |
| peak1_dSL | 0.0514 | 0.0505 | -0.0009 |
| ktrP_k63 | 0.0386 | 0.0402 | +0.0016 |
| restretchSlopeStart | 0.0287 | 0.0313 | +0.0026 |
| **peak2** | 0.0967 | 0.1259 | **+0.0291** |
| **vall_t** | 0.1161 | 0.1455 | **+0.0293** |
| TOTAL | 0.7988 | 0.6950 | -0.1037 |
| stairs (held out) | 0.0390 | 0.0384 | -0.0006 |

**What worked.** The re-floored force block is exactly where the gain came from: `vall_y`
-56 %, `steady` -44 %, `peak1_y` -20 %. Sec. 7's diagnosis was right and the fix does what
it was built to do -- these terms were invisible before and the search fixed them as soon
as they could be seen.

**What did not.** The ktr ladder, the entire point of the round:

| source | seg1 | seg2 | seg3 | seg4 | seg5 | slope | span |
|---|---|---|---|---|---|---|---|
| pooled data | 51.61 | 49.50 | 47.87 | 46.32 | 43.96 | -1.848 | -14.8 % |
| refit5_seed | 47.23 | 47.18 | 47.20 | 48.20 | 48.32 | +0.319 | +2.3 % |
| refit5_opt | 47.76 | 47.63 | 47.63 | 48.17 | 48.52 | +0.207 | +1.6 % |

**5 % of the gap closed.** `P_bound_max` did not move at all (1.2082 -> 1.2082), precisely
as sec. 9 predicted: the search will not pay 27 units of F-V for a term worth 0.03.

**And two terms went backwards.** `peak2` +30 % and `vall_t` +25 %, together +0.058 -- more
than half of what `steady` and `peak1_y` gained. That is a genuine tradeoff inside the
re-stretch, not a scoring artefact.

**Held out, ktr protocol** (`validateKtrProtocol('refit5_opt')`): k63 ratio 1.11 / 1.31 /
1.42, mean **1.26** (refit4_opt 1.27); rmse/noise 1.45 / 1.96 / 2.94 against refit4_opt's
1.37 / 1.85 / 2.83 -- **slightly worse on all three days**. The staircase is marginally
better (0.0390 -> 0.0384) and still overshoots every step.

**Read.** `ktr_shape` at 4 % of the objective is too weak to buy a slope that costs real
money elsewhere. The mechanisms exist (sec. 4: three cheap levers at -1.4 to -3.8 per +15 %)
and three of them were compulsory every round, so this is a pricing problem, not a
capability problem. Options, in order of how much they claim:
  * raise `ktr_shape` to `s := 1 x s_shape` (w 0.70 -> 2.78, fl 0.0204 -> 0.0051), worth
    ~0.14 or **21 %** of the objective. Defensible only if we accept the n=3 shape spread
    at face value, and `peak2`/`vall_t` are already regressing, so expect to trade.
  * accept that the ladder slope is out of reach of this parametrisation and report it as
    a known discrepancy -- the model has no length-dependent *rate* mechanism whose price
    the F-V curve will tolerate.
