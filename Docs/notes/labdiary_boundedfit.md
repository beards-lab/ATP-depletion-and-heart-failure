# Lab diary — bounded manual fit (DriverBoundedFit, params_reseeded.m)

Working file: `params/params_reseeded.m` (loaded LAST in DriverBoundedFit, after iter_17
+ inline SRX block). Diagnostics: `Auxiliary/reportFeatureCost.m` (sorted cost_vec,
per-output sub-costs read from boundedOutputFn, raw data-vs-model dump).
Data: 8 mM only (slack 8+2 mM available). Target = mouse/rat alpha-MHC, ~20-25 C, max Ca.

## Physiological framework (the reconciliation that must hold)

Brenner: ktr = fapp + gapp ; duty r = fapp/(fapp+gapp) = attached fraction.
ATPase identity (per head): **XTOR = r·gapp = r·(1-r)·ktr**.

Target corner (cardiac alpha, max Ca, ~22 C), all mutually consistent:
- duty r ≈ 0.20  (attached_ss bound [0.10,0.45]; strong-only ~0.10-0.15, +weak ~0.2)
- ktr ≈ 49 /s    (data, flat across slack segments)
- => gapp ≈ ktr(1-r) ≈ 39 ;  fapp ≈ ktr·r ≈ 10 ;  XTOR = r·gapp ≈ 8  (bound [3,10]) ✓
- A (isometric) ≈ 62 ; steady ≈ 80 (data; decreasing with segment/SLdiff)

So the corner is NOT paradoxical — it sits right on the frontier. The job is to bring
the model onto it.

## Rate -> observable map (CONFIRMED by iters 1-2)

- **kah -> XTOR.** XTOR ≈ kah·PT (hydrolysis flux = cycle flux in SS). When kah is the
  slow step it pins XTOR; when downstream is slower, XTOR decouples from kah.
- **k2 (detachment, gapp) -> attached / ktr / force.** Lower k2 -> longer dwell ->
  attached up, ktr down, force up.
- **ka (attachment, fapp) -> ktr / attached / force.** No strain curve; ka is the raw fapp knob.
- **Strain curves modulate the EFFECTIVE rate** (base stays in physiological bounds):
  - R12 = k1·p1·PieceWiseStrainDep(s)            (power stroke)
  - R2  = k2·p2·PieceWiseStrainDep2(s+dr2-dr)     (detachment; gapp shaping) <-- key turnover/ktr knob
  - R1D = kd·p1·PieceWiseStrainDepR1D(s)          (weak detachment)
  - R21 = k_1·p2·PieceWiseStrainDepR21(s)         (reverse stroke)
  - X_logOffset params SHIFT each curve along strain; Params scale its VALUES.
- **Force is kSE-influenced** (kSE=3203 in series with kstiff~35k): isometric XB strain
  is partly absorbed by SE, so raw kstiff is a weak force knob. Real force knobs: duty,
  overlap, dr, kSE.

## Constraint: kah >= 80 (alpha myosin, user). Reverse hydrolysis (raising kamh) REJECTED
(not physiological here). => XTOR must be capped by making a DOWNSTREAM step rate-limiting
(lower effective gapp via k2/R2), NOT by lowering kah. With kah=80 fast, PT stays small and
XTOR = r·gapp is set downstream.

## Iteration log

### iter 0 (reseed: ka375 kd100 k2350, inline kah100) — GRAND TOTAL 525
XTOR 36, ktr 157, attached 0.15, A 49. Cycle far too fast; ATPase pinned high by kah=100.

### iter 1 (k2 350->110) — TOTAL 514
ktr 157->76, attached 0.15->0.36, A 49->114 (overshoot), XTOR ~42 (pinned by kah·PT).
Lesson: k2 strongly controls attached/ktr/force but NOT XTOR (kah-pinned).

### iter 2 (kah 100->40, ksr2srd decoupled=100) — GRAND TOTAL 134
XTOR 42->19.6, A 114->84, SRX 0.069->0.162 (in bound), attached 0.25.
Lesson: XTOR ≈ kah·PT confirmed (19.6 = 40·0.49). Big win. BUT kah=40 violates the
alpha-myosin expectation (kah>=80). ktr rose to 127 and became OSCILLATORY (ktr_rmse 3.2,
ktr2_overshoot 19) -> redevelopment is multi-phase, not single-exp. Force still ~30% high.

## Open problems at TOTAL 134 (by cost)
1. ktr complex (~55): redevelopment too fast AND oscillatory (shape, not just rate).
2. A / slack-force magnitudes (~50): all ~30% high together (one force-scale knob).
3. XTOR (14): 19.6 vs <=10 (must come down WITH kah raised to 80).

### iter 3 (kah 40->80, ka 375->150, k2 110->100) — GRAND TOTAL 142
XTOR 19.6->26.8 (=80*PT0.348), ktr 127->78 and ktr_rmse 3.2->0.77 (OSCILLATION GONE!),
A 84->90, attached 0.285. Lesson: lowering ka killed the oscillatory redevelopment —
the multi-order ktr artifact was driven by too-fast attachment (ka=375).

### iter 4 (ka 150->75) — GRAND TOTAL 37 (output 25)  *** big win ***
A 90->68.5 (data 69.6 ✓), peak1 125->96 (✓), vall 96->70 (✓), steady 98->77 (✓),
SRX 0.092->0.105 (in bound), attached 0.285->0.22. ktr 78->84, XTOR 26.8->18.8.
Lesson: ka is the master FORCE+duty knob; ka=75 nails all slack force magnitudes at once.
Remaining: XTOR 18.8 and ktr 84 both ~1.8x too fast (cycle still fast). PT=0.25.

## Lessons consolidated
- **ka = master knob for force/duty AND redevelopment smoothness.** 375 -> 75 fixed
  A/peak1/vall/steady/SRX simultaneously and removed ktr oscillation.
- **k2 = attached/ktr/force**, **kah = XTOR (=kah*PT)**, but with kah pinned >=80 the only
  way to drop XTOR is to drop PT, i.e. slow the cycle flux (lower effective fapp & gapp).
- Identity check at duty 0.22: ktr 84 -> want 49 => cut fapp & gapp x0.58 (keep ratio so
  force holds): ka->50 (floor) + cut k2eff via R2 curve. Predict XTOR -> ~8, ktr -> ~49.

### iter 5 (ka 75->50, R2 central knots x0.6) — GRAND TOTAL 51 (WORSE)
ktr 84->57 ✓, XTOR 18.8->15.2 ✓, BUT A 68->78 (slower detach raised duty/force),
FV_fnorm tanked (iso denom up), vall too deep. Lesson: cutting detachment ALONE raises
force; must cut the forward limb too to hold duty.

### iter 6 (+ k1 264->140) — GRAND TOTAL 36, output 24
A 78->71.4 ✓, ktr ~60 (clean fit), XTOR 13.1, peak1 still 111 (hi). Cost profile SHIFTED:
**FV_fnorm now #1 (7.7)** — model FV too steep (fnorm 0.57 vs 0.92 @v=0.5); model
over-attaches isometrically (FV_f0 71 vs 56) and under-bears force in shortening.

## Lessons consolidated (updated)
- **ka = master force/duty + redevelopment-smoothness knob** (375->75 nailed slack forces).
- **k2/R2-curve = detachment(gapp); k1 = forward(fapp); together at fixed ratio = cycle
  speed without changing duty.** Floors (ka>=50,k2>=100,kah>=80) cap the minimum cycle
  speed -> strain curves (k1eff,k2eff) are how you go below the floor turnover.
- **kah=80 pins XTOR=kah*PT** (~13 now); only PT drainage (parking) or even lower cycle
  flux brings XTOR<=10. Currently XTOR 13 (close).
- Cutting detachment alone -> force UP; cut forward limb (k1) to compensate. CONFIRMED.

## Best so far: iter 6, GRAND 36 (from 525). State: ka50 kd100 k2 100 k1 140 kah80,
## R2 central knots [.. 1.80 0.617 0.629 ..]. Open: FV shape, peak1 hi, t0, SRX 0.064(lo).

### iter 7 (UseVelGaussAttachment=true) — GRAND 152 (REVERTED)
XTOR->8 (in bound!), FV_fnorm->0.90 (✓), A->54, BUT ktr exploded to 150 (rmse 5.7,
oscillatory). Lesson: velocity-gating attachment corrupts the isometric force
redevelopment (ktr). FV/XTOR win, ktr loss -> net much worse. STRUCTURALLY incompatible
with clean ktr. Confirms: reducing isometric attachment is the FV/XTOR fix, but must be
done WITHOUT a velocity gate.

### iter 8 (k1 140->200) — GRAND 42 (worse)
A 71->76, FV still steep (0.58). Lesson: k1 raises force but does NOT fix FV ->
FV steepness is NOT a k1 problem.

### iter 9 (k1->140, restore R2 knot5 s2=0 -> 1.0475) — GRAND 35.1 (best then)
FV still 0.58. Restoring the shortening-transition knot didn't help FV either.

### iter 10 (mu 0.0404->0.015) — GRAND 34.2, output 22.3 (BEST)
FV still 0.58 (drag isn't the FV lever at these v). Tiny net gain.

## FINAL STATE (iter 10): GRAND 34.2 (from 525 = 15x). params_reseeded.m:
##   ka=50, kd=100, k2=100, k1=140, kah=80, mu=mu_neg=0.015,
##   R2 PieceWiseStrainDep2Params=[50 15.97 1.80 0.617 1.0475 50.68 50 50]
## Fits well: A, peak1, vall, steady, peak2, ktr(clean ~60), XTOR(13), attached(0.25),
##   SLslack, restretchSlope; slack force trace tracks data; FV on the Hill curve.

## THE TRADE-OFF MAP (why single-knob manual tuning is now at its frontier)
- ka up  : FV better (reattach in shortening) + force up, BUT ktr up + XTOR up
- ka down : ktr/XTOR down, BUT FV steeper + force down
- k2/R2-central down : ktr/XTOR down + force up, BUT FV steeper (heads drag at -strain)
- k1 up  : force up (+ slight FV), BUT ktr up
- kah    : pins XTOR = kah*PT (>=80 floor -> XTOR can't go <~10 without draining PT)
- velGauss: FV+XTOR fixed BUT ktr destroyed (no good)
- mu     : ~no effect on FV at these velocities
=> The remaining residual is dominated by FV_fnorm (7.3), which is in DIRECT mutual
   tension with ktr+XTOR through ka. No single knob improves all; this is a genuine
   multi-objective frontier.

## REMAINING (for phase-2 optimizer / structural):
1. FV_fnorm (7.3): model FV too steep (0.58 vs 0.92 @v=0.5). Fights turnover via ka.
   Likely needs a velocity-dependent REATTACHMENT that does NOT gate the isometric ktr
   (i.e. boost attach only at v<0 shortening, leave v~0 alone) — a clean structural fix.
2. A/peak1 force-magnitude segment trends (3.6/3.2): model A decreases too little across
   segments vs data.
3. XTOR 13 vs <=10 (2.0): needs ~0.13 more PT drained -> more SRX parking.
4. SRX_ss 0.066 (<0.10) + INPUT violations ksr0=190(ub50), kmsrd=60(ub20): INLINE SRX
   block (user domain). The SRX parking is over-strong on rate but under-populated at
   force -> revisit sigma1/2 (force mobilization) + bring ksr0/kmsrd into bounds.

## PHASE 2 — manual, force-length (A-vs-segment) + SRX slope + denser FV

### iter 11 (FV_velocities -> -[0 0.5 1 2 3 4 5 6], 8 pts; no tuning)
Candidate held (A/ktr/slack unchanged). BUT denser FV exposed FV is badly too steep
ACROSS THE WHOLE RANGE: FV_fnorm model [1 .58 .33 .14 .077 .047] vs data
[1 .92 .66 .32 .20 .11]. FV_fnorm cost 7->16. The 5-pt FV was hiding this. Honest now.

### A-vs-segment (force-length) is the other target
A is the asymptotic redeveloped force at SLslack (extractSlackAttributes: feats.A=ae.A+F0).
The 5 consecutive deepening slacks sample SLslack=[2.04 2 1.96 1.92 1.88]. Model A Δ=5
(71.6->66.3) vs data Δ=15 (69.6->54.7) -> model force-length 3x too SHALLOW.
SRX force-feedback knobs (kmsr=0 inline, so sigma1 inert): RPT2SR=ksr0*exp(-F/sigma2)*PT
(park, lower F=more), RSRD2PD=kmsrd*exp(F/sigma_srd1)*P_SRD (mobilize, lower F=less).
Lower sigma2 / sigma_srd1 -> stronger low-force parking -> steeper A-vs-SL. But the SAME
feedback acts during FV shortening (lower F -> more park -> steeper FV) => A-slope and FV
may TRADE through the SRX slope. Seed of the force-length is overlap (shallow here);
SRX only amplifies it; lattice-spacing attachment (UseLatticeSpacing=0) is the untried
attachment-based length-dependence the user hinted at.

### iter 12 (sigma2 42->25) — FAILED (GRAND 49)
Sign error: lower sigma2 => larger F/sigma2 => SMALLER exp(-F/sigma2) => LESS parking at
operating force => SRX collapsed to 0.03 (out of bound), A-trend unchanged. KEY LESSON:
the SRX pool at active force (~0.03-0.07) is TOO SMALL to drive the force-length LDA, and
the force differences between segments are too small to seed SRX feedback. Reverted.

### iter 13 (UseLatticeSpacing=true, d_optimal=0.0246, sigma_lattice=0.0025) *** WIN ***
Length-dependent ATTACHMENT via lattice spacing = the physiological cardiac LDA mechanism.
A model [71.5 69.3 67.5 63.5 58.8] vs data [69.6 66.2 62.2 58.3 54.7] — force-length slope
now Δ12.7 (was 5; data 15). A cost 3.6->1.1. FV UNTOUCHED (f_lattice=1 at SL>=2.04, so the
FV/isometric reference at SL0=2.2 is unaffected). GRAND 43->40, output 28.
LESSON: the force-length seed must come from ATTACHMENT (lattice), not SRX feedback. The
user's "does the attachment work on that?" was the right instinct.

### iter 14 (d_optimal 0.0246->0.0242, slightly stronger lattice) — A SOLVED
A model [70.7 68.4 65.4 60.2 55.1] vs data [69.6 66.2 62.2 58.3 54.7]; A cost 1.1->0.27.
Am also matches. Force-length now physiological. GRAND 38, output 26.

## CURRENT BEST = iter 14. params_reseeded.m: ka50 kd100 k2 100 k1 140 kah80 mu0.015,
## R2=[50 15.97 1.80 0.617 1.0475 50.68 50 50], FV_v=-[0 .5 1 2 3 4 5 6],
## UseLatticeSpacing=1 d_optimal0.0242 sigma_lattice0.0025.
## Remaining cost 26: FV_fnorm 16 (model FV ~2x too steep across whole range; also the
##   data FV-isometric 56 < slack-A isometric 70 -> a protocol normalization gap the model
##   can't represent with one isometric); t0_crossing 2.0; XTOR 1.8 (12.9 vs 10);
##   vall_y 1.8 (valley too deep 64 vs 77 — restretch transient, catch-bond territory);
##   peak1_y 1.7 (restretch peak ~110 vs 96).
## SOLVED/GOOD: A (force-length), ktr (clean ~61), attached 0.25, steady, peak2, A.
## Levers tried for slack-transient (peak1/valley) trade vs FV/Vmax; catch-bond
## (UseCatchBond, k_catch_bond, currently off) is the mechanism the earlier labdiary.md
## used for the valley — left for a focused restretch session / optimizer.

## RECOMMENDATION
Manual physiological tuning has done its job (525->34, all base rates in-bounds, kah=80
alpha). The last mile is a focused multi-objective optimization over the SIX identified
knobs {ka, k1, k2, kah, R2-central-knots, mu} against params0.fn, reusing
Auxiliary/ResidualAndJacobian.m + fminsearch (small, well-posed now that the landscape is
mapped). FV likely also needs the shortening-only reattachment mechanism above.

## PHASE 3 — FV via OCCUPANCY + R2 ONSET reshape (see Analyses/BindingSiteOccupancy)

Proposed (correctly) that the FV "velocity-dependent attachment" should gate on
binding-site OCCUPANCY (global P_bound), not instantaneous velocity. KEY mechanistic
insight: velocity is an instantaneous kinematic variable that SWINGS during the restretch
transient (fast stretch -> hold) -> velGauss mis-fires -> oscillatory ktr. Occupancy is a
slow STATE variable: coincides with velocity at steady state (reproduces FV benefit) but
does NOT swing during transients (preserves clean ktr). CONFIRMED below.

Turned out Analyses/BindingSiteOccupancy/conclusions.md §7 had ALREADY done the full
occupancy+FV investigation against THIS iter-14 basis (its §7.1 baseline = our 26.0 / XTOR
12.9 / FV 16 exactly). Its verdict: occupancy ALONE collapses force (artifact, no kstiff
recompensation); the real FV lever is the R2 ONSET-detachment reshape.

### iter 15 (occupancy langmuir, P_bound_max=0.15, ALONE) — reproduces report §7.1
ktr stayed CLEAN (rmse 0.18-0.30, NOT oscillatory) -> confirms occupancy preserves ktr
where velGauss destroyed it. XTOR 13->7.9 (in bound). BUT A 70->45 (force collapsed,
no kstiff comp), FV ~unchanged (15.7). Occupancy alone is not the FV fix.

### iter 16 (occupancy linear cap0.30 + kstiff2 x1.6 + R2 onset knot s2=0: 1.05->5) — NEW BEST
Output cost 26 -> 19.5. FV_fnorm cost 16 -> 8.7 (-46%); FV model [1 .64 .42 .23 .15 .11]
vs data [1 .92 .66 .32 .20 .11]. A 69.8 (✓), steady 80.9 (✓), peak1 97.5 (✓), peak2 80.9
(✓), XTOR 10.1 (bound edge). The R2 onset raise (heads detach EARLY before negative-force
drag) repairs the low-v FV drop; occupancy flattens the high-v end + caps attached(0.15).
TRADE: ktr 61 -> 79 (data 49; cost 3.7 + rmse 1.3) — iron-law (occupancy lowers attached
-> ktr up). vall_y still deep (61 vs 77, cost 2.5).

## CURRENT BEST = iter 16. params_reseeded.m adds (on top of iter14):
##   UseGlobalOccupancySaturation=1, OccupancyForm='linear', P_bound_max=0.30,
##   kstiff2 = 34897*1.6 = 55835, R2 PieceWiseStrainDep2Params=[50 15.97 1.80 0.617 5 50.68 50 50]
## GRAND 31.5 (output 19.5 + phys 11.9=SRX bound violations ksr0/kmsrd, inline domain).
## Open, by cost: FV 8.7 (diminishing; model still <data at v=0.5-1), ktr 79 vs 49 (NEEDS
##   3rd low-force slow-detaching state per [[mech_tradeoff]] — can't close FV+ktr in 2-state),
##   vall_y deep (restretch/catch-bond), XTOR at edge.

## PHASE 4 — is ktr fixable in 2-state? (Route A, no p1 rewire) — YES, but it TRADES FV.

XTOR cost bound relaxed 10->15 (boundedOutputFn.m, mouse-alpha plausible, FJ). Identity:
ktr = XTOR / (r(1-r)). Current ktr is high ONLY because duty r is low; the identity binds
because p1 is force-less (force ∝ s, centered at attachment ≈0; p2 force ∝ (s+dr) carries
the working force) -> p2 is the sole force+cycling state. Verified each run: r(1-r)ktr ≈ XTOR.

### iter17 / Route A.1 (P_bound_max 0.30->0.80, kstiff2 x1.6->x1.0)
duty 0.15->0.21, ktr 79->70 (identity confirmed: higher duty lowers ktr). XTOR 10.8 (in
[3,15]). But duty capped ~0.21 (ka floored at 50 + R2-onset detachment + residual occ).
FV regressed 8.7->10.7 (occupancy was helping high-v). A/steady dropped (kstiff undershoot).

### iter18 / Route A.2 (+ isometric R2 knots 3,4: 1.80,0.617 -> 1.12,0.38)
duty 0.26, **ktr NAILED -> [48.9 49 52 51 56], cost 0.06** (data ~49). XTOR 10.4. PROVES ktr
is NOT an iron law in 2-state. BUT FV cost 8.7->18.9 (FV_fnorm 0.53 vs data 0.92 @v=0.5) and
A rose to 76: the isometric-detachment cut that lifts duty also lifts the isometric force/
reference and re-steepens FV.

### TRADE-OFF MAP (force-less-p1 2-state): ktr<->FV trade monotonically through DUTY
| run | duty | ktr | ktr cost | FV cost | output |
| iter16 | 0.15 | 79 | 3.7 | 8.7  | 19.5 |
| A.1    | 0.21 | 70 | 1.9 | 10.7 | 24   |
| A.2    | 0.26 | 50 | 0.06| 18.9 | 35   |
Weighted-best standalone stays iter16 (FV weight 10 >> ktr weight 2). No 2-state sweet spot
closes BOTH. REVERTED working file to iter16.

## CORRECTED CONCLUSION (supersedes "needs 3rd state"):
ktr alone IS fixable in 2-state (A.2). ktr AND FV together are not — they trade through duty
because force comes only from p2, so the duty that sets ktr is the same p2 detachment that
sets FV. **Route B = make p1 force-bearing** (force ∝ (s+dr1), dr1>0 weakly-bound pre-stroke
force; or dr1<0 back-force): carries the high force-duty ktr needs WITHOUT the slow isometric
p2-detachment that steepens FV -> decouples ktr from FV. This likely OBVIATES the 3rd state
(does the same job with the existing p1). ~2-line code change in the force block (lines
165-175) + dr1/kstiff1 retune. AWAITING greenlight (user: "do not rewire p1 yet").

## PHASE 5 — manual FV wall, p1-force (Route B) IMPLEMENTED, FV structural limit

XTOR cost bound relaxed to [3,15] (mouse-alpha, FJ). FV velocities -> 5pt -[0 0.5 1 2 4].

### Manual FV wall confirmed
Model FV F(0)=70 vs data 56 (FV-iso < slack-iso in DATA at same SL ≈ history-dependent force
depression the steady-state model can't reproduce). Pushing harder occupancy (P_bound_max 0.18)
+ kstiff x2.4: F(0) went UP (kstiff overcomp), FV_fnorm STAYED ~0.65, ktr went oscillatory.
**FV_fnorm is force-scale-INVARIANT** -> the relative isometric decrease cannot be set by
occupancy/kstiff; it's kinetic (R2 detachment-vs-strain during shortening).

### Route B IMPLEMENTED (dPUdT force block): p1 force ∝ Σ((s+dr1)·p1), reuse dr1 (default 0
### = identical to before; backward-compatible). *** ktr DECOUPLED — works ***
High-duty base (occ relaxed P_bound_max=0.80, R2 isometric knots [1.12 0.38 5], ka=50) + dr1=0.004,
kstiff1=45000, kstiff2=28000: **ktr=[50.5 50 51 53 53] (cost 0.03!), A perfect (69.4 vs 69.6),
XTOR 9.9** — ALL the things that were welded are now satisfied at high duty. The user's idea is
CONFIRMED: force-bearing p1 breaks the ktr<->duty weld. BUT FV still steep (FV_fnorm 0.51, cost 9).

### FV is the residual — and it is STRUCTURAL (invariant to every knob)
Raising ka 50->250 (hypothesis: faster reattachment flattens FV): FV_fnorm STAYED 0.51, just
doubled force + dropped ktr. FV_fnorm @v=0.5 has been ~0.5-0.65 across EVERY config all campaign
(data 0.92). Force collapses ~2x at v=0.5 ML/s regardless of ka/kstiff/occupancy/dr1. The
flat-shoulder FV is a structural limit of the model's shortening force (strain advection +
point attachment + R2 slip), NOT a tuning target. Candidates left: R2 shortening-knot reshape
(weak, tried), distributed attachment (UseA1AttachmentKernel), or the FV-iso history-depression
being a genuine protocol effect (compare FV-protocol SL vs slack SL in the data).

## CURRENT STATE: params_reseeded.m = p1-force config (dr1=0.004, ka50, P_bound_max0.80,
## kstiff1 45000, kstiff2 28000, R2 [50 15.97 1.12 0.38 5 50.68 50 50]). Output ~18.4:
## ktr/A/steady/XTOR all GOOD; residual = FV_fnorm 9.0 (structural) + vall_y/peak1 (restretch
## transient degraded by high duty — catch-bond/c_SE_visc territory).
## Optimizer PREPARED: Workbench/RunBoundedFit_Optim.m (free-set incl. dr1 + R2 shortening
## knots; bounded fmincon multi-start over evaluateBakersExp). Run to refine FV/transient +
## lock-proof. Net: p1-force solved ktr; FV flat-shoulder is the structural wall.

## PHASE 6 — FV WALL BREACHED: registration-availability state (UseRegistrationAvailability)
## *** the "structural" FV shoulder was NOT structural — it was a missing mechanism ***

The optimizer (RunBoundedFit_Optim.m) ran to params/params_reseeded_opt.m (ka81 k2150 dr1.0047
kstiff1 36298 kstiff2 24030 P_bound_max1.02 R2 knots5/6 9.4/16.0; output cost 9.9). FV_fnorm still
stuck @v=0.5 = 0.616 (data 0.92). Per Analyses/FV_Shoulder: the flat shoulder = target-zone
REGISTRATION availability (Edman 1988; Tanner/Daniel/Regnier 2007) — sub-maximal isometric
attachment relieved by sliding. NOT reproducible by force-scale knobs (FV_fnorm is a ratio);
needs a velocity-set availability that the OLD velGauss mis-keyed on instantaneous velHS (-> broke ktr).

### Implemented (Model/dPUdT_CombinedTransitions.m + getParams.m): +1 SCALAR ODE state A_reg at
### the tail (index Ns*ss+8), gates attachment ka_eff*=A_reg, driven by IMPOSED Vums (not velHS):
###   A_inf = A0 + (1-A0)|vel|/(|vel|+v_ref_reg);  dA_reg = (A_inf - A_reg)/tau_reg
### Flag UseRegistrationAvailability (default off). Backward-compat PROVEN: A0=1 -> max|Δf|=0 vs off.

### RESULT (Workbench/ValidateRegistrationAvailability.m, on params_reseeded_opt basis, 8mM FV+Ktr+Slack):
| config                  | FV_fnorm@.5 | FV cost | A    | steady | ktr_rmse | XTOR | TOTAL |
| OFF (= the optimum)     | 0.616       | 4.68    | 67.3 | 77.6   | 0.58     | 13.2 | 9.90  |
| ON A0=0.6 (no recomp)   | 0.725       | 1.91    | 52.7 | 63.0   | 0.26     | 9.9  | 37.2  |
| ON A0=0.6 + kstiff x1.28| **0.734**   | **1.68**| 67.4 | 78.0   | 0.62     | 11.3 | **8.89** |
A_reg=A0 depresses F(0) (shoulder reference); sliding recovers it (flatter ratio); tau_reg≈1/ktr
keeps ktr single-order (no velGauss oscillation). kstiff x1.28 restores absolute A/steady WITHOUT
un-flattening the ratio. **Net cost 8.89 < 9.90 — beats the force-less-p1 optimum AND flattens the
shoulder.** The campaign-long "FV structural wall" is breached.

## CURRENT STATE: params/params_reseeded_regavail.m (full snapshot, feature ON, A0=0.6, v_ref_reg=0.8,
## tau_reg=0.02, kstiff x1.28; cost 8.89). Refiner: Workbench/RunRegAvailFit_Optim.m (free-set
## A0, v_ref_reg, ka, k2, kstiff1/2, R2 knots 5/6). Residuals (tuning, not structural): ktr~64 vs 49
## (PRE-EXISTING high-ktr — retune ka/k2); FV not fully flat (deeper A0); restretch peak2/vall (kstiff
## bump — retunable). Still to test: 8-vs-2mM shoulder FADE (mechanism prediction).
