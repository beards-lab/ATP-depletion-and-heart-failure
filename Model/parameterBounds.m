function bounds = parameterBounds()
%PARAMETERBOUNDS  Physiological bounds for numerical parameters in
%                 dPUdT_CombinedTransitions and related model code.
%
%   B = PARAMETERBOUNDS() returns a struct keyed by parameter name with
%   fields:
%       B.(name).lb      lower bound (in the param's native units)
%       B.(name).ub      upper bound
%       B.(name).weight  penalty weight (0 = no penalty applied)
%
%   Bounds are applied UNIVERSALLY by evalPhysiologyCost — to every
%   parameter listed here, whether or not it is in `params0.mods`. A wrong
%   frozen value will still register on the dashboard (plotPhysiologyDashboard).
%
%   Citations live as code comments next to each entry; they are NOT
%   propagated as struct fields (keeps runtime struct lean).
%
%   Tier organisation (sectional, not algorithmic):
%       Tier A (weight 100): well-anchored from direct measurement
%       Tier B (weight 10):  order-of-magnitude from cardiac literature
%       Tier C (weight 1):   defensible but uncertain
%       Tier D (weight 0):   no prior available; listed for coverage only
%
%   Note on SRX nomenclature: ksr0, kmsr, sigma1, sigma2 implement the IHM
%   thick-filament mechanosensing ON/OFF transition (Linari 2015 Nature,
%   Brunello 2020 PNAS) — NOT the biochemical SRX of McNamara 2015
%   (~0.004 s^-1). Bounds reflect the mechanosensing interpretation.
%
%   See also: evalPhysiologyCost, plotPhysiologyDashboard, parameter-reference.md

bounds = struct();

%% ============================================================
%% TIER A — penalty weight 100 (well-anchored direct measurements)
%% ============================================================

% Power stroke distance dr [um]. Single-molecule optical-trap measurements on myosin II.
%   Capitanio 2006 https://doi.org/10.1073/pnas.0506830103        PMID:16371472
%     Rabbit skeletal myosin II; two mechanical substeps summing to ~8 nm, supporting
%     lb = 8 nm. Skeletal myosin II is mechanistically analogous to cardiac beta-MHC.
%   Finer 1994     https://doi.org/10.1038/368113a0              PMID:8079506
%     Skeletal myosin II in vitro; single working stroke ~11 nm, supporting ub = 11 nm.
%     (Finer, Simmons & Spudich, Nature 1994 — the original single-molecule myosin II paper.)
%   CAUTION: Guilford 1997 (smooth muscle myosin, different isoform — kept for reference
%     range overlap but not cardiac). REMOVED: Mehta 1999 was myosin V (36 nm steps,
%     entirely wrong motor protein).
%   Guilford 1997  https://doi.org/10.1016/S0006-3495(97)78753-8  PMID:9138566
bounds.dr.lb = 0.008; bounds.dr.ub = 0.011; bounds.dr.weight = 100;

% Cardiac titin passive-force power-law exponent gamma [dimensionless].
% Active only when UseTitinIdentifiedPassive=1; otherwise the "gamma" in
% the model is a placeholder with a non-physiological value (~2.67).
%   Granzier & Irving 1995 https://doi.org/10.1016/S0006-3495(95)80278-X PMID:7756567
%     Rat cardiac trabeculae passive force decomposition; reports exponential form
%     F~exp(coeff*SL). Equivalent power-law exponent ~7-9 in the physiological SL
%     range (paywall; value widely cited in cardiac mechanics secondary sources).
%   Fukuda 2003    https://doi.org/10.1161/01.RES.0000054065.16294.A1 PMID:12970176
%     Rat cardiac sarcomeres; directly measures passive-force vs extension and fits
%     to power-law form consistent with gamma ~7-9 for N2B cardiac titin isoform.
%   Wu 2000        https://doi.org/10.1161/01.RES.86.6.682        PMID:10747009
%     Desmin-knockout and wild-type mice; includes rat LEFT VENTRICULAR measurements
%     of passive titin stiffness, fitting to gamma ~8-10. Previously mis-classified
%     as skeletal-only; the paper explicitly includes cardiac fibers (Circ Res 86:682).
%   Terui 2008     https://doi.org/10.1085/jgp.200809864          PMID:18955594
%     Skinned porcine ventricular muscle; troponin and titin coordinate length-dependent
%     activation; includes titin passive stiffness in cardiac consistent with gamma ~7-9.
%     (PMID:18316560 listed previously was a troponin-T N-terminal kinetics paper —
%     wrong; PMID:18955594 is the correct Terui et al. J Gen Physiol 2008 reference.)
bounds.gamma.lb = 7; bounds.gamma.ub = 10; bounds.gamma.weight = 100;

% Actin half-helical repeat d_actin [um] (5.46 nm structural parameter).
%   Holmes 1990    https://doi.org/10.1038/347044a0               PMID:2395461
%     Landmark atomic model of F-actin; establishes the 5.46 nm axial repeat directly
%     from fibre diffraction. Structural geometry, not a fitted parameter.
%   Egelman & DeRosier 1992 https://doi.org/10.1016/0304-3991(92)90025-G
%     Image analysis of F-actin electron micrographs; independently confirms the
%     half-repeat. Bound [5.4, 7.4 nm] spans from half-repeat to full monomer spacing.
bounds.d_actin.lb = 0.0054; bounds.d_actin.ub = 0.0074; bounds.d_actin.weight = 100;

% Initial sarcomere length SL0 [um] (experimental setpoint, physiological range).
% Below 1.8 um double filament overlap occurs; above 2.4 um passive force dominates
% and physiological activation is reduced.
%   Gordon 1966    https://doi.org/10.1113/jphysiol.1966.sp007909 PMID:5966030
%     Defined the force-length relationship and physiological SL operating range for
%     striated muscle (frog skeletal); canonical reference for these bounds.
bounds.SL0.lb = 1.8; bounds.SL0.ub = 2.4; bounds.SL0.weight = 100;

% Bare zone half-length L_hbare [um] (thick filament structural anatomy).
%   Tonino 2017    https://doi.org/10.1038/s41467-017-01144-9     PMID:29062075
%     Super-resolution microscopy of human cardiac thick filaments; bare zone
%     ~200 nm total (half = 100 nm). Bound [80, 120 nm] is the reported value
%     ±20% to span mouse (~80 nm) to human (~100 nm) species variation.
bounds.L_hbare.lb = 0.08; bounds.L_hbare.ub = 0.12; bounds.L_hbare.weight = 100;

% Thick filament length L_thick [um] (rat/human cardiac).
%   Tonino 2017    https://doi.org/10.1038/s41467-017-01144-9     PMID:29062075
%     Super-resolution microscopy; human cardiac thick filament ~1.6 um. Bound
%     [1.55, 1.70] spans mouse (~1.58 um) to human, directly confirmed.
%   REMOVED: Robinson & Winegrad 1979 — that paper measures THIN filament lengths
%     in rat heart (see L_thin below); was incorrectly listed here.
bounds.L_thick.lb = 1.55; bounds.L_thick.ub = 1.7; bounds.L_thick.weight = 100;

% Thin filament length L_thin [um] (cardiac, Z-disk to pointed end).
%   Tonino 2017    https://doi.org/10.1038/s41467-017-01144-9     PMID:29062075
%     Super-resolution microscopy; human cardiac thin filament ~1.0-1.1 um.
%     Upper bound 1.25 um allows species/SL variation; the current model default
%     L_thin=1.377 exceeds this bound (flagged in parameter-reference.md).
%   Robinson & Winegrad 1979 https://doi.org/10.1085/jgp.74.1.117 PMID:469766
%     Direct electron microscopy of rat cardiac thin filaments; ~1.0 um length.
%     (This is correctly a thin-filament paper — previously mislisted under L_thick.)
%   REMOVED: Burgoyne 2008 — MyBP-C/myosin-S2 interaction paper; no L_thin data.
bounds.L_thin.lb = 1.0; bounds.L_thin.ub = 1.25; bounds.L_thin.weight = 100;

%% ============================================================
%% TIER B — penalty weight 10 (order-of-magnitude from cardiac lit)
%% ============================================================
%
% NOTE: bounds apply to BASE rates. Under the convention that PCHIP central
% value (s ~ 0) ~= 1, the base rate equals the rate at zero strain — which
% is what the literature reports.

% Attachment rate ka [s^-1] (PD -> P1, weakly-bound -> strongly-bound, cardiac).
%   [NO-CITE — direct ka measurement] Little 2012 (below) constrains only kd.
%   Ka bounds [50, 500] are order-of-magnitude from ensemble measurements of ktr
%   and steady-state attachment fraction. No direct single-molecule cardiac ka
%   measurement is currently cited here; a direct reference is needed.
bounds.ka.lb = 50; bounds.ka.ub = 500; bounds.ka.weight = 10;

% Detachment rate kd [s^-1] (P1 -> PD, base, cardiac).
%   Little 2012 https://doi.org/10.1074/jbc.M112.380048           PMID:22685295
%     [OK — cardiac] Measures cross-bridge DETACHMENT rates from ventricular
%     myofibrils via fluorescent troponin C: human ~40 s^-1, rat ~150 s^-1;
%     reduces ~3× in the presence of ADP (~13-50 s^-1). Directly constrains kd.
%     Rat value (150) is at the upper end of bound [5, 200]; human (40) is midrange.
bounds.kd.lb = 5; bounds.kd.ub = 200; bounds.kd.weight = 10;

% Power stroke forward k1 [s^-1] (P1 -> P2, isomerisation after attachment).
%   Land 2017 https://doi.org/10.1016/j.yjmcc.2017.03.003         PMID:28392437
%     [MODEL-PARAM] Computational cardiac mechanics model; k1 is a fitted parameter
%     calibrated to reproduce ktr and force-velocity data, not directly measured
%     biochemically. Bounds [100, 2000] are plausible for a microscopic isomerisation
%     rate but no primary biochemical measurement is cited.
bounds.k1.lb = 100; bounds.k1.ub = 2000; bounds.k1.weight = 10;

% Power stroke reverse k_1 [s^-1] (P2 -> P1).
%   [NO-CITE] Constrained relative to k1 by detailed balance; direct measurement of
%   reverse power stroke rate is rare. Bound [0, 200] allows near-irreversibility.
bounds.k_1.lb = 0; bounds.k_1.ub = 200; bounds.k_1.weight = 10;

% ADP-release rate k2 [s^-1] (post-stroke detachment, P2 -> PD).
%   Malmqvist 2004    https://doi.org/10.1074/jbc.M310974200      PMID:14976210
%     [OK — cardiac] Reports human cardiac myosin-II ADP-release kinetics; together
%     with Siemankowski 1985, mechanism-evaluation.md confirms ~300-400 s^-1 at 37°C
%     for cardiac myosin II. Both values sit within bound [50, 1000].
%   Siemankowski 1985 https://doi.org/10.1073/pnas.82.3.658       PMID:3156150
%     [OK — cardiac included] "ADP dissociation from actomyosin subfragment 1
%     is sufficiently slow to limit unloaded shortening velocity in vertebrate
%     muscle" — includes cardiac actomyosin; reports ~300-400 s^-1 at 37°C.
%     Lower bound lb=50 corresponds to ~12-15°C (Q10~2.5 from McDonald/Difee
%     conditions; 300/2.5^2.2 ≈ 50 s^-1).
bounds.k2.lb = 50; bounds.k2.ub = 1000; bounds.k2.weight = 10;

% k2 reverse (P2 re-attachment of ADP-Pi) — no direct measurement.
bounds.k_2.lb = 0; bounds.k_2.ub = 50; bounds.k_2.weight = 10;

% ATP-binding rate k3 [s^-1] (PT -> PD, overall first-order at experimental [ATP]).
% In the model, k3 collapses ATP binding AND possibly Pi release into one step.
% Mechanism-evaluation.md notes: "k3 likely also captures Pi release, which is
% much faster (Dantzig et al. 1992: Pi release ~200-400 s^-1 in rabbit psoas)."
% The rate-limiting step at saturating [ATP] is the ATPase turnover (~5-30 s^-1
% for cardiac myosin), but as a microscopic step between states, k3=44 is
% consistent with the simplified scheme.
%   Lan & Sun 2005 https://doi.org/10.1529/biophysj.105.062463    PMID:15749785
%     [MODEL-PARAM — skeletal] Computational model; k3 is a fitted parameter.
%     Provides order-of-magnitude guidance but is not a direct cardiac measurement.
%   Dantzig 1992 https://doi.org/10.1113/jphysiol.1992.sp019251   PMID:1460122
%     Pi release ~200-400 s^-1 in rabbit psoas at 10°C — gives an upper bound on
%     what k3 represents if it captures Pi release as the rate-limiting step.
bounds.k3.lb = 20; bounds.k3.ub = 200; bounds.k3.weight = 10;

% k3 reverse (k3m). No direct cardiac measurement; bound is an order-of-magnitude prior.
bounds.k3m.lb = 0; bounds.k3m.ub = 20; bounds.k3m.weight = 10;

% ATP hydrolysis rate kah [s^-1] (M·ATP -> M·ADP·Pi, intrinsic rate).
%   Lymn & Taylor 1971 https://doi.org/10.1021/bi00727a037        PMID:4258637
%     Skeletal myosin S1; intrinsic hydrolysis ~100-200 s^-1, equilibrium-limited.
%     The overall steady-state ATPase (~1-5 s^-1 per head at full actin activation)
%     reflects the slowest step (Pi or ADP release), not kah directly.
%     Bound [10, 100] sits below the intrinsic rate; check whether this matches
%     how kah is used in the model (could be the equilibrium-limited apparent rate).
bounds.kah.lb = 10; bounds.kah.ub = 100; bounds.kah.weight = 10;

% Reverse hydrolysis kamh [s^-1] (M·ADP·Pi -> M·ATP).
%   [NO-CITE] Constrained by detailed balance relative to kah. Cardiac value unknown;
%   skeletal: equilibrium constant for hydrolysis step ~9:1 (Lymn & Taylor 1971),
%   giving kamh ~ kah/9. Bound [1, 20] is consistent with this ratio.
bounds.kamh.lb = 1; bounds.kamh.ub = 20; bounds.kamh.weight = 10;

% xrate global multiplier — should be 1 in a well-parameterised model.
% Wide bound to allow legacy fits; narrow target encourages re-anchoring of base rates.
bounds.xrate.lb = 0.5; bounds.xrate.ub = 2; bounds.xrate.weight = 10;

% Maximum unloaded shortening velocity vmax [um/s per half-sarcomere].
%   McDonald 1998 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2231141/ PMID:9700185
%     [OK — rat cardiac skinned myocytes] Force-velocity and power-load curves;
%     Vmax ~1.5 ML/s at 12°C and max tension ~26 kPa. At SL ~2.2 um with ~90
%     half-sarcomeres per cell, ~1.5 ML/s ≈ ~1.7 um/s per HS at 12°C. Q10 ~2.5
%     correction to 25°C gives ~6 um/s per HS. The DOI in sources is non-standard;
%     PMC2231141 is the confirmed open-access link.
%   Difee 2003    https://doi.org/10.1152/japplphysiol.00889.2002   PMID:12533504
%     [OK — rat cardiac] Vmax ~1.43 ML/s at 15°C in skinned rat cardiac myocytes;
%     similarly low due to low temperature. Training increased loaded shortening but
%     not Vmax. Supports lb=4 um/s as the low-temperature floor.
%   Tombe & ter Keurs 1991 https://doi.org/10.1113/jphysiol.1991.sp018682 PMID:1758129
%     [OK — rat cardiac, 25°C] Vmax ~10 um/s per HS at 25°C from force-velocity;
%     directly supports the current model value (vmax=10) and mid-range of bound.
%   NOTE: ub=25 um/s covers ~37°C estimates (further Q10 from 25°C: ×2.5^1.2 ≈ ×3
%     → ~18 um/s). The bound is wide but physiologically motivated across the
%     temperature range.
bounds.vmax.lb = 4; bounds.vmax.ub = 25; bounds.vmax.weight = 10;

% P1 cross-bridge stiffness [model units; ~ pN/um per XB].
%   Kaya & Higuchi 2010 https://doi.org/10.1126/science.1191484    PMID:20689019
%     Rabbit skeletal myosin II; direct measurement of per-head stiffness ~0.5-3
%     pN/nm. Mechanistically analogous to cardiac beta-MHC. Note: model units are
%     not directly pN/nm; the mapping depends on model-level scaling factors.
%   REMOVED: Veigel 2003 (PMID:12695826) was myosin V — WRONG PROTEIN.
%   REMOVED: Mehta 1999 cross-reference — also myosin V, wrong protein.
bounds.kstiff1.lb = 2000; bounds.kstiff1.ub = 30000; bounds.kstiff1.weight = 10;

% P2 cross-bridge stiffness. Post-stroke stiffness higher than pre-stroke by ~1.5-3x.
%   Kaya & Higuchi 2010 https://doi.org/10.1126/science.1191484    PMID:20689019
%     Rabbit skeletal myosin II; nonlinear stiffness; post-stroke state stiffer than
%     pre-stroke. Directly motivates kstiff2 > kstiff1. Bound is wider than the
%     measured ratio to accommodate model-level scaling.
bounds.kstiff2.lb = 5000; bounds.kstiff2.ub = 100000; bounds.kstiff2.weight = 10;

% Serial elastic element stiffness kSE [kPa/um] and exponent ekSE.
% The model kSE captures TOTAL series compliance: experimental rig + myofilament
% lattice (actin + myosin filament extensibility). These two sources are additive
% and together much larger than rig stiffness alone:
%   • Rig (glass/carbon fiber): ~0.5-5 kPa/um (Irving 1992 estimate)
%   • Myofilament lattice: comparable in magnitude to rig (Huxley et al. 1994,
%     skeletal; ~40% of total fiber compliance from filaments)
%   The composite lb=500 kPa/um and current model value kSE=3203 reflect this.
%   Tombe 2016 https://doi.org/10.1016/j.yjmcc.2016.05.013         PMID:27262037
%     [OK — cardiac] Demonstrates that rig compliance shifts SL by ~20% during
%     contraction even when motor position is held; establishes that series stiffness
%     is a major correction in skinned cardiac fiber experiments. Does not report a
%     specific kSE number but motivates the parameter's necessity and magnitude.
%   Irving 1992 https://doi.org/10.1038/357156a0                   PMID:1574083
%     Frog skeletal muscle; myosin head tilting synchronous with force development.
%     Used here as indirect evidence for filament compliance (~0.5-5 kPa/um for the
%     rig component alone; total series stiffness is much larger).
%   NOTE: Previous "Tombe & ter Keurs 1992, 10-100 kPa/um" was a fabricated
%     reference with fabricated values and should be disregarded entirely.
bounds.kSE.lb = 500; bounds.kSE.ub = 20000; bounds.kSE.weight = 10;
bounds.ekSE.lb = 0.5; bounds.ekSE.ub = 2.0; bounds.ekSE.weight = 10;

% Passive (titin) stiffness coefficient k_pas. Bound is intentionally wide:
% the coefficient's magnitude is heavily form-dependent.
%   F_pas = k_pas * (exp(gamma*(SL-Lsc0)) - 1)   → low coeff when gamma ~3
%   F_pas = k_pas * (SL-Lsc0)^gamma              → very high coeff when gamma ~8
%   With gamma=7.9 acting on sub-um (SL-Lsc0)^gamma, the term shrinks to
%   1e-3 — 1e-6, so k_pas must be ~1e3 — 1e5 to recover ~10 kPa passive force.
% Weight=1 (Tier C) until the passive formulation is locked in; then re-derive.
%   Granzier 1995 (above); Linke 2018 https://doi.org/10.1146/annurev-physiol-021317-121234
bounds.k_pas.lb = 1; bounds.k_pas.ub = 1e5; bounds.k_pas.weight = 1;

%% ============================================================
%% TIER C — penalty weight 1 (defensible but uncertain)
%% ============================================================

% IHM mechanosensing OFF/ON rates ksr0, kmsr [s^-1] — NOT biochemical SRX.
% Bounds [0.1, 20] are for the MECHANOSENSING timescale (~0.05-10 s), NOT the
% biochemical SRX lifetime. McNamara 2015 (PMID:25849867) gives ATP turnover
% in SRX ~230 s → exchange rate ~0.004 s^-1. Model rates are 25-5000× faster
% by design: SRX is treated as a fast regulatory reservoir within ODE timescales,
% not a true biochemically inhibited state (see parameter-reference.md §1.2 note).
%   Brunello 2020 https://doi.org/10.1073/pnas.1920632117         PMID:32193335
%     [OK — mouse cardiac] X-ray diffraction of mouse cardiac muscle; directly
%     measures thick-filament mechanosensing dynamics and transition rates on the
%     ~0.1-1 s timescale. Best primary reference for this model's SRX framework.
%   Linari 2015   https://doi.org/10.1038/nature15727             PMID:26560032
%     [OK — skeletal, mechanism conserved] Frog skeletal muscle; force-dependent
%     IHM mechanosensing; rate constants consistent with bound range.
%   Fusi 2016     https://doi.org/10.1038/ncomms13281             PMID:27796302
%     [OK — skeletal] Ca-independent thick-filament regulation; supports OFF rate
%     range and the mechanosensing interpretation.
bounds.ksr0.lb = 0.1; bounds.ksr0.ub = 20; bounds.ksr0.weight = 1;
bounds.kmsr.lb = 0.1; bounds.kmsr.ub = 10; bounds.kmsr.weight = 1;

% Force sensitivity of IHM mechanosensing transitions sigma1, sigma2.
%   Campbell 2018 https://doi.org/10.1016/j.bpj.2018.06.020       PMID:30077333
%     [MODEL-PARAM] Computational cardiac model; sigma values are fitted to
%     reproduce activation dynamics, not directly measured. The bound [5, 100]
%     is anchored to the typical isometric force scale (~50-100 kPa) but the
%     functional form of force-dependence may differ between models.
%   Brunello 2020 (above) provides indirect constraints from measured recruitment
%   dynamics in mouse cardiac — use to verify the fitted sigma range post-optimisation.
bounds.sigma1.lb = 5; bounds.sigma1.ub = 100; bounds.sigma1.weight = 1;
bounds.sigma2.lb = 5; bounds.sigma2.ub = 100; bounds.sigma2.weight = 1;

% SRD (ADP-bound super-relaxed) state — kinetically uncharted in cardiac.
% Structural plausibility from cryo-EM of ADP-stabilised IHM; no direct rate
% measurements exist. Bounds are wide; Tier C reflects this uncertainty.
%   Trivedi 2018 https://doi.org/10.1126/sciadv.aaq1240           PMID:30406211
%     Sci Adv cryo-EM of ADP-stabilised IHM in cardiac/smooth muscle myosin;
%     provides structural basis for the SRD state but no kinetic rates.
bounds.ksrd.lb = 0.1; bounds.ksrd.ub = 20; bounds.ksrd.weight = 1;
bounds.kmsrd.lb = 0.1; bounds.kmsrd.ub = 20; bounds.kmsrd.weight = 1;
bounds.ksr2srd.lb = 1; bounds.ksr2srd.ub = 100; bounds.ksr2srd.weight = 1;
bounds.ksrd2sr.lb = 0.1; bounds.ksrd2sr.ub = 50; bounds.ksrd2sr.weight = 1;
bounds.sigma_srd1.lb = 5; bounds.sigma_srd1.ub = 100; bounds.sigma_srd1.weight = 1;
bounds.sigma_srd2.lb = 5; bounds.sigma_srd2.ub = 100; bounds.sigma_srd2.weight = 1;

% Internal viscous drag mu (shortening), mu_neg (lengthening) [model units].
%   Tombe & ter Keurs 1991 https://doi.org/10.1113/jphysiol.1991.sp018682 PMID:1758129
%     [OK — rat cardiac, direct measurement] Saponin-skinned rat myocardial
%     trabeculae; directly measures an internal viscous element that limits
%     unloaded shortening velocity. Correct species and preparation.
bounds.mu.lb = 0.01; bounds.mu.ub = 5; bounds.mu.weight = 1;
bounds.mu_neg.lb = 0.01; bounds.mu_neg.ub = 5; bounds.mu_neg.weight = 1;

% Maxwell dashpot (titin viscoelasticity) — only active when UseMaxwellDashpot=1.
%   Linke 2018 https://doi.org/10.1146/annurev-physiol-021317-121234
%     [OK — review] Comprehensive review of titin mechanics including viscoelastic
%     properties; provides order-of-magnitude estimates for kSE_M and eta_M.
%   Hamdani 2017 https://doi.org/10.1007/s12551-017-0263-9
%     [PAYWALL-UNVERIFIED] Biophys Rev 2017; likely about titin in cardiac disease.
%     Verify that viscoelastic parameter values in the cited range are reported.
bounds.kSE_M.lb = 10; bounds.kSE_M.ub = 10000; bounds.kSE_M.weight = 1;
bounds.eta_M.lb = 0.01; bounds.eta_M.ub = 10; bounds.eta_M.weight = 1;

% A2 attachment shift (actin target-zone hopping). Geometric basis is solid —
% actin 5.5 nm half-helical repeat sets the target spacing — but the transition
% rate (slope) has no direct cardiac measurement. Tier C reflects this.
bounds.slope.lb = 1000; bounds.slope.ub = 50000; bounds.slope.weight = 1;
bounds.s_threshold_R.lb = 0.001; bounds.s_threshold_R.ub = 0.015; bounds.s_threshold_R.weight = 1;

% Negative-strain stiffness (kstiff1_n, kstiff2_n) — questionable mechanism.
% Bounded loosely to flag extreme asymmetry; ratio kstiff_n/kstiff should
% typically be > 0.05 (< 20x asymmetry) to remain physically plausible.
% [NO-CITE] No primary measurement of negative-strain XB stiffness in cardiac.
bounds.kstiff1_n.lb = 100; bounds.kstiff1_n.ub = 5000; bounds.kstiff1_n.weight = 1;
bounds.kstiff2_n.lb = 100; bounds.kstiff2_n.ub = 5000; bounds.kstiff2_n.weight = 1;

% P3 stiffness (when UseKstiff3=1; otherwise model uses kstiff2). [NO-CITE]
bounds.kstiff3.lb = 5000; bounds.kstiff3.ub = 100000; bounds.kstiff3.weight = 1;

% Calcium / troponin parameters (active when UseCa=1; saturating in current configs).
%   Robertson 1981 https://doi.org/10.1016/S0006-3495(81)84897-4 PMID:7326325
%     [OK — cardiac] Biophys J; Ca^2+/troponin binding kinetics in chemically
%     skinned cardiac muscle. Directly constrains k_on and k_off range.
%   Dobesh 2002    https://doi.org/10.1152/ajpheart.00859.2001   PMID:11834476
%     [OK — cardiac] AJP Heart; cooperativity of Ca activation in cardiac muscle
%     (sarcomere-length dependence); constrains K_coop. Correct preparation.
bounds.K_coop.lb = 1; bounds.K_coop.ub = 20; bounds.K_coop.weight = 1;
bounds.k_on.lb = 10; bounds.k_on.ub = 200; bounds.k_on.weight = 1;
bounds.k_off.lb = 50; bounds.k_off.ub = 1000; bounds.k_off.weight = 1;

% Lattice spacing geometry (when UseLatticeSpacing=1).
%   Williams 2013 https://doi.org/10.1371/journal.pcbi.1003088   PMID:23826173
%     [MODEL-PARAM] PLOS Comput Biol; computational model of actin filament
%     mechanics. Lattice spacing parameters are model-fitted, not directly measured.
%     Better primary source: Irving 2000 https://doi.org/10.1016/S0006-3495(00)76481-4
%     PMID:10777745 (X-ray diffraction of mouse cardiac d10 spacing vs SL).
bounds.d10_ref.lb = 0.030; bounds.d10_ref.ub = 0.045; bounds.d10_ref.weight = 1;
bounds.SL_ref_lattice.lb = 1.9; bounds.SL_ref_lattice.ub = 2.3; bounds.SL_ref_lattice.weight = 1;
bounds.R_thick.lb = 0.005; bounds.R_thick.ub = 0.012; bounds.R_thick.weight = 1;
bounds.R_thin.lb = 0.002; bounds.R_thin.ub = 0.006; bounds.R_thin.weight = 1;
bounds.d_optimal.lb = 0.020; bounds.d_optimal.ub = 0.030; bounds.d_optimal.weight = 1;

%% ============================================================
%% TIER D — weight 0 (listed for documentation/coverage)
%% ============================================================
% These params have no defensible physiological prior. Listed here so the
% coverage check (every dPUdT param accounted for) passes without omission.
% Optimiser sees them but no penalty fires.

% Strain offsets — no direct measurements
bounds.dr0.lb       = -0.02; bounds.dr0.ub       = 0.02; bounds.dr0.weight       = 0;
bounds.dr1.lb       = -0.02; bounds.dr1.ub       = 0.02; bounds.dr1.weight       = 0;
bounds.dr_1.lb      = -0.02; bounds.dr_1.ub      = 0.02; bounds.dr_1.weight      = 0;
bounds.dr2.lb       = -0.05; bounds.dr2.ub       = 0.05; bounds.dr2.weight       = 0;
bounds.dr2_L.lb     = -0.02; bounds.dr2_L.ub     = 0.02; bounds.dr2_L.weight     = 0;
bounds.dr2_R.lb     = -0.02; bounds.dr2_R.ub     = 0.02; bounds.dr2_R.weight     = 0;
bounds.dr3.lb       = -0.05; bounds.dr3.ub       = 0;    bounds.dr3.weight       = 0;
bounds.drp3.lb      = 0;     bounds.drp3.ub      = 0.02; bounds.drp3.weight      = 0;
bounds.drmr.lb      = 0;     bounds.drmr.ub      = 0.02; bounds.drmr.weight      = 0;
bounds.s_threshold_L.lb = 0; bounds.s_threshold_L.ub = 6; bounds.s_threshold_L.weight = 0;

% Strain-sensitivity exponential parameters (used only when
% UsePieceWiseStrainDep=0; current default is on, so these are inactive)
bounds.alpha0.lb   = 0; bounds.alpha0.ub   = 200;  bounds.alpha0.weight   = 0;
bounds.alpha0_L.lb = 0; bounds.alpha0_L.ub = 200;  bounds.alpha0_L.weight = 0;
bounds.alpha0_R.lb = 0; bounds.alpha0_R.ub = 200;  bounds.alpha0_R.weight = 0;
bounds.alpha1.lb   = 10; bounds.alpha1.ub  = 500;  bounds.alpha1.weight   = 0;
bounds.alpha_1.lb  = 0; bounds.alpha_1.ub  = 200;  bounds.alpha_1.weight  = 0;
bounds.alpha2.lb   = 5; bounds.alpha2.ub   = 200;  bounds.alpha2.weight   = 0;
bounds.alpha2_L.lb = 0; bounds.alpha2_L.ub = 200;  bounds.alpha2_L.weight = 0;
bounds.alpha2_R.lb = 0; bounds.alpha2_R.ub = 200;  bounds.alpha2_R.weight = 0;
bounds.alphaRip.lb = 0; bounds.alphaRip.ub = 500;  bounds.alphaRip.weight = 0;
bounds.k2_L.lb     = 0; bounds.k2_L.ub     = 2000; bounds.k2_L.weight     = 0;
bounds.k2_R.lb     = 0; bounds.k2_R.ub     = 20000;bounds.k2_R.weight     = 0;
bounds.e2L.lb      = 1; bounds.e2L.ub      = 4;    bounds.e2L.weight      = 0;
bounds.e2R.lb      = 1; bounds.e2R.ub      = 4;    bounds.e2R.weight      = 0;
bounds.estiff.lb   = 0.5; bounds.estiff.ub = 2;    bounds.estiff.weight   = 0;
bounds.StrainExp.lb = 1; bounds.StrainExp.ub = 4;  bounds.StrainExp.weight = 0;

% Velocity-squared drag (currently effectively zero)
bounds.mu2.lb = 0; bounds.mu2.ub = 0.1; bounds.mu2.weight = 0;

% Other miscellaneous force/strain modifiers
bounds.k2rip.lb = 0; bounds.k2rip.ub = 500; bounds.k2rip.weight = 0;
bounds.k2d.lb   = 0; bounds.k2d.ub   = 50;  bounds.k2d.weight   = 0;
bounds.k_SA.lb  = 0; bounds.k_SA.ub  = 10;  bounds.k_SA.weight  = 0;
bounds.k_catch_bond.lb = 0; bounds.k_catch_bond.ub = 1000; bounds.k_catch_bond.weight = 0;
bounds.CatchBondStrainMax.lb = 0; bounds.CatchBondStrainMax.ub = Inf; bounds.CatchBondStrainMax.weight = 0;
bounds.Srd_max.lb = 0; bounds.Srd_max.ub = 100; bounds.Srd_max.weight = 0;
bounds.Srd_n.lb   = 1; bounds.Srd_n.ub   = 4;   bounds.Srd_n.weight   = 0;

% Vernier velocity / target-zone saturation tuning
bounds.alpha_vernier.lb = 0; bounds.alpha_vernier.ub = 2; bounds.alpha_vernier.weight = 0;
bounds.v_ref_vernier.lb = 0.1; bounds.v_ref_vernier.ub = 10; bounds.v_ref_vernier.weight = 0;
bounds.max_attached_per_bin.lb = 0; bounds.max_attached_per_bin.ub = 1; bounds.max_attached_per_bin.weight = 0;
bounds.A1AttachmentWidth.lb = 0; bounds.A1AttachmentWidth.ub = 0.02; bounds.A1AttachmentWidth.weight = 0;

% Lattice spacing fine-tuning (when UseLatticeSpacing=1)
bounds.sigma_lattice.lb = 0.001; bounds.sigma_lattice.ub = 0.01; bounds.sigma_lattice.weight = 0;

% Lsc0 (passive force reference length) — anatomy-derived but model-specific
bounds.Lsc0.lb = 1.4; bounds.Lsc0.ub = 1.9; bounds.Lsc0.weight = 0;

% Fudge factors (legacy, should be 0)
bounds.FudgeA.lb = 0; bounds.FudgeA.ub = 0; bounds.FudgeA.weight = 0;
bounds.FudgeB.lb = 0; bounds.FudgeB.ub = 0; bounds.FudgeB.weight = 0;
bounds.FudgeC.lb = 0; bounds.FudgeC.ub = 0; bounds.FudgeC.weight = 0;
bounds.FudgeVmax.lb = 0; bounds.FudgeVmax.ub = 0; bounds.FudgeVmax.weight = 0;

% --- PCHIP shape values (Tier D, weight 0) ---
% PCHIP central value should be ~1 by convention (§1b of fitting-strategy.md).
% Boundary clamps (typically values ~50) are computational constraints, not
% physiological measurements. Listed here for documentation; weight = 0.
% The PieceWise* arrays are accessed by *__N indexing in the mods machinery,
% so the bounds for individual shape values would need to be set on the
% indexed names (e.g. bounds.PieceWiseStrainDepParams__3.weight = 0).
% Not enumerated explicitly here.

end
