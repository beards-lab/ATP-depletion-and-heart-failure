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
%   ================ TARGET SPECIES / ISOFORM (read this first) ================
%   The modelling target is MOUSE (and adult rat) CARDIAC muscle. Adult mouse
%   ventricle is ~100% ALPHA-myosin heavy chain (α-MHC / MYH6); adult rat
%   ventricle is ≥90% α-MHC. α-MHC is the FAST cardiac isoform.
%
%   This matters enormously for the RATE constants. Most "classic" cardiac
%   myosin biochemistry was done on the SLOW β-MHC isoform (human ventricle,
%   bovine, porcine, rabbit) or on skeletal myosin. β-MHC and skeletal-slow
%   rates are ~2–10× SLOWER than α-MHC for the kinetically distinguishing
%   steps (ADP release, ATP hydrolysis, attachment turnover, Vmax). Bounds
%   below were re-anchored (2026-06) onto α-MHC data where it exists; entries
%   that still rely on a β/skeletal proxy are flagged [β/SKEL PROXY] and given
%   a lower weight to reflect the isoform mismatch.
%
%   Reference temperature: the model's base rate constants appear calibrated
%   near ~20–25 °C (e.g. default kah=80 ≈ α-MHC forward hydrolysis at 20 °C;
%   default ka=373 ≈ Beard 2022 25 °C fit). Most α-MHC transient-kinetics data
%   are at 20 °C. Bounds are therefore stated for ~20–25 °C, with deliberate
%   upper headroom for Q10 extrapolation toward 37 °C (Q10 ≈ 2–2.5 for the
%   ATPase steps; Q10 ≈ 4–5 for Vmax in cardiac trabeculae, de Tombe 1990).
%
%   NOTE on BASE rates vs xrate: bounds apply to BASE rates (the rate at zero
%   strain, where the PCHIP central value ≈ 1 — which is what the literature
%   reports). updateRates.m multiplies several turnover rates by `xrate`; if
%   xrate<1 the effective in-ODE rates can fall below these bounds (flagged in
%   parameter-reference.md). That is a modelling concern, not something the
%   bounds here are meant to absorb.
%
%   Tier organisation (weights now encode BOTH evidence quality AND isoform
%   match for our mouse-α target):
%       Tier A (weight 100): direct mouse/rat α-MHC measurement, or a hard
%                            structural constant. High trust.
%       Tier B (weight 10):  order-of-magnitude / model-fit / α-MHC proxy.
%       Tier C (weight 1):   defensible but uncertain, or wrong-isoform proxy,
%                            or model-state mapping ambiguous.
%       Tier D (weight 0):   no prior available; listed for coverage only.
%
%   Note on SRX nomenclature: ksr0, kmsr, sigma1, sigma2 implement the IHM
%   thick-filament mechanosensing ON/OFF transition (Linari 2015 Nature,
%   Brunello 2020 PNAS [mouse/rat cardiac]) — NOT the slow biochemical SRX of
%   McNamara 2015 (~0.004 s^-1). Bounds reflect the mechanosensing
%   interpretation (rate constants ~10–450 s^-1 in Marcucci 2017), but are
%   left wide because the model's SRX timescale is itself a modelling choice.
%
%   Passive/titin parameters (gamma, k_pas, kSE_M, eta_M, Lsc0) are OWNED BY
%   the separate titin identification and were intentionally NOT revised here.
%
%   See also: evalPhysiologyCost, plotPhysiologyDashboard, parameter-reference.md

bounds = struct();

%% ============================================================
%% TIER A — penalty weight 100 (direct mouse/rat α-MHC or hard structural)
%% ============================================================

% Power stroke distance dr [um]. Single-molecule optical-trap, CARDIAC myosin.
%   Re-anchored 2026-06 from skeletal (Capitanio/Finer) to direct CARDIAC
%   measurements — the working stroke does NOT differ measurably between α and
%   β cardiac isoforms (only the kinetics do), so cardiac single-molecule data
%   apply directly to our mouse-α target. Consensus unloaded stroke 7–10 nm.
%   Sugiura 1998   https://doi.org/10.1161/01.RES.82.10.1029     PMID:9622155
%     Rat cardiac V1(α) and V3(β) by optical trap: ~9–10 nm; no isoform
%     difference in unitary displacement. Direct rat-α datum.
%   Tyska 2000     https://doi.org/10.1161/01.RES.86.7.737       PMID:10764406
%     Mouse α-cardiac single-molecule (R403Q HCM model vs WT); confirms fast
%     mouse α (2.3x ATPase, 2.2x force vs WT); working stroke in the 6–10 nm range.
%   Caremani 2016  https://doi.org/10.1073/pnas.1525057113       PMID:26984499
%     Rat cardiac IN SITU: working stroke 3 nm/hs at high load → ~8 nm unloaded
%     (apparent reduction under load = XB elasticity, not a smaller intrinsic stroke).
%   Palmiter 1999  https://doi.org/10.1111/j.1469-7793.1999.0669n.x PMID:10457082
%     Rabbit cardiac V1/V3: ~7 nm; confirms no α-vs-β step-size difference.
%   Bound widened to 6–11 nm to span mouse-α single-molecule (6 nm) to the
%   upper cardiac estimates (10–11 nm). Model default dr=10 nm sits in range.
bounds.dr.lb = 0.006; bounds.dr.ub = 0.011; bounds.dr.weight = 100;

% --- ATP HYDROLYSIS rate kah [s^-1] (M·ATP -> M·ADP·Pi) — PROMOTED to Tier A ---
% This is the step the user specifically asked to anchor. DIRECT mouse/human
% α-MHC measurements now exist; the hydrolysis step is ~10× FASTER in α than β,
% so old skeletal/Lymn-Taylor framing is replaced by α-MHC numbers.
%   Deacon 2012   https://doi.org/10.1007/s00018-012-0927-3      PMID:22349210
%     Tissue-purified MOUSE cardiac α-S1, 20 °C: hydrolysis k+3+k-3 = 150 s^-1
%     (single fast phase). Human α-S1 = 168 s^-1; human β-S1 = 17 s^-1 (10× slower).
%   Mijailovich/Walklate 2017 https://doi.org/10.1016/j.bpj.2017.01.021 PMID:28297657
%     Full-cycle model fit: human α-S1 FORWARD hydrolysis kH = 77 s^-1 (β = 12.5),
%     20 °C. Model default kah=80 ≈ this α forward value.
%   Lymn & Taylor 1971 https://doi.org/10.1021/bi00801a004       PMID:4258719
%     Original skeletal burst ~100–150 s^-1, 20 °C — α-cardiac behaves like fast
%     skeletal, confirming the α value (NOT the slow β-cardiac value).
%   Bound [40,200] brackets the α-MHC forward step at 20–25 °C with Q10 headroom.
%   NOTE: kah is the intrinsic HYDROLYSIS step rate, NOT the steady-state
%   actin-activated ATPase per head (~4–8 s^-1 for mouse α; that is the XTOR
%   model output, rate-limited downstream by Pi/ADP release — do not conflate).
%   Promoted weight 10 -> 100: direct mouse-α measurement, high trust.
bounds.kah.lb = 40; bounds.kah.ub = 200; bounds.kah.weight = 100;

% --- DETACHMENT rate kd [s^-1] (P1 -> PD, base) — PROMOTED to Tier A ---
% Direct rat/mouse α-MHC myofibril/fibre detachment now available.
%   Little 2012  https://doi.org/10.1074/jbc.M111.337295         PMID:22718768
%     [OK — α-MHC] Ventricular myofibril XB detachment via fluorescent cTnC:
%     RAT (α) = 159 s^-1 (0 ADP) / 42 s^-1 (+2 mM ADP) at 15 °C; human (β) = 41 s^-1.
%     (Corrects the previously-listed DOI 10.1074/jbc.M112.380048 / PMID:22685295,
%      which was wrong.) Rat α detaches ~4× faster than human β.
%   Wang 2013    https://doi.org/10.1016/j.yjmcc.2012.10.010     PMID:23123290
%     [OK — α-MHC] Mouse α-MyHC skinned strips, sinusoidal analysis: detachment
%     ~111 s^-1 at 17 °C; rat α ~65 s^-1. Direct mouse datum.
%   Model default kd≈103 sits mid-range. Bound [20,250]: lower end covers the
%   ADP-slowed regime (~40) and rat-at-physiological-ADP; upper covers 0-ADP rat
%   (159) plus modest Q10 headroom. Promoted weight 10 -> 100.
bounds.kd.lb = 20; bounds.kd.ub = 250; bounds.kd.weight = 100;

% Cardiac titin passive-force power-law exponent gamma [dimensionless].
% [PASSIVE/TITIN — owned by separate titin identification; NOT revised 2026-06]
% Active only when UseTitinIdentifiedPassive=1; otherwise the "gamma" in
% the model is a placeholder with a non-physiological value (~2.67).
%   Granzier & Irving 1995 https://doi.org/10.1016/S0006-3495(95)80278-X PMID:7756567
%   Fukuda 2003    https://doi.org/10.1161/01.RES.0000054065.16294.A1 PMID:12970176
%   Wu 2000        https://doi.org/10.1161/01.RES.86.6.682        PMID:10747009
%   Terui 2008     https://doi.org/10.1085/jgp.200809864          PMID:18955594
bounds.gamma.lb = 7; bounds.gamma.ub = 10; bounds.gamma.weight = 100;

% Actin half-helical repeat d_actin [um] (5.46 nm structural constant).
% Species-independent; sets the myosin target-site periodicity for A2 hopping.
%   Holmes 1990    https://doi.org/10.1038/347044a0               PMID:2395461
%     Atomic F-actin model; 5.46 nm axial half-repeat from fibre diffraction.
%   Egelman & DeRosier 1992 https://doi.org/10.1016/0304-3991(92)90025-G
%     EM image analysis; independently confirms the half-repeat.
%   Bound [5.4, 7.4 nm] spans half-repeat to full monomer spacing.
bounds.d_actin.lb = 0.0054; bounds.d_actin.ub = 0.0074; bounds.d_actin.weight = 100;

% Initial sarcomere length SL0 [um] (experimental setpoint, mouse cardiac).
% Re-anchored to MOUSE cardiac skinned-preparation operating range.
%   Cazorla 2001   https://doi.org/10.1161/hh1001.090876         PMID:11375272
%     [OK — mouse] Skinned MOUSE cardiac myocytes tested over SL 1.9–2.3 um;
%     maximal Ca-activated force toward the long end of that range.
%   Ait-Mou 2016   https://doi.org/10.1073/pnas.1516732113       PMID:26862168
%     [OK — mouse] Mouse cardiac structural rearrangements measured over SL 1.9–2.3 um.
%   Gordon 1966    https://doi.org/10.1113/jphysiol.1966.sp007909 PMID:5966030
%     Canonical force-length operating range (frog skeletal); historical anchor.
%   Bound [1.8, 2.4] retained (slack ~1.75–1.85; experimental optimum ~2.1).
bounds.SL0.lb = 1.8; bounds.SL0.ub = 2.4; bounds.SL0.weight = 100;

% Bare zone half-length L_hbare [um] (thick filament anatomy; structural).
% Total central bare zone ~160–170 nm -> half ~80–85 nm. Mouse-confirmed.
%   Woodhead & Craig 2023 https://doi.org/10.1038/s41586-023-06690-5  PMID:37914933
%     Cryo-EM of NATIVE cardiac thick filament: 170 nm bare zone (half = 85 nm).
%   Zoghbi 2008    https://doi.org/10.1073/pnas.0708912105        PMID:18272487
%     3D structure of vertebrate cardiac myosin filament; bare-zone dimensions.
%   Bound tightened to [0.075, 0.105] (half-zone 75–105 nm). Model default
%   L_hbare=0.10 sits at the upper edge (the measured value is ~0.08–0.085).
bounds.L_hbare.lb = 0.075; bounds.L_hbare.ub = 0.105; bounds.L_hbare.weight = 100;

% Thick filament length L_thick [um]. Robust vertebrate constant ~1.60 um,
% directly confirmed in native cardiac. (Interpreted here as full bipolar length.)
%   Woodhead & Craig 2023 https://doi.org/10.1038/s41586-023-06690-5  PMID:37914933
%     Cryo-EM native cardiac thick filament = 1.6 um.
%   Al-Khayat 2013 https://doi.org/10.1073/pnas.1212708110       PMID:23169664
%     Atomic model of human cardiac thick filament; 1600 nm.
%   Bound [1.55, 1.70] retained to admit the model default L_thick=1.67, but the
%   MEASURED value is 1.60 um. FLAG: some saved param sets use L_thick=1.14 um,
%   which is far outside this range and geometrically wrong (see parameter-reference.md).
bounds.L_thick.lb = 1.55; bounds.L_thick.ub = 1.7; bounds.L_thick.weight = 100;

% Thin filament length L_thin [um] (cardiac, Z-disk to pointed end).
% Re-anchored to DIRECT mouse/rat cardiac measurement — strongest geometric prior.
%   Burgoyne 2008  https://doi.org/10.1093/cvr/cvm117            PMID:18178575
%     [OK — mouse/rat] Electron tomography of MOUSE and RAT cardiac: mean thin
%     filament length = 1.04 um (range 0.94–1.10).
%   Tonino 2017    https://doi.org/10.1016/j.yjmcc.2017.04.010   PMID:27139341
%     [OK — mouse] Mouse cardiac confocal/STORM: ~1.03 um at rest; mildly
%     SL-dependent; independent of titin/nebulin isoform.
%   Robinson & Winegrad 1979 https://doi.org/10.1038/280099a0    PMID:
%     Classic EM ~1.0 um in cardiac thin filaments.
%   Bound tightened to [0.95, 1.15] (best 1.03–1.05). FLAG: model default
%   L_thin=1.20 (getParams) and the 1.377 seen in some param sets BOTH exceed
%   this — they should be revised down toward ~1.04 um.
bounds.L_thin.lb = 0.95; bounds.L_thin.ub = 1.15; bounds.L_thin.weight = 100;

%% ============================================================
%% TIER B — penalty weight 10 (order-of-magnitude / model-fit / α-MHC proxy)
%% ============================================================
%
% NOTE: bounds apply to BASE rates. Under the convention that PCHIP central
% value (s ~ 0) ~= 1, the base rate equals the rate at zero strain — which
% is what the literature reports.

% Attachment rate ka [s^-1] (PD -> P1, weakly-bound -> strongly-bound).
% No direct single-molecule "ka"; constrained by ktr (tension redevelopment)
% and by the parent model's own fit.
%   Beard 2022   https://doi.org/10.1016/j.bpj.2022.07.029       PMID:35918899
%     [OK — same architecture] The parent model (this codebase's ancestor) fits
%     ka = 373 s^-1 to rat cardiac FV/ktr data at 25 °C. Model default ka=373.23.
%   de Tombe 2007 https://doi.org/10.1113/jphysiol.2007.138693   PMID:17717017
%     [OK — rat α] Temperature dependence of ktr in skinned rat myocardium:
%     ktr,max rises with temperature (Q10 ~2.4–3.5); at 22–30 °C ktr ~15–30 s^-1,
%     implying an attachment attempt rate of a few hundred s^-1 given the duty ratio.
%     (Previously mis-cited here as "Stelzer 2007 / jgp.200709815" — wrong paper.)
%   Bound widened to [50, 600] (was 50–500): ktr + duty-ratio arguments place
%   the attachment attempt rate in the few-hundred range for α at 20–25 °C.
bounds.ka.lb = 50; bounds.ka.ub = 600; bounds.ka.weight = 10;

% Power stroke forward k1 [s^-1] (P1 -> P2, isomerisation) — DEMOTED to Tier C.
% Model-state mapping is genuinely ambiguous: the RAW mechanical power-stroke
% rate is very fast, whereas ODE-lumped "k1" is much slower. Low weight + wide
% bound reflects that we cannot pin which quantity this parameter represents.
%   Caremani 2016 https://doi.org/10.1073/pnas.1525057113        PMID:26984499
%     [OK — rat α] Rat cardiac IN SITU working-stroke SPEED 1,000–6,000 s^-1
%     (load-dependent) — the raw mechanical transition.
%   Land 2017    https://doi.org/10.1016/j.yjmcc.2017.03.008      PMID:28392437
%     [MODEL-PARAM, β/human] Lumped weak->strong rate; ODE models use 40–200 s^-1.
%   Model default k1≈40 reflects the LUMPED convention; bound [40,3000] admits
%   both that and the raw mechanical rate. Weight 10 -> 1 (ambiguous mapping).
bounds.k1.lb = 40; bounds.k1.ub = 3000; bounds.k1.weight = 1;

% Power stroke reverse k_1 [s^-1] (P2 -> P1) — DEMOTED to Tier C.
% Constrained only via the stroke equilibrium K1 = k1/k_1 ≈ 3–10 (Caremani 2016
% implies K2~3–5 for cardiac). No direct α reverse-rate measurement.
%   Caremani 2016 (above) — load-dependent stroke implies reversibility, K~3–5.
%   Bound widened to [0,500] to track k1's upper range while keeping K1 plausible.
bounds.k_1.lb = 0; bounds.k_1.ub = 500; bounds.k_1.weight = 1;

% ADP-release / post-stroke detachment rate k2 [s^-1] (P2 -> P3).
% RE-ANCHORED to α-MHC: ADP release is the kinetic step that most separates α
% (fast) from β (slow); the old citations (Malmqvist/Siemankowski) were β-MHC
% and ~3× too slow for our target.
%   Lowey 2013   https://doi.org/10.1074/jbc.M113.462846         PMID:23564459
%     [OK — mouse α] Transgenic MOUSE α-S1 ADP release ~350 s^-1 at 20 °C
%     (β-S1 = 130). ADP affinity is ~4× weaker in α (K_AD ~100 vs ~25 uM).
%   Deacon 2012  https://doi.org/10.1007/s00018-012-0927-3       PMID:22349210
%     [OK — mouse/human α] ADP release >1000 s^-1 at 20 °C for α-S1 — so fast it
%     does NOT limit detachment (unlike β). Biochemical upper anchor.
%   Wang 2013    https://doi.org/10.1016/j.yjmcc.2012.10.010     PMID:23123290
%     [OK — α] Fibre-level post-stroke detachment ~65–160 s^-1 at 17 °C (lower
%     than isolated-S1 biochemistry due to ensemble/load effects).
%   Bound shifted UP to [100,1500] (was 50–1000): brackets fibre-level α (~100)
%   to isolated-S1 α biochemistry (~350 to >1000). Model default k2≈419.
%   Weight kept at 10 (Tier B) ONLY because the P2->P3 state ↔ measurable mapping
%   is ambiguous (biochemical ADP release vs effective fibre detachment).
bounds.k2.lb = 100; bounds.k2.ub = 1500; bounds.k2.weight = 10;

% k2 reverse k_2 [s^-1] (P2 re-attachment / ADP rebinding to post-stroke state).
% Second-order ADP on-rate k+AD ≈ 2.8–3.5 uM^-1 s^-1 for mouse/rat α (Lowey
% 2013 PMID:23564459; Pereira 2001 PMID:11076938). The first-order rate scales
% with [MgADP]; the default protocol runs at MgADP=0, so the effective reverse
% rate ~ 0. Bound [0,50] retained for low-but-nonzero ADP conditions.
bounds.k_2.lb = 0; bounds.k_2.ub = 50; bounds.k_2.weight = 10;

% ATP-binding + detachment rate k3 [s^-1] (P3 -> PT).
% In the model k3 collapses ATP binding AND (largely) Pi release into one step.
% At saturating [MgATP] (5–8 mM) the ATP-binding sub-step is effectively
% instantaneous for α (K1k+2 ≈ 1.6–2.0 uM^-1 s^-1 -> >10,000 s^-1 apparent;
% mouse α, Lowey 2013 PMID:23564459 / Deacon 2012 PMID:22349210). The rate-
% LIMITING content of k3 is therefore Pi release.
%   Mijailovich/Walklate 2017 https://doi.org/10.1016/j.bpj.2017.01.021 PMID:28297657
%     Modeled α-cardiac Pi release ≈ 32 s^-1 at 20 °C (β ≈ 16; rabbit skel ≈ 45).
%   Dantzig 1992 https://doi.org/10.1113/jphysiol.1992.sp019251  PMID:1460122
%     Pi release ~200–400 s^-1 in rabbit psoas at 10 °C — upper conceptual bound.
%   Model default k3≈44 matches the α Pi-release estimate at 20 °C. Bound
%   [15,200] (was 20–200): lower covers Pi release at cool temps, upper covers
%   Q10 toward 37 °C. Weight kept 10 (α Pi release is modeled, not directly measured).
bounds.k3.lb = 15; bounds.k3.ub = 200; bounds.k3.weight = 10;

% k3 reverse k3m [s^-1] — DEMOTED to Tier C. No direct cardiac measurement and
% no detailed-balance anchor in the current scheme; model default k3m=0.
bounds.k3m.lb = 0; bounds.k3m.ub = 20; bounds.k3m.weight = 1;

% Reverse hydrolysis kamh [s^-1] (M·ADP·Pi -> M·ATP).
% Set by the hydrolysis equilibrium K_hyd = kah/kamh. Classic skeletal/β value
% K_hyd ≈ 9–10 (Lymn & Taylor 1971; Bagshaw & Trentham 1974 PMID:4154547), but
% the α-cardiac fit gives a LOWER K_hyd ≈ 4 (Mijailovich/Walklate 2017
% PMID:28297657), i.e. relatively more M·ATP. No direct α measurement exists.
%   With kah ~80 and K_hyd 4–10 -> kamh ~8–20 s^-1. Model uses kamh='=0.1*kah'
%   (K_hyd=10). Bound widened to [1,30] to admit the α K_hyd≈4 case (kamh up to
%   ~20–30), which can matter for ATP-depletion behaviour. Weight kept 10.
bounds.kamh.lb = 1; bounds.kamh.ub = 30; bounds.kamh.weight = 10;

% xrate global multiplier — should be 1 in a well-parameterised model.
% Wide bound to allow legacy fits; narrow target encourages re-anchoring of base rates.
bounds.xrate.lb = 0.5; bounds.xrate.ub = 2; bounds.xrate.weight = 0;

% Maximum unloaded shortening velocity vmax [um/s per half-sarcomere].
% Mouse/rat are α-MHC and FAST; Vmax also has an unusually strong temperature
% dependence in cardiac (Q10 ≈ 4–5), so the plausible band is wide.
%   de Tombe & ter Keurs 1990 https://doi.org/10.1161/01.RES.66.5.1239  PMID:2335024
%     [OK — rat α, intact] Vmax ~10 um/s/HS at 20 °C rising to ~25 um/s/HS at
%     25 °C; Q10 ≈ 4.6. Directly supports model default vmax=10 at the cool end.
%   McDonald 1998 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2231141/ PMID:9706028
%     [OK — rat α, skinned] Vmax ~1.5 ML/s at 12 °C — sets the low-temperature floor.
%   de Tombe & ter Keurs 1992 https://doi.org/10.1113/jphysiol.1992.sp019283 PMID:1474506
%     [OK — rat α, 25 °C] FV-derived Vmax ~10 um/s/HS; mid-range anchor.
%   Bound [4,30] (was 4–25): low-T floor 4; upper raised to 30 for warm rat-α
%   (~25–30 um/s/HS at 25 °C). NB unit convention (per-HS vs ML/s) is the main
%   source of spread — confirm against the Baker-lab rig temperature before
%   tightening. Weight kept 10.
bounds.vmax.lb = 4; bounds.vmax.ub = 30; bounds.vmax.weight = 10;

% P1 cross-bridge stiffness [model units; ~ pN/um per XB].
% Re-cited to CARDIAC: per-head stiffness of cardiac (α) myosin is ~1 pN/nm,
% i.e. 2–3× SOFTER than fast skeletal (~2–3 pN/nm). Model units are not pN/nm,
% so the mapping is indirect -> weight stays 1.
%   Pinzauti 2018 https://doi.org/10.1113/JP275579   PMC6023834 (J Physiol 596:2581)
%     [OK — rat α] Rat cardiac trabecula: per-head stiffness ~1 pN/nm; force
%     ~6 pN/head at Tmax. Cardiac-specific.
%   Kaya & Higuchi 2010 https://doi.org/10.1126/science.1191484  PMID:20689019
%     Skeletal nonlinear stiffness; cross-reference only (softer in cardiac).
bounds.kstiff1.lb = 2000; bounds.kstiff1.ub = 100000; bounds.kstiff1.weight = 1;

% P2 cross-bridge stiffness. Post-stroke stiffness > pre-stroke (state-dependent).
%   Pinzauti 2018 (above) / Caremani 2016 PMID:26984499 — rat α cardiac stiffness.
%   Kaya & Higuchi 2010 PMID:20689019 — post-stroke state stiffer than pre-stroke
%     (motivates kstiff2 > kstiff1). Model units indirect -> weight 1.
bounds.kstiff2.lb = 5000; bounds.kstiff2.ub = 100000; bounds.kstiff2.weight = 1;

% Serial elastic element stiffness kSE [kPa/um] and exponent ekSE.
% Captures TOTAL series compliance (experimental rig + myofilament lattice),
% not titin — so retained here (NOT part of the skipped passive/titin set).
%   Tombe 2016 https://doi.org/10.1016/j.yjmcc.2016.05.013       PMID:27262037
%     [OK — cardiac] Rig compliance shifts SL ~20% during contraction; establishes
%     series stiffness as a major correction in skinned cardiac fibre experiments.
%   Irving 1992 https://doi.org/10.1038/357156a0                 PMID:1574083
%     Indirect evidence for filament/rig compliance (~0.5–5 kPa/um rig component;
%     total series stiffness is much larger). Composite lb=500; model kSE=3203.
bounds.kSE.lb = 500; bounds.kSE.ub = 20000; bounds.kSE.weight = 10;
bounds.ekSE.lb = 0.5; bounds.ekSE.ub = 2.0; bounds.ekSE.weight = 10;

% Passive (titin) stiffness coefficient k_pas.
% [PASSIVE/TITIN — owned by separate titin identification; NOT revised 2026-06]
% Wide, form-dependent bound retained at weight 1 until the passive formulation
% is locked in by the titin identification, then re-derive.
%   Granzier 1995 (above); Linke 2018 https://doi.org/10.1146/annurev-physiol-021317-121234
bounds.k_pas.lb = 1; bounds.k_pas.ub = 1e5; bounds.k_pas.weight = 1;

%% ============================================================
%% TIER C — penalty weight 1 (defensible but uncertain)
%% ============================================================

% IHM mechanosensing OFF/ON rates ksr0, kmsr [s^-1] — NOT biochemical SRX.
% RE-ANCHORED 2026-06. The header note (and Brunello 2020) frame these as the
% fast THICK-FILAMENT MECHANOSENSING transition, for which fitted cardiac rate
% constants are ~10–450 s^-1 — far faster than the old bounds implied, and the
% old kmsr bound [0.1,10] actually EXCLUDED the model's own default kmsr=250.
%   Marcucci 2017 https://doi.org/10.1038/s41598-017-05999-2     PMID:28717219
%     [OK — rat cardiac] Mechanosensing model fit to rat trabecula (27 °C):
%     OFF->ON kmin = 10.2 s^-1, kmax = 442 s^-1; ON->OFF ≈ 262 s^-1
%     (force-independent). Model defaults ksr0=9 (≈kmin) and kmsr=250 (≈kON->OFF)
%     ARE these mechanosensing rates.
%   Brunello 2020 https://doi.org/10.1073/pnas.1920632117        PMID:32193335
%     [OK — mouse/rat cardiac] X-ray diffraction; OFF<->ON structural transitions
%     on the ~10–100 ms timescale; rapid OFF recovery after force drops (kmsr>>ksr0).
%   McNamara 2015 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5425749/ PMID:25849867
%     Context only: the SLOW biochemical SRX is ~0.004 s^-1 (~230 s) — a different
%     phenomenon, NOT what these rates model.
%   Bounds left WIDE because the model's chosen SRX timescale is itself a design
%   choice: ksr0 [0.1,50] and kmsr [0.1,500] now admit BOTH the slow-reservoir
%   regime used by some optimized param sets (kmsr~1) and the fast mechanosensing
%   regime of the getParams default (kmsr~250). Weight kept 1 (timescale genuinely
%   uncertain — exactly the case for a low weight).
bounds.ksr0.lb = 0.1; bounds.ksr0.ub = 50;  bounds.ksr0.weight = 1;
bounds.kmsr.lb = 0.1; bounds.kmsr.ub = 500; bounds.kmsr.weight = 1;

% Force sensitivity of IHM mechanosensing transitions sigma1, sigma2 [kPa].
% Force scale over which recruitment occurs is ~10–40 kPa (threshold ~17 kPa,
% half-max ~20–40 kPa) for rat/mouse cardiac.
%   Marcucci 2017 (above) — threshold ~17 kPa, saturating ~12 kPa above threshold.
%   Park-Holohan 2021 https://doi.org/10.1073/pnas.2023706118    PMID:33850019
%     [OK — rat cardiac] Folded-motor sensitivity to active stress; refines the
%     force-dependence of OFF<->ON.
%   Campbell 2018 https://doi.org/10.1016/j.bpj.2018.06.020      PMID:30077333
%     [MODEL-PARAM] sigma values are fitted, not directly measured.
%   Bound [5,100] retained (model sigma1=33 in range). CAUTION: getParams default
%   sigma2=1e6 is a "force-independent" sentinel and will read as a large penalty
%   here — optimized param sets use sigma2~40; treat the sentinel as out-of-scope.
bounds.sigma1.lb = 5; bounds.sigma1.ub = 100; bounds.sigma1.weight = 1;
bounds.sigma2.lb = 5; bounds.sigma2.ub = 100; bounds.sigma2.weight = 1;

% SRD (ADP-bound super-relaxed) state — kinetically uncharted in cardiac.
% Structural plausibility from cryo-EM of ADP-stabilised IHM; no direct rate
% measurements exist. Bounds are wide; Tier C reflects this uncertainty.
%   Trivedi 2018 https://doi.org/10.1126/sciadv.aaq1240           PMID:30406211
%     Sci Adv cryo-EM of ADP-stabilised IHM in cardiac/smooth muscle myosin;
%     structural basis for the SRD state but no kinetic rates.
bounds.ksrd.lb = 0.1; bounds.ksrd.ub = 20; bounds.ksrd.weight = 1;
bounds.kmsrd.lb = 0.1; bounds.kmsrd.ub = 20; bounds.kmsrd.weight = 1;
bounds.ksr2srd.lb = 1; bounds.ksr2srd.ub = 100; bounds.ksr2srd.weight = 1;
bounds.ksrd2sr.lb = 0.1; bounds.ksrd2sr.ub = 50; bounds.ksrd2sr.weight = 1;
bounds.sigma_srd1.lb = 5; bounds.sigma_srd1.ub = 1000; bounds.sigma_srd1.weight = 1;
bounds.sigma_srd2.lb = 5; bounds.sigma_srd2.ub = 1000; bounds.sigma_srd2.weight = 1;

% Internal viscous drag mu (shortening), mu_neg (lengthening) [model units].
%   de Tombe & ter Keurs 1992 https://doi.org/10.1113/jphysiol.1992.sp019283 PMID:1474506
%     [OK — rat cardiac, direct] Saponin-skinned rat trabeculae; an internal
%     viscous element limits unloaded shortening velocity. Correct species/prep.
bounds.mu.lb = 0.01; bounds.mu.ub = 5; bounds.mu.weight = 1;
bounds.mu_neg.lb = 0.01; bounds.mu_neg.ub = 5; bounds.mu_neg.weight = 1;

% Maxwell dashpot (titin viscoelasticity) — only active when UseMaxwellDashpot=1.
% [PASSIVE/TITIN — owned by separate titin identification; NOT revised 2026-06]
%   Linke 2018 https://doi.org/10.1146/annurev-physiol-021317-121234
%   Hamdani 2017 https://doi.org/10.1007/s12551-017-0263-9   [PAYWALL-UNVERIFIED]
bounds.kSE_M.lb = 10; bounds.kSE_M.ub = 10000; bounds.kSE_M.weight = 1;
bounds.eta_M.lb = 0.01; bounds.eta_M.ub = 10; bounds.eta_M.weight = 1;

% A2 attachment shift (actin target-zone hopping). Geometric basis is solid —
% actin 5.46 nm half-helical repeat sets the target spacing (see d_actin) — but
% the transition rate (slope) has no direct cardiac measurement. Tier C.
bounds.slope.lb = 1000; bounds.slope.ub = 50000; bounds.slope.weight = 1;
bounds.s_threshold_R.lb = 0.001; bounds.s_threshold_R.ub = 0.015; bounds.s_threshold_R.weight = 1;

% Negative-strain stiffness (kstiff1_n, kstiff2_n) — questionable mechanism.
% Single-molecule data suggest myosin stiffness under compression is COMPARABLE
% to tension (Kaya & Higuchi 2010), so the ~100× softening used in some param
% sets is unphysical (mechanism-evaluation.md §1.8); the physically-correct knob
% for limiting negative-strain force is fast detachment, not reduced stiffness.
% Bounded loosely to flag extreme asymmetry; ratio kstiff_n/kstiff should
% typically be > 0.05 (< 20× asymmetry). [NO-CITE — no primary cardiac measurement.]
bounds.kstiff1_n.lb = 100; bounds.kstiff1_n.ub = 50000; bounds.kstiff1_n.weight = 1;
bounds.kstiff2_n.lb = 100; bounds.kstiff2_n.ub = 50000; bounds.kstiff2_n.weight = 1;

% P3 stiffness (when UseKstiff3=1; otherwise model uses kstiff2). [NO-CITE]
bounds.kstiff3.lb = 5000; bounds.kstiff3.ub = 100000; bounds.kstiff3.weight = 1;

% Calcium / troponin parameters (active when UseCa=1; SATURATING in current
% configs, so not rate-limiting -> low weight despite good data).
% NOTE: cTnC Ca-off rate is strongly context-dependent — isolated cTnC ~700–1000
% s^-1, full thin filament + cycling S1 ~13 s^-1 (all at 15 °C). Model k_off=307
% sits in the isolated/complex range. Re-cited to rat cardiac cTnC.
%   Davis 2007 https://doi.org/10.1529/biophysj.106.095406       PMID:17293397
%     [OK — rat/human cardiac cTnC] koff: cTn complex ~42, thin filament ~105,
%     thin filament + S1 ~13 s^-1 (15 °C). Definitive context-dependence.
%   Tikunova & Davis 2004 https://doi.org/10.1074/jbc.M405413200 PMID:15205455
%     [OK — cardiac cTnC] WT site-II koff 700–800 s^-1 at 4 °C; kon ~200–400 uM^-1 s^-1.
%   Dobesh 2002 https://doi.org/10.1152/ajpheart.00859.2001      PMID:11834476
%     [OK — cardiac] Cooperativity / SL-dependence; nH ~3–4 (constrains K_coop).
bounds.K_coop.lb = 1; bounds.K_coop.ub = 20; bounds.K_coop.weight = 1;
bounds.k_on.lb = 10; bounds.k_on.ub = 200; bounds.k_on.weight = 1;
bounds.k_off.lb = 50; bounds.k_off.ub = 1000; bounds.k_off.weight = 1;

% Lattice spacing geometry (when UseLatticeSpacing=1).
%   Irving 2000 https://doi.org/10.1016/S0006-3495(00)76481-4   PMID:10777745
%     [OK — mouse cardiac] X-ray diffraction of mouse cardiac d10 spacing vs SL —
%     better primary source than the model-fitted Williams 2013.
%   Williams 2013 https://doi.org/10.1371/journal.pcbi.1003088  PMID:23826173
%     [MODEL-PARAM] Lattice parameters are model-fitted, not directly measured.
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
bounds.max_attached_per_bin.lb = 0; bounds.max_attached_per_bin.ub = 1; bounds.max_attached_per_bin.weight = 0; % [DEPRECATED, dS-dependent]
% Global occupancy ceiling: max total bound fraction. Upper bound = rat-cardiac
% max-Ca isometric bound fraction (~0.12-0.15; range to ~0.20). Lower bound keeps
% it from collapsing attachment entirely.
bounds.P_bound_max.lb = 0.02; bounds.P_bound_max.ub = 0.40; bounds.P_bound_max.weight = 0;
% Per-bin density cap (dS-invariant) for the deprecated per-strain-bin form.
bounds.rho_attach_max.lb = 1; bounds.rho_attach_max.ub = 100; bounds.rho_attach_max.weight = 0;
bounds.A1AttachmentWidth.lb = 0; bounds.A1AttachmentWidth.ub = 0.02; bounds.A1AttachmentWidth.weight = 0;

% Lattice spacing fine-tuning (when UseLatticeSpacing=1)
bounds.sigma_lattice.lb = 0.001; bounds.sigma_lattice.ub = 0.01; bounds.sigma_lattice.weight = 0;

% Lsc0 (passive force reference length) — anatomy-derived but model-specific
% [PASSIVE/TITIN — owned by separate titin identification; NOT revised 2026-06]
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
