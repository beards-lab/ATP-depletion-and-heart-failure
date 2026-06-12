% params_reseeded.m
% =========================================================================
% Manual rate-reseeding overlay for DriverBoundedFit (USE_RESEEDED=true).
% Runs LAST in the driver — AFTER iter_17 and AFTER the inline SRX/parking
% block — so values set here are the final word for the tuned rates.
%
% Scope: the base kinetic rates that iter_17 left BELOW their physiological
% floors (parameterBounds.m, mouse/rat alpha-MHC, ~20-25 C). kah is
% deliberately NOT set here — the inline block sets kah=100 and couples it to
% the SRX mixing (ksr2srd='=kah'); tune kah there.
%
% Sources (see Model/parameterBounds.m, Docs/sources.md):
%   ka  = 375  Beard 2022 parent-model fit to rat cardiac FV/ktr at 25C;
%              matches getParams default. (was 37.7, floor 50)
%   kd  = 100  Wang 2013 mouse alpha skinned ~111/s @17C; Little 2012 rat
%              alpha ~159/s @0-ADP. (was 10.0, floor 20)
%   k2  = 350  Lowey 2013 mouse alpha-S1 ADP release ~350/s @20C.
%              (was 69.5, floor 100)  <-- primary XTOR lever; lower toward
%              the floor (100) to slow cycling / reduce isometric ATPase.
%
% MANUAL FITTING NOTES (update as we iterate):
%   - XTOR ~ kah*PT (hydrolysis flux) = k2eff*attached (detach flux) in SS.
%     With kah=100 (inline) and k2=350, XTOR runs ~36 (target [3,10]).
%     Levers to lower XTOR: lower k2 (slows cycle, raises attached, lowers
%     ktr), or lower inline kah, or increase SRX parking.
% =========================================================================

params0.ka = 50;      % attachment (floor 50) -- iter5: 75->50 (cut fapp further; XTOR/ktr down)
params0.kd = 100;     % P1 detachment    (floor 20)

% iter5: cut k2eff via the R2 detachment strain curve, keeping k2 base at floor (100)
% and preserving the +/- strain wings (positive=stretch detach, negative=slip->Vmax).
% Central knots (small positive strain, where isometric p2 sits) scaled ~x0.6 to slow
% working-stroke detachment -> lowers gapp (ktr) and cycle flux (XTOR) without dropping
% force. Physiological: heads bearing force at the working stroke detach slowly; fast
% only when strained away. iter_17 knots: [50 15.97 3.006 1.029 1.0475 50.68 50 50].
% iter9: cut ONLY the isometric working-stroke knots (3,4 @ s2~0.026,0.013) to slow
% isometric detachment (ktr/XTOR down), but RESTORE knot5 (s2=0, the shortening-transition
% region) to original 1.0475 so heads still detach during shortening -> FV not dragged.
% iter_17: [50 15.97 3.006 1.029 1.0475 50.68 50 50]; iter5-8 over-cut knot5 to 0.629.
params0.PieceWiseStrainDep2Params = [50 15.9726 1.80 0.617 1.0475 50.6847 50 50];

% iter7 REVERTED: UseVelGaussAttachment=true fixed FV+XTOR (reduced isometric attach)
% but DESTROYED ktr (velocity-gating corrupts the isometric force redevelopment ->
% ktr 150, rmse 5.7, oscillatory). Structurally incompatible with clean ktr. Keep OFF.
params0.UseVelGaussAttachment = false;

% iter10: lower viscous drag (mu) to flatten FV — drag subtracts force during shortening,
% so less drag -> more force at velocity -> flatter FV_fnorm, without touching kinetics.
params0.mu     = 0.015;   % was 0.0404
params0.mu_neg = 0.015;

% iter11: more FV datapoints for robustness (match Data_ATP velocities 0..6 ML/s,
% was only [0 0.5 1 2 4]). Denser FV constraint guards against overfitting 5 points.
params0.FV_velocities = -[0 0.5 1 2 3 4 5 6];

% iter12 REVERTED: sigma2 42->25 had the sign backwards — larger F/sigma2 => SMALLER
% exp(-F/sigma2) => LESS parking at operating force => SRX collapsed to 0.03 (out of
% bound) and A-trend did not steepen. SRX pool is too small at active force to drive LDA.
% params0.sigma2 = 25;

% iter13: length-dependent ATTACHMENT via lattice spacing (the physiological basis of
% cardiac LDA / Frank-Starling). d10 widens at shorter SL -> f_lattice cuts attachment ->
% lower redeveloped A at the deeper slacks. Retune so the cut spans the slack SL range
% (2.04->1.88): d_optimal=0.0246 (f=1 at SL>=2.04), sigma_lattice=0.0025 -> f_lattice ~
% [1 .99 .95 .90 .82] -> ~18% force drop (data A drops 21%). Was too weak at default
% d_optimal=0.025 (only ~6%). At SL>=2.04 (incl. FV isometric at SL0=2.2) f_lattice=1.
params0.UseLatticeSpacing = true;
params0.d_optimal     = 0.0242;   % iter14: 0.0246->0.0242 (slightly stronger LDA to
                                  % center the A curve, which sat ~3 kPa high everywhere)
params0.sigma_lattice = 0.0025;


% --- kah -> XTOR (CONFIRMED, manual fit iter 2) ---
% XTOR ~= kah*PT (hydrolysis-flux-limited): at kah=100, PT~0.42 => XTOR~42;
% dropping kah to its Tier-A floor (40) gave XTOR~19.6 = 40*0.49, and pulled
% SRX_ss and attached_ss into bounds. Grand total 525 -> 134.
% SRX mixing held fixed at 100 (decoupled from the inline '=kah') so only the
% hydrolysis flux changed. NOTE the structural tension: XTOR target is [3,10]
% but kah floor is 40, so XTOR can't reach <=10 unless PT is drained below
% ~0.25 (more SRX parking) OR hydrolysis is made near-equilibrium (raise kamh
% toward Khyd~1-4 so XTOR becomes downstream-limited instead of = kah*PT).
params0.kah     = 80;    % iter3: 40->80 (alpha-myosin physiological; turnover capped downstream)
params0.ksr2srd = 100;   % hold SRX mixing (decoupled from inline '=kah')
% k2 lowered 350 -> 110 (manual fit iter 1): at k2=350 model ktr~157 (data ~49,
% 3x too fast) and XTOR~36 (target [3,10]). ktr scales ~with k2eff, so
% k2 ~= 350*(49/157) ~= 110 should pull ktr->~50 and XTOR->~10, and lengthen
% attached dwell (raising attached_ss / A / steady). Near floor 100; defensible
% (Wang 2013 fibre-level ADP-release ~65-160/s).
params0.k2 = 100;     % P2 ADP release (floor 100) -- iter3: 110->100 (min gapp base)
params0.k1 = 140;     % power stroke P1->P2 -- iter9: back to 140 (iter6 best; k1=200 raised A,
                      % didn't fix FV -> FV steepness is from R2 cut, not k1)
