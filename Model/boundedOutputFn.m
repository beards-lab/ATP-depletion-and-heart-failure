function fn = boundedOutputFn(assertScale)
%BOUNDEDOUTPUTFN  Reusable feature/boundary list (params0.fn) for bounded fitting.
%
%   FN = BOUNDEDOUTPUTFN() returns the cell array assigned to params0.fn,
%   combining three layers:
%     (a) DATA-FIT features   — scored vs Baker-lab experimental curves
%                               (FV, ktr, slack transient). NOT literature bounds.
%     (b) OUTPUT boundaries   — 'field|lb-ub|weight' literature ranges on model
%                               OUTPUTS (per-head ATPase, state occupancies),
%                               grouped into one joint-cost tile. Mouse/rat
%                               alpha-MHC, ~21-25 C, maximal Ca, isometric
%                               pre-slack steady state unless noted.
%     (c) INPUT guard         — 'AssertParams|XX' auto-expands (in RunBakersExp)
%                               to bound every ACTIVE parameter against
%                               parameterBounds.m. XX scales the penalty weight.
%
%   FN = BOUNDEDOUTPUTFN(ASSERTSCALE) sets the AssertParams weight scale
%   (default 0.001 = soft input-bound hint). Pass 0 (or []) to OMIT the
%   input guard entirely (e.g. when input bounds are enforced separately via
%   evalPhysiologyCost).
%
%   Consumed by: evalFeatureCost / plotFeatures (see those for the '|lb-ub|w'
%   grammar and grouped-tile semantics).
%
%   ----------------------------------------------------------------------
%   OUTPUT-BOUND PROVENANCE (all DOIs/PMIDs verified 2026-06 to resolve)
%   ----------------------------------------------------------------------
%   TARGET: mouse / adult-rat ventricle = ALPHA-MHC (MYH6), the fast isoform.
%   Skeletal/beta sources are flagged and alpha-corrected.
%
%   XTOR  — per-head cross-bridge ATP turnover (= kah*PT flux), ISOMETRIC,
%           maximal Ca. Whole-population average (incl. parked heads), NOT the
%           solution actin-activated Vmax. Best: rat alpha ~4-5 s^-1; mouse
%           alpha ~6-8 s^-1 at 20-25 C. Bound [3,10].
%     Ebus & Stienen 1996  https://doi.org/10.1006/jmcc.1996.0164  PMID:8877784
%       Rat skinned cardiac, 20 C, max Ca: isometric ATPase 0.48 mM/s, force
%       53 kN/m^2 -> ~4-5 s^-1 per head.
%     Rossmanith 1995      https://doi.org/10.1111/j.1440-1681.1995.tb02034.x
%       Rat V1 (100% alpha), 24 C: tension-cost slope; alpha ~2.5-3x beta.
%     de Tombe 2007        https://doi.org/10.1113/jphysiol.2007.138693  PMID:17717017
%       Rat myocardium; Q10 ~2.4-3.5 for ktr/ATPase (temperature scaling).
%     CAUTION: kah (hydrolysis step ~80-150 s^-1) is NOT XTOR; XTOR is the
%       downstream-limited steady-state turnover (~order 5 s^-1). See parameterBounds.m.
%
%   XTOR_vmax — per-head ATP turnover at UNLOADED / max shortening velocity
%           (Fenn effect). ATPase rises above isometric during shortening.
%           Fenn factor ~1.3-2.5x in cardiac. Bound [5,15].
%     Ebus & Stienen 1996  PMID:8877784  — Fenn peak ~1.7x isometric ATPase
%       during length changes in rat cardiac (same prep as XTOR above).
%     Stienen 2004         https://doi.org/10.1161/01.RES.0000080536.06932.E3
%       Creatine-phosphate consumption during cardiac shortening; ATPase-ktr link.
%     Fusi 2017            https://doi.org/10.1113/JP273299       PMID:27763660
%       [SKEL] Duty ratio at V0 ~0.05 (vs ~0.22 isometric); per-attached-head
%       turnover elevated at Vmax. Skeletal benchmark for the shortening regime.
%
%   SRX_ss — super-relaxed / thick-filament OFF fraction (model SR+SRD),
%           ISOMETRIC maximal activation (NOT resting ~40-60%). Biochemical
%           SRX ~0.10-0.30; structural OFF (X-ray) higher (0.50-0.80). Model SR
%           pool is treated as the biochemical/regulatory OFF reservoir, so
%           bound [0.10,0.40]. (At Vmax the OFF fraction RISES as active stress
%           falls — mechanosensing; add an SRX_vmax feature if needed.)
%     Hooijman 2011        https://doi.org/10.1016/j.bpj.2011.02.061  PMID:21504733
%       Cardiac SRX ~40-48% at rest; residual SRX persists at activation.
%     Reconditi 2017       https://doi.org/10.1073/pnas.1619484114
%       Rat trabecula X-ray; ON fraction proportional to force; substantial
%       structural OFF remains at peak force.
%     Brunello 2020        https://doi.org/10.1073/pnas.1920632117  PMID:32193335
%       Mouse/rat cardiac; ~10% of motors bear peak force; OFF recovers rapidly.
%     Ma 2022              https://doi.org/10.1085/jgp.202213213
%       Synthetic cardiac thick filaments: SRX ~15% -> ~3% from pCa8 -> pCa4
%       (isolated filaments; intact-fibre residual is higher).
%     Park-Holohan 2021    https://doi.org/10.1073/pnas.2023706118
%       Active stress ~8x more effective than passive at switching filament ON;
%       predicts re-OFF during low-force (shortening) phases.
%
%   attached_ss — attached fraction (model p1_0+p2_0 = weak P1 + strong P2,
%           i.e. TOTAL attached incl. weakly bound, NOT the strong-only duty
%           ratio), ISOMETRIC maximal activation. Strong-only ~0.10-0.15;
%           total incl. weak ~0.15-0.45. Bound [0.10,0.45]. (At Vmax the
%           attached fraction COLLAPSES to ~0.01-0.06.)
%     Pinzauti 2018        https://doi.org/10.1113/JP275579   PMC6023834 (J Physiol 596:2581)
%       Rat cardiac trabecula: per-head stiffness 1.07 pN/nm, force ~6 pN/head;
%       strongly-attached fraction ~0.10-0.15 at maximal Ca.
%     Reconditi 2017       https://doi.org/10.1073/pnas.1619484114
%       ~10% of motors bear peak systolic force (rat trabecula X-ray).
%     Fusi 2017            https://doi.org/10.1113/JP273299   PMID:27763660
%       [SKEL] Attached fraction falls ~4x from isometric to V0 — basis for
%       the steep drop at Vmax.
%
%   PT_ss  — disordered-relaxed (DRX), detached-but-available fraction.
%           Bound [0.20,0.70] by conservation: SRX + attached + DRX + (PuR+NP)
%           = 1; DRX is the dominant pool at maximal isometric activation
%           (~0.40-0.55). Hooijman 2011 (PMID:21504733): fast (DRX) phase
%           dominates in activated cardiac muscle.
%
%   DATA-FIT features (FV_fnorm|FV_v, ktr, A, peak*, vall*, steady, t0_crossing,
%   restretchSlope*, ovrsht*) are scored against the Baker-lab experimental
%   curves (params0.FV_dataset, default 'Baker2022'), NOT against literature
%   constants. 'ktr_rmse|0-0.2' is a fit-quality (residual) bound.
%
%   See also: parameterBounds, evalFeatureCost, plotFeatures, assertActiveParams,
%             extractSlackAttributes, DriverBoundedFit

if nargin < 1 || isempty(assertScale)
    assertScale = 0.001;            % soft input-bound hint by default
end

% (a) data-fit features + (b) grouped literature OUTPUT bounds.
% [1] selects ATP/segment index 1 (8 mM). Weights: XTOR 15, SRX_ss 5,
% attached_ss 10 (others default 1) — carried from the working parametrization.
fn = {
    'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', ...
    'ktr_rmse|0-0.2|.1', ...
    ['XTOR[1]|3-10|15,' ...          % per-head ATPase, isometric  (Ebus&Stienen 1996; Rossmanith 1995)
     'XTOR_vmax[1]|5-15,' ...        % per-head ATPase at Vmax      (Fenn ~1.7x; Ebus&Stienen 1996; Fusi 2017)
     'SRX_ss[1]|0.10-0.40|5,' ...    % SRX/OFF fraction, active     (Reconditi 2017; Brunello 2020; Ma 2022)
     'attached_ss[1]|0.10-0.45|10,' ...% total attached (weak+strong)(Pinzauti 2018; Reconditi 2017)
     'PT_ss[1]|0.20-0.70'], ...      % DRX available fraction       (conservation; Hooijman 2011)
    't0_crossing|SLdiff|2', ...
    'restretchSlopeStart|0.1', ...
    'peak1_y|10', 'peak1_dSL|0.2', ...
    'vall_y|10', 'vall_t|0.2', 'peak2|5', ...
    'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1' ...
};

% (c) input-parameter guard (auto-expanded in RunBakersExp via assertActiveParams)
if assertScale > 0
    fn{end+1} = sprintf('AssertParams|%g', assertScale);
end

end
