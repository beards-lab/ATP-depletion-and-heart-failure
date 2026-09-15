% BuildLean5Seed.m -- the `lean5` campaign seed: slack cycle 5 back in the
% objective, the ablation's pool cut applied, and the free mechanisms actually
% deleted rather than merely un-drawn.
%
% SUPERSEDES BuildFlooredSeed.m. Four changes from `floored_seed`, in the order
% they are applied:
%
% 1. SLACK CYCLE 5 IS SCORED AGAIN.  `rsK80[1:4]` -> `rsK80`.
%
%    The restriction came from Analyses/RecoveryOscillation sec. 3: injecting
%    each run's own measured noise put rsK80's measurement noise at 32.7 %
%    against a between-preparation spread of 16.9 %, which is impossible, and
%    the per-cycle noise of cycle 5 (29-245 %) was the cause. That reasoning is
%    sound about the NOISE and wrong about the CONCLUSION, for two reasons:
%
%      a) Cycle 5 is not a noisy repeat of cycles 1-4, it is a DIFFERENT
%         MANOEUVRE -- the only one in the protocol that re-stretches to
%         ~2.0 um SL instead of 2.2, by the smallest amplitude (rsA 0.107 vs
%         ~0.18) at the slowest ramp, to a plateau of 0.761 vs ~0.99 of the
%         reference. Dropping it removes the protocol's ONLY probe of how the
%         re-stretch recovery depends on the length it recovers at.
%
%      b) It carries a real ATP signal. Measured here, cycle 5 alone has
%         between-preparation spread 19.5 % (vs 15.8 % for cycles 1-4) and
%         ATP z = 3.01, and its 2 mM/8 mM ratio is 0.390 against 0.404 /
%         0.371 / 0.369 / 0.374 for cycles 1-4. A quantity that was 29-245 %
%         transducer noise cannot reproduce a x0.39 ATP effect across three
%         preparations to within 5 % of the other four cycles. The noise
%         estimate and this observation cannot both be complete; the ATP
%         ratio is the harder evidence.
%
%    Scored over 1:5 the spread is 16.9 % (weight 35.0) against 15.8 % (40.1)
%    for 1:4 -- so putting cycle 5 back costs cycles 1-4 about 13 % of their
%    weight. That is the honest price and it is paid deliberately.
%
%    Model/pooledTarget8mM.m MUST be regenerated with cycle 5 in scope first;
%    a chi-2 is only a chi-2 if the residual and the spread come from the same
%    manoeuvre. RunLean5Campaign.m checks this and refuses to start otherwise.
%
% 2. MECHANISMS DELETED (MechanismAblation sec. 0a.C, measured free):
%      UseA2Reattaching      -> 0   : frees kA2re, kA2hop, sA2hop  (+0.0001)
%      UseMaxwellTensionOnly -> 0   : frees x_M_slack (+0.037, rsK80 0.06 BETTER)
%      UseDynamicPassive     -> 0   : exactly inert (c_titin_visc = 0)
%    UseSerialStiffness is NOT touched: it is dead code, so setting it changes
%    nothing. Deleting that flag is a source edit, not a seed edit.
%
% 3. LATTICE CONSTANTS FIXED (MechanismAblation sec. 7.1). The six lattice
%    parameters carry exactly two identifiable combinations. The four that are
%    OBSERVATIONS (an X-ray d10 and the length it was taken at; two EM filament
%    radii) are fixed at the incumbent's own values -- numerically a no-op,
%    which is the point: it changes nothing except that the optimiser stops
%    spending draws on a ridge. d_optimal and sigma_lattice, the two that are
%    genuinely unmeasured, carry the whole mechanism and stay in the pool.
%
% 4. OBJECTIVE: coolDownLS retired, rsK80 x2, F-V pin 0.40 -> 0.25.
%    coolDownLS is an unbounded raw MSE (extractSlackAttributes.m:456), not a
%    preparation-width chi-2 like every other data term, and it dominated the
%    movers of nearly every ablation arm -- swinging by up to 700 in a 0.9
%    objective (sec. 5, sec. 8). rsK80 measures the same window on a proper
%    floored chi-2 and takes double weight to inherit the role. The F-V pin was
%    always arbitrary and the campaign gave no reason to keep it high; the
%    floor already stops F-V being exploitable, so this is belt and braces.
%
% NOT COMPARABLE to 0.9081. The objective changed in three places. The printout
% below reports the base under both objectives so the move is auditable.
%
%   matlab -batch "run('Analyses/RecoveryOscillation/BuildLean5Seed.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('BASE','var');      BASE      = 'rigaware_a_opt'; end
if ~exist('BASETAG','var');   BASETAG   = 'lean5_base';     end
if ~exist('SEEDTAG','var');   SEEDTAG   = 'lean5_seed';     end
if ~exist('NWORKERS','var');  NWORKERS  = 5;                end

%% ---- 0. the target must already carry cycle 5 --------------------------
PT = pooledTarget8mM();
kR = find(strcmp(PT.name,'rsK80'), 1);
if isempty(kR) || ~isequal(sort(PT.seg{kR}(:)'), 1:5)
    error('BuildLean5Seed:targetSegments', ...
        ['Model/pooledTarget8mM.m has rsK80 spread over %s, not 1:5.\n' ...
         'Regenerate it first:\n' ...
         '    clear OPTS; run(''Analyses/DataPooling/WritePooledTarget.m'')'], ...
        mat2str(PT.seg{kR}));
end
fprintf('pooled target OK: rsK80 s = %.1f %% over cycles %s\n', ...
        100*PT.s(kR), mat2str(PT.seg{kR}));

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

%% ---- 1. the reduced model, written as its own base snapshot ------------
p = getParams(); s = loadParams(BASE);
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunKtr=false; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
p.UsePooledTarget = true;
p.mods = {}; p.g = ones(1,0);

fprintf('\n=== %s, before deletion ===\n', BASE);
c_before = evaluateBakersExp(1, p, true, false);
fprintf('  %.4f\n', c_before);

p.UseA2Reattaching      = 0;  p.kA2re = 0; p.kA2hop = 0;
p.UseMaxwellTensionOnly = 0;  p.x_M_slack = 0;
p.UseDynamicPassive     = 0;  p.c_titin_visc = 0;

% The four lattice observations, at the incumbent's own values (see header 3).
p.d10_ref        = 0.0365153;
p.SL_ref_lattice = 2.00579;
p.R_thick        = 0.00825173;
p.R_thin         = 0.00389595;

fprintf('=== after deleting the three free mechanisms ===\n');
c_after = evaluateBakersExp(1, p, true, false);
% MEASURED +0.1145, against the ablation's +0.0001 and +0.037. Not a
% contradiction: those were REOPT figures -- each deletion was followed by a
% budget-matched round with its compensators made compulsory. This is the raw,
% un-refitted cost of applying all three at once, and recovering it is the
% campaign's first job. If 15 h does not recover it, the deletions were not
% free after all and that is a result worth having.
fprintf('  %.4f   (%+.4f raw, no refit; ablation reopt figures were +0.0001 and +0.037)\n', ...
        c_after, c_after - c_before);

outBase = fullfile(root, 'params', [BASETAG '.m']);
writeParams17g(outBase, p, sprintf(['%s with UseA2Reattaching/UseMaxwellTensionOnly/' ...
    'UseDynamicPassive deleted and the four lattice observations fixed; ' ...
    'cost %.6f (was %.6f)'], BASE, c_after, c_before));
fprintf('written: %s\n', outBase);

%% ---- 2. re-derive every weight from the regenerated target -------------
% Everything below is BuildWeightedObjective's job: it re-reads the pooled
% spreads, recomputes w = lambda*R/s^2 for every surviving term, re-emits the
% one-preparation-width floors against the NEW s, and picks lambda so the data
% subtotal is preserved. Hand-patching weights here would desynchronise them
% from the target, which is the failure campaignSpec.m exists to prevent.
BASE          = BASETAG;
SPREAD_SOURCE = 'pooled';
POOLED_TARGET = true;
FLOOR_CHI2    = 1.0;                        % one preparation-width of deadband
FN_RENAME     = {'rsK80[1:4]', 'rsK80'};    % cycle 5 back (header 1)
FN_DROP       = {'coolDownLS'};             % retired   (header 4)
W_BOOST       = {'rsK80', 2.0};             % inherits coolDownLS's role
ADD_FEATS     = {};
FV_SHARE      = 0.25;
SHARE_CAP     = 0.30;
ANCHOR_GUARD  = true;
ANCHOR_ENTRY  = 'steady[1]|52.5-87.7|5';

fprintf('\nBuildLean5Seed: %s -> params/%s.m\n\n', BASE, SEEDTAG);
run(fullfile(root, 'Analyses', 'FeatureWeighting', 'BuildWeightedObjective.m'));
%% ---- note --------------------------------------------------------------
% The base snapshot is written by Auxiliary/writeParams17g -- %.17g, not
% writeParamsToMFile's %g. Six significant digits perturbs every parameter by
% ~1e-6 relative, and on this objective %.6g serialisation alone is worth 0.60
% of cost (Docs objective_noise_floor). No local functions live in this file:
% it `run`s BuildWeightedObjective, which brings its own.
