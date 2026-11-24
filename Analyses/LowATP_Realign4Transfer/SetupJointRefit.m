% SetupJointRefit.m -- build and VALIDATE the joint 8+2 mM objective before
% committing hours of optimiser time to it.
%
% Checks, in order:
%   1. both legs are fully scorable (checkFeatureCoverage) -- an unscorable fn
%      entry costs a flat MISSING_FEATURE_COST=100 per eval and silently
%      dominates. This is what produced the +1.200 artifact in the last round.
%   2. the joint seed cost, split into its 8 mM and 2 mM legs.
%   3. the 8 mM leg is INVARIANT to the ATP gates (it must be, by construction).
%   4. what the structural lever (UseAtpDetachR2 + UseStrainDepAtpGate, replacing
%      the strain-uniform pAtp_k2) does at the seed, over a small (kAtpStrain,
%      sAtpStrain0) grid -- before letting the optimiser loose on it.
%
% Writes results/jointsetup.mat and params/jointrefit_seed.m.
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/SetupJointRefit.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

VT2 = 'protocol_03_27_2026_2mM_slack.mat';

%% ---- 1. the seed: current best of both worlds ---------------------------
% realign4 core + the refitted gates. The gate block already carries MgATP=2 etc,
% so start from realign4 and copy ONLY the gate exponents across.
p = getParams(loadParams('realign4_opt'), [], true, false);
gsrc = getParams(loadParams('lowATP_realign4_gates_opt'), [], true, false);

GATES = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd', ...
         'ksrd2sr','tau_reg','A0','eta_M','kstiff2','kSE','kSE_M','mu_neg', ...
         'dr','dr2','P_bound_max'};
p.UseAtpPowerGates = 1; p.MgATP_ref = 8;
nz = {};
for i = 1:numel(GATES)
    f = ['pAtp_' GATES{i}];
    v = gsrc.(f);
    if abs(v) < 1e-3; v = 0; end        % optimiser dust -> exactly zero
    p.(f) = v;
    if v ~= 0; nz{end+1} = GATES{i}; end %#ok<SAGROW>
end
fprintf('active gates (%d): %s\n', numel(nz), strjoin(nz,', '));

% 8 mM condition + full battery
p.MgATP = 8; p.MgADP = 0; p.Pi = 0;
p.RunForceVelocity = true; p.RunSlack = true; p.RunSlackPassive = true;
p.RunKtr = false; p.RunStairs = false; p.RunForceVelocityTime = false;
p.EvalFeatures = true; p.BreakOnODEUnstable = false; p.OptimizeOn = 'Feats';
p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 300;
p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';

%% ---- 2. the 2 mM leg's feature list -------------------------------------
% Drop FV*/PS_* (no 2 mM data), de-duplicate, then drop anything the 2 mM leg
% cannot actually produce -- verified below, not assumed.
fn2 = p.fn; [~,iu] = unique(fn2,'stable'); fn2 = fn2(iu);
prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn2 = fn2(~(startsWith(prim(fn2),'FV') | startsWith(prim(fn2),'PS_')));

pc2 = p; pc2.MgATP = 2; pc2.Pi = 0.8; pc2.velocitytableonfile = VT2;
pc2.RunForceVelocity = false; pc2.RunSlackPassive = false;
pc2.fn = fn2; pc2.mods = {}; pc2.g = ones(1,0);
r2 = checkFeatureCoverage(pc2);
if ~r2.ok
    bad = unique([r2.missingModel(:); r2.missingData(:)]');
    fprintf(2,'2 mM leg cannot score: %s  -> dropping\n', strjoin(bad,', '));
    fn2 = fn2(~ismember(prim(fn2), bad));
    pc2.fn = fn2; r2 = checkFeatureCoverage(pc2);
end
fprintf('2 mM leg: %d features, coverage ok = %d\n', numel(fn2), r2.ok);

pc8 = p; pc8.mods = {}; pc8.g = ones(1,0);
r8 = checkFeatureCoverage(pc8);
if ~r8.ok
    bad8 = unique([r8.missingModel(:); r8.missingData(:)]');
    fprintf(2,'8 mM leg cannot score: %s  -> dropping\n', strjoin(bad8,', '));
    p.fn = p.fn(~ismember(prim(p.fn), bad8));
    pc8 = p; pc8.mods = {}; pc8.g = ones(1,0);
    r8 = checkFeatureCoverage(pc8);
end
fprintf('8 mM leg: %d features, coverage ok = %d\n', numel(p.fn), r8.ok);

%% ---- 3. wire the joint objective ----------------------------------------
p.UseAtpJointFit = true;
p.AtpJoint_w8 = 1; p.AtpJoint_w2 = 1;
p.AtpJoint_MgATP2 = 2; p.AtpJoint_Pi2 = 0.8;
p.AtpJoint_vt2 = VT2; p.AtpJoint_fn2 = fn2;
p.AtpJoint_RunFV2 = false; p.AtpJoint_RunPS2 = false;
p.mods = {}; p.g = [];

fprintf('\n--- seed, joint ---\n');
t0 = tic; [Et, ~, fc, fnm] = evaluateBakersExp(1, p, true, false);
i8 = startsWith(fnm,'hi:'); i2 = startsWith(fnm,'lo:');
c8 = sum(fc(i8)); c2 = sum(fc(i2));
fprintf('JOINT SEED  Et = %.4f   [ 8 mM %.4f | 2 mM %.4f ]   (%.0f s)\n', Et, c8, c2, toc(t0));

%% ---- 4. the 8 mM leg MUST be gate-invariant ------------------------------
pz = p; for i=1:numel(GATES); pz.(['pAtp_' GATES{i}]) = 0; end
[~,~,fcz,fnz] = evaluateBakersExp(1, pz, true, false);
c8z = sum(fcz(startsWith(fnz,'hi:'))); c2z = sum(fcz(startsWith(fnz,'lo:')));
fprintf('gates OFF   ->  [ 8 mM %.4f | 2 mM %.4f ]\n', c8z, c2z);
fprintf('8 mM invariance |d| = %.3e  (MUST be 0)\n', abs(c8z - c8));
assert(abs(c8z-c8) < 1e-9, 'ATP gates are leaking into the 8 mM leg');

%% ---- 5. the structural lever, at the seed --------------------------------
% Replace the strain-UNIFORM pAtp_k2 with the strain-DEPENDENT Michaelis gate:
%   Km_eff(s) = K_T2 * (1 + kAtpStrain * max(0, s - sAtpStrain0))
% Both are 1 at MgATP_ref, so 8 mM stays invariant either way.
fprintf('\n--- structural lever scan (UseAtpDetachR2 + UseStrainDepAtpGate) ---\n');
fprintf('%10s %10s %10s %10s %10s\n','kAtpStrain','sAtpStrain0','joint','8 mM','2 mM');
KA = [0 20 50 100 200]; S0 = [0 0.005 0.01];
best = struct('c',Inf);
G = [];
for ka = KA
    for s0 = S0
        if ka == 0 && s0 ~= 0; continue; end
        q = p; q.pAtp_k2 = 0;                 % hand k2's ATP dependence to the gate
        q.UseAtpDetachR2 = 1; q.K_T2 = 2.5;
        q.UseStrainDepAtpGate = (ka > 0);
        q.kAtpStrain = ka; q.sAtpStrain0 = s0;
        try
            [~,~,fcg,fng] = evaluateBakersExp(1, q, true, false);
            a8 = sum(fcg(startsWith(fng,'hi:'))); a2 = sum(fcg(startsWith(fng,'lo:')));
            fprintf('%10g %10g %10.4f %10.4f %10.4f\n', ka, s0, a8+a2, a8, a2);
            G(end+1,:) = [ka s0 a8+a2 a8 a2]; %#ok<SAGROW>
            if a8+a2 < best.c; best = struct('c',a8+a2,'ka',ka,'s0',s0,'c8',a8,'c2',a2); end
        catch ME
            fprintf('%10g %10g   FAILED: %s\n', ka, s0, ME.message);
        end
    end
end
fprintf('\nbest structural point: kAtpStrain=%g sAtpStrain0=%g -> %.4f [%.4f | %.4f]\n', ...
        best.ka, best.s0, best.c, best.c8, best.c2);
fprintf('vs plain-gate seed %.4f [%.4f | %.4f]\n', Et, c8, c2);

save(fullfile(resdir,'jointsetup.mat'), 'Et','c8','c2','G','best','fn2','GATES','nz');

%% ---- 6. write the seed snapshot the optimiser will start from -----------
seed = p;
if best.c < Et - 0.25          % only adopt the lever if it is a real improvement
    seed.pAtp_k2 = 0; seed.UseAtpDetachR2 = 1; seed.K_T2 = 2.5;
    seed.UseStrainDepAtpGate = best.ka > 0;
    seed.kAtpStrain = best.ka; seed.sAtpStrain0 = best.s0;
    fprintf('\nADOPTING the structural lever into the seed.\n');
else
    fprintf('\nStructural lever NOT adopted at the seed (no clear gain); optimiser can still find it.\n');
end
writeParamsToMFile(fullfile(root,'params','jointrefit_seed.m'), seed, [], ...
    sprintf('joint 8+2 mM refit seed: realign4 + refitted gates%s', ...
            ternary(seed.UseStrainDepAtpGate,' + strain-dep ATP gate','')));
fprintf('\nSETUP DONE\n');

function out = ternary(c,a,b); if c; out=a; else; out=b; end; end
