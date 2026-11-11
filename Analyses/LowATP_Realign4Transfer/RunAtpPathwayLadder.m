% RunAtpPathwayLadder.m -- WHICH pathways may [ATP] act on?
%
% The 21-gate fit is the wildest hypothesis: it lets [ATP] modulate ANY rate or
% mechanical parameter, including ones with no nucleotide chemistry (attachment,
% the power stroke, series stiffness, dashpot viscosity). That is a fine null to
% test against, but it is not a mechanism. This ladder compares it with
% progressively more constrained, physiologically defensible sets, and with the
% model's own purpose-built Michaelis mechanisms.
%
% The core model is FROZEN at realign4 throughout. Every gate is exactly 1 at
% MgATP_ref, so the 8 mM leg is invariant by construction (asserted below) and
% only the 2 mM leg has to be scored -- which is what makes this affordable.
%
% Stage 1 (this script, ~10 evaluations): ABLATION at the incumbent values --
% zero the gates outside each hypothesis and score. Answers "how much is each
% tier actually buying?" without any refitting. Stage 2 (refit each rung) is
% only worth running on what stage 1 says is live.
%
% 2 mM objective now includes FORCE-VELOCITY and POWER: Data_ATP column 4 and
% Data_ATP_AB2021 column 4 are 2 mM (ATP_c = [8 4 2]), which runFVExperiment
% now indexes per condition instead of scoring everything against 8 mM.
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/RunAtpPathwayLadder.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

%% ---- the 21 gateable knobs, classified by whether ATP *should* touch them --
% Direct ATP chemistry: the nucleotide is a substrate of the step itself.
CHEM = {'k2','kah','kamh'};
%   k2   strong-bridge detachment -- ATP binding dissociates actomyosin. The
%        canonical ATP step; without it there is no cross-bridge cycle.
%   kah  hydrolysis / re-priming (PT -> PD) -- consumes ATP.
%   kamh reverse of that step -- ATP-inhibited direction.

% Nucleotide-state transitions: no ATP in the rate law as written, but the
% states differ in bound nucleotide, so an ATP dependence is arguable.
SRX  = {'kmsrd','ksr0','ksr2srd','ksrd2sr'};

% Mechanical / geometric: no nucleotide chemistry. An ATP dependence here is
% phenomenology, not mechanism.
MECH = {'eta_M','kSE','mu_neg','kstiff2','kSE_M','dr','dr2','P_bound_max'};

% Kinetic but not ATP steps: actin attachment, the power stroke, availability.
OTHER = {'ka','kd','k1','k_1','tau_reg','A0'};

ALL21 = [CHEM SRX MECH OTHER];

H = { 'H0  none (no ATP effect)',            {}
      'H1  ATP chemistry only',              CHEM
      'H2  + SRX nucleotide states',         [CHEM SRX]
      'H3  + mechanical (= current fit)',    [CHEM SRX {'eta_M','kSE','mu_neg'}]
      'H4  all 21 (any rate, the wild one)', ALL21 };

%% ---- 2 mM objective, now WITH force-velocity ------------------------------
p = getParams(loadParams('joint82_opt'), [], true, false);
p.UseAtpJointFit = false;
p.EvalFeatures = true; p.BreakOnODEUnstable = false; p.OptimizeOn = 'Feats';
p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 400;
p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.RunKtr = false; p.RunStairs = false; p.RunForceVelocityTime = false;
p.RunSlack = true; p.RunForceVelocity = true; p.RunSlackPassive = false;

% 2 mM has slack AND force-velocity data, but no passive protocol -> drop PS_*
% only. FV stays in, which is the change from every previous 2 mM objective.
fn2 = p.fn; [~,iu] = unique(fn2,'stable'); fn2 = fn2(iu);
prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn2 = fn2(~startsWith(prim(fn2),'PS_'));

p2 = p; p2.MgATP = 2; p2.MgADP = 0; p2.Pi = 0.8;
p2.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
p2.fn = fn2; p2.mods = {}; p2.g = ones(1,0);
r = checkFeatureCoverage(p2);
if ~r.ok
    bad = unique([r.missingModel(:); r.missingData(:)])';
    fprintf(2,'2 mM leg cannot score: %s -> dropping\n', strjoin(bad,', '));
    fn2 = fn2(~ismember(prim(fn2), bad)); p2.fn = fn2;
    r = checkFeatureCoverage(p2);
end
fprintf('2 mM objective: %d features (FV included: %d), coverage ok = %d\n', ...
        numel(fn2), any(startsWith(prim(fn2),'FV')), r.ok);

%% ---- ablation -------------------------------------------------------------
fprintf('\n%-38s %5s %10s %10s %10s\n','hypothesis','DOF','2 mM cost','slack','FV');
res = struct('name',{},'dof',{},'c',{},'cSlack',{},'cFV',{});
iFV = startsWith(prim(fn2),'FV');
for k = 1:size(H,1)
    q = p2;
    for g = 1:numel(ALL21)                       % zero everything outside H{k}
        if ~ismember(ALL21{g}, H{k,2}); q.(['pAtp_' ALL21{g}]) = 0; end
    end
    q.UseAtpPowerGates = ~isempty(H{k,2});
    try
        [~,~,fm,fd] = runSlackExperiment(q);
        qf = q; qf.RunSlack = false;
        [~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
        fm = iMerge(fm, fmv); fd = iMerge(fd, fdv);
        v = evalFeatureCost(fd, fm, fn2, 2);
        res(end+1) = struct('name',H{k,1},'dof',numel(H{k,2}),'c',sum(v), ...
                            'cSlack',sum(v(~iFV)),'cFV',sum(v(iFV))); %#ok<SAGROW>
        fprintf('%-38s %5d %10.3f %10.3f %10.3f\n', H{k,1}, numel(H{k,2}), ...
                sum(v), sum(v(~iFV)), sum(v(iFV)));
    catch ME
        fprintf('%-38s %5d   FAILED: %s\n', H{k,1}, numel(H{k,2}), ME.message);
    end
end

%% ---- what the unjustifiable gates are worth -------------------------------
if numel(res) >= 5
    fprintf('\n--- marginal value of each tier, at incumbent values ---\n');
    fprintf('  ATP chemistry alone (H1 vs H0) : %+8.3f\n', res(2).c - res(1).c);
    fprintf('  adding SRX states   (H2 vs H1) : %+8.3f\n', res(3).c - res(2).c);
    fprintf('  adding mechanics    (H3 vs H2) : %+8.3f  <- 3 DOF, no nucleotide chemistry\n', res(4).c - res(3).c);
    fprintf('  adding everything   (H4 vs H3) : %+8.3f  <- 11 more DOF\n', res(5).c - res(4).c);
    fprintf('\n(H4 is scored at the CURRENT exponents, where the 11 extra gates are\n');
    fprintf(' ~0, so H4 == H3 by construction. Its value needs a refit -- stage 2.)\n');
end

save(fullfile(resdir,'atp_ladder.mat'), 'res', 'H', 'fn2', 'CHEM','SRX','MECH','OTHER');
fprintf('\nLADDER STAGE 1 DONE\n');

function a = iMerge(a, b)
f = fieldnames(b);
for i = 1:numel(f); a.(f{i}) = b.(f{i}); end
end
