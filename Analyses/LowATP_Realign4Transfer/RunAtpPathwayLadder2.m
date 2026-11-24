% RunAtpPathwayLadder2.m -- stage 2: REFIT each rung of the pathway ladder.
%
% Stage 1 ablated at the incumbent exponents, which is biased AGAINST the
% constrained hypotheses: those exponents were fitted with all 10 gates active
% and with NO 2 mM force-velocity data in the objective, so zeroing a subset
% leaves the survivors at values optimised for a different problem. Every
% constrained rung can therefore only improve when refitted. This does that, so
% the comparison is fair.
%
% Core frozen at realign4; only the gate exponents move. 8 mM is invariant by
% construction (asserted once at the end), so only the 2 mM leg is scored --
% now including force-velocity and power at 2 mM.
%
% Coordinate descent on ABSOLUTE exponent values (not multipliers) so a gate can
% move off exactly zero and can change sign. Accept threshold 0.25 = the measured
% peak-to-peak scatter; anything finer is fitting the solver.
%
% Run: matlab -batch "BUDGET_H=2.5; run('Analyses/LowATP_Realign4Transfer/RunAtpPathwayLadder2.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

if ~exist('BUDGET_H','var'); BUDGET_H = 2.5; end
ACCEPT = 0.25;

L = load(fullfile(resdir,'atp_ladder.mat'));
CHEM = L.CHEM; SRX = L.SRX; MECH = L.MECH; OTHER = L.OTHER;
ALL21 = [CHEM SRX MECH OTHER];
fn2 = L.fn2;

H = { 'H1  ATP chemistry only',              CHEM
      'H2  + SRX nucleotide states',         [CHEM SRX]
      'H3  + mechanical (= current fit)',    [CHEM SRX {'eta_M','kSE','mu_neg'}]
      'H4  all 21 (any rate, the wild one)', ALL21 };

%% ---- the 2 mM objective ---------------------------------------------------
p = getParams(loadParams('joint82_opt'), [], true, false);
p.UseAtpJointFit = false;
p.EvalFeatures = true; p.BreakOnODEUnstable = false; p.OptimizeOn = 'Feats';
p.RunSlackSegments = 'AllPar'; p.MaxRunTime = 400;
p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.RunKtr = false; p.RunStairs = false; p.RunForceVelocityTime = false;
p.RunSlack = true; p.RunForceVelocity = true; p.RunSlackPassive = false;
p.MgATP = 2; p.MgADP = 0; p.Pi = 0.8;
p.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
p.fn = fn2; p.mods = {}; p.g = ones(1,0);

cost = @(q) iCost(q, fn2, ATP_c, Data_ATP);
nEval = 0; tAll = tic;

%% ---- refit each rung ------------------------------------------------------
OUT = struct('name',{},'dof',{},'c0',{},'c',{},'p',{},'evals',{});
budgetPerRung = BUDGET_H*3600 / size(H,1);

for k = 1:size(H,1)
    gates = H{k,2};
    q = p;
    for g = 1:numel(ALL21)
        if ~ismember(ALL21{g}, gates); q.(['pAtp_' ALL21{g}]) = 0; end
    end
    q.UseAtpPowerGates = true;
    c = cost(q); nEval = nEval + 1; c0 = c;
    fprintf('\n=== %s  (%d gates) seed %.4f ===\n', H{k,1}, numel(gates), c);

    tR = tic; improved = true; sweep = 0;
    while improved && toc(tR) < budgetPerRung
        improved = false; sweep = sweep + 1;
        for gi = 1:numel(gates)
            if toc(tR) > budgetPerRung; break; end
            f = ['pAtp_' gates{gi}]; base = q.(f);
            steps = [0.15 -0.15 0.30 -0.30];
            if sweep > 1; steps = [0.07 -0.07]; end
            bestC = c; bestV = base;
            for s = steps
                cand = q; cand.(f) = base + s;
                cc = cost(cand); nEval = nEval + 1;
                if cc < bestC - ACCEPT; bestC = cc; bestV = cand.(f); end
                if toc(tR) > budgetPerRung; break; end
            end
            if bestV ~= base
                fprintf('  [%5.1f min|%3d] %-10s %+.3f -> %+.3f   %.4f (%+.4f)\n', ...
                        toc(tAll)/60, nEval, gates{gi}, base, bestV, bestC, bestC-c);
                q.(f) = bestV; c = bestC; improved = true;
            end
        end
    end
    fprintf('  %s: %.4f -> %.4f  (%d sweeps, %.1f min)\n', H{k,1}, c0, c, sweep, toc(tR)/60);
    OUT(end+1) = struct('name',H{k,1},'dof',numel(gates),'c0',c0,'c',c,'p',q,'evals',nEval); %#ok<SAGROW>
    save(fullfile(resdir,'atp_ladder2.mat'), 'OUT', 'H');
end

%% ---- verdict --------------------------------------------------------------
fprintf('\n\n================ PATHWAY LADDER, REFITTED ================\n');
fprintf('%-38s %5s %10s %10s %10s\n','hypothesis','DOF','ablated','refit','per DOF');
for k = 1:numel(OUT)
    fprintf('%-38s %5d %10.3f %10.3f %10.3f\n', OUT(k).name, OUT(k).dof, ...
            OUT(k).c0, OUT(k).c, OUT(k).c/max(OUT(k).dof,1));
end
if numel(OUT) >= 4
    fprintf('\nmarginal cost bought by the INDEFENSIBLE gates:\n');
    fprintf('  H3 - H2 (3 mechanical DOF) : %+.3f\n', OUT(3).c - OUT(2).c);
    fprintf('  H4 - H2 (14 extra DOF)     : %+.3f\n', OUT(4).c - OUT(2).c);
    fprintf('  H4 - H3 (11 extra DOF)     : %+.3f\n', OUT(4).c - OUT(3).c);
end

%% ---- 8 mM invariance, once, for all rungs ---------------------------------
fprintf('\n8 mM invariance check:\n');
b8 = p; b8.MgATP = 8; b8.Pi = 0; b8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
c8ref = cost(b8);
for k = 1:numel(OUT)
    q8 = OUT(k).p; q8.MgATP = 8; q8.Pi = 0;
    q8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
    fprintf('  %-38s |d| = %.3e\n', OUT(k).name, abs(cost(q8) - c8ref));
end

fprintf('\ntotal %d evaluations in %.2f h\n', nEval, toc(tAll)/3600);
fprintf('LADDER STAGE 2 DONE\n');

%% ==========================================================================
function c = iCost(q, fn2, ATP_c, Data_ATP)
[~,~,fm,fd] = runSlackExperiment(q);
qf = q; qf.RunSlack = false;
[~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
f = fieldnames(fmv); for i=1:numel(f); fm.(f{i}) = fmv.(f{i}); end
f = fieldnames(fdv); for i=1:numel(f); fd.(f{i}) = fdv.(f{i}); end
c = sum(evalFeatureCost(fd, fm, fn2, 2));
end
