% RunGateRefit.m -- refit the 8 mM-inert ATP gates on the realign4 baseline.
%
% Round 3. Round 2 (Analyses/LowATP_AbsoluteFit/RunPowerGateOptim2.m) fitted these same
% 21 exponents on the OLD baseline (rskR2_w025_opt) against the OLD feature list. The
% baseline has since moved to realign4_opt and the objective with it, so the exponents
% are seeded from there but have never been optimised here.
%
% Three changes from round 2, all forced by Analyses/AtpGateHillForm/conclusions.md 3a:
%
%   1. ACCEPT THRESHOLD 0.25, not 1e-4. The objective carries numerical scatter of
%      sigma ~0.076 / peak-to-peak 0.25 under parameter nudges of order 1e-9. Round 2's
%      1e-4 is three orders of magnitude below that, so its late fine-tuning was fitting
%      the solver. A move must now beat the seed by more than the scatter to be taken.
%   2. STEPS ARE COARSE for the same reason -- a step whose true effect is under ~0.1 is
%      not measurable, so there is no point spending evaluations on one.
%   3. SERIALISE AT %.17g. Round 2 wrote %.6g, which cost 0.60 of the objective all by
%      itself -- a saved optimum was NOT the optimum that was found.
%
% Staged, cheapest first, because a 2 mM slack evaluation costs ~100 s:
%   stage 1  the 9 gates that are actually non-zero in the seed, coarse steps
%   stage 2  the 12 near-zero gates -- do any of them earn a place on THIS baseline?
%   stage 3  refine stage-1 winners
%
% Every gate stays (MgATP/MgATP_ref)^p, so 8 mM is bit-identical for ANY exponent and
% none of this can touch the realign4 fit. That is asserted at the end, not assumed.
%
% Run: matlab -batch "BUDGET_H=3; run('Analyses/LowATP_Realign4Transfer/RunGateRefit.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',5); end

if ~exist('BUDGET_H','var'); BUDGET_H = 3.0; end
if ~exist('ACCEPT', 'var'); ACCEPT  = 0.25; end   % measured peak-to-peak of the scatter
T0 = tic; BUDGET_S = BUDGET_H*3600;
STATE = fullfile(resdir,'gaterefit_state.mat');

TIERA = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M'};
TIERB = {'kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};
P  = [TIERA TIERB];
LB = -1.5*ones(1,numel(P)); UB = 1.5*ones(1,numel(P));
ib = numel(TIERA)+1:numel(P); LB(ib) = -0.8; UB(ib) = 0.8;   % shape levers: tighter box

%% ---- base configuration: realign4 at 2 mM, slack only ----------------------
b = getParams(loadParams('lowATP_realign4_gates'), [], true, false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false;  b.RunSlack=true;  b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=600;
b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.UseAtpDetachR2=0; b.UseAtpKmsrd=0; b.UseStrainDepAtpGate=0;
b.UseR2Ceiling=0; b.UseD0State=0;
b.MgATP_ref=8; b.UseAtpPowerGates=1;
b.MgATP=2; b.MgADP=0; b.Pi=0.8;
b.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';

% This baseline's OWN feature list, slack subset (no 2 mM FV or passive data exists).
fn = b.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
prim = cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), fn, 'uni', 0);
FN = fn(~(startsWith(prim,'FV') | startsWith(prim,'PS_'))); b.fn = FN;
fprintf('objective: %d features, accept threshold %.3f\n%s\n\n', numel(FN), ACCEPT, strjoin(FN,' | '));

%% ---- the 8 mM floor on the identical list ----------------------------------
q8 = b; q8.MgATP=8; q8.Pi=0; q8.UseAtpPowerGates=0;
q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
[~,~,fm8,fd8] = runSlackExperiment(q8);
FLOOR = sum(evalFeatureCost(fd8,fm8,FN,2));
fprintf('8 mM floor on the same list: %.4f\n', FLOOR);
fprintf('  (of which %.3f is the k2/ksrd/kmsrd bound term a slack-only run cannot produce)\n', 1.2);

%% ---- seed = the transferred exponents, full precision ----------------------
x0 = zeros(1,numel(P));
for i = 1:numel(P)
    f = ['pAtp_' P{i}];
    if isfield(b,f); x0(i) = b.(f); end
end
seed.x = x0; seed.c = objGR(x0,b,P,LB,UB); nEval = 1;
best = seed;
fprintf('seed cost (transferred, unrefitted): %.4f\n\n', seed.c);
save(STATE,'best','seed','P','FLOOR','FN','nEval','ACCEPT');

ACTIVE = find(abs(x0) > 1e-3);
DUST   = find(abs(x0) <= 1e-3);
fprintf('stage 1: %d active gates (%s)\n', numel(ACTIVE), strjoin(P(ACTIVE),','));
fprintf('stage 2: %d near-zero gates\n\n', numel(DUST));

%% ---- staged coordinate descent ---------------------------------------------
STAGES = { ACTIVE, [-0.30 -0.15 0.15 0.30]
           DUST,   [-0.20 0.20]
           ACTIVE, [-0.08 0.08] };
LOG = struct('stage',{},'gate',{},'from',{},'to',{},'cost',{});
for s = 1:size(STAGES,1)
    idx = STAGES{s,1}; st = STAGES{s,2};
    fprintf('--- stage %d: %d gates x %d steps ---\n', s, numel(idx), numel(st));
    for i = idx(:)'
        if toc(T0) > BUDGET_S; fprintf('  (budget)\n'); break; end
        bl = best.c; bv = best.x(i);
        for d = st
            xt = best.x; xt(i) = best.x(i) + d;
            if xt(i) < LB(i) || xt(i) > UB(i); continue; end
            c = objGR(xt,b,P,LB,UB); nEval = nEval + 1;
            if c < bl - ACCEPT; bl = c; bv = xt(i); end
        end
        if bl < best.c - ACCEPT
            fprintf('  [%5.1f min|%4d] %-12s %+.3f -> %+.3f   cost %.4f (-%.3f)\n', ...
                    toc(T0)/60, nEval, P{i}, best.x(i), bv, bl, best.c-bl);
            LOG(end+1) = struct('stage',s,'gate',P{i},'from',best.x(i),'to',bv,'cost',bl); %#ok<SAGROW>
            best.x(i) = bv; best.c = bl;
            save(STATE,'best','seed','P','FLOOR','FN','nEval','ACCEPT','LOG');
        end
    end
    if toc(T0) > BUDGET_S; break; end
end
save(STATE,'best','seed','P','FLOOR','FN','nEval','ACCEPT','LOG');

%% ---- report -----------------------------------------------------------------
fprintf('\n=== REFIT: %.4f -> %.4f  (8 mM floor %.4f, %d evals, %.1f min) ===\n', ...
        seed.c, best.c, FLOOR, nEval, toc(T0)/60);
if best.c >= seed.c - ACCEPT
    fprintf('NO MOVE cleared the %.2f threshold -- the transferred gates are already at a\n', ACCEPT);
    fprintf('local optimum on this baseline, to within the objective''s own resolution.\n');
end
fprintf('\n%-14s %9s %9s   %s\n','gate','seed','refit','x at 2 mM (seed -> refit)');
for i = 1:numel(P)
    if abs(best.x(i)) > 1e-3 || abs(seed.x(i)) > 1e-3
        tier = 'A'; if i > numel(TIERA); tier = 'B'; end
        fprintf('[%s] pAtp_%-9s %+9.4f %+9.4f   x%.3f -> x%.3f\n', tier, P{i}, ...
                seed.x(i), best.x(i), (2/8)^seed.x(i), (2/8)^best.x(i));
    end
end

qb = applyGR(b,P,best.x);
[~,out,fm,fd] = runSlackExperiment(qb);
cS = evalFeatureCost(fd8,fm8,FN,2);           % 8 mM floor, per feature
cB = evalFeatureCost(fd,fm,FN,2);             % refit, per feature
qs = applyGR(b,P,seed.x); [~,~,fms,fds] = runSlackExperiment(qs);
cA = evalFeatureCost(fds,fms,FN,2);           % seed, per feature
fprintf('\nper-feature 2 mM cost:\n%-46s %9s %9s %9s\n','feature','seed','refit','8mM floor');
for k = 1:numel(FN)
    fprintf('%-46s %9.3f %9.3f %9.3f\n', FN{k}(1:min(46,end)), cA(k), cB(k), cS(k));
end
fprintf('%-46s %9.3f %9.3f %9.3f\n','TOTAL',sum(cA),sum(cB),sum(cS));

fprintf('\nATP effect (2 mM / 8 mM), refit vs data:\n');
fprintf('%-10s %10s %10s %8s | %8s %8s\n','feature','model 8mM','model 2mM','ratio','data','m/d');
for f = {'steady','A','peak1_y','vall_y','peak2','ktr','rsK'}
    n = f{1}; if ~isfield(fm,n) || ~isfield(fm8,n); continue; end
    a8 = mean(double(fm8.(n)(:)),'omitnan'); a2 = mean(double(fm.(n)(:)),'omitnan');
    b8 = mean(double(fd8.(n)(:)),'omitnan'); b2 = mean(double(fd.(n)(:)),'omitnan');
    fprintf('%-10s %10.4g %10.4g %8.3f | %8.3f %8.3f\n', n, a8, a2, a2/a8, b2/b8, (a2/a8)/(b2/b8));
end

%% ---- 8 mM invariance MUST be exactly zero -----------------------------------
g8 = qb; g8.MgATP=8; g8.Pi=0; g8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
n8 = g8; n8.UseAtpPowerGates=0;
[~,~,h1,~] = runSlackExperiment(g8); [~,~,h0,~] = runSlackExperiment(n8);
d1 = max(abs(double(h1.steady(:))-double(h0.steady(:))));
d2 = max(abs(double(h1.rsK(:))   -double(h0.rsK(:))));
fprintf('\n8 mM invariance: max|dsteady| %.3g  max|drsK| %.3g   (MUST be 0)\n', d1, d2);
assert(d1==0 && d2==0, 'gates leaked into the 8 mM condition');

%% ---- write the snapshot at FULL precision ------------------------------------
if best.c < seed.c - ACCEPT
    dst = fullfile(root,'params','lowATP_realign4_gates_opt.m');
    copyfile(fullfile(root,'params','realign4_opt.m'), dst);
    fid = fopen(dst,'a');
    fprintf(fid,'\n%%%% ==== Low-ATP gate refit on the realign4 baseline (%s) ====\n', datestr(now,'yyyy-mm-dd'));
    fprintf(fid,'%%%% 2 mM cost %.6f (seed %.6f, 8 mM floor %.6f on the same %d-term list).\n', best.c, seed.c, FLOOR, numel(FN));
    fprintf(fid,'%%%% Accept threshold %.2f = the objective''s measured peak-to-peak scatter.\n', ACCEPT);
    fprintf(fid,'%%%% Serialised at %%.17g: %%.6g rounding alone is worth ~0.60 of this cost.\n');
    fprintf(fid,'params0.UseAtpPowerGates = 1;\nparams0.MgATP_ref = 8;\n');
    fprintf(fid,'params0.UseAtpDetachR2 = 0;\nparams0.UseAtpKmsrd = 0;\n');
    for i=1:numel(P); fprintf(fid,'params0.pAtp_%s = %.17g;\n', P{i}, best.x(i)); end
    fprintf(fid,'params0.MgATP = 2;\nparams0.MgADP = 0;\nparams0.Pi = 0.8;\n');
    fprintf(fid,'params0.velocitytableonfile = ''protocol_03_27_2026_2mM_slack.mat'';\n');
    fclose(fid);
    fprintf('\nwrote %s\n', dst);
else
    fprintf('\nno snapshot written -- nothing cleared the threshold.\n');
end
save(STATE,'best','seed','P','FLOOR','FN','nEval','ACCEPT','LOG','cA','cB','cS');
fprintf('\nGATE REFIT DONE\n');

%% =============================================================================
function q = applyGR(q,P,x)
    q.UseAtpPowerGates = 1;
    for i = 1:numel(P); q.(['pAtp_' P{i}]) = x(i); end
end
function c = objGR(x,b,P,LB,UB)
    if any(x<LB) || any(x>UB); c = 1e6; return; end
    q = applyGR(b,P,x);
    try
        [~,~,fm,fd] = runSlackExperiment(q);
        c = sum(evalFeatureCost(fd,fm,b.fn,2));
        if ~isfinite(c); c = 1e6; end
    catch; c = 1e6; end
end
