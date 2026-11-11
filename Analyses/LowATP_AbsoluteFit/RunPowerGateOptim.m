% RunPowerGateOptim.m
% =========================================================================
% Fit the 2 mM data HARD, using 14 general ATP power gates.
%   rate_eff = rate * (MgATP/MgATP_ref)^p     (exactly 1 at 8 mM for ANY p)
% => the 8 mM fit is BIT-IDENTICAL by construction (verified max|d| = 0), so every
% exponent is a free low-ATP degree of freedom that costs the high-ATP fit nothing.
% At the only two concentrations we ever evaluate (2 and 8 mM) this form is strictly
% MORE general than the Michaelis gates it replaces: any K_T2/K_srx maps onto an
% equivalent exponent, so the previous best is inside the search space by construction.
%
% Objective: ABSOLUTE 2 mM slack feature cost, rsK down-weighted to 0.1
% (rsK is compliance-suspect, unreachable and numerically discontinuous -- see
% ../LowATP_RestretchRecovery/mechanism_hunt.md).
%
% Phase A  seed = power-gate equivalent of the current best (K_T2=1.80, K_srx=9.00)
% Phase B  coordinate line-search sweeps (robust, monotone)
% Phase C  Nelder-Mead polish
% State is saved on EVERY improvement, so the run can be inspected or killed anytime.
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

if ~exist('BUDGET_H','var'); BUDGET_H = 3.0; end
T0 = tic; BUDGET_S = BUDGET_H*3600;
STATE = fullfile(resdir,'powergate_state.mat');

P = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M'};
LB = -1.5*ones(1,numel(P));  UB = 1.5*ones(1,numel(P));

b=getParams(loadParams('rskR2_w025_opt'),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=300;
b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.UseAtpDetachR2=0; b.UseAtpKmsrd=0; b.UseStrainDepAtpGate=0;
b.UseR2Ceiling=0; b.UseD0State=0; b.MgATP_ref=8; b.UseAtpPowerGates=1;
b.MgATP=2; b.MgADP=0; b.Pi=0.8;
b.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';

fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn=cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput',false);
FN=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
b.fn=FN;
fprintf('objective: ABSOLUTE 2 mM, %d feature terms, rsK weight 0.1\n', numel(FN));

% ---- Phase A: seeds ----
% Michaelis gate K -> equivalent exponent at 2 mM:  p = ln(g)/ln(2/8),
%   g = [2/(2+K)] / [8/(8+K)]
mich2p = @(K) log( (2/(2+K)) / (8/(8+K)) ) / log(2/8);
x0 = zeros(1,numel(P));
x0(strcmp(P,'k2'))    = mich2p(1.80);
x0(strcmp(P,'kmsrd')) = mich2p(9.00);
fprintf('seed from current best: pAtp_k2 = %.4f, pAtp_kmsrd = %.4f\n', ...
    x0(strcmp(P,'k2')), x0(strcmp(P,'kmsrd')));

cZero = objPG(zeros(1,numel(P)), b, P, LB, UB);
cSeed = objPG(x0,                b, P, LB, UB);
fprintf('no ATP mechanism at all : %.4f\ncurrent-best equivalent  : %.4f\n\n', cZero, cSeed);

best.x = x0; best.c = cSeed; nEval = 2;
save(STATE,'best','P','cZero','cSeed','nEval');

% ---- Phase B: coordinate line-search sweeps ----
STEPS = {[-0.6 -0.3 -0.15 0.15 0.3 0.6], [-0.25 -0.12 0.12 0.25], [-0.1 -0.05 0.05 0.1]};
for sweep = 1:numel(STEPS)
    st = STEPS{sweep};
    fprintf('--- sweep %d (steps %s) ---\n', sweep, num2str(st,'%.2f '));
    improvedAny = false;
    for i = 1:numel(P)
        if toc(T0) > BUDGET_S; break; end
        bestLocal = best.c; bestVal = best.x(i);
        for d = st
            xt = best.x; xt(i) = best.x(i) + d;
            if xt(i) < LB(i) || xt(i) > UB(i); continue; end
            c = objPG(xt, b, P, LB, UB); nEval = nEval + 1;
            if c < bestLocal - 1e-4
                bestLocal = c; bestVal = xt(i);
            end
        end
        if bestLocal < best.c - 1e-4
            best.x(i) = bestVal; best.c = bestLocal; improvedAny = true;
            fprintf('  [%5.1f min | %4d ev] %-9s -> %+.3f   cost %.4f\n', ...
                toc(T0)/60, nEval, P{i}, bestVal, best.c);
            save(STATE,'best','P','cZero','cSeed','nEval');
        end
    end
    if toc(T0) > BUDGET_S; break; end
    if ~improvedAny; fprintf('  (no improvement this sweep)\n'); end
end
fprintf('\nafter coordinate sweeps: %.4f (%d evals, %.1f min)\n', best.c, nEval, toc(T0)/60);

% ---- Phase C: Nelder-Mead polish ----
if toc(T0) < BUDGET_S
    opts = optimset('Display','iter','MaxFunEvals',4000,'MaxIter',4000, ...
                    'TolFun',1e-4,'TolX',1e-4, ...
                    'OutputFcn', @(x,ov,st) toc(T0) > BUDGET_S);
    [xo, co] = fminsearch(@(x) objPG(x,b,P,LB,UB), best.x, opts);
    if co < best.c
        best.x = xo; best.c = co;
        save(STATE,'best','P','cZero','cSeed','nEval');
    end
end

fprintf('\n=== BEST 2 mM absolute cost: %.4f  (seed %.4f, no-mechanism %.4f) ===\n', ...
    best.c, cSeed, cZero);
for i=1:numel(P)
    if abs(best.x(i))>1e-3
        fprintf('  pAtp_%-9s %+7.3f   -> rate x%.3f at 2 mM\n', P{i}, best.x(i), (2/8)^best.x(i));
    end
end

% ---- report ----
q = applyP(b,P,best.x);
[~,~,fm,fd] = runSlackExperiment(q);
[cc,~,~] = evalFeatureCost(fd,fm,FN,2);
[~,o]=sort(cc,'descend');
fprintf('\nper-feature (2 mM absolute):\n');
for k=1:min(12,numel(o)); if cc(o(k))<1e-3; break; end
    fprintf('   %-24s %8.3f\n', FN{o(k)}, cc(o(k))); end
sh={'steady','A','Am','peak1_y','vall_y','peak2','ktr','rsK','t0'};
g=@(x)num2str(double(x(:))','%7.3g');
fprintf('\n%-10s %-40s | %s\n','feature','model 2 mM','data 2 mM');
for i=1:numel(sh); f=sh{i}; if ~isfield(fm,f); continue; end
    fprintf('%-10s %-40s | %s\n', f, g(fm.(f)), g(fd.(f))); end

% ---- verify 8 mM untouched, then bake ----
q8=q; q8.MgATP=8; q8.Pi=0; q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
r8=q8; r8.UseAtpPowerGates=0;
[~,~,h1,~]=runSlackExperiment(q8); [~,~,h0,~]=runSlackExperiment(r8);
fprintf('\n8 mM invariance: max|dsteady| = %.3g, max|drsK| = %.3g (MUST be 0)\n', ...
    max(abs(double(h1.steady(:))-double(h0.steady(:)))), ...
    max(abs(double(h1.rsK(:))-double(h0.rsK(:)))));

dst=fullfile(root,'params','lowATP_powergate_opt.m');
copyfile(fullfile(root,'params','rskR2_w025_opt.m'), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ==== Low-ATP fit via 8 mM-inert ATP power gates (%s) ====\r\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid,'%%%% rate_eff = rate*(MgATP/MgATP_ref)^p; = 1 at 8 mM for any p => 8 mM bit-identical.\r\n');
fprintf(fid,'%%%% 2 mM absolute feature cost %.4f (was %.4f with the Michaelis gates).\r\n', best.c, cSeed);
fprintf(fid,'params0.UseAtpPowerGates = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'params0.UseAtpDetachR2 = 0;\r\nparams0.UseAtpKmsrd = 0;\r\n');
for i=1:numel(P); fprintf(fid,'params0.pAtp_%s = %.6g;\r\n', P{i}, best.x(i)); end
fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n');
fclose(fid);
save(STATE,'best','P','cZero','cSeed','nEval','cc','FN');
fprintf('\nSaved %s\nPOWERGATE OPTIM DONE\n', dst);

% ================= local functions =================
function q = applyP(q,P,x)
    q.UseAtpPowerGates = 1;
    for i=1:numel(P); q.(['pAtp_' P{i}]) = x(i); end
end
function c = objPG(x, b, P, LB, UB)
    if any(x < LB) || any(x > UB); c = 1e6; return; end
    q = applyP(b,P,x);
    try
        [~,~,fm,fd] = runSlackExperiment(q);
        c = sum(evalFeatureCost(fd,fm,b.fn,2));
        if ~isfinite(c); c = 1e6; end
    catch
        c = 1e6;
    end
end
