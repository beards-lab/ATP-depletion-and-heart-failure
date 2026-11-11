% RunOtherOptions.m -- options for a better low-ATP fit BEYOND power-law gates.
%
% Constraint that eliminates most candidates immediately: the framework relies on
% the 8 mM fit being invariant, which requires every ATP mechanism to be exactly
% 1 at MgATP_ref. That rules out UseD0State (apo/rigor pool) and UseR2Ceiling --
% both change 8 mM when switched on, so adopting either costs a full 8 mM refit
% and forfeits the clean separation. They are also already falsified for the
% 8 mM restretch problem (Analyses/ApoPoolDetachment), though on the old core.
%
% What remains, all ref-normalised and therefore 8 mM-invariant:
%
%   O1  MICHAELIS instead of power law on the detachment step.
%       UseAtpDetachR2: R2 *= [ATP/(ATP+K_T2)] / [8/(8+K_T2)]
%       One parameter (K_T2) with a physical meaning -- an apparent Km -- rather
%       than an exponent. The interesting outcome is not the cost but the FITTED
%       Km: if it lands near the solution value for ATP-induced actomyosin
%       dissociation the mechanism is corroborated; if it needs several mM it is
%       falsified the way kmsrd was in AtpGateHillForm.
%
%   O2  STRAIN-DEPENDENT Michaelis (UseStrainDepAtpGate, nested in O1):
%       Km_eff(s) = K_T2 * (1 + kAtpStrain * max(0, s - sAtpStrain0))
%       The only mechanism in the model that can slow a HIGH-STRAIN process more
%       than an isometric one, i.e. separate rsK from ktr.
%
%   O3  SATURABLE (Hill) gates instead of unbounded power law, on the best
%       plausible set. Same 2 mM factor by construction, so the cost is nearly
%       unchanged -- the point is the implied chemistry and the extrapolation.
%
%   O4  IMPLIED-Km AUDIT (no simulation). For every fitted exponent p, the 2 mM
%       factor is f = 0.25^p, and the Michaelis Km that reproduces the same
%       factor is K = 16(1-f)/(8f-2). A gate needing an implausible Km is
%       phenomenology wearing a mechanism's name.
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunOtherOptions.m')"

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));

if ~exist('CORE','var');       CORE       = 'spreadw_iter_018'; end
if ~exist('BUDGET_MIN','var'); BUDGET_MIN = 10;   end
if ~exist('ACCEPT','var');     ACCEPT     = 0.32; end
if ~exist('NWORKERS','var');   NWORKERS   = 5;    end
resDir = fullfile(here,'results');

%% ---- O4 first: costs nothing, and it frames the rest ---------------------
PS = fullfile(resDir,'plaus_state.mat');
fprintf('=== O4  implied Michaelis Km for every fitted gate ===\n');
fprintf('   f = 0.25^p is the 2 mM factor; K = 16(1-f)/(8f-2) reproduces it\n');
fprintf('   as a Michaelis gate. Only valid for f > 0.25 (K > 0).\n\n');
fprintf('%-6s %-12s %8s %8s %10s\n','set','gate','p','f(2mM)','Km (mM)');
if isfile(PS)
    L = load(PS);
    for k = 1:numel(L.OUT)
        o = L.OUT(k); f = fieldnames(o.exps);
        for i = 1:numel(f)
            p = o.exps.(f{i}); if abs(p) < 1e-9; continue; end
            fa = 0.25^p; Km = 16*(1-fa)/(8*fa-2);
            if fa <= 0.25 || Km < 0
                fprintf('%-6s %-12s %8.3f %8.3f %10s\n', o.id, f{i}, p, fa, 'n/a');
            else
                fprintf('%-6s %-12s %8.3f %8.3f %10.2f\n', o.id, f{i}, p, fa, Km);
            end
        end
    end
else
    fprintf('   (no plaus_state.mat yet -- run RunPlausibleSets.m first)\n');
end
fprintf(['\n   Reference: apparent Km(MgATP) for ATP-induced actomyosin\n' ...
         '   dissociation in fibres is sub-millimolar. A gate needing several\n' ...
         '   mM is not simple ATP binding.\n\n']);

%% ---- setup for the simulated options -------------------------------------
LoadData;
if exist('parpool','file')==2 && isempty(gcp('nocreate')); parpool('local',NWORKERS); end

p0 = getParams(loadParams(CORE), [], true, false);
p0.UseAtpJointFit=false; p0.UseAtpPowerGates=true; p0.MgATP_ref=8;
p0.UseAtpDetachR2=0; p0.UseAtpKmsrd=0; p0.UseStrainDepAtpGate=0;
p0.UseR2Ceiling=0; p0.UseD0State=0;
p0.EvalFeatures=true; p0.BreakOnODEUnstable=false; p0.OptimizeOn='Feats';
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=400;
p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0; p0.PlotFullscreen=0;
p0.RunKtr=false; p0.RunStairs=false; p0.RunForceVelocityTime=false;
p0.RunSlackPassive=false; p0.RunSlack=true; p0.RunForceVelocity=true;
p0.mods={}; p0.g=ones(1,0);

prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn2 = p0.fn; [~,iu]=unique(fn2,'stable'); fn2=fn2(iu);
fn2 = fn2(~startsWith(prim(fn2),'PS_'));
GALL = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr', ...
        'tau_reg','A0','eta_M','kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};

p2 = p0; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8;
p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; p2.fn=fn2;
for i=1:numel(GALL); p2.(['pAtp_' GALL{i}])=0; end
rep = checkFeatureCoverage(p2);
if ~rep.ok
    bad = unique([rep.missingModel(:);rep.missingData(:)])';
    fn2 = fn2(~ismember(prim(fn2),bad)); p2.fn=fn2;
end
cost = @(q) iCost(q, fn2, ATP_c, Data_ATP);

% best plausible set so far, as the reference point
best = [];
if isfile(PS)
    L = load(PS); [~,bi] = min([L.OUT.c]); best = L.OUT(bi);
    fprintf('reference: %s = %.4f (%s)\n\n', best.id, best.c, strjoin(best.gates,','));
end

RES = struct('id',{},'what',{},'dof',{},'c',{},'note',{});

%% ---- O1: Michaelis on detachment, K_T2 fitted -----------------------------
fprintf('=== O1  Michaelis detachment gate, K_T2 fitted ===\n');
qM = p2; qM.UseAtpDetachR2 = 1; qM.pAtp_k2 = 0;   % gate replaces the power law
bestK = NaN; bestC = Inf;
for K = [0.1 0.3 0.6 1.0 1.6 2.5 4.0 6.0 10.0]
    q = qM; q.K_T2 = K; c = cost(q);
    fprintf('   K_T2 = %5.2f mM -> %.4f\n', K, c);
    if c < bestC; bestC = c; bestK = K; end
end
fprintf('   best K_T2 = %.2f mM -> %.4f\n\n', bestK, bestC);
RES(end+1) = struct('id','O1','what','Michaelis on k2 (K_T2 fitted)','dof',1, ...
                    'c',bestC,'note',sprintf('K_T2 = %.2f mM', bestK));

%% ---- O2: strain-dependent Michaelis ---------------------------------------
fprintf('=== O2  strain-dependent Km (separates rsK from ktr) ===\n');
b2 = struct('c',Inf,'ka',NaN,'s0',NaN);
for ka = [0 25 50 100 200]
    for s0 = [0 0.005 0.010]
        if ka==0 && s0~=0; continue; end
        q = qM; q.K_T2 = bestK;
        q.UseStrainDepAtpGate = ka>0; q.kAtpStrain = ka; q.sAtpStrain0 = s0;
        c = cost(q);
        fprintf('   kAtpStrain %5g  s0 %.3f -> %.4f\n', ka, s0, c);
        if c < b2.c; b2 = struct('c',c,'ka',ka,'s0',s0); end
    end
end
fprintf('   best: kAtpStrain %g, s0 %.3f -> %.4f\n\n', b2.ka, b2.s0, b2.c);
RES(end+1) = struct('id','O2','what','+ strain-dependent Km','dof',3, ...
                    'c',b2.c,'note',sprintf('kAtpStrain %g, s0 %.3f', b2.ka, b2.s0));

%% ---- O3: saturable (Hill) form of the best plausible set ------------------
if ~isempty(best)
    fprintf('=== O3  Hill (saturable) form of %s ===\n', best.id);
    qH = p2; f = fieldnames(best.exps);
    for i=1:numel(f); qH.(['pAtp_' f{i}]) = best.exps.(f{i}); end
    cP = cost(qH);
    qH.AtpGateForm = 'hill';
    bH = struct('c',Inf,'K',NaN);
    for K = [1 2 4 8 16 32]
        q = qH; for i=1:numel(f); q.(['KAtp_' f{i}]) = K; end
        c = cost(q);
        fprintf('   KAtp = %4g mM -> %.4f\n', K, c);
        if c < bH.c; bH = struct('c',c,'K',K); end
    end
    fprintf('   power law %.4f | best Hill %.4f at KAtp %g\n\n', cP, bH.c, bH.K);
    RES(end+1) = struct('id','O3','what',sprintf('Hill form of %s',best.id), ...
                        'dof',best.dof+1,'c',bH.c,'note',sprintf('KAtp %g mM',bH.K));
end

%% ---- summary --------------------------------------------------------------
fprintf('\n============ OTHER OPTIONS, RANKED ============\n');
if ~isempty(best)
    fprintf('%-4s %-34s %4s %10s  %s\n','id','mechanism','DOF','cost','note');
    fprintf('%-4s %-34s %4d %10.4f  %s\n','ref', ...
        sprintf('best power-law set (%s)',best.id), best.dof, best.c, strjoin(best.gates,','));
end
[~,o] = sort([RES.c]);
for r = o(:)'
    fprintf('%-4s %-34s %4d %10.4f  %s\n', RES(r).id, RES(r).what, RES(r).dof, RES(r).c, RES(r).note);
end
save(fullfile(resDir,'other_options.mat'),'RES','bestK','b2');
fprintf('\nOTHER OPTIONS DONE\n');

%% ==========================================================================
function c = iCost(q, fn2, ATP_c, Data_ATP)
try
    [~,~,fm,fd] = runSlackExperiment(q);
    qf=q; qf.RunSlack=false;
    [~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
    f=fieldnames(fmv); for i=1:numel(f); fm.(f{i})=fmv.(f{i}); end
    f=fieldnames(fdv); for i=1:numel(f); fd.(f{i})=fdv.(f{i}); end
    c = sum(evalFeatureCost(fd, fm, fn2, 2));
catch ME
    fprintf(2,'   [eval failed: %s]\n', ME.message); c = 1e6;
end
end
