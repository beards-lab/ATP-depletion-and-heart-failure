% RunJointFitRsk010.m — best joint fit with rsK down-weighted to 0.1.
% Objective: Et = E(8 mM, absolute) + E(2 mM, absolute), rsK at weight 0.1 in BOTH legs.
% Free: shared {k2, kstiff2, ka, eta_M} + low-ATP-only gates {K_T2, K_srx}.
% eta_M is included because it is the CLEAN rsK level lever (peak1/steady untouched);
% k2/kstiff2/ka are the standard compensating trio.
% Bounded fminsearch in log-multiplier space; wall-clock stop via OutputFcn.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(root,'Analyses','LowATP_AbsoluteFit','results');
if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

WRSK     = 0.1;
BUDGET_S = 70*60;
t0       = tic;

base = getParams(loadParams('rskR2_w025_opt'), [], true, false);
base.RunForceVelocity=true; base.RunForceVelocityTime=false; base.RunKtr=false;
base.RunStairs=false; base.RunSlackPassive=true; base.RunSlack=true;
base.EvalFeatures=true; base.BreakOnODEUnstable=false; base.OptimizeOn='Feats';
base.RunSlackSegments='AllPar'; base.MaxRunTime=300;
base.PlotEachSeparately=0; base.PlotFeatureFitting=0;
base.UseAtpDetachR2=1; base.UseAtpKmsrd=1; base.MgATP_ref=8;
base.K_T2=1.80; base.K_srx=9.0; base.UseR2Ceiling=0; base.UseD0State=0;
base.MgATP=8; base.MgADP=0; base.Pi=0;
base.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';

% ---- feature lists with rsK re-weighted to WRSK ----
fn = base.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn = cellfun(@(s) reweight(s,'rsK',WRSK), fn, 'UniformOutput', false);
fn2 = fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
base.fn = fn;
base.UseAtpJointFit=true; base.AtpJoint_w8=1; base.AtpJoint_w2=1;
base.AtpJoint_fn2=fn2; base.AtpJoint_vt2='protocol_03_27_2026_2mM_slack.mat';
base.AtpJoint_MgATP2=2; base.AtpJoint_Pi2=0.8;
base.AtpJoint_RunFV2=false; base.AtpJoint_RunPS2=false;

P    = {'k2','kstiff2','ka','eta_M','K_T2','K_srx'};
LB   = [0.4 0.6 0.5 0.15 0.4 0.4];
UB   = [1.6 1.6 2.0 3.0  4.0 3.5];
seed = [base.k2 base.kstiff2 base.ka base.eta_M base.K_T2 base.K_srx];

objfun = @(x) jointCost(x, base, P, seed, LB, UB);
c0 = objfun(zeros(1,numel(P)));
fprintf('\nJOINT baseline (rsK w=%.2f) = %.4f\nfree: %s\n\n', WRSK, c0, strjoin(P,', '));

opts = optimset('Display','iter','MaxFunEvals',400,'MaxIter',300, ...
                'TolFun',1e-3,'TolX',1e-3, ...
                'OutputFcn', @(x,ov,st) toc(t0) > BUDGET_S);
[xopt, bestC] = fminsearch(objfun, zeros(1,numel(P)), opts);
bestX = exp(xopt);

fprintf('\n=== BEST joint (rsK w=%.2f): %.4f (from %.4f) in %.1f min ===\n', ...
    WRSK, bestC, c0, toc(t0)/60);
for k=1:numel(P)
    fprintf('  %-9s x%.3f  ->  %.4g\n', P{k}, bestX(k), seed(k)*bestX(k));
end

% ---- verify + per-leg report at the optimum ----
q = base; for k=1:numel(P); q.(P{k}) = seed(k)*bestX(k); end
q = getParams(q, [], true, false);
[Et,~,fc,fnm] = evaluateBakersExp(1, q, true, false);
fprintf('\nverify Et=%.4f\n', Et);
[~,o]=sort(fc,'descend');
for i=1:min(14,numel(o)); fprintf('   %-26s %8.3f\n', fnm{o(i)}, fc(o(i))); end

% ---- rsK / force outcome at both ATP ----
r=q; r.UseAtpJointFit=false; r.RunForceVelocity=false; r.RunSlackPassive=false;
a=r; a.MgATP=8; a.Pi=0; a.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
[~,~,m8,d8]=runSlackExperiment(a);
c=r; c.MgATP=2; c.Pi=0.8; c.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
[~,~,m2,d2]=runSlackExperiment(c);
f=@(x)mean(double(x),'omitnan');
fprintf('\n            rsK      steady   peak1    ktr\n');
fprintf('8mM model  %6.1f  %7.1f %7.1f %6.1f\n', f(m8.rsK),f(m8.steady),f(m8.peak1_y),f(m8.ktr));
fprintf('8mM DATA   %6.1f  %7.1f %7.1f %6.1f\n', f(d8.rsK),f(d8.steady),f(d8.peak1_y),f(d8.ktr));
fprintf('2mM model  %6.1f  %7.1f %7.1f %6.1f\n', f(m2.rsK),f(m2.steady),f(m2.peak1_y),f(m2.ktr));
fprintf('2mM DATA   %6.1f  %7.1f %7.1f %6.1f\n', f(d2.rsK),f(d2.steady),f(d2.peak1_y),f(d2.ktr));
fprintf('rsK ratio  model %.2f | data %.2f\n', f(m2.rsK)/f(m8.rsK), f(d2.rsK)/f(d8.rsK));

% ---- bake ----
dst=fullfile(root,'params','joint_rsk010_opt.m');
copyfile(fullfile(root,'params','rskR2_w025_opt.m'), dst);
fid=fopen(dst,'a');
fprintf(fid,'\r\n%%%% ==== JOINT 8+2 mM absolute fit, rsK weight %.2f (%s) ====\r\n', WRSK, datestr(now,'yyyy-mm-dd'));
for k=1:numel(P); fprintf(fid,'params0.%s = %.6g;\r\n', P{k}, seed(k)*bestX(k)); end
fprintf(fid,'params0.UseAtpDetachR2 = 1;\r\nparams0.UseAtpKmsrd = 1;\r\nparams0.MgATP_ref = 8;\r\n');
fprintf(fid,'%%%% joint cost %.4f (from %.4f); low-ATP condition: MgATP=2, Pi=0.8\r\n', bestC, c0);
fclose(fid);
save(fullfile(resdir,'joint_rsk010.mat'),'bestC','bestX','P','seed','c0','fc','fnm');
fprintf('\nSaved %s\nJOINT rsK010 DONE\n', dst);

% ================= local functions =================
function c = jointCost(x, base, P, seed, LB, UB)
    m = exp(x);
    if any(m < LB) || any(m > UB); c = 1e6; return; end
    q = base;
    for k = 1:numel(P); q.(P{k}) = seed(k)*m(k); end
    try
        q = getParams(q, [], true, false);
        c = evaluateBakersExp(1, q, true, false);
        if ~isfinite(c); c = 1e6; end
    catch
        c = 1e6;
    end
end

function s = reweight(s, feat, w)
    if startsWith(s, [feat '|']) || strcmp(s, feat); s = sprintf('%s|%g', feat, w); end
end
