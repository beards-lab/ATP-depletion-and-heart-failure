% RunHillGateConversion.m -- replace the power-law ATP gates with saturable (Hill)
% ones WITHOUT refitting, and prove the two forms are indistinguishable on the data.
%
% The incumbent low-ATP fit gates rates as       rate_eff = rate * (ATP/ATP_ref)^p.
% The requested classic form, renormalised so 8 mM stays untouched, is
%       rate_eff = rate * Corr_ref * ATP^n/(K^n + ATP^n),  Corr_ref = (K^n+ATP_ref^n)/ATP_ref^n.
%
% KEY FACT this script exploits: the data live at exactly TWO [ATP] (8 and 2 mM).
% Both forms are pinned to 1 at 8 mM by construction, so a gate is fully described by
% the single number g2 = factor at 2 mM. g2 is what the optimiser actually identified;
% p is only one parameterisation of it. Therefore K can be solved in CLOSED FORM
% (solveHillK) to reproduce g2 exactly -- no refit, no cost change, at either condition.
%
% What the script does:
%   1. reads the incumbent gates, converts p -> (n = +-1, K) and reports feasibility
%   2. reports the shared-Km alternative (one Km for all steps, Hill coefficient free)
%   3. writes params/lowATP_hillgate_opt.m
%   4. VERIFIES: analytic factor identity, 2 mM feature/cost identity, 8 mM invariance
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunHillGateConversion.m')"
% Knobs: SKIP_SIM = true  -> conversion + snapshot only (no MATLAB-side simulation)

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if ~exist('SKIP_SIM','var'); SKIP_SIM = false; end

BASE    = 'lowATP_powergate2_opt';  % incumbent power-gate fit (Analyses/LowATP_AbsoluteFit)
TAG     = 'lowATP_hillgate_opt';    % snapshot written by this script
ATP2    = 2;                        % the low-ATP condition (mM)
PTHRESH = 1e-3;                     % |p| below this is optimiser noise -> gate dropped
KSHARED = [1 2.5 5 10];             % shared-Km scenarios reported in section 2

GATED = {'k2','kmsrd','kah','kamh','ka','kd','k1','k_1', ...
         'ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M', ...
         'kstiff2','kSE','kSE_M','mu_neg','dr','dr2','P_bound_max'};
TIERA = 14;   % first 14 are rate constants; the rest are shape levers

b0   = getParams(loadParams(BASE), [], true, false);
Aref = b0.MgATP_ref;  r = ATP2/Aref;
fprintf('\n=== ATP gate re-parameterisation: power law -> saturable Hill ===\n');
fprintf('base snapshot : %s\n', BASE);
fprintf('ATP_ref = %g mM, low-ATP condition = %g mM  (ratio r = %.4g)\n\n', Aref, ATP2, r);

%% ---- 1. classic Michaelis conversion (Hill coefficient n = +-1) ---------------
G = struct('name',{},'p',{},'g2',{},'n',{},'K',{},'ok',{},'tier',{});
for i = 1:numel(GATED)
    p = b0.(['pAtp_' GATED{i}]);
    if abs(p) < PTHRESH, continue; end
    g2 = r^p;                            % the ONLY thing the 2-condition data identify
    [K, n, ok] = solveHillK(g2, ATP2, Aref);
    G(end+1) = struct('name',GATED{i},'p',p,'g2',g2,'n',n,'K',K,'ok',ok, ...
                      'tier', char('A' + (i > TIERA)));  %#ok<SAGROW>
end

fprintf('--- 1. per-gate conversion, classic Michaelis (n = +1, or -1 = ATP-inhibited) ---\n');
fprintf('%-14s %-4s %9s %9s %5s %11s   %s\n','rate','tier','p (power)','g2 (2 mM)','n','K (mM)','feasible?');
for i = 1:numel(G)
    if G(i).ok, fb = 'yes'; else, fb = sprintf('NO  (needs |n| > %.2f)', abs(log(G(i).g2)/log(r))); end
    fprintf('%-14s   %s  %+9.4f %9.4f %+5d %11.4g   %s\n', ...
        G(i).name, G(i).tier, G(i).p, G(i).g2, G(i).n, G(i).K, fb);
end
nBad = sum(~[G.ok]);
fprintf('\n%d of %d active gates convert to a classic Michaelis term; %d do not.\n', ...
    sum([G.ok]), numel(G), nBad);

% analytic identity check: does the Hill factor reproduce the power factor at 2 mM?
dmax = 0;
for i = 1:numel(G)
    if ~G(i).ok, continue; end
    dmax = max(dmax, abs(atpGateFactor(ATP2, Aref, G(i).K, G(i).n) - G(i).g2)/G(i).g2);
end
fprintf('max relative error of the converted factors at %g mM: %.3g (machine precision => exact)\n\n', ATP2, dmax);

%% ---- 2. shared-Km alternative (one Km for every ATP-binding step) -------------
% Physiologically the gated steps all bind the SAME ATP, so a single Km is the more
% honest model; the Hill coefficient then absorbs the per-step difference. Solve
% g2 = f(K_shared, n) for n. n is what a shared Km would have to claim about
% cooperativity -- values far from 1 are the cost of that parsimony.
fprintf('--- 2. shared-Km scenario: one Km for all steps, Hill coefficient solved per gate ---\n');
fprintf('%-14s %9s', 'rate', 'g2');
for k = KSHARED, fprintf('   n(K=%.3g)', k); end, fprintf('\n');
nShared = nan(numel(G), numel(KSHARED));
for i = 1:numel(G)
    fprintf('%-14s %9.4f', G(i).name, G(i).g2);
    for j = 1:numel(KSHARED)
        Kj = KSHARED(j);
        fn = @(n) log(atpGateFactor(ATP2, Aref, Kj, n)) - log(G(i).g2);
        try
            nShared(i,j) = fzero(fn, [-40 40]);
            fprintf('   %8.3f', nShared(i,j));
        catch
            fprintf('   %8s', 'n/a');
        end
    end
    fprintf('\n');
end
fprintf(['(a shared Km always has a solution, but note how the required coefficient blows up as\n' ...
         ' Km falls below the measured concentrations -- a saturated step cannot carry an ATP effect)\n\n']);

%% ---- 3. write the converted snapshot -----------------------------------------
dst = fullfile(root,'params',[TAG '.m']);
copyfile(fullfile(root,'params',[BASE '.m']), dst);
fid = fopen(dst,'a');
fprintf(fid,'\r\n%%%% ==== Saturable (Hill/Michaelis) ATP gates, converted from %s (%s) ====\r\n', BASE, datestr(now,'yyyy-mm-dd'));
fprintf(fid,'%%%% rate_eff = rate * Corr_ref * ATP^n/(KAtp^n + ATP^n), Corr_ref = (KAtp^n+ATP_ref^n)/ATP_ref^n\r\n');
fprintf(fid,'%%%% Km solved in closed form (Analyses/AtpGateHillForm) so the %g mM factor is IDENTICAL to\r\n', ATP2);
fprintf(fid,'%%%% the power-law fit it replaces; %g mM is untouched by construction. No refit involved.\r\n', Aref);
fprintf(fid,'params0.UseAtpPowerGates = 1;\r\nparams0.AtpGateForm = ''hill'';\r\nparams0.MgATP_ref = %g;\r\n', Aref);
for i = 1:numel(GATED)                       % start from a clean slate
    fprintf(fid,'params0.pAtp_%s = 0;\r\nparams0.KAtp_%s = Inf;\r\n', GATED{i}, GATED{i});
end
for i = 1:numel(G)
    if ~G(i).ok, continue; end
    fprintf(fid,'params0.pAtp_%s = %d;\r\nparams0.KAtp_%s = %.10g;\r\n', G(i).name, G(i).n, G(i).name, G(i).K);
end
for i = 1:numel(G)                           % anything that would not convert stays power-law
    if G(i).ok, continue; end
    fprintf(fid,'params0.pAtp_%s = %.10g;\r\nparams0.KAtp_%s = Inf;  %% not Michaelis-representable\r\n', ...
        G(i).name, G(i).p, G(i).name);
end
fprintf(fid,'params0.MgATP = %g;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n', ATP2);
fclose(fid);
fprintf('--- 3. wrote %s ---\n\n', dst);

save(fullfile(resdir,'conversion.mat'),'G','KSHARED','nShared','BASE','TAG','ATP2','Aref');
if SKIP_SIM, fprintf('SKIP_SIM -> stopping before the simulations.\n'); return; end

%% ---- 4. verification by simulation -------------------------------------------
if isempty(gcp('nocreate')); parpool('local',6); end

pw = iCfg(getParams(loadParams(BASE), [], true, false));   % incumbent, power form
hl = iCfg(getParams(loadParams(TAG),  [], true, false));   % converted,  hill  form
FN  = pw.fn;

fprintf('--- 4a. LOW ATP (%g mM): converted gates must reproduce the power-law fit exactly ---\n', ATP2);
[~,~,fm_pw,fd2] = runSlackExperiment(pw);
[~,~,fm_hl,~  ] = runSlackExperiment(hl);
c_pw = sum(evalFeatureCost(fd2,fm_pw,FN,2));
c_hl = sum(evalFeatureCost(fd2,fm_hl,FN,2));
fprintf('2 mM feature cost:  power %.6f   hill %.6f   |diff| %.3g\n', c_pw, c_hl, abs(c_pw-c_hl));
iCmp(fm_pw, fm_hl, 'power vs hill @ 2 mM');

fprintf('\n--- 4b. HIGH ATP (%g mM): the gate must be inert under BOTH forms ---\n', Aref);
p8 = pw; p8.MgATP = Aref; p8.Pi = 0; p8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
h8 = hl; h8.MgATP = Aref; h8.Pi = 0; h8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
o8 = p8; o8.UseAtpPowerGates = 0;                       % gates switched off entirely
[~,~,fm_p8,fd8] = runSlackExperiment(p8);
[~,~,fm_h8,~  ] = runSlackExperiment(h8);
[~,~,fm_o8,~  ] = runSlackExperiment(o8);
fprintf('8 mM feature cost:  gates off %.6f   power %.6f   hill %.6f\n', ...
    sum(evalFeatureCost(fd8,fm_o8,FN,2)), sum(evalFeatureCost(fd8,fm_p8,FN,2)), sum(evalFeatureCost(fd8,fm_h8,FN,2)));
iCmp(fm_o8, fm_h8, 'gates-off vs hill @ 8 mM (MUST be 0)');
iCmp(fm_o8, fm_p8, 'gates-off vs power @ 8 mM (MUST be 0)');

save(fullfile(resdir,'verification.mat'),'c_pw','c_hl','fm_pw','fm_hl','fm_p8','fm_h8','fm_o8','FN');
fprintf('\nHILLGATE CONVERSION DONE\n');

%% ---- helpers ------------------------------------------------------------------
function q = iCfg(q)
% same evaluation configuration the power-gate fit was scored under
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false; q.RunStairs=false;
q.RunSlackPassive=false;  q.RunSlack=true; q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=300;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0; q.UseR2Ceiling=0; q.UseD0State=0;
q.MgATP=2; q.MgADP=0; q.Pi=0.8;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
fn = q.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
fn = cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput', false);
q.fn = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
end

function iCmp(fa, fb, label)
% largest relative difference over every numeric feature the two runs share
fns = intersect(fieldnames(fa), fieldnames(fb));
worst = 0; wname = '-';
for i = 1:numel(fns)
    x = fa.(fns{i}); y = fb.(fns{i});
    if ~isnumeric(x) || ~isnumeric(y) || ~isequal(size(x),size(y)) || isempty(x), continue; end
    d = max(abs(double(x(:)) - double(y(:))) ./ max(abs(double(x(:))), 1e-12));
    if isfinite(d) && d > worst, worst = d; wname = fns{i}; end
end
fprintf('  %-42s max relative feature difference %.3g  (worst: %s)\n', label, worst, wname);
end
