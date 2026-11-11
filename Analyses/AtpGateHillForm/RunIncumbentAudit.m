% RunIncumbentAudit.m -- reconcile the incumbent 2 mM cost with what its header claims.
%
% params/lowATP_powergate2_opt.m records "2 mM absolute cost 15.9849" (2026-08-10), but
% re-evaluating the same snapshot today gives ~16.59. Nothing in Model/ or Auxiliary/
% has been committed since that fit and the 8 mM floor reproduces exactly (4.685491 vs
% the recorded 4.6855), so this is not a code regression -- but the number matters,
% because anyone comparing the Hill conversion (16.52) against the header (15.98) would
% read a real improvement as a regression, or vice versa.
%
% Splits the gap into:
%   (a) SERIALISATION -- the snapshot stores the exponents at %.6g, the optimiser held
%       them at full precision. Re-evaluate the exact saved best.x.
%   (b) DRIFT -- anything left after (a) is the evaluation itself having moved.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunIncumbentAudit.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',6); end

STATE = fullfile(root,'Analyses','LowATP_AbsoluteFit','results','powergate2_state.mat');
BASE  = 'lowATP_powergate2_opt';
S = load(STATE);
fprintf('\nsaved optimiser state: best.c = %.6f  (%d gates, %d evals)\n', S.best.c, numel(S.P), S.nEval);
if isfield(S,'FLOOR'); fprintf('saved 8 mM floor: %.6f\n', S.FLOOR); end

% (1) the snapshot as committed (exponents rounded to %.6g)
q = iCfg(getParams(loadParams(BASE), [], true, false));
[~,~,fm,fd] = runSlackExperiment(q);
c_snap = sum(evalFeatureCost(fd,fm,q.fn,2));

% (2) the optimiser's own vector, at full precision, on the same base it used
b = iCfg(getParams(loadParams('rskR2_w025_opt'), [], true, false));
b.UseAtpPowerGates = 1;
for i = 1:numel(S.P); b.(['pAtp_' S.P{i}]) = S.best.x(i); end
[~,~,fmX,fdX] = runSlackExperiment(b);
c_exact = sum(evalFeatureCost(fdX,fmX,b.fn,2));

fprintf('\n%-42s %12s\n','configuration','2 mM cost');
fprintf('%-42s %12.6f\n','optimiser''s recorded best.c', S.best.c);
fprintf('%-42s %12.6f\n','best.x re-evaluated at full precision', c_exact);
fprintf('%-42s %12.6f\n','the committed snapshot (%.6g rounding)', c_snap);
fprintf('\nsplit of the gap:\n');
fprintf('  drift since the fit (best.c -> re-evaluated) : %+.4f\n', c_exact - S.best.c);
fprintf('  serialisation rounding (exact -> snapshot)   : %+.4f\n', c_snap - c_exact);
fprintf('  total                                        : %+.4f\n', c_snap - S.best.c);

% how far apart are the two parameter vectors, really?
d = 0;
for i = 1:numel(S.P)
    if isfield(q, ['pAtp_' S.P{i}])
        d = max(d, abs(q.(['pAtp_' S.P{i}]) - S.best.x(i)));
    end
end
fprintf('\nlargest |p_snapshot - p_best|: %.3g  (pure %%.6g serialisation error)\n', d);
fprintf('\nINCUMBENT AUDIT DONE\n');

function q = iCfg(q)
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
