% CheckJointAtpBaseline.m
% Verify (a) the merged tree + new baseline params/rskR2_w025_opt.m evaluate here,
% (b) UseAtpJointFit is inert when OFF, (c) what the ABSOLUTE low-ATP cost is at the
% seed ATP levers, and which features carry it.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root));
refreshPool(6);
LoadData;

params0 = getParams(); run(fullfile(root,'params','rskR2_w025_opt.m'));
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity=true; params0.RunKtr=false; params0.RunSlack=true;
params0.RunStairs=false; params0.RunForceVelocityTime=false; params0.RunSlackPassive=true;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.OptimizeOn='Feats';
params0.RunSlackSegments='AllPar'; params0.MaxRunTime=600;
params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0;
params0.mods={}; params0.g=[];

% ---- (1) high-ATP baseline, joint flag OFF (regression) ----
[Et8, ~, fc8, fn8] = evaluateBakersExp(1, params0, true, false);
fprintf('\n=== HIGH-ATP baseline (rskR2_w025_opt), joint OFF: Et = %.3f ===\n', Et8);
[~,ord]=sort(fc8,'descend');
fprintf('  top terms: ');
for i=1:min(8,numel(ord)); fprintf('%s %.2f | ', fn8{ord(i)}, fc8(ord(i))); end; fprintf('\n');

% ---- (2) LOW-ATP absolute, slack only, ATP levers ON at seed ----
fn2 = params0.fn;
drop = startsWith(fn2,'FV') | startsWith(fn2,'PS_');
fn2 = fn2(~drop);
fn2 = fn2(~cellfun(@(s) any(strcmp(s, fn2(1:find(strcmp(fn2,s),1)-1)))  , fn2)); % de-duplicate
fprintf('\nlow-ATP feature list: %d terms (dropped %d FV/PS + dupes)\n', numel(fn2), sum(drop));

p2 = params0;
p2.UseAtpDetachR2=1; p2.K_T2=1.97; p2.UseAtpKmsrd=1; p2.K_srx=13; p2.MgATP_ref=8;
p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8;
p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
p2.RunForceVelocity=false; p2.RunForceVelocityTime=false; p2.RunSlackPassive=false;
p2.fn = fn2;
[Et2, ~, fc2, fn2n] = evaluateBakersExp(1, p2, true, false);
fprintf('\n=== LOW-ATP absolute (2 mM), levers K_T2=1.97 K_srx=13: Et = %.3f ===\n', Et2);
[~,o2]=sort(fc2,'descend');
fprintf('  per-feature (desc):\n');
for i=1:numel(o2)
    if fc2(o2(i))<1e-4; break; end
    fprintf('   %-22s %8.3f\n', fn2n{o2(i)}, fc2(o2(i)));
end

% ---- (3) same feature list at 8 mM, for an apples-to-apples floor ----
p8s = params0; p8s.fn = fn2;
p8s.RunForceVelocity=false; p8s.RunForceVelocityTime=false; p8s.RunSlackPassive=false;
p8s.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; p8s.MgATP=8; p8s.Pi=0;
[Et8s,~,fc8s,~] = evaluateBakersExp(1, p8s, true, false);
fprintf('\n=== SAME list at 8 mM (the floor set by the baseline): Et = %.3f ===\n', Et8s);
fprintf('  -> low-ATP absolute excess over the 8 mM floor: %.3f\n', Et2-Et8s);

% ---- (4) joint-fit wiring check ----
pj = params0; pj.UseAtpJointFit=true; pj.AtpJoint_w8=1; pj.AtpJoint_w2=1;
pj.AtpJoint_fn2=fn2; pj.AtpJoint_vt2='protocol_03_27_2026_2mM_slack.mat';
pj.AtpJoint_MgATP2=2; pj.AtpJoint_Pi2=0.8;
pj.UseAtpDetachR2=1; pj.K_T2=1.97; pj.UseAtpKmsrd=1; pj.K_srx=13;
Etj = evaluateBakersExp(1, pj, true, false);
fprintf('\n=== JOINT (w8*hi + w2*lo): Et = %.3f ===\n', Etj);
save(fullfile(root,'params','joint_baseline_check.mat'),'Et8','Et2','Et8s','Etj','fn2','fc2','fc8','fn8');
disp('CHECK DONE');
