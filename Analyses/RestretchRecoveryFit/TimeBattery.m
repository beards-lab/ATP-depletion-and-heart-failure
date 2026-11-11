% TimeBattery.m — where does one objective evaluation spend its time, and
% what MaxRunTime does this snapshot actually need?
%
% The first overnight launch returned a baseline cost of 1016441 (the
% 1e6 + 1e3*unsimulated-time ODE penalty) because cfg.MaxRunTime = 60 s was
% below what the slack init run needs at params/optfull7_opt_mov.m. An
% optimiser seeded from a failed baseline is worthless, so measure first.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if isempty(gcp('nocreate')); parpool('local', 5); end

p = getParams(loadParams('params/optfull7_opt_mov.m'), [], true, false);
p.MaxRunTime = 1200;
p.EvalFeatures = 1; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;
p.RunKtr = 0; p.RunStairs = 0; p.RunForceVelocityTime = 0;

modelFcn = str2func(p.modelFcn);
ds = load(fullfile('data', p.velocitytableonfile));
vt = ds.velocitytable; vt(1,1) = -20;

% --- 1. the slack INIT run: this is the one that timed out ---------------
pi_ = p; pi_.Velocity = 0;
if isfield(pi_,'PU0'); pi_ = rmfield(pi_,'PU0'); end
pi_ = getParams(pi_, pi_.g, true);
t1 = tic; [~, o_init] = evaluateModel(modelFcn, [-10 vt(3,1)], pi_); T_init = toc(t1);
fprintf('slack INIT run  ([-10 %.2f], %.1f s simulated) : %6.1f s wall\n', ...
        vt(3,1), vt(3,1)+10, T_init);

% --- 2. one slack CHUNK (5 run in parallel) ------------------------------
pc = p; pc.PU0 = o_init.PU(end,:); pc.Velocity = vt(2:6,2);
pc = getParams(pc, pc.g, true); pc.PU0 = o_init.PU(end,:);
t2 = tic; evaluateModel(modelFcn, vt(2:6,1), pc); T_chunk = toc(t2);
fprintf('one slack CHUNK                                : %6.1f s wall\n', T_chunk);

% --- 3. protocol blocks --------------------------------------------------
cfgs = {'FV only',      struct('RunForceVelocity',1,'RunSlack',0,'RunSlackPassive',0); ...
        'Slack only',   struct('RunForceVelocity',0,'RunSlack',1,'RunSlackPassive',0); ...
        'Passive only', struct('RunForceVelocity',0,'RunSlack',0,'RunSlackPassive',1); ...
        'FULL',         struct('RunForceVelocity',1,'RunSlack',1,'RunSlackPassive',1)};

ATP_c = p.MgATP; FV_velocities = p.FV_velocities; %#ok<NASGU>
% NOTE: RunBakersExp is a SCRIPT and reuses `i`/`j` internally, so the loop
% counter here must be a name it cannot clobber.
fprintf('\n%-14s %10s\n', 'block', 'wall (s)');
for i_cfg_ = 1:size(cfgs,1)
    nm_ = cfgs{i_cfg_,1};
    params0 = p;
    f_ = fieldnames(cfgs{i_cfg_,2});
    for j_ = 1:numel(f_); params0.(f_{j_}) = cfgs{i_cfg_,2}.(f_{j_}); end
    t3 = tic;
    try
        RunBakersExp;
        fprintf('%-14s %10.1f\n', nm_, toc(t3));
    catch ME
        fprintf('%-14s   FAILED after %.1f s: %s\n', nm_, toc(t3), ME.message);
    end
end

fprintf('\nRECOMMENDATION: MaxRunTime must exceed the INIT run (%.0f s) with margin.\n', T_init);
fprintf('   -> use MaxRunTime >= %d\n', max(120, ceil(T_init*3/10)*10));
