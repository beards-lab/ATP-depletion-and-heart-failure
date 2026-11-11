% RunPhase2.m
%
% Run Phase 2 of the slack identification: ViscoSE dashpot optimisation.
% Loads Phase 1 result from IdentifyLastSlack_state.mat and continues.
%
% Run this interactively in MATLAB to see live fminsearch progress:
%   >> RunPhase2
%
% Or re-run to continue from last checkpoint (150 evals per run).

cd('C:\home\git\ATP-depletion-and-heart-failure');
addpath(genpath('.'));

STATE_FILE  = 'IdentifyLastSlack_state.mat';
MAX_EVALS   = 300;   % increase for longer runs
fn_target   = {'ktr|SLslack','A|SLslack','peak1_y','peak1_dSL', ...
               'peak2','vall_y','vall2_dy','steady'};
mods_phase1 = {'kSE','mu','ka','kd','k1','k2','slope','k_pas'};

%% Load state
load(STATE_FILE,'state');
fprintf('Phase 1 best: FeatTot=%.4f\n', sum(state.costs1));
fprintf('Phase 1 g:    %s\n', num2str(state.g_opt1,'%.3f '));

%% Setup params
% clear params0; params0=getParams(); 
ModelParamsInitManualLastSlack; LoadData;
params0.RunSlackSegments='Last'; params0.RunForceVelocity=false;
params0.RunForceLengthEstim=false; params0.PlotEachSeparately=false;
params0.ShowResidualPlots=false; params0.ghostLoad='oj'; params0.ghostSave='';
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.fn=fn_target;
% params0.g=ones(1,numel(mods_phase1));

params0.justPlotStateTransitionsFlag = false;
params0.PlotEachSeparately = true;

params0.mods = [mods_phase1, {'c_SE_visc', 'kstiff2'}];

% Start from Phase1 result + c_SE_visc_g=1.0 (=6) + kstiff2_g=0.867 (=13000)
params0.g = [state.g_opt1, 1.0, 0.867];

params0.RunSlackSegments = 'AllPar';
% params0.RunForceVelocity = true;
RunBakersExp
%% Phase 2 setup: kstiff2 (−13%) + ViscoSE dashpot
% Sweep analysis showed:
%  - kstiff2=13000 (g=0.867) + c_SE_visc=6 reduces peak1 from 1.556→0.416
%  - s_threshold_R locked at 0.0046 (catastrophic cliff if lowered)
%  - Other params need re-tuning to restore pk2 (currently 0.921→target 0.745)
params0.UseViscoelasticSE = true;
params0.c_SE_visc         = 6;     % base; g multiplies
params0.kstiff2           = 15000; % base; g multiplies (g=0.867 → 13000)
params0.kstiff1           = '=kstiff2';
params0.mods = [mods_phase1, {'c_SE_visc', 'kstiff2'}];

% Start from Phase1 result + c_SE_visc_g=1.0 (=6) + kstiff2_g=0.867 (=13000)
g0 = [state.g_opt1, 1.0, 0.867];
g_opt2 = g0;

opts = optimset('Display','iter','TolFun',1e-4,'MaxFunEvals',MAX_EVALS,'MaxIter',MAX_EVALS);

fprintf('\n=== Phase 2: kstiff2+ViscoSE (10 params, max %d evals) ===\n', MAX_EVALS);
fprintf('Starting: c_SE_visc=%.1f  kstiff2=%.0f\n', ...
    params0.c_SE_visc*g0(end-1), params0.kstiff2*g0(end));
g_opt2 = fminsearch(@(g) slack_cost_p2(g, params0, fn_target), g0, opts);

%% Evaluate result
clf
p2 = params0; p2.g = g_opt2;
% p2.PlotEachSeparately = true;
% p2.c_SE_visc = 6;
% p2.kstiff2 = 15000*0.65;
% p2.kstiff1 = "=kstiff2";

[~,out,fm2,fd2] = runSlackExperiment(p2);
c2 = evalFeatureCost(fd2,fm2,fn_target,1);
plotFeatures(fd2, fm2, [], fn_target)

fprintf('\n%-20s  %8s  %8s  %8s\n','Feature','Baseline','Phase1','ViscoSE');
fprintf('%s\n',repmat('-',1,58));
for i=1:numel(fn_target)
    fprintf('%-20s  %8.3f  %8.3f  %8.3f\n', fn_target{i}, ...
        state.costs_base(i), state.costs1(i), c2(i));
end
fprintf('%s\n',repmat('-',1,58));
fprintf('%-20s  %8.3f  %8.3f  %8.3f\n','TOTAL', ...
    sum(state.costs_base), sum(state.costs1), sum(c2));
fprintf('\nc_SE_visc = %.2f  kstiff2 = %.0f\n', ...
    params0.c_SE_visc*g_opt2(end-1), params0.kstiff2*g_opt2(end));
fprintf('g: %s\n', num2str(g_opt2,'%.3f '));

%% Save
state.phase   = 2;
state.g_opt2  = g_opt2;
state.costs2  = c2;
state.best_feat = sum(c2);
save(STATE_FILE,'state');
fprintf('\nSaved to %s\n', STATE_FILE);

%% Local function
function E = slack_cost_p2(g, params0, fn_target)
    p = params0; p.g = g;
    try
        [~,~,fm,fd] = runSlackExperiment(p);
        E = sum(evalFeatureCost(fd,fm,fn_target,1));
    catch
        E = 1e6;
    end
end
