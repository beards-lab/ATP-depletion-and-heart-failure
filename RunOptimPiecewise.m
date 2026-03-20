% RunOptimPiecewise.m
%
% Identification using PieceWise strain-dependent rate spline knots directly.
%
% Key insight from sweep analysis:
%   - Phase1 reduced k2 to 0.52x, which halved R2 at ALL strains including
%     the high-strain detachment knots that drive the valley.
%   - Instead of a single k2 scalar, optimise the PieceWise knot VALUES
%     independently — this decouples isometric R2 (affects steady force)
%     from high-strain R2 (drives valley depth).
%   - UseCatchBond = false (piecewise already encodes all strain dependence)
%
% Free parameters (14 total):
%   Kinetics:  kSE, mu, kd, k1, k_pas  (5)
%   R2 spline: PieceWiseStrainDep2Params indices 6,7,8,9  (4 knots, s=0..0.023)
%   R12 spline: PieceWiseStrainDepParams indices 4,5  (2 knots, s=0.0046..0.02)
%   R1D spline: PieceWiseStrainDepR1DParams indices 5,6  (2 knots, s=0.004..0.01)
%   kstiff2, c_SE_visc  (2)
%
% NOTE: ka and k2 are left FREE only as global rate scalars to prevent
%       degeneracy with the spline knots.

cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));

STATE_FILE  = 'RunOptimPiecewise_state.mat';
MAX_EVALS   = 300;

fn_target = {'ktr|SLslack','A|SLslack','peak1_y','peak1_dSL', ...
             'peak2','vall_y','vall2_dy','steady'};

fn_target = {'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart', 'vall2_dy'};

%% Setup
clear params0; params0=getParams(); ModelParamsInitManualLastSlack; LoadData;
params0.RunSlackSegments='Last'; params0.RunForceVelocity=false;
params0.RunForceLengthEstim=false; params0.PlotEachSeparately=false;
params0.ShowResidualPlots=false; params0.ghostLoad=''; params0.ghostSave='';
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.fn = fn_target;
params0.UseCatchBond = false;   % piecewise splines handle all strain dependence

% ViscoSE: c_SE_visc=6 base (helps peak1 overshoot)
params0.UseViscoelasticSE = true;
params0.c_SE_visc         = 6;
params0.kstiff2           = 15000;
params0.kstiff1           = '=kstiff2';

params0.RunSlackSegments = 'Fourth'
% Define mods: mix of global rate params + individual spline knots
% Naming convention for array elements: 'FieldName__index'
mods = { ...
    'kSE', 'mu', 'kd', 'k1', 'k_pas', ...      % 5 mechanical/kinetic
    'c_SE_visc', 'kstiff2', ...                  % 2 new mechanisms
    'PieceWiseStrainDep2Params__6',  ...  % R2 at s=0       (val=0.5)
    'PieceWiseStrainDep2Params__7',  ...  % R2 at s=0.006   (val=1.0)
    'PieceWiseStrainDep2Params__8',  ...  % R2 at s=0.010   (val=20)
    'PieceWiseStrainDep2Params__9',  ...  % R2 at s=0.0233  (val=50)
    'PieceWiseStrainDepParams__4',   ...  % R12 at s=0.0046 (val=1.155)
    'PieceWiseStrainDepParams__5',   ...  % R12 at s=0.02   (val=0.015)
    'PieceWiseStrainDepR1DParams__5',...  % R1D at s=0.004  (val=8)
    'PieceWiseStrainDepR1DParams__6' ...  % R1D at s=0.01   (val=120)
};
params0.mods = mods;
Np = numel(mods);

%% Load or init state
if exist(STATE_FILE,'file')
    load(STATE_FILE,'optstate');
    fprintf('Resuming: evals=%d, best=%.4f\n', optstate.total_evals, optstate.best_cost);
    g0 = optstate.g_opt;
else
    % Start from ones (all at baseline values) except:
    % kstiff2: start at 0.867 (=13000/15000, best from sweep)
    % c_SE_visc: start at 1.0 (base=6)
    g0 = ones(1, Np);
    g0(strcmp(mods,'kstiff2'))  = 0.867;
    optstate = struct('total_evals',0,'best_cost',Inf,'g_opt',g0,'costs',[]);
end

opts = optimset('Display','iter','TolFun',1e-4,'MaxFunEvals',MAX_EVALS,'MaxIter',MAX_EVALS);
%% Init and test
p = params0;
p.PlotEachSeparately = true;
tic
init_costs = pw_cost(g0, p, fn_target)
toc

%% Print current state
fprintf('\n=== RunOptimPiecewise: %d free params, max %d evals ===\n', Np, MAX_EVALS);
fprintf('Params: %s\n\n', strjoin(mods,', '));

%% Run optimisation
g_new = fminsearch(@(g) pw_cost(g, params0, fn_target), g0, opts);
% writeParamsToMFile('ModelOptParamsLastSlack_PieceWiseFunc', p_opt)
%% Evaluate result
figure(232);clf;
p_opt = params0; p_opt.g = g_new;
p_opt.PlotEachSeparately = true;
p_opt.RunSlackSegments = 'AllPar';
p_opt.RunForceVelocity = true;
params0 = p_opt;

RunBakersExp;
% [~,~,features_model,features_data] = runSlackExperiment(p_opt);

c_new = evalFeatureCost(features_data, features_model, fn_target, 2);
optstate.total_evals = optstate.total_evals + MAX_EVALS;
optstate.costs = c_new;

fprintf('\n%-20s  %8s  %8s\n','Feature','Baseline','Optimised');
fprintf('%s\n',repmat('-',1,42));
for i=1:numel(fn_target)
    fprintf('%-20s  %8.3f  %8.3f\n', fn_target{i}, ...
        optstate.costs(i) * (numel(optstate.costs)>=i), c_new(i));
end
fprintf('%s\n',repmat('-',1,42));
fprintf('%-20s  %8.3f  %8.3f\n','TOTAL',sum(optstate.costs),sum(c_new));

% Print resolved spline knot values
p_res = getParams(p_opt, g_new, false, true);
fprintf('\nResolved R2 spline (PieceWiseStrainDep2Params at positive strain):\n');
fprintf('  s=0:      %.3f (g=%.3f)\n', p_res.PieceWiseStrainDep2Params(6),  g_new(strcmp(mods,'PieceWiseStrainDep2Params__6')));
fprintf('  s=0.006:  %.3f (g=%.3f)\n', p_res.PieceWiseStrainDep2Params(7),  g_new(strcmp(mods,'PieceWiseStrainDep2Params__7')));
fprintf('  s=0.010:  %.3f (g=%.3f)\n', p_res.PieceWiseStrainDep2Params(8),  g_new(strcmp(mods,'PieceWiseStrainDep2Params__8')));
fprintf('  s=0.023:  %.3f (g=%.3f)\n', p_res.PieceWiseStrainDep2Params(9),  g_new(strcmp(mods,'PieceWiseStrainDep2Params__9')));
fprintf('  kstiff2 = %.0f  c_SE_visc = %.2f\n', p_res.kstiff2, p_res.c_SE_visc);
%%
fn = {'FV_f|FV_v', 'ktr|SLslack', 'A|SLslack', 't0|SLslack', 'peak1_y', 'peak1_dSL', 'peak2', 'steady', 'XTOR|0.1', 'vall_y', 'restretchSlopeStart'};
figure(233);plotFeatures(features_data, features_model, [], fn)

%% Save
optstate.g_opt      = g_new;
optstate.best_cost  = sum(c_new);
optstate.costs      = c_new;
save(STATE_FILE,'optstate');
fprintf('\nSaved to %s\n', STATE_FILE);

%% Local function
function E = pw_cost(g, params0, fn_target)
    p = params0; p.g = g;
    try
        [~,~,fm,fd] = runSlackExperiment(p);
        E = sum(evalFeatureCost(fd, fm, fn_target, 1));
    catch
        E = 1e6;
    end
end
