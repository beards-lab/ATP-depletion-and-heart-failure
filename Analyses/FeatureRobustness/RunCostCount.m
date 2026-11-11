%% RunCostCount.m -- what each solver setting costs, measured in RHS calls (contention-proof).
% Wall time on a shared machine is not a measurement, so the cost is counted:
% every call of the RHS, split into SOLVER calls (ode15s, incl. its numerical
% Jacobians) and OUTPUT calls (evaluateModel's storeOutputs, one per stored
% point). Two runs per setting, both serial (no pool):
%   slack : runSlackExperiment, AllPar semantics; calls are binned into the
%           init hold and the five cycles, so the parallel CRITICAL PATH
%           (init + the slowest cycle, what a 5-worker campaign waits for) is known
%   full  : evaluateBakersExp, the whole objective; total calls and CPU time
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root)); res = fullfile(here, 'results');
fx = @(varargin) struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, varargin{:});
if ~exist('CFG', 'var')
CFG = { 'A old: 1e-2, step output, hard rectifier', struct('SmoothRectifierEps', 0, 'ClampStrainDepExtrapolation', 0, 'odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 0)
        'B SiLU 1e-3, step output',                 fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 0)
        'C SiLU 1e-3, dense in cycles only',        fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'D SiLU 1e-2, dense in cycles only',        fx('odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'E SiLU 1e-4, dense in cycles only',        fx('odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'F SiLU 1e-3, dense everywhere (campaign)', fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4) };
end
P = getParams(loadParams(PARAMSET), [], true, false);
P.MaxRunTime = 3000; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0; P.BreakOnODEUnstable = false;
P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
P.modelFcn = '@dPUdT_CountedTyped';
D = load(fullfile(root, 'data', P.velocitytableonfile)); vt = D.velocitytable;
cyc = find(vt(:,2) < -1)'; edges = [-inf, vt(cyc,1)', inf];     % init | cycle 1 .. 5
global RHS_T RHS_K RHS_N %#ok<GVMIS>
delete(gcp('nocreate'));
Q = struct();
for c = 1:size(CFG, 1)
    Pc = P; ov = CFG{c,2}; f = fieldnames(ov); for k = 1:numel(f), Pc.(f{k}) = ov.(f{k}); end
    Pc = getParams(Pc, Pc.g, true, false);
    % -- slack alone, binned by protocol window
    Ps = Pc; Ps.RunSlackSegments = 'AllPar';
    RHS_T = zeros(1e6,1); RHS_K = false(1e6,1); RHS_N = 0;
    c0 = cputime; t0 = tic; [~, ~, fm] = runSlackExperiment(Ps); ws = toc(t0); cs = cputime - c0;
    tt = RHS_T(1:RHS_N); kk = RHS_K(1:RHS_N);
    nS = histcounts(tt(~kk), edges); nO = histcounts(tt(kk), edges);
    crit = nS(1) + nO(1) + max(nS(2:end) + nO(2:end));
    % -- full objective, serial
    RHS_T = zeros(1e6,1); RHS_K = false(1e6,1); RHS_N = 0;
    c0 = cputime; t0 = tic; E = evaluateBakersExp(1, Pc, true, false); wf = toc(t0); cf = cputime - c0;
    kk = RHS_K(1:RHS_N);
    Q(c).name = CFG{c,1}; Q(c).slackSolver = nS; Q(c).slackOutput = nO; Q(c).crit = crit;
    Q(c).slackWall = ws; Q(c).slackCpu = cs; Q(c).fullSolver = nnz(~kk); Q(c).fullOutput = nnz(kk);
    Q(c).fullWall = wf; Q(c).fullCpu = cf; Q(c).cost = E;
    fprintf('%-44s slack: init %6d+%6d  worst cycle %6d  crit-path %7d calls | full: solver %8d output %8d  cpu %6.0f s  wall %6.0f s  cost %.4f\n', ...
        CFG{c,1}, nS(1), nO(1), max(nS(2:end) + nO(2:end)), crit, Q(c).fullSolver, Q(c).fullOutput, cf, wf, E);
    save(fullfile(res, ['costcount_' PARAMSET '.mat']), 'Q', 'CFG');
end
fprintf('DONE\n');
