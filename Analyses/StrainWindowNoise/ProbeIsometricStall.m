%% ProbeIsometricStall.m -- which non-smooth term makes ode15s grind at isometric rest?
%
% Isometric only (one zero-velocity segment, 0..T s from the getParams initial
% state), RHS calls counted with dPUdT_Counted, at a tight tolerance. Each row of
% VARS switches one suspect off or smooths it. Diagnostic: the variants change
% the model, the question is only which one removes the stall.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
if ~exist('TEND', 'var'), TEND = 5; end
if ~exist('BUDGET', 'var'), BUDGET = 90; end
if ~exist('TOL', 'var'), TOL = [1e-3 1e-6]; end
if ~exist('VARS', 'var')
    VARS = { 'as fitted',                 struct()
             'mu_neg = mu',               struct('mu_neg', NaN)
             'no A_reg state',            struct('UseRegistrationAvailability', 0)
             'no load-realign state',     struct('UseLoadRealign', 0) };
end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
global RHS_T RHS_N %#ok<GVMIS>
for v = 1:size(VARS, 1)
    P = getParams(loadParams(PARAMSET), [], true, false);
    ov = VARS{v, 2}; f = fieldnames(ov);
    for k = 1:numel(f)
        if strcmp(f{k}, 'mu_neg') && isnan(ov.mu_neg), P.mu_neg = P.mu; else, P.(f{k}) = ov.(f{k}); end
    end
    P.Velocity = 0; P.MaxRunTime = BUDGET; P.odeRelTol = TOL(1); P.odeAbsTol = TOL(2);
    P.BreakOnODEUnstable = false;
    P = getParams(P, [], true, false);
    RHS_T = zeros(1e6, 1); RHS_N = 0; tic; st = 'finished';
    try
        [F, out] = evaluateModel(@dPUdT_Counted, [0 TEND], P);
        Fend = F(end);
    catch err
        st = ['stopped: ' err.message]; Fend = NaN;
    end
    tt = RHS_T(1:RHS_N);
    fprintf('%-24s %6.1f s %9d calls  reached t=%.3f  F_end %.3f  (%s)\n', VARS{v,1}, toc, RHS_N, max([tt; 0]), Fend, st);
end
