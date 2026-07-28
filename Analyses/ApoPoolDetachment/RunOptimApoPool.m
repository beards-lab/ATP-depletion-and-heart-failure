% RunOptimApoPool.m
% =========================================================================
% Refit driver for the ATP-ceiling / apo-pool mechanism (UseR2Ceiling +
% UseD0State). See conclusions.md — READ THE PHYSIOLOGY REVIEW FIRST; three of
% the four original justifications do not survive, and this driver is set up so
% that the surviving question gets a clean answer:
%
%   Can the post-restretch recovery be brought to ~49 /s WITHOUT
%   (a) pushing the detachment ceiling below the alpha-MHC ADP-release range,
%   (b) pushing k2rip past its own bound, and
%   (c) wrecking slack / ktr / FV?
%
% If the answer is no, the mechanism is falsified — and that is a result, not a
% failed run. `cfg.w.phys` is what makes the falsification visible: with the
% alpha-MHC bounds in the objective, a fit that only works out of bounds pays
% for it. Set cfg.w.phys = 0 only to see how far out of bounds it wants to go.
%
% Usage:
%   cd(root); addpath(genpath('.')); refreshPool(5);   % see CLAUDE.md
%   RunOptimApoPool
% Resume from a saved state:
%   RESUME = 'apoPool_v1'; RunOptimApoPool
% =========================================================================

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','ApoPoolDetachment','results');
if ~exist(resdir,'dir'); mkdir(resdir); end
TAG = 'apoPool_v1';

%% ---------------- seed ----------------
% TNF + the two fixes established in RestretchVsKtrRecovery (defects 1 and 2).
base = getParams(loadParams('ThursdayNightFever'), [], true, false);
base.UseMaxwellTensionOnly = 1;      % defect 1: titin dashpot must discharge on release
base.tau_reg               = 0.002;  % defect 2: A_reg must not remember the release
base.MaxRunTime            = 400;
base.BreakOnODEUnstable    = false;

% the mechanism under test
base.UseR2Ceiling = true;
base.k2max = 800;      % SEEDED INSIDE bounds.k2 = [50 1500] (alpha-MHC ADP release).
base.K_T2c = 0.15;     % kept only for the low-ATP arm; inert at 8 mM by construction.
base.UseD0State = true;
base.k2rip   = 400;    % SEEDED INSIDE bounds.k2rip.ub = 500 (the scan used 3000-10000).
base.alphaRip = 300;
base.dr3      = -0.025;
base.k_D0T    = 60;    % NOT ATP binding (that is >10^4 /s) — a post-rupture refractory rate.
base.K_D0T    = 0.15;
base.D0FromR2 = 0;
base = getParams(base, base.g, true);

%% ---------------- free parameters ----------------
% Multiplicative modifiers on the seed values. Ordered roughly by expected leverage.
cfg.mods = { 'k2max'    % detachment ceiling            <- the contested one
             'k2rip'    % mechanical rupture rate       <- the other contested one
             'alphaRip' % rupture strain sharpness
             'dr3'      % rupture strain knee
             'k_D0T'    % post-rupture refractory rate
             'k2'       % base post-stroke detachment (compensates the ceiling)
             'ka'       % attachment (compensates lost duty)
             'kstiff2'  % force scale (compensates lost attached fraction)
             'tau_reg'  % re-check: the ceiling changes what A_reg has to do
             };
cfg.seed = base;

% Search box on the modifiers (multiples of the seed). Deliberately generous on
% the compensating parameters and tight on the two contested ones — the point is
% to find out whether compensation can rescue an in-bounds mechanism.
lo = [0.20 0.25 0.35 0.50 0.20 0.40 0.50 0.60 0.25];
hi = [1.90 1.25 3.00 2.00 5.00 2.50 2.00 1.70 4.00];

%% ---------------- targets (03/27, 8 mM) ----------------
S = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
K = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
T.svt = S.velocitytable;  T.kvt = K.velocitytable;
T.iS  = find(T.svt(:,2) < -1);  T.iR = find(T.svt(:,2) > 1);
FisoD = median(S.datatable(S.datatable(:,1) > T.svt(T.iS(1),1)-0.30 & ...
                           S.datatable(:,1) < T.svt(T.iS(1),1)-0.01, 3));
Wd  = recoveryWindows(S.datatable(:,1), S.datatable(:,3), T.svt, 'slack', FisoD);
Wdk = recoveryWindows(K.datatable(:,1), K.datatable(:,3), T.kvt, 'ktr',   1.0);
gd  = @(ty,c,f) getfieldOf(Wd, ty, c, f);
T.kS     = arrayfun(@(c) gd('postSlack',c,'k63'),     1:4);
T.kR     = arrayfun(@(c) gd('postRestretch',c,'k63'), 1:4);
T.kK     = Wdk.k63;
T.peak1  = mean(arrayfun(@(c) max(S.datatable(S.datatable(:,1) >= T.svt(T.iR(c),1) & ...
                        S.datatable(:,1) <= T.svt(T.iR(c)+1,1)+0.005, 3))/FisoD, 1:4));
T.valley = mean(arrayfun(@(c) gd('postRestretch',c,'B'), 1:4));

% Trace costs are normalised by the CURRENT best so they enter the objective at
% O(1) and cannot silently dominate the rate terms.
T.E_slack_ref = 180;  T.E_ktr_ref = 0.06;  T.E_fv_ref = 0.02;
if exist('LoadData','file') == 2 && ~evalin('base','exist(''Data_ATP'',''var'')')
    try
        LoadData
    catch
        warning('RunOptimApoPool:data','LoadData failed — force-velocity term disabled.');
    end
end
T.ATP_c    = evalinSafe('ATP_c',    []);
T.Data_ATP = evalinSafe('Data_ATP', []);
cfg.targets = T;

fprintf('TARGETS (03/27, 8 mM): post-slack %s | post-restretch %s | ktr %.1f\n', ...
        mat2str(round(T.kS,1)), mat2str(round(T.kR,1)), T.kK);
fprintf('              peak1 %.3f  valley %.3f  F_iso %.1f kPa\n', T.peak1, T.valley, FisoD);

%% ---------------- weights ----------------
cfg.w = struct( ...
    'kS',   1.0, ...   % post-slack rate      (already fits — hold it)
    'kR',   4.0, ...   % post-restretch rate  <- THE TARGET, weighted up
    'kK',   1.0, ...   % ktr rate             (already fits — hold it)
    'peak', 2.0, ...   % restretch peak       <- the ceiling's known failure mode
    'vall', 2.0, ...   % post-restretch valley
    'slack',1.0, 'ktr', 1.0, 'fv', 1.0, ...
    'phys', 1.0);      % alpha-MHC bounds ACTIVE. Set 0 only to probe how far out it wants to go.
cfg.bounds = parameterBounds();

%% ---------------- baseline evaluation ----------------
g0 = ones(numel(cfg.mods),1);
lo = lo(:); hi = hi(:);
if any(g0 <= lo) || any(g0 >= hi)
    error('RunOptimApoPool:box','the seed (g=1) must be strictly inside [lo,hi].');
end
fprintf('\nEvaluating seed...\n'); tic;
[C0, D0] = apoPoolCost(g0, base, cfg); t1 = toc;
reportRow('SEED', C0, D0);
fprintf('one evaluation = %.0f s -> ~%.0f min for 300 iterations\n', t1, 300*t1/60);
if ~isfinite(C0) || C0 >= 1e8
    error('RunOptimApoPool:seed','the seed does not evaluate — fix it before optimising.');
end

%% ---------------- optimise ----------------
% Bounded by a logistic transform so fminsearch cannot leave the box.
%   g = lo + (hi-lo)./(1+exp(-x));  x0 is the inverse map of the seed g0 = 1.
x0 = log((g0 - lo) ./ (hi - g0));
apoPoolObjective('reset', resdir, TAG);

if exist('RESUME','var') && ~isempty(RESUME)
    f = fullfile(resdir, ['state_' RESUME '.mat']);
    if exist(f,'file')
        L = load(f);
        x0 = log((L.state.gbest - lo) ./ (hi - L.state.gbest));
        fprintf('resumed %s at C = %.4f\n', RESUME, L.state.best);
    else
        warning('RunOptimApoPool:resume','no saved state "%s" — starting from the seed.', RESUME);
    end
end

opts = optimset('Display','iter','MaxFunEvals',400,'MaxIter',400,'TolX',1e-3,'TolFun',1e-3);
xopt = fminsearch(@(x) apoPoolObjective(x, base, cfg, lo, hi), x0, opts);

st   = apoPoolObjective('state');
gopt = st.gbest;                      % best seen, not necessarily the last simplex point

%% ---------------- report ----------------
[Cf, Df] = apoPoolCost(gopt, base, cfg);
reportRow('FINAL', Cf, Df);
fprintf('\n%-10s %12s %12s %12s   %s\n','param','seed','fitted','x','in bounds?');
pf = base; for i=1:numel(cfg.mods); pf.(cfg.mods{i}) = cfg.seed.(cfg.mods{i})*gopt(i); end
pf = getParams(pf, pf.g, true);
for i = 1:numel(cfg.mods)
    nm = cfg.mods{i}; v = pf.(nm); ok = 'n/a';
    if isfield(cfg.bounds, nm)
        b = cfg.bounds.(nm);
        if v >= b.lb && v <= b.ub; ok = sprintf('yes [%g %g]', b.lb, b.ub);
        else; ok = sprintf('*** OUT [%g %g]', b.lb, b.ub); end
    end
    fprintf('%-10s %12.4g %12.4g %12.3f   %s\n', nm, cfg.seed.(nm), v, gopt(i), ok);
end
[~, viol] = evalPhysiologyCost(pf, cfg.bounds);
fprintf('\nphysiology violations: %d\n', numel(viol));
fprintf('\nVERDICT: if k2max landed below ~350 (alpha-MHC ADP release) or k2rip above 500,\n');
fprintf('the mechanism only works outside species-appropriate kinetics — see conclusions.md.\n');
try
    writeParamsToMFile(pf, fullfile(root,'params',[TAG '_opt.m']));
    fprintf('saved params/%s_opt.m\n', TAG);
catch e
    warning('RunOptimApoPool:save','writeParamsToMFile failed (%s); state .mat still written.', e.message);
end
save(fullfile(resdir, ['final_' TAG '.mat']), 'pf', 'gopt', 'Cf', 'Df', 'cfg', '-v7.3');
fprintf('saved results/final_%s.mat and results/state_%s.mat\n', TAG, TAG);

%% ---------------- helpers ----------------
function reportRow(tag, C, D)
    if ~isfield(D,'kR'); fprintf('%-14s cost %.4f (no detail)\n', tag, C); return; end
    fprintf(['%-14s C=%8.4f | rate %6.3f shape %6.3f trace %6.3f phys %6.3f' ...
             ' | kS %5.1f kR %5.1f kK %5.1f | peak %.3f vall %.3f' ...
             ' | E slack %7.1f ktr %6.3f fv %6.3f\n'], ...
        tag, C, D.c_rate, D.c_shape, D.c_trace, D.c_phys, ...
        D.kS, D.kR, D.kK, D.peak1, D.valley, D.E_slack, D.E_ktr, D.E_fv);
end

function v = getfieldOf(W, ty, c, f)
    w = W(strcmp({W.type},ty) & [W.cyc]==c);
    if isempty(w); v = NaN; else; v = w.(f); end
end

function v = evalinSafe(nm, dflt)
    try
        v = evalin('base', nm);
    catch
        v = dflt;
    end
end
