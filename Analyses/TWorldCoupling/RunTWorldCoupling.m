%% RunTWorldCoupling.m -- T-World (human ventricular myocyte) as the Ca source of the XB model
%
% Compares the three couplings of CaInput = 'tworld' (tworldCaLayer.m) on a paced
% isometric twitch, against the standalone T-World beat they all start from:
%   C0  standalone T-World, one beat from the 1 Hz 'TW_endo' steady state (ThirdParty/TWorld)
%   C1  'signal'  C0's cytosolic Ca as a periodic table drives the thin filament (no extra state)
%   C2  'land'    T-World's 93 states inside the XB ODE; Land's troponin buffers the Ca (open loop)
%   C3  'xb'      as C2, but THIS model's thin filament buffers the troponin Ca (closed loop)
% each with two thin-filament sets:
%   'human'   the EHT-twitch fit of params/human_twitch_t1_isoT.m (tf_k_on 9.2, k_off 5.1, K_coop 20)
%   'marzban' the Marzban 2020 rat defaults (101.2, 723.9, 9.68), the fast/weak thin filament
% Checks: C2 reproduces C0's Ca (same equations, different integrator/tolerance), C1 = C2 up
% to the table interpolation, and C3 - C2 is the feedback of the XB thin filament on Ca.
%
% Protocol: SL 2.2 um fixed ends (isometric), 8 mM MgATP as in the snapshot, PRE s of rest
% (T-World frozen at its steady state, XB model relaxing), then NB beats at 1 Hz.
% Results: results/tworld_coupling.mat, results/tworld_coupling.csv, results/*.png.
% Runtime ~5 min (R2023a, one worker); every case is cached in results/case_*.mat, so
% delete those to rerun.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');  if ~exist(res, 'dir'), mkdir(res); end
addpath(genpath(root));

PRE = 0.3;  NB = 3;  T = 1;  tEnd = PRE + NB*T;
TFs = struct('human',   struct(), ...
             'marzban', struct('tf_k_on', 101.185, 'tf_k_off', 723.852, 'tf_K_coop', 9.6846));
couplings = {'signal', 'land', 'xb'};

base = loadParams('human_twitch_t1_isoT');
base.CaModel = 'twoState';  base.CaInput = 'tworld';
base.CaPeriod = T;  base.CaFirstBeat = PRE;
base.odeRelTol = 1e-3;                      % T-World's own default (evaluateModel's 1e-2 is loose for an AP)
base.Velocity = 0;

R = struct();
tfNames = fieldnames(TFs);
for a = 1:numel(tfNames)
    tfn = tfNames{a};
    for b = 1:numel(couplings)
        cp = couplings{b};  key = sprintf('%s_%s', tfn, cp);
        f = fullfile(res, ['case_' key '.mat']);
        if exist(f, 'file'), S = load(f);  R.(key) = S.r;  fprintf('%-16s cached\n', key);  continue; end
        P = base;
        for g = fieldnames(TFs.(tfn))', P.(g{1}) = TFs.(tfn).(g{1}); end
        P.CaTworld = struct('coupling', cp);
        P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
        tic;  [~, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 tEnd], P);  wall = toc;
        r = struct('key', key, 'tf', tfn, 'coupling', cp, 'wall', wall, 'nSteps', numel(out.t), ...
            't', out.t(:), 'F', out.Force(:), 'Ca', out.Ca(:), 'gate', out.CaGate(:), ...
            'N', out.CaState(:, 1), 'bound', out.p1_0(:) + out.p2_0(:), 'SR', out.SR(:), ...
            'ref', P.ca.tw.ref, 'tf_k_on', P.tf_k_on, 'tf_k_off', P.tf_k_off, 'tf_K_coop', P.tf_K_coop);
        if P.ca.tw.coupling > 1
            Y = out.CaState(:, 2:end);                       % the 93 T-World states
            r.Ta = landTa(Y);  r.CaTRPN = Y(:, 71);  r.TnC = 1e3*Y(:, 13);   % kPa, -, uM
            r.CaSR = Y(:, 10);  r.V = Y(:, 1);
        end
        save(f, 'r');  R.(key) = r;
        fprintf('%-16s %6.1f s wall, %5d points, Fmax %.2f kPa, Ca %.3f..%.3f uM\n', key, wall, numel(out.t), max(r.F), min(r.Ca), max(r.Ca));
    end
end

%% self-check of the Land tension formula on the standalone beat
ref = R.human_land.ref;
fprintf('landTa check: standalone beat, max |Ta(getCurrentsStructure) - Ta(formula)| = %.2e kPa (Ta max %.2f kPa)\n', ...
    max(abs(ref.Ta(:) - landTa(ref.X))), max(ref.Ta));

%% metrics per case on the last beat (and beat 1 vs the standalone beat)
rows = {};
for a = 1:numel(tfNames)
    for b = 1:numel(couplings)
        key = sprintf('%s_%s', tfNames{a}, couplings{b});  r = R.(key);
        tb = PRE + (NB-1)*T;  m = r.t >= tb & r.t <= tb + T;  tt = r.t(m) - tb;
        [CaPk, i] = max(r.Ca(m));  CaTTP = tt(i);  Cadia = min(r.Ca(m));
        Fd = r.F(m) - r.F(find(m, 1));  [FPk, j] = max(Fd);  FTTP = tt(j);
        k = find(Fd(j:end) <= 0.5*FPk, 1);  RT50 = NaN;  if ~isempty(k), RT50 = tt(j + k - 1) - FTTP; end
        % beat 1 vs the standalone beat C0 (same starting state)
        m1 = r.t >= PRE & r.t <= PRE + T;  Ca1 = interp1(r.t(m1) + 1e-9*(0:nnz(m1)-1)', r.Ca(m1), PRE + r.ref.t, 'linear');
        dCa0 = max(abs(Ca1 - r.ref.Ca));
        % the same beat in the 'land' run of the same thin filament
        rl = R.(sprintf('%s_land', tfNames{a}));
        Cal = interp1(rl.t + 1e-9*(0:numel(rl.t)-1)', rl.Ca, r.t, 'linear');
        dCaLand = max(abs(r.Ca - Cal));
        rows(end+1, :) = {key, tfNames{a}, couplings{b}, r.wall, r.nSteps, Cadia, CaPk, 1e3*CaTTP, FPk, 1e3*FTTP, 1e3*RT50, max(r.gate(m)), dCa0, dCaLand}; %#ok<SAGROW>
    end
end
Tb = cell2table(rows, 'VariableNames', {'case', 'thinFilament', 'coupling', 'wall_s', 'nPoints', 'Ca_dia_uM', 'Ca_peak_uM', 'Ca_TTP_ms', ...
    'Fdev_peak_kPa', 'F_TTP_ms', 'F_RT50_ms', 'gate_max', 'maxdCa_beat1_vs_standalone_uM', 'maxdCa_vs_land_uM'});
disp(Tb);
writetable(Tb, fullfile(res, 'tworld_coupling.csv'));
save(fullfile(res, 'tworld_coupling.mat'), 'R', 'Tb', 'PRE', 'NB', 'T');

%% figure: Ca, force, thin filament, per thin-filament set
col = struct('signal', [0.3 0.3 0.3], 'land', [0 0.45 0.74], 'xb', [0.85 0.33 0.1]);
for a = 1:numel(tfNames)
    tfn = tfNames{a};
    fh = figure('Position', [50 50 1300 850], 'Color', 'w');
    tl = tiledlayout(fh, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('T-World as the Ca source, thin filament ''%s'' (tf\\_k\\_on %.3g, k\\_off %.3g, K\\_coop %.3g)', ...
        tfn, R.([tfn '_land']).tf_k_on, R.([tfn '_land']).tf_k_off, R.([tfn '_land']).tf_K_coop), 'Interpreter', 'tex');
    tb = PRE + (NB-1)*T;
    % (1) Ca, whole run
    ax = nexttile(1);  hold(ax, 'on');
    ref = R.([tfn '_land']).ref;
    for k = 0:NB-1, plot(ax, PRE + k*T + ref.t, ref.Ca, ':', 'Color', [0 0 0], 'LineWidth', 1.2); end
    for b = 1:numel(couplings), r = R.([tfn '_' couplings{b}]);  plot(ax, r.t, r.Ca, 'Color', col.(couplings{b}), 'LineWidth', 1); end
    ylabel(ax, '[Ca^{2+}]_i (\muM)');  title(ax, 'cytosolic Ca (dotted: standalone T-World beat, repeated)');
    legend(ax, [{'standalone'}, couplings], 'Location', 'northeast');  xlim(ax, [0 tEnd]);
    % (2) Ca, last beat zoom + difference
    ax = nexttile(2);  hold(ax, 'on');
    plot(ax, 1e3*ref.t, ref.Ca, ':k', 'LineWidth', 1.2);
    for b = 1:numel(couplings), r = R.([tfn '_' couplings{b}]);  m = r.t >= tb & r.t <= tb + T;  plot(ax, 1e3*(r.t(m) - tb), r.Ca(m), 'Color', col.(couplings{b}), 'LineWidth', 1); end
    xlabel(ax, 'ms after stimulus');  ylabel(ax, '[Ca^{2+}]_i (\muM)');  title(ax, sprintf('last beat (%d)', NB));  xlim(ax, [0 600]);
    % (3) force
    ax = nexttile(3);  hold(ax, 'on');
    for b = 1:numel(couplings), r = R.([tfn '_' couplings{b}]);  plot(ax, r.t, r.F, 'Color', col.(couplings{b}), 'LineWidth', 1); end
    ylabel(ax, 'force (kPa)');  title(ax, 'XB model force');  xlim(ax, [0 tEnd]);
    % (4) force last beat, with Land Ta of the embedded T-World (its own contraction model, scaled)
    ax = nexttile(4);  hold(ax, 'on');
    for b = 1:numel(couplings)
        r = R.([tfn '_' couplings{b}]);  m = r.t >= tb & r.t <= tb + T;
        plot(ax, 1e3*(r.t(m) - tb), r.F(m) - r.F(find(m, 1)), 'Color', col.(couplings{b}), 'LineWidth', 1);
    end
    r = R.([tfn '_land']);  m = r.t >= tb & r.t <= tb + T;
    plot(ax, 1e3*(r.t(m) - tb), (r.Ta(m) - r.Ta(find(m, 1)))*max(r.F(m) - r.F(find(m, 1)))/max(r.Ta(m) - r.Ta(find(m, 1))), '--', 'Color', [0.5 0.5 0.5]);
    xlabel(ax, 'ms after stimulus');  ylabel(ax, 'developed force (kPa)');  title(ax, 'last beat; dashed: T-World''s Land tension, scaled to the same peak');
    legend(ax, [couplings, {'Land Ta (scaled)'}], 'Location', 'northeast');  xlim(ax, [0 1000]);
    % (5) thin filament: permissive fraction vs Land Ca_TRPN
    ax = nexttile(5);  hold(ax, 'on');
    for b = 1:numel(couplings), r = R.([tfn '_' couplings{b}]);  m = r.t >= tb & r.t <= tb + T;  plot(ax, 1e3*(r.t(m) - tb), 1 - r.N(m), 'Color', col.(couplings{b}), 'LineWidth', 1); end
    r = R.([tfn '_land']);  m = r.t >= tb & r.t <= tb + T;  plot(ax, 1e3*(r.t(m) - tb), r.CaTRPN(m), '--', 'Color', [0.5 0.5 0.5]);
    xlabel(ax, 'ms after stimulus');  ylabel(ax, 'fraction');  title(ax, 'permissive fraction 1-N (XB thin filament); dashed: Land Ca\_TRPN');
    legend(ax, [couplings, {'Land Ca_{TRPN}'}], 'Location', 'northeast');  xlim(ax, [0 1000]);
    % (6) troponin-bound Ca in T-World's bookkeeping state (buffers_TnClow): Land's vs ours
    ax = nexttile(6);  hold(ax, 'on');
    for b = 2:numel(couplings), r = R.([tfn '_' couplings{b}]);  m = r.t >= tb & r.t <= tb + T;  plot(ax, 1e3*(r.t(m) - tb), r.TnC(m), 'Color', col.(couplings{b}), 'LineWidth', 1); end
    xlabel(ax, 'ms after stimulus');  ylabel(ax, 'Ca bound to TnC (\muM)');  title(ax, 'troponin Ca (T-World state 13): Land (land) vs XB thin filament (xb)');
    legend(ax, couplings(2:end), 'Location', 'northeast');  xlim(ax, [0 1000]);
    exportgraphics(fh, fullfile(res, sprintf('tworld_coupling_%s.png', tfn)), 'Resolution', 110);
end
fprintf('done: %s\n', res);

function Ta = landTa(Y)
% Land 2017 active tension as in model_TWorld.m/getContractionLand (intact mode, lambda = 1,
% no PKA): Ta = Lfac*PKAForceMultiplier*(Tref/dr)*((ZETAS+1)*XS + ZETAW*XW), Lfac = 1, Tref = 80, dr = 0.25.
XS = max(0, Y(:, 69));  XW = max(0, Y(:, 70));  ZS = Y(:, 73);  ZW = Y(:, 74);
Ta = (80/0.25)*((ZS + 1).*XS + ZW.*XW);
end
