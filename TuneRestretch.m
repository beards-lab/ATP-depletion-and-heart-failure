% TuneRestretch.m
% Investigate and fix the restretch double-peak discrepancy.
% See Docs/restretch-discrepancy.md for full analysis.
%
% Three hypotheses tested:
%   1. kSE too low  — serial elastic element absorbs rapid restretch,
%                     cross-bridges never see the strain increase
%   2. R1 too fast at near-zero strain — bridges lost during restretch
%                     even before they can be yanked to high positive strain
%   3. ka too low   — post-valley reattachment too slow, second peak blunted
%
% Strategy: start from the RunMechanismEvaluation baseline, apply each
% fix in isolation and in combination, compare to data.

%% Baseline params (identical to RunMechanismEvaluation §2 block)
clear params0;
params0 = getParams();
ModelParamsSlackKtrOpt;

params0.UseTitinIdentifiedPassive = 0;
params0.k_pas   = 77*2*2*2;
params0.gamma   = 7.9;
params0.UseOverlap       = true;
params0.UseOverlapFactor = true;
params0.L_thick = 1.65;
params0.L_thin  = 1.15;
params0.UseNegativeKstiff = 0;
params0.kstiff2 = 15000;
params0.kstiff1 = "=kstiff2";
params0.MaxSlackNegativeForce = -1;
params0.kSE     = 2000;

params0.ka  = 67*2*2;
params0.kd  = 9.55*2;
params0.k1  = 270*2;
params0.k2  = 133*0.85;
params0.kah = 30*1.5;

params0.A1AttachmentWidth      = 0.004;
params0.UseTargetZoneSaturation = false;
params0.max_attached_per_bin   = 0.05;

% R1D
params0.PieceWiseStrainDepR1DX     = [-0.05, -0.015, -0.004, 0, 0.004, 0.01, 0.1];
params0.PieceWiseStrainDepR1DParams = [1000, 1000, 10.147, 1, 2, 30, 500]*4;
params0.dr2    = "=dr";
params0.mu     = 0.633*0.5;
params0.mu_neg = "=mu";

% R1 (detachment from strongly attached state)
params0.PieceWiseStrainDepX      = [-1,    -0.0083,  -0.0023,  0.0046,  0.02];
params0.PieceWiseStrainDepParams = [10.0,   10.8001,   2.8824,  1.1554,  0.015];

% R2
params0.PieceWiseStrainDep2X      = [-0.01, -0.0105, -0.0085, -0.0055, -0.004, 0.000, 0.006, 0.01, 0.0233, 0.1];
params0.PieceWiseStrainDep2Params = [50,     50,      50.6847,  1,       1,     0.5,   1,     20,   50,     50];

% R21
params0.PieceWiseStrainDepR21X      = [-.1  .1];
params0.PieceWiseStrainDepR21Params = [ 1    1];

params0.slope                 = 11000;
params0.UseA2AttachmentShift  = true;
params0.k2d                   = 80000;
params0.UseA2MechanicalRecocking = false;
params0.BreakOnODEUnstable    = false;
params0.dS                    = 0.0006;
params0.MaxStrainArraySize    = 40;
params0.UseLatticeSpacing     = false;
params0.sigma_lattice         = 0.0006;
params0.UseSuperRelaxed       = false;
params0.UseSuperRelaxedADP    = false;
params0.xrate                 = 0.5;

params0.RunSlack          = true;
params0.RunSlackSegments  = 'Last';
params0.RunForceVelocity  = false;
params0.RunForceLengthEstim = false;
params0.PlotEachSeparately  = false;
params0.ShowResidualPlots   = false;
params0.ghostLoad = '';
params0.ghostSave = '';

LoadData;

%% Define variants
% Each row: {label, struct of param overrides}
%
% Key insight (Docs/restretch-discrepancy.md):
%   UseA2AttachmentShift IS the reattachment mechanism for restretch — p2 heads
%   at s > s_threshold_R hop by -d_actin toward zero strain, producing the rapid
%   force recovery seen in data. But with kSE=2000, bridges never reach the
%   threshold. Primary fix: increase kSE so bridges see the restretch strain.
%   A2Shift then fires automatically. Secondary: tune slope / s_threshold_R.

% kstiff2 is FIXED — it sets steady-state force. Cannot reduce it.
% Spike amplitude = attached_fraction × kstiff2 × bridge_strain.
% To limit peak without touching kstiff2: lower s_threshold_R so A2Shift
% fires earlier, capping bridge strain before it builds up too high.
% kSE still needed to transmit the restretch to bridges at all.

ka_base = 67*2*2;
variants = {
    'Baseline',                         struct()
    'kSE x3 + early A2Shift',           struct('kSE', 6000, 'slope', 22000, 's_threshold_R', 0.002)
    'ViscoSE c=1000',                   struct('UseViscoelasticSE', true, 'c_SE_visc', 1000)
    'ViscoSE c=5000',                   struct('UseViscoelasticSE', true, 'c_SE_visc', 5000)
    'ViscoSE c=5000 + kSE x2',          struct('UseViscoelasticSE', true, 'c_SE_visc', 5000, 'kSE', 4000)
    'Catch k=0.05',                     struct('UseCatchBond', true, 'k_catch_bond', 0.05)
    'Catch k=0.05 + kSE=3000',          struct('UseCatchBond', true, 'k_catch_bond', 0.05, 'kSE', 3000)
    'Catch k=0.05 + ViscoSE c=1000',    struct('UseCatchBond', true, 'k_catch_bond', 0.05, ...
                                                'UseViscoelasticSE', true, 'c_SE_visc', 1000)
    'Catch k=0.05 + ViscoSE c=5000',    struct('UseCatchBond', true, 'k_catch_bond', 0.05, ...
                                                'UseViscoelasticSE', true, 'c_SE_visc', 5000)
    'Catch k=0.05 + kSE=3000 + Visco',  struct('UseCatchBond', true, 'k_catch_bond', 0.05, ...
                                                'kSE', 3000, 'UseViscoelasticSE', true, 'c_SE_visc', 1000)
};

N = size(variants, 1);
results = cell(N, 1);
costs   = zeros(N, 1);

fprintf('\n%-30s  E_slack\n', 'Variant');
fprintf('%s\n', repmat('-', 1, 42));
for i = 1:N
    p = params0;
    mods = variants{i, 2};
    fnames = fieldnames(mods);
    for k = 1:numel(fnames)
        p.(fnames{k}) = mods.(fnames{k});
    end
    [E, out, ~, ~] = runSlackExperiment(p);
    results{i} = out;
    costs(i)   = sum(E);
    fprintf('%-30s  %.4f\n', variants{i,1}, costs(i));
end

%% Comparison plots
t_zoom  = [2.87, 3.05];
colors  = lines(N);
data_t  = results{1}.datatable(:, 1);
data_F  = results{1}.datatable(:, 3);

figure(210); clf; tiledlayout(3, 1, 'TileSpacing', 'compact');

% --- Force traces ---
ax1 = nexttile; hold on;
plot(data_t, data_F, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Data');
for i = 1:N
    out_i = results{i};
    nr = makeMonotonous(out_i.t);
    plot(out_i.t(nr), out_i.Force(nr), '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s  (E=%.3f)', variants{i,1}, costs(i)));
end
xlim(t_zoom); grid on;
ylabel('Force (rel.)');
title('Restretch force comparison');
legend('Location', 'southeast', 'FontSize', 8);

% --- Residuals ---
ax2 = nexttile; hold on;
yline(0, 'k--', 'LineWidth', 0.8);
for i = 1:N
    out_i = results{i};
    nr = makeMonotonous(out_i.t);
    Fi  = interp1(out_i.t(nr), out_i.Force(nr), data_t, 'linear', NaN);
    err = data_F - Fi;
    plot(data_t, err, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', variants{i,1});
end
xlim(t_zoom); grid on;
ylabel('Error (data - model)');
title('Residuals');
legend('Location', 'best', 'FontSize', 8);

% --- SL vs LXB for each variant (shows how kSE change affects LXB tracking) ---
ax3 = nexttile; hold on;
for i = 1:N
    out_i = results{i};
    nr = makeMonotonous(out_i.t);
    gap = out_i.SL(nr) - out_i.LXB(nr);
    plot(out_i.t(nr), gap, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', variants{i,1});
end
xlim(t_zoom); grid on;
ylabel('SL - LXB (\mum)');
xlabel('t (s)');
title('Serial element stretch (SL - LXB)  — smaller = kSE stiffer');
legend('Location', 'best', 'FontSize', 8);

linkaxes([ax1 ax2 ax3], 'x');

%% Print best variant
[~, ibest] = min(costs);
fprintf('\nBest variant: %s  (E = %.4f)\n', variants{ibest,1}, costs(ibest));
