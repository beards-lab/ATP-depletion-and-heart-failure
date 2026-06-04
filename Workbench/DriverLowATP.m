% DriverLowATP.m  — Heart-failure ATP/Pi/ADP metabolite effects on cross-bridge mechanics
%
% Parameterisation (two-Km scheme):
%   UseAtpK2  + K_T1   = 1.5 mM  : k2_eff  = k2  * [ATP]/([ATP]+1.5)  — sets Vmax sensitivity
%   UseAtpKah + K_T_kah= 4.0 mM  : kah_eff = kah * [ATP]/([ATP]+4.0)  — sets Ktr sensitivity
%
% With both active, at 2 mM ATP vs 8 mM baseline:
%   k2  ratio = 0.68  -> Vmax  -32%
%   kah ratio = 0.50  -> Ktr   ~-50% (product of both terms)
%   F0  effect: k2 slower -> heads stay in P2 longer -> +force
%               kah slower -> fewer PD heads -> -force  (net: +~10-20%)
%
% Metabolite coupling (heart failure literature):
%   Normal: ATP=8, Pi=1, ADP=0.02 mM
%   Mild:   ATP=4, Pi=2, ADP=0.05 mM
%   Mod:    ATP=2, Pi=3, ADP=0.10 mM
%   Severe: ATP=1, Pi=5, ADP=0.20 mM

clear functions;
cd('C:\home\git\ATP-depletion-and-heart-failure');
addpath(genpath('.'));

%% Load latest optimised parameters
run('params/ModelOptimParam_TL5_iter_81.m');
params0.ghostLoad  = [];
params0.ghostSave  = [];
params0.drawPlots  = false;
params0.drawPlots2 = false;

%% Enable two-Km ATP dependence
params0.UseAtpK2  = true;
params0.K_T1      = 1.5;    % mM — Km for k2 (detachment), sets Vmax sensitivity
params0.UseAtpKah = true;
params0.K_T_kah   = 4.0;    % mM — Km for kah (hydrolysis), sets Ktr sensitivity

%% Metabolite scenarios
scenarios = {
    struct('label','Normal',      'MgATP',8, 'Pi',1,  'MgADP',0.02),
    struct('label','Mild HF',     'MgATP',4, 'Pi',2,  'MgADP',0.05),
    struct('label','Moderate HF', 'MgATP',2, 'Pi',3,  'MgADP',0.10),
    struct('label','Severe HF',   'MgATP',1, 'Pi',5,  'MgADP',0.20),
};
N = numel(scenarios);

%% Helper: Hill fit -> Vmax
hillFit = @(fv_f, fv_v) fitHill(fv_f, fv_v);

%% Run each scenario
FV_curves = cell(N,1);
F0_vals   = zeros(N,1);
Vmax_vals = zeros(N,1);
Pmax_vals = zeros(N,1);

fprintf('%-14s  ATP  k2fac  kahfac   F0(kPa) Vmax(ML/s) Pmax  FVerr\n','Scenario');
for kk = 1:N
    sc = scenarios{kk};
    params0.MgATP = sc.MgATP;  params0.Pi = sc.Pi;  params0.MgADP = sc.MgADP;
    RunBakersExp;

    fv_f = fm_fv(1,1).FV_f;
    fv_v = fm_fv(1,1).FV_v;
    FV_curves{kk} = [fv_f, fv_v];
    F0_vals(kk)   = fv_f(1);               % isometric (v=0 point)
    Vmax_vals(kk) = hillFit(fv_f, fv_v);   % Hill extrapolation
    Pmax_vals(kk) = max(fv_f .* fv_v);

    % Report Michaelis-Menten factors at this [ATP]
    k2fac  = sc.MgATP/(sc.MgATP + params0.K_T1);
    kahfac = sc.MgATP/(sc.MgATP + params0.K_T_kah);
    fprintf('%-14s  %3g   %.3f   %.3f   %6.1f   %6.2f     %5.2f  %.4f\n', ...
        sc.label, sc.MgATP, k2fac, kahfac, F0_vals(kk), Vmax_vals(kk), Pmax_vals(kk), E_fv);
end

%% Normalised changes vs Normal
fprintf('\n--- Normalised to Normal ---\n');
fprintf('%-14s  dF0     dVmax   dPmax\n','Scenario');
for kk = 1:N
    fprintf('%-14s  %+.0f%%    %+.0f%%    %+.0f%%\n', scenarios{kk}.label, ...
        (F0_vals(kk)/F0_vals(1)-1)*100, ...
        (Vmax_vals(kk)/Vmax_vals(1)-1)*100, ...
        (Pmax_vals(kk)/Pmax_vals(1)-1)*100);
end

%% Plot 1: FV curves
colors  = {'k','b','r','m'};
lstyles = {'-','--','-.', ':'};
figure('Name','FV: HF progression','Position',[100 80 680 480]); hold on;
for kk = 1:N
    fv = FV_curves{kk};
    % Add Hill fit curve
    Vm = Vmax_vals(kk); F0 = F0_vals(kk);
    pfit = fitHillParams(fv(:,1), fv(:,2));
    F_fit = linspace(0, pfit(3)*0.98, 80);
    v_fit = pfit(2)*(pfit(3)-F_fit)./(F_fit+pfit(1));
    lbl = sprintf('%s  (F0=%d kPa, Vmax=%.1f ML/s, Pmax=%.1f)', ...
        scenarios{kk}.label, round(F0_vals(kk)), Vmax_vals(kk), Pmax_vals(kk));
    plot(F_fit, v_fit, [colors{kk} lstyles{kk}], 'LineWidth',2.2, 'DisplayName', lbl);
    plot(fv(:,1), fv(:,2), [colors{kk} 'o'], 'MarkerSize',7, 'HandleVisibility','off');
end
plot(fd_fv(1,1).FV_f, fd_fv(1,1).FV_v, 'ks', 'MarkerSize',9, 'LineWidth',1.5, 'DisplayName','data (8mM)');
xlabel('Force (kPa)'); ylabel('Velocity (ML/s)');
legend('Location','northeast','FontSize',8);
title('FV: ATP/Pi/ADP heart failure progression  (K_{Tm,k2}=1.5 mM, K_{Tm,kah}=4.0 mM)');
grid on;

%% Plot 2: Summary bar (normalised)
figure('Name','Summary','Position',[800 80 500 420]);
lbls = cellfun(@(s) s.label, scenarios, 'UniformOutput', false);
metrics = 100*[F0_vals/F0_vals(1), Vmax_vals/Vmax_vals(1), Pmax_vals/Pmax_vals(1)];
b = bar(metrics);
b(1).FaceColor=[0.2 0.4 0.8]; b(2).FaceColor=[0.8 0.2 0.2]; b(3).FaceColor=[0.2 0.7 0.2];
set(gca,'XTickLabel',lbls,'YLim',[0 125]);
yline(100,'k--','LineWidth',1.2);
ylabel('% of Normal'); title('F_0, V_{max}, P_{max}  (K_{T1}=1.5 mM, K_{T,kah}=4.0 mM)');
legend({'F_0','V_{max}','P_{max}'},'Location','southwest'); grid on;

%% --- local helpers ---
function Vmax = fitHill(fv_f, fv_v)
    p = fitHillParams(fv_f, fv_v);
    Vmax = p(2)*p(3)/p(1);
end

function pfit = fitHillParams(fv_f, fv_v)
    hill = @(p,F) p(2)*(p(3)-F)./(F+p(1));
    p0   = [fv_f(1)/4, max(fv_v)/2, fv_f(1)*1.1];
    opts = optimoptions('lsqcurvefit','Display','off');
    try
        pfit = lsqcurvefit(hill, p0, fv_f, fv_v, [0.1 0.1 fv_f(1)], [], opts);
    catch
        pfit = p0;
    end
end
