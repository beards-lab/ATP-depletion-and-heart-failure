% AnalyseRestretchMechanisms.m — Systematic sweep of restretch mechanisms
cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));
load('IdentifyLastSlack_state.mat','state');

fn_target = {'ktr|SLslack','A|SLslack','peak1_y','peak1_dSL', ...
             'peak2','vall_y','vall2_dy','steady'};
mods_p1   = {'kSE','mu','ka','kd','k1','k2','slope','k_pas'};

%% Base params
clear params0; params0=getParams(); ModelParamsInitManualLastSlack; LoadData;
params0.RunSlackSegments='Last'; params0.RunForceVelocity=false;
params0.RunForceLengthEstim=false; params0.PlotEachSeparately=false;
params0.ShowResidualPlots=false; params0.ghostLoad=''; params0.ghostSave='';
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false;
params0.fn=fn_target; params0.mods=mods_p1; params0.g=state.g_opt1;
p_base = params0;

p_resolved = getParams(p_base, state.g_opt1, false, true);

fprintf('\n========== RESTRETCH MECHANISM ANALYSIS ==========\n');
fprintf('Phase 1 g: %s\n\n', num2str(state.g_opt1,'%.3f '));

[c0,~,~] = evalP(p_base, fn_target);
fprintf('E0 Phase1 baseline:\n');
printRow('Phase1 opt', c0, c0);
fprintf('\n');

%% E1: c_SE_visc sweep
fprintf('--- E1: c_SE_visc sweep (kstiff2 unchanged) ---\n');
c_SE_vals = [0, 3, 6, 10, 15, 20, 30];
res_E1 = nan(numel(c_SE_vals), numel(fn_target));
for i = 1:numel(c_SE_vals)
    p = p_base;
    p.UseViscoelasticSE = (c_SE_vals(i) > 0);
    p.c_SE_visc = c_SE_vals(i);
    [c,~,~] = evalP(p, fn_target);
    res_E1(i,:) = c;
    printRow(sprintf('c_SE_visc=%g', c_SE_vals(i)), c, c0);
end
fprintf('\n');

%% E2: kstiff2 sweep at c_SE_visc=6
fprintf('--- E2: kstiff2 sweep (c_SE_visc=6) ---\n');
kst2_vals = [15000, 13000, 11000, 9750, 8000];
res_E2 = nan(numel(kst2_vals), numel(fn_target));
for i = 1:numel(kst2_vals)
    p = p_base;
    p.UseViscoelasticSE = true; p.c_SE_visc = 6;
    p.kstiff2 = kst2_vals(i); p.kstiff1 = '=kstiff2';
    [c,~,~] = evalP(p, fn_target);
    res_E2(i,:) = c;
    printRow(sprintf('kstiff2=%g + c_SE=6', kst2_vals(i)), c, c0);
end
fprintf('\n');

%% E3: s_threshold_R sweep at c_SE_visc=6
fprintf('--- E3: s_threshold_R sweep (c_SE_visc=6) ---\n');
str_vals = [0.002, 0.003, 0.0046, 0.006, 0.008, 0.012, 0.020];
res_E3 = nan(numel(str_vals), numel(fn_target));
for i = 1:numel(str_vals)
    p = p_base;
    p.UseViscoelasticSE = true; p.c_SE_visc = 6;
    p.s_threshold_R = str_vals(i);
    [c,~,~] = evalP(p, fn_target);
    res_E3(i,:) = c;
    printRow(sprintf('s_thresh_R=%.4f + c_SE=6', str_vals(i)), c, c0);
end
fprintf('\n');

%% E4: c_SE_visc x s_threshold_R grid
fprintf('--- E4: c_SE_visc × s_threshold_R grid ---\n');
c_se_grid = [6, 15, 30];
str_grid  = [0.002, 0.003, 0.0046, 0.008];
best_E4 = Inf; best_cfg = struct('c_SE_visc',6,'s_threshold_R',0.0046,'costs',c0);
row = 0;
res_E4 = nan(numel(c_se_grid)*numel(str_grid), numel(fn_target));
for i = 1:numel(c_se_grid)
    for j = 1:numel(str_grid)
        row = row+1;
        p = p_base;
        p.UseViscoelasticSE = true; p.c_SE_visc = c_se_grid(i);
        p.s_threshold_R = str_grid(j);
        [c,~,~] = evalP(p, fn_target);
        res_E4(row,:) = c;
        printRow(sprintf('c_SE=%g  str=%.4f', c_se_grid(i), str_grid(j)), c, c0);
        if sum(c) < best_E4
            best_E4 = sum(c);
            best_cfg = struct('c_SE_visc',c_se_grid(i),'s_threshold_R',str_grid(j),'costs',c);
        end
    end
end
fprintf('\n  Best E4: c_SE_visc=%g  s_threshold_R=%.4f  tot=%.3f\n\n', ...
    best_cfg.c_SE_visc, best_cfg.s_threshold_R, best_E4);

%% E5: slope (A2Shift rate) sweep at best E4
fprintf('--- E5: slope multiplier sweep at best E4 ---\n');
slope_mults = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0];
res_E5 = nan(numel(slope_mults), numel(fn_target));
for i = 1:numel(slope_mults)
    p = p_base;
    p.UseViscoelasticSE = true;
    p.c_SE_visc     = best_cfg.c_SE_visc;
    p.s_threshold_R = best_cfg.s_threshold_R;
    p.slope = p_resolved.slope * slope_mults(i);
    [c,~,~] = evalP(p, fn_target);
    res_E5(i,:) = c;
    printRow(sprintf('slope×%.2f + best_E4', slope_mults(i)), c, c0);
end
fprintf('\n');

%% Summary
fprintf('=================== SUMMARY ===================\n');
fprintf('  %-36s  %6s  %6s  %6s  %6s  %6s\n','Config','Total','pk1','valley','pk2','ktr');
fprintf('  %s\n', repmat('-',1,62));
psr = @(lbl,c) fprintf('  %-36s  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n', ...
    lbl, sum(c), c(3), c(6), c(5), c(1));
psr('E0 Phase1 baseline', c0);
[~,bi]=min(sum(res_E1,2)); psr(sprintf('E1 best c_SE=%g',     c_SE_vals(bi)),  res_E1(bi,:));
[~,bi]=min(sum(res_E2,2)); psr(sprintf('E2 best kst2=%g',     kst2_vals(bi)),  res_E2(bi,:));
[~,bi]=min(sum(res_E3,2)); psr(sprintf('E3 best str=%.4f',    str_vals(bi)),   res_E3(bi,:));
psr(sprintf('E4 best c_SE=%g str=%.4f', best_cfg.c_SE_visc, best_cfg.s_threshold_R), best_cfg.costs);
[~,bi]=min(sum(res_E5,2)); psr(sprintf('E5 best slope×%.2f', slope_mults(bi)), res_E5(bi,:));
fprintf('\n  pk1=peak1_y  vly=vall_y  pk2=peak2  ktr=ktr|SLslack\n');

%% Save results
save('AnalyseRestretchMechanisms_results.mat', ...
    'c0','res_E1','res_E2','res_E3','res_E4','res_E5', ...
    'c_SE_vals','kst2_vals','str_vals','c_se_grid','str_grid','slope_mults', ...
    'best_cfg','fn_target');
fprintf('\nResults saved to AnalyseRestretchMechanisms_results.mat\n');
fprintf('Analysis complete.\n');

%% ── Local functions ─────────────────────────────────────────────────────
function [E, fm, fd] = evalP(p, fn_target)
    try
        [~,~,fm,fd] = runSlackExperiment(p);
        E = evalFeatureCost(fd, fm, fn_target, 1);
    catch ex
        fprintf('  ERROR: %s\n', ex.message);
        E = nan(1, numel(fn_target)); fm=[]; fd=[];
    end
end

function printRow(label, costs, ref_costs)
    delta = sum(costs) - sum(ref_costs);
    fprintf('  %-36s  tot=%6.3f  d=%+6.3f  [pk1=%5.3f  vly=%5.3f  pk2=%5.3f  ktr=%5.3f]\n', ...
        label, sum(costs), delta, costs(3), costs(6), costs(5), costs(1));
end
