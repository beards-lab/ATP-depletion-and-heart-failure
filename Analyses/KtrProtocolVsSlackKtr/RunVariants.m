% RunVariants.m -- is the ktr/post-slack ratio a parametrisation contradiction
% or a missing mechanism?  Scan the two slow attachment gates and report BOTH
% the recovery rates and the full current feature objective, so any "fix" is
% priced.
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

base = getParams(loadParams('estiff125_opt'), [], true, false);
base.PlotEachSeparately = 0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures = 1; base.BreakOnODEUnstable = 0;

svt = load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
kvt = load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;
iS  = find(svt(:,2) < -1);

% data reference (03/27, the file the objective is fitted to)
Sd = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
Kd = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
FisoD = median(Sd.datatable(Sd.datatable(:,1) > svt(iS(1),1)-0.30 & Sd.datatable(:,1) < svt(iS(1),1)-0.01, 3));
wdS = recoveryWindows(Sd.datatable(:,1), Sd.datatable(:,3), svt, 'slack', FisoD);
wdS = wdS(strcmp({wdS.type},'postSlack'));
wdK = recoveryWindows(Kd.datatable(:,1), Kd.datatable(:,3), kvt, 'ktr', 1.0);
kD_S2 = wdS(2).k63; kD_K = wdK(1).k63;

sf = @(p,f,v) setfield(p,f,v); %#ok
V = { 'V0  baseline (tau_reg 26.6ms)',  @(p) p
      'V1  tau_reg = 10 ms',            @(p) sf(p,'tau_reg',0.010)
      'V2  tau_reg = 5 ms',             @(p) sf(p,'tau_reg',0.005)
      'V3  tau_reg = 2 ms',             @(p) sf(p,'tau_reg',0.002)
      'V4  tau_reg = 1 ms',             @(p) sf(p,'tau_reg',0.001)
      'V5  tau_reg = 0.2 ms',           @(p) sf(p,'tau_reg',2e-4)
      'V6  A_reg OFF (diagnostic)',     @(p) sf(p,'UseRegistrationAvailability',0)
      'V7  A_def OFF (diagnostic)',     @(p) sf(p,'UseLoadRealign',0)
      'V8  tau_reg=2ms + A_def OFF',    @(p) sf(sf(p,'tau_reg',0.002),'UseLoadRealign',0) };

R = struct('name',{},'kS',{},'kK',{},'ratio',{},'Et',{},'featCost',{},'Ektr',{},'Fiso',{});
for v = 1:size(V,1)
    p = V{v,2}(base);
    p = getParams(p, p.g, true);

    % --- slack ---
    [~, os] = runSlackExperiment(p);
    FisoM = median(os.Force(os.t >= svt(iS(2),1)-0.05 & os.t <= svt(iS(2),1)-0.005));
    ws = recoveryWindows(os.t, os.Force, svt, 'slack', FisoM);
    ws = ws(strcmp({ws.type},'postSlack'));

    % --- ktr protocol (held out of the objective) ---
    pk = p; pk.SL0 = 2*kvt(1,4); pk = getParams(pk, pk.g, true); pk.Velocity = kvt(:,2);
    [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk);
    FisoK = interp1(ok.t, ok.Force, kvt(2,1), 'linear','extrap');
    wk = recoveryWindows(ok.t, ok.Force, kvt, 'ktr', FisoK);
    % ktr trace MSE, the cost runKtrExperiment would charge
    Fi  = interp1(ok.t, ok.Force/FisoK, Kd.datatable(:,1), 'linear','extrap');
    Ektr = mean((Kd.datatable(:,3) - Fi).^2) * 10;

    % --- full current objective ---
    [Et, ~, fc, fn] = evaluateBakersExp([], p, false, false);

    R(v).name=V{v,1}; R(v).kS=[ws.k63]; R(v).kK=wk(1).k63;
    R(v).ratio=wk(1).k63/ws(2).k63; R(v).Et=Et; R(v).featCost=fc;
    R(v).Ektr=Ektr; R(v).Fiso=[FisoM FisoK];
    fprintf('%-32s kS2=%5.1f kK=%6.1f ratio=%4.2f  Et=%7.3f  Ektr=%6.3f  Fiso=%5.1f/%5.1f\n', ...
        V{v,1}, ws(2).k63, wk(1).k63, R(v).ratio, Et, Ektr, FisoM, FisoK);
end

fprintf('\n================== SUMMARY (8 mM, 03/27) ==================\n');
fprintf('%-32s %7s %7s %7s %9s %8s\n','variant','kS2','kKtr','ratio','Et(feat)','E_ktr');
fprintf('%-32s %7.1f %7.1f %7.2f %9s %8s\n','DATA 03/27',kD_S2,kD_K,kD_K/kD_S2,'--','--');
for v=1:numel(R)
    fprintf('%-32s %7.1f %7.1f %7.2f %9.3f %8.3f\n', R(v).name, R(v).kS(2), R(v).kK, R(v).ratio, R(v).Et, R(v).Ektr);
end

fprintf('\n---- per-feature cost, baseline vs tau_reg=2ms ----\n');
fn = fn(:)'; d = R(4).featCost - R(1).featCost;
[~,ord] = sort(abs(d),'descend');
fprintf('%-28s %9s %9s %9s\n','feature','V0','V3','delta');
for i = ord(1:min(12,numel(ord)))
    fprintf('%-28s %9.3f %9.3f %+9.3f\n', fn{i}, R(1).featCost(i), R(4).featCost(i), d(i));
end
save(fullfile(resdir,'variants.mat'),'R','fn','kD_S2','kD_K');
