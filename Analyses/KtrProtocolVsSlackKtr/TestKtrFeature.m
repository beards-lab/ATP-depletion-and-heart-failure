% TestKtrFeature.m -- verify the new ktrP_* feature hook.
%  1. it reproduces the analysis-folder estimator (recoveryWindows) exactly
%  2. it does NOT fire on the FV2 files that share the same runner
%  3. with RunKtr=1 the objective can now score the ktr protocol
%  4. with RunKtr=0 nothing changes (Et identical to before the hook)
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
base=getParams(loadParams('estiff125_opt'),[],true,false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad=''; base.EvalFeatures=1;
base.BreakOnODEUnstable=0;
kvt=load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;

fprintf('\n--- 1. agreement with recoveryWindows (the analysis estimator) ---\n');
for tau=[base.tau_reg 2e-4]
    p=base; p.tau_reg=tau; p.RunKtr=1; p=getParams(p,p.g,true);
    [~,ok,fm,fd]=runKtrExperiment(p,[]);
    FisoK=interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
    w=recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK);
    fprintf('tau_reg=%7.4g ms | ktrP_k63 model %6.2f  (recoveryWindows %6.2f, diff %+.3f)  data %6.2f  F20 %.3f/%.3f\n', ...
        tau*1e3, fm.ktrP_k63, w(1).k63, fm.ktrP_k63-w(1).k63, fd.ktrP_k63, fm.ktrP_F20, fd.ktrP_F20);
end

fprintf('\n--- 2. FV2 files must emit NOTHING (they end on a lengthening ramp) ---\n');
for f={'protocol_04_03_2026_velocitytable_FV2.mat','protocol_04_10_2026_velocitytable_ktr.mat'}
    p=base; p.RunKtr=1; p.ktr_velocitytableonfile=f{1}; p=getParams(p,p.g,true);
    [~,~,fm,~]=runKtrExperiment(p,[]);
    fprintf('%-46s -> features: %s\n', f{1}, strjoin([{'(none)'} fieldnames(fm)'],' '));
end

fprintf('\n--- 3./4. objective with the hook ---\n');
p=base; p.RunKtr=0; p=getParams(p,p.g,true);
Et0=evaluateBakersExp([],p,false,false);
fprintf('RunKtr=0, fn unchanged            : Et = %.4f   (was 0.9559 before the hook)\n', Et0);
p=base; p.RunKtr=1; p=getParams(p,p.g,true);
Et1=evaluateBakersExp([],p,false,false);
fprintf('RunKtr=1, fn unchanged            : Et = %.4f   (feature present but unscored)\n', Et1);
for tau=[base.tau_reg 2e-4]
    p=base; p.RunKtr=1; p.tau_reg=tau; p.fn=[base.fn, {'ktrP_k63|1'}]; p=getParams(p,p.g,true);
    [Et,~,fc,fn]=evaluateBakersExp([],p,false,false);
    i=find(strcmp(fn,'ktrP_k63'));
    fprintf('RunKtr=1 + ktrP_k63 scored, tau=%7.4g ms : Et = %8.4f   (ktrP_k63 term %.4f)\n', tau*1e3, Et, fc(i));
end
