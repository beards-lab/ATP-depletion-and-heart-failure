% RunTauScan.m -- identifiability of tau_reg.
% The ktr protocol is the only experiment that sees it; it is switched off in
% the objective. Scan tau_reg over 3 decades and put the fitted objective next
% to the held-out ktr protocol on the same axis.
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

base = getParams(loadParams('estiff125_opt'), [], true, false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=1; base.BreakOnODEUnstable=0;

svt = load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
kvt = load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;
iS  = find(svt(:,2)<-1);
Kd  = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));

TAU = [26.6487 20 15 10 7 5 3 2 1.5 1 0.7 0.5 0.3 0.2 0.1]*1e-3;
S = struct('tau',{},'kS',{},'kK',{},'Et',{},'Ektr',{},'fc',{},'Fiso',{},'secs',{});
for i=1:numel(TAU)
    tic;
    p = base; p.tau_reg = TAU(i); p = getParams(p,p.g,true);
    [~, os] = runSlackExperiment(p);
    FisoM = median(os.Force(os.t>=svt(iS(2),1)-0.05 & os.t<=svt(iS(2),1)-0.005));
    ws = recoveryWindows(os.t,os.Force,svt,'slack',FisoM); ws=ws(strcmp({ws.type},'postSlack'));
    pk = p; pk.SL0=2*kvt(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=kvt(:,2);
    [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk);
    FisoK = interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
    wk = recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK);
    Fi = interp1(ok.t,ok.Force/FisoK,Kd.datatable(:,1),'linear','extrap');
    Ektr = mean((Kd.datatable(:,3)-Fi).^2)*10;
    [Et,~,fc,fn] = evaluateBakersExp([],p,false,false);
    S(i)=struct('tau',TAU(i),'kS',[ws.k63],'kK',wk(1).k63,'Et',Et,'Ektr',Ektr,'fc',fc,'Fiso',[FisoM FisoK],'secs',toc);
    fprintf('tau=%6.2f ms  kS2=%5.1f  kKtr=%6.1f  ratio=%4.2f  Et=%6.3f  Ektr=%6.4f  (%4.1f s)\n', ...
        TAU(i)*1e3, ws(2).k63, wk(1).k63, wk(1).k63/ws(2).k63, Et, Ektr, S(i).secs);
end
fprintf('\nEt over the whole scan: min %.3f  max %.3f  spread %.3f\n', min([S.Et]), max([S.Et]), max([S.Et])-min([S.Et]));
fprintf('Ektr  over the whole scan: min %.4f  max %.4f  ratio %.1fx\n', min([S.Ektr]), max([S.Ektr]), max([S.Ektr])/min([S.Ektr]));
save(fullfile(resdir,'tauscan.mat'),'S','fn','TAU');
