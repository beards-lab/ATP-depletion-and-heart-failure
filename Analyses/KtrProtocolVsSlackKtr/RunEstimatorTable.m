% RunEstimatorTable.m
% (A) The same three estimators applied to BOTH protocols for data / incumbent
%     / fast-gate, so the model-vs-data verdict cannot be an estimator artefact.
% (B) Does a compensating ka lift the fast gate's absolute rate back without
%     re-opening the protocol split?
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

base=getParams(loadParams('estiff125_opt'),[],true,false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=1; base.BreakOnODEUnstable=0;
svt=load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
kvt=load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;
iS=find(svt(:,2)<-1);
Sd=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
Kd=load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));

%% ---------- (A) estimator x protocol ----------
FisoD=median(Sd.datatable(Sd.datatable(:,1)>svt(iS(1),1)-0.30 & Sd.datatable(:,1)<svt(iS(1),1)-0.01,3));
rowD_s=recoveryWindows(Sd.datatable(:,1),Sd.datatable(:,3),svt,'slack',FisoD);
rowD_s=rowD_s(strcmp({rowD_s.type},'postSlack')); rowD_s=rowD_s(2);
rowD_k=recoveryWindows(Kd.datatable(:,1),Kd.datatable(:,3),kvt,'ktr',1.0);

M=struct('name',{},'s',{},'k',{});
for v=1:2
    p=base; if v==2; p.tau_reg=2e-4; end; p=getParams(p,p.g,true);
    [~,os]=runSlackExperiment(p);
    FisoM=median(os.Force(os.t>=svt(iS(2),1)-0.05 & os.t<=svt(iS(2),1)-0.005));
    w=recoveryWindows(os.t,os.Force,svt,'slack',FisoM); w=w(strcmp({w.type},'postSlack'));
    pk=p; pk.SL0=2*kvt(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=kvt(:,2);
    [~,ok]=evaluateModel(str2func(pk.modelFcn),kvt(:,1),pk);
    FisoK=interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
    wkk=recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK);
    M(v).name=sprintf('model tau_reg=%.4g ms',p.tau_reg*1e3);
    M(v).s=w(2); M(v).k=wkk(1); M(v).os=os; M(v).ok=ok; M(v).FisoM=FisoM; M(v).FisoK=FisoK;
end

fprintf('\n============ ESTIMATOR x PROTOCOL (8 mM, 03/27) ============\n');
fprintf('postSlack cycle 2 (SL 2.00 um) and the ktr protocol (SL 2.00 um)\n');
fprintf('%-26s %8s %8s %8s %8s\n','','k63','kfixA','kfree','kfree58');
E={'k63','kfixA','kfit','kfit58'};
pr=@(lab,w) fprintf('%-26s %8.1f %8.1f %8.1f %8.1f\n',lab,w.k63,w.kfixA,w.kfit,w.kfit58);
pr('DATA  postSlack2',rowD_s);      pr('DATA  ktr protocol',rowD_k(1));
fprintf('%-26s %8.2f %8.2f %8.2f %8.2f  <= DATA ratio\n','',rowD_k(1).k63/rowD_s.k63, ...
    rowD_k(1).kfixA/rowD_s.kfixA, rowD_k(1).kfit/rowD_s.kfit, rowD_k(1).kfit58/rowD_s.kfit58);
for v=1:2
    pr([M(v).name ' | slack2'],M(v).s); pr([M(v).name ' | ktr'],M(v).k);
    fprintf('%-26s %8.2f %8.2f %8.2f %8.2f  <= ratio\n','',M(v).k.k63/M(v).s.k63, ...
        M(v).k.kfixA/M(v).s.kfixA, M(v).k.kfit/M(v).s.kfit, M(v).k.kfit58/M(v).s.kfit58);
end
fprintf('\nsingle-exponential r2 (how well ONE rate describes the window):\n');
fprintf('%-26s slack2 %.3f   ktr %.3f\n','DATA',rowD_s.r2,rowD_k(1).r2);
for v=1:2; fprintf('%-26s slack2 %.3f   ktr %.3f\n',M(v).name,M(v).s.r2,M(v).k.r2); end

%% ---------- (B) compensating ka at the fast gate ----------
fprintf('\n============ COMPENSATING ka AT THE FAST GATE ============\n');
KA=[1 1.10 1.20 1.30 1.45];
C=struct('xka',{},'ktrFeat',{},'kS2',{},'kK',{},'Et',{},'Ektr',{},'Fiso',{},'fc',{});
for i=1:numel(KA)
    p=base; p.tau_reg=2e-4; p.ka=base.ka*KA(i); p=getParams(p,p.g,true);
    [Et,~,fc,fn]=evaluateBakersExp([],p,false,false);
    [~,os,fm,fd]=runSlackExperiment(p);
    FisoM=median(os.Force(os.t>=svt(iS(2),1)-0.05 & os.t<=svt(iS(2),1)-0.005));
    w=recoveryWindows(os.t,os.Force,svt,'slack',FisoM); w=w(strcmp({w.type},'postSlack'));
    pk=p; pk.SL0=2*kvt(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=kvt(:,2);
    [~,ok]=evaluateModel(str2func(pk.modelFcn),kvt(:,1),pk);
    FisoK=interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
    wk=recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK);
    Fi=interp1(ok.t,ok.Force/FisoK,Kd.datatable(:,1),'linear','extrap');
    Ektr=mean((Kd.datatable(:,3)-Fi).^2)*10;
    C(i)=struct('xka',KA(i),'ktrFeat',fm.ktr(2),'kS2',w(2).k63,'kK',wk(1).k63,'Et',Et, ...
                'Ektr',Ektr,'Fiso',[FisoM FisoK],'fc',fc);
    fprintf('ka x%.2f  slackKtrFeat=%5.1f (data %5.1f)  kS2=%5.1f  kKtr=%5.1f  ratio=%4.2f  Fiso=%5.1f  Et=%6.3f  Ektr=%6.4f\n', ...
        KA(i), fm.ktr(2), fd.ktr(2), w(2).k63, wk(1).k63, wk(1).k63/w(2).k63, FisoM, Et, Ektr);
end
save(fullfile(resdir,'estimator_table.mat'),'M','rowD_s','rowD_k','C','fn','-v7.3');
