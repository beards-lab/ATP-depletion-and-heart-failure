% RunNoiseAndProbe.m
% (A) the data's own noise floor in each window, so the trace RMSEs can be read
% (B) which feature the naive ka compensation broke
% (C) a 2-D directional probe (ka, kd) at the fast gate: can the slack
%     redevelopment rate be restored without re-opening the protocol split?
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
base=getParams(loadParams('estiff125_opt'),[],true,false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=1; base.BreakOnODEUnstable=0;
svt=load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
kvt=load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;
iS=find(svt(:,2)<-1); iR=find(svt(:,2)>1);
Sd=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
Kd=load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
FisoD=median(Sd.datatable(Sd.datatable(:,1)>svt(iS(1),1)-0.30 & Sd.datatable(:,1)<svt(iS(1),1)-0.01,3));

%% (A) noise floor = RMS residual of the DATA about its own best single exponential
fprintf('\n=========== DATA NOISE FLOOR IN EACH WINDOW ===========\n');
nf = @(t,y) localNoise(t,y);
mS=Sd.datatable(:,1)>=svt(iS(2)+1,1) & Sd.datatable(:,1)<=svt(iR(2),1);
mK=Kd.datatable(:,1)>=kvt(end-1,1)   & Kd.datatable(:,1)<=kvt(end,1);
nS=nf(Sd.datatable(mS,1), Sd.datatable(mS,3)/FisoD);
nK=nf(Kd.datatable(mK,1), Kd.datatable(mK,3));
LR=load(fullfile(resdir,'tracerms.mat')); RM=LR.RM;
fprintf('%-30s %12s %12s\n','','postSlack2','ktr');
fprintf('%-30s %12.4f %12.4f\n','data noise floor (RMS resid)',nS,nK);
fprintf('%-30s %12.4f %12.4f\n','model incumbent trace RMSE',RM(1,1),RM(1,2));
fprintf('%-30s %12.4f %12.4f\n','model fast-gate trace RMSE',RM(2,1),RM(2,2));
fprintf('%-30s %11.2fx %11.2fx\n','incumbent / noise',RM(1,1)/nS,RM(1,2)/nK);
fprintf('%-30s %11.2fx %11.2fx\n','fast gate / noise',RM(2,1)/nS,RM(2,2)/nK);

%% (B) what the naive ka x1.10 broke
Le=load(fullfile(resdir,'estimator_table.mat')); C=Le.C; fn=Le.fn(:)';
fprintf('\n=========== WHAT ka x1.10 BROKE (fast gate) ===========\n');
d=C(2).fc-C(1).fc; [~,o]=sort(abs(d),'descend');
for i=o(1:6); fprintf('%-26s %9.4f -> %9.4f  (%+.3f)\n',fn{i},C(1).fc(i),C(2).fc(i),d(i)); end

%% (C) 2-D probe: ka up, kd up (rate up at ~constant force)
fprintf('\n=========== (ka, kd) PROBE AT THE FAST GATE ===========\n');
XKA=[1.00 1.15 1.30]; XKD=[1.00 1.30 1.70];
P=struct('xka',{},'xkd',{},'ktrFeat',{},'kS2',{},'kK',{},'Et',{},'Ektr',{},'Fiso',{});
n=0;
for a=XKA
  for b=XKD
    n=n+1;
    p=base; p.tau_reg=2e-4; p.ka=base.ka*a; p.kd=base.kd*b; p=getParams(p,p.g,true);
    try
      [Et,~,fc,~]=evaluateBakersExp([],p,false,false);
      [~,os,fm,fd]=runSlackExperiment(p);
      FisoM=median(os.Force(os.t>=svt(iS(2),1)-0.05 & os.t<=svt(iS(2),1)-0.005));
      w=recoveryWindows(os.t,os.Force,svt,'slack',FisoM); w=w(strcmp({w.type},'postSlack'));
      pk=p; pk.SL0=2*kvt(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=kvt(:,2);
      [~,ok]=evaluateModel(str2func(pk.modelFcn),kvt(:,1),pk);
      FisoK=interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
      wk=recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK);
      Fi=interp1(ok.t,ok.Force/FisoK,Kd.datatable(:,1),'linear','extrap');
      Ektr=mean((Kd.datatable(:,3)-Fi).^2)*10;
      P(n)=struct('xka',a,'xkd',b,'ktrFeat',fm.ktr(2),'kS2',w(2).k63,'kK',wk(1).k63, ...
                  'Et',Et,'Ektr',Ektr,'Fiso',FisoM);
      fprintf('ka x%.2f kd x%.2f | slackKtrFeat=%5.1f (tgt %4.1f) kS2=%5.1f kKtr=%5.1f ratio=%4.2f Fiso=%5.1f Et=%8.3f Ektr=%6.4f\n', ...
          a,b,fm.ktr(2),fd.ktr(2),w(2).k63,wk(1).k63,wk(1).k63/w(2).k63,FisoM,Et,Ektr);
    catch ME
      fprintf('ka x%.2f kd x%.2f | FAILED: %s\n',a,b,ME.message);
    end
  end
end
save(fullfile(resdir,'noise_probe.mat'),'nS','nK','P');

function r = localNoise(t,y)
% RMS residual of the trace about its own best B + A(1-exp(-k(t-t0)))
    t=t(:); y=y(:); t=t-t(1); B=min(y); A0=max(y)-B;
    obj=@(q) sum((y-(B+abs(q(1)).*(1-exp(-abs(q(2)).*max(t-max(q(3),0),0))))).^2);
    q=fminsearch(obj,[A0,45,0.005],optimset('Display','off','MaxFunEvals',8e3,'MaxIter',8e3));
    pred=B+abs(q(1)).*(1-exp(-abs(q(2)).*max(t-max(q(3),0),0)));
    r=sqrt(mean((y-pred).^2));
end
