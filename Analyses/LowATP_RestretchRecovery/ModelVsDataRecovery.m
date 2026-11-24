% ModelVsDataRecovery.m
% Same post-restretch window analysis as AnalyzeRecoveryData.m, but MODEL vs DATA at both ATP.
% Question: does the model exaggerate the low-ATP post-restretch valley depth and slow rise?
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(root,'Analyses','LowATP_RestretchRecovery','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

base=getParams(loadParams('lowATP_TNF_best'),[],true,false);
base.RunForceVelocity=false; base.RunKtr=false; base.RunSlack=true; base.RunStairs=false;
base.RunForceVelocityTime=false; base.EvalFeatures=true; base.BreakOnODEUnstable=false;
base.RunSlackSegments='All'; base.MaxRunTime=150;

cond(1)=struct('nm','2 mM','vt','protocol_03_27_2026_2mM_slack.mat','atp',2,'pi',0.8);
cond(2)=struct('nm','8 mM','vt','protocol_03_27_2026_8mM_slack.mat','atp',8,'pi',0);

Res=struct();
for c=1:2
  p=base; p.velocitytableonfile=cond(c).vt; p.MgATP=cond(c).atp; p.MgADP=0; p.Pi=cond(c).pi;
  [~,out,~,~]=runSlackExperiment(p);
  vt=load(fullfile(root,'data',cond(c).vt)).velocitytable;
  Res(c).nm=cond(c).nm;
  Res(c).data=windows(out.datatable(:,1), out.datatable(:,3), vt);
  Res(c).model=windows(out.t(:),          out.Force(:),       vt);
  Res(c).out=out;
end

fprintf('\n============ POST-RESTRETCH RECOVERY: MODEL vs DATA ============\n');
fprintf('%-6s %-6s | %-30s | %-22s\n','cond','src','valley (abs / %%steady / drop)','tau_ms  (post-slack tau)');
for c=1:2
  for s={'data','model'}
    W=Res(c).(s{1});
    v=mean([W.cyc(1:4).valley]); st=mean([W.cyc(1:4).steady]);
    tr=1000*mean([W.cyc(1:4).tauR]); ts=1000*mean([W.cyc(1:4).tauS]);
    fprintf('%-6s %-6s |  %6.1f   %5.1f%%   drop %5.1f     |  %6.1f   (%6.1f)   ratio %.2fx\n', ...
      Res(c).nm, s{1}, v, 100*v/st, st-v, tr, ts, tr/ts);
  end
end

fprintf('\n--- per-cycle valley depth (drop below steady) ---\n');
fprintf('%-6s %-6s | %s\n','cond','src','cycles 1-5');
for c=1:2; for s={'data','model'}
  W=Res(c).(s{1}); d=[W.cyc.steady]-[W.cyc.valley];
  fprintf('%-6s %-6s | %s\n', Res(c).nm, s{1}, num2str(d,'%7.1f'));
end; end

fprintf('\n--- per-cycle post-restretch tau (ms) ---\n');
for c=1:2; for s={'data','model'}
  W=Res(c).(s{1});
  fprintf('%-6s %-6s | %s\n', Res(c).nm, s{1}, num2str(1000*[W.cyc.tauR],'%7.1f'));
end; end

% ================= FIGURE =================
fig=figure('Position',[40 40 1500 850],'Color','w');
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for c=1:2
  nexttile(tl); hold on; box on; grid on;
  for k=1:4
    d=Res(c).data.cyc(k); m=Res(c).model.cyc(k);
    plot(d.t-d.t(1), d.F,'k-','LineWidth',1.0,'HandleVisibility',tern(k==1,'on','off'),'DisplayName','data');
    plot(m.t-m.t(1), m.F,'r-','LineWidth',1.4,'HandleVisibility',tern(k==1,'on','off'),'DisplayName','model');
  end
  xlabel('time since valley (s)'); ylabel('Force (rel.)');
  title(sprintf('%s — post-restretch recovery (cycles 1-4)', Res(c).nm)); legend('Location','southeast');
end
% valley depth + tau bars
nexttile(tl); hold on; box on; grid on;
vd=zeros(2,2); for c=1:2; vd(c,1)=mean([Res(c).data.cyc(1:4).steady]-[Res(c).data.cyc(1:4).valley]);
                          vd(c,2)=mean([Res(c).model.cyc(1:4).steady]-[Res(c).model.cyc(1:4).valley]); end
b=bar(vd); b(1).FaceColor=[0.2 0.2 0.2]; b(2).FaceColor=[0.85 0.2 0.2];
set(gca,'XTickLabel',{Res.nm}); ylabel('valley drop below steady'); legend({'data','model'},'Location','northwest');
title('Valley depth after restretch');
nexttile(tl); hold on; box on; grid on;
td=zeros(2,2); for c=1:2; td(c,1)=1000*mean([Res(c).data.cyc(1:4).tauR]);
                          td(c,2)=1000*mean([Res(c).model.cyc(1:4).tauR]); end
b=bar(td); b(1).FaceColor=[0.2 0.2 0.2]; b(2).FaceColor=[0.85 0.2 0.2];
set(gca,'XTickLabel',{Res.nm}); ylabel('\tau post-restretch (ms)'); legend({'data','model'},'Location','northwest');
title('Post-restretch recovery timescale');
sgtitle('Post-restretch recovery: model vs data');
exportgraphics(fig, fullfile(resdir,'recovery_model_vs_data.png'),'Resolution',140);
save(fullfile(resdir,'recovery_model_vs_data.mat'),'Res');
disp('MODEL-VS-DATA RECOVERY DONE');

% ================= helpers =================
function W=windows(t,F,vt)
  iS=find(vt(:,2)<-1); iR=find(vt(:,2)>1); n=min(numel(iS),numel(iR));
  for k=1:n
    ta=vt(iS(k)+1,1); tb=vt(iR(k),1); m=t>=ta&t<=tb;
    [tauS,~,~,~]=fitExpRise(t(m),F(m));
    tc=vt(iR(k)+1,1); if k<numel(iS); td=vt(iS(k+1),1); else; td=max(t); end
    m2=t>=tc&t<=td; tt=t(m2); FF=F(m2);
    mv=tt<=tt(1)+0.08; [vmin,iv]=min(FF(mv)); tv=tt(iv);
    mf=tt>=tv; [tauR,Finf,F0,r2]=fitExpRise(tt(mf),FF(mf));
    st=median(FF(tt>=td-0.05));           % steady just before the next event
    W.cyc(k)=struct('tauS',tauS,'tauR',tauR,'r2',r2,'valley',vmin,'steady',st, ...
                    't',tt(mf),'F',FF(mf),'Finf',Finf,'F0',F0);
  end
end
function [tau,Finf,F0,r2]=fitExpRise(t,F)
  t=t(:)-t(1); F=F(:); ok=isfinite(t)&isfinite(F); t=t(ok); F=F(ok);
  if numel(t)<10; tau=NaN;Finf=NaN;F0=NaN;r2=NaN; return; end
  obj=@(p) sum((F-(p(2)-(p(2)-p(3)).*exp(-t./max(p(1),1e-4)))).^2);
  p=fminsearch(obj,[max(t)/3,max(F),F(1)],optimset('Display','off','MaxFunEvals',4000,'MaxIter',4000));
  tau=max(p(1),1e-4); Finf=p(2); F0=p(3);
  pred=Finf-(Finf-F0).*exp(-t./tau); r2=1-sum((F-pred).^2)/max(sum((F-mean(F)).^2),eps);
end
function o=tern(c,a,b); if c; o=a; else; o=b; end; end
