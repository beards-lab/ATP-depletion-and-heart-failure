% AnalyzeRecoveryData.m
% =========================================================================
% What does the DATA actually do after the restretch, at 2 mM vs 8 mM ATP?
% Protocol (velocitytable): [slack -> hold_short -> restretch ramp -> hold_long] x5.
% The post-restretch HOLD window (~0.28 s; 0.84 s on the last cycle) is the
% "slow rise" in question. Compare its timescale to the post-slack force
% redevelopment (ktr) measured in the SAME trace.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_RestretchRecovery/AnalyzeRecoveryData
% =========================================================================
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root));
here=fullfile(root,'Analyses','LowATP_RestretchRecovery'); resdir=fullfile(here,'results');
if ~exist(resdir,'dir'); mkdir(resdir); end

D{1}=load('data/protocol_03_27_2026_2mM_slack.mat'); nm{1}='2 mM';
D{2}=load('data/protocol_03_27_2026_8mM_slack.mat'); nm{2}='8 mM';

R=struct();
for c=1:2
  dt=D{c}.datatable; vt=D{c}.velocitytable;
  t=dt(:,1); ML=dt(:,2); F=dt(:,3);
  iS=find(vt(:,2)<-1); iR=find(vt(:,2)>1);          % slack rows / restretch rows
  n=min(numel(iS),numel(iR));
  fprintf('\n================ %s ================\n', nm{c});
  fprintf('%3s | %-28s | %-28s | %s\n','cyc','POST-SLACK redevelopment','POST-RESTRETCH recovery','windows');
  fprintf('%3s | %7s %7s %7s | %7s %7s %7s %6s | %s\n','', 'tau_ms','ktr_1/s','r2','tau_ms','k_1/s','r2','%done','');
  for k=1:n
    % --- post-slack redevelopment window: end of slack ramp -> start of restretch ramp
    ta=vt(iS(k)+1,1); tb=vt(iR(k),1);
    m=t>=ta & t<=tb; [tauS,~,~,r2S]=fitExpRise(t(m),F(m));
    % --- post-restretch recovery window: end of restretch ramp -> next slack (or trace end)
    tc=vt(iR(k)+1,1);
    if k<numel(iS); td=vt(iS(k+1),1); else; td=max(t); end
    m2=t>=tc & t<=td; tt=t(m2); FF=F(m2);
    % valley = min within the first 80 ms of the hold (after the elastic spike decays)
    mv=tt<=tt(1)+0.08; [~,iv]=min(FF(mv)); tv=tt(iv);
    mf=tt>=tv; [tauR,FinfR,F0R,r2R]=fitExpRise(tt(mf),FF(mf));
    pct=100*(FF(end)-F0R)/max(FinfR-F0R,eps);      % how much of the fitted rise is complete
    R(c).cyc(k)=struct('tauS',tauS,'tauR',tauR,'r2S',r2S,'r2R',r2R,'pct',pct, ...
      'win',[tc td],'tv',tv,'t',tt(mf),'F',FF(mf),'Finf',FinfR,'F0',F0R, ...
      'dur',td-tc,'MLhold',vt(iR(k)+1,4));
    fprintf('%3d | %7.1f %7.1f %7.3f | %7.1f %7.1f %7.3f %6.0f | rec %.3f s @ML %.4f\n', ...
      k, 1000*tauS, 1/tauS, r2S, 1000*tauR, 1/tauR, r2R, pct, td-tc, vt(iR(k)+1,4));
  end
end

% ---- headline comparison ----
fprintf('\n=== TIMESCALE COMPARISON (mean over cycles 1-4, equal-length windows) ===\n');
for c=1:2
  tS=arrayfun(@(s)s.tauS,R(c).cyc(1:4)); tR=arrayfun(@(s)s.tauR,R(c).cyc(1:4));
  fprintf('%s : post-slack tau %5.1f ms (ktr %4.1f/s) | post-restretch tau %6.1f ms (k %4.1f/s) | ratio %.2fx SLOWER\n', ...
    nm{c}, 1000*mean(tS), 1/mean(tS), 1000*mean(tR), 1/mean(tR), mean(tR)/mean(tS));
end
fprintf('\nLast cycle (long 0.84 s window, returns to SHORTER length):\n');
for c=1:2
  s=R(c).cyc(5);
  fprintf('%s : post-restretch tau %6.1f ms, %.0f%% complete at window end (F %.1f -> fitted Finf %.1f)\n', ...
    nm{c}, 1000*s.tauR, s.pct, s.F(end), s.Finf);
end

% ================= FIGURE =================
fig=figure('Position',[40 40 1500 900],'Color','w');
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
cols={[0.85 0.1 0.1],[0.1 0.3 0.85]};
% (a) post-restretch recoveries, absolute
nexttile(tl,[1 2]); hold on; box on; grid on;
for c=1:2; for k=1:5; s=R(c).cyc(k);
  plot(s.t-s.t(1), s.F, '-', 'Color',[cols{c} 0.55], 'LineWidth',1.1, ...
    'HandleVisibility', tern(k==1,'on','off'), 'DisplayName',nm{c});
end; end
xlabel('time since valley (s)'); ylabel('Force (rel.)');
title('Post-restretch recovery, all cycles — absolute'); legend('Location','southeast');
% (b) normalised to the fitted rise
nexttile(tl); hold on; box on; grid on;
for c=1:2; for k=1:5; s=R(c).cyc(k);
  y=(s.F-s.F0)/max(s.Finf-s.F0,eps);
  plot(s.t-s.t(1), y, '-','Color',[cols{c} 0.55],'LineWidth',1.1, ...
    'HandleVisibility', tern(k==1,'on','off'), 'DisplayName',nm{c});
end; end
yline(1,'k--','HandleVisibility','off'); ylim([-0.1 1.4]);
xlabel('time since valley (s)'); ylabel('fraction of fitted rise');
title('Normalised — is it still climbing?'); legend('Location','southeast');
% (c) tau comparison bars
nexttile(tl); hold on; box on; grid on;
tS=zeros(2,5); tR=zeros(2,5);
for c=1:2; tS(c,:)=1000*arrayfun(@(s)s.tauS,R(c).cyc); tR(c,:)=1000*arrayfun(@(s)s.tauR,R(c).cyc); end
b=bar([mean(tS(:,1:4),2) mean(tR(:,1:4),2)]); b(1).FaceColor=[0.6 0.6 0.6]; b(2).FaceColor=[0.85 0.3 0.1];
set(gca,'XTickLabel',nm); ylabel('\tau (ms)'); legend({'post-slack (ktr)','post-restretch'},'Location','northwest');
title('Recovery timescale: restretch vs slack');
% (d) per-cycle tau
nexttile(tl); hold on; box on; grid on;
for c=1:2
  plot(1:5, tS(c,:), 'o--','Color',cols{c},'LineWidth',1.3,'DisplayName',[nm{c} ' post-slack']);
  plot(1:5, tR(c,:), 's-','Color',cols{c},'LineWidth',1.6,'DisplayName',[nm{c} ' post-restretch']);
end
xlabel('cycle (deeper slack \rightarrow)'); ylabel('\tau (ms)'); legend('Location','northwest','FontSize',7);
title('Per-cycle timescales');
% (e) the long last-cycle window
nexttile(tl); hold on; box on; grid on;
for c=1:2; s=R(c).cyc(5);
  plot(s.t-s.t(1), s.F,'-','Color',cols{c},'LineWidth',1.5,'DisplayName',nm{c});
  plot(s.t-s.t(1), s.Finf-(s.Finf-s.F0)*exp(-(s.t-s.t(1))/s.tauR),'k:','LineWidth',1.1,'HandleVisibility','off');
end
xlabel('time since valley (s)'); ylabel('Force (rel.)');
title('Last cycle: 0.84 s window (dotted = single-exp fit)'); legend('Location','southeast');
sgtitle('Post-restretch recovery in the DATA: 2 mM vs 8 mM');
exportgraphics(fig, fullfile(resdir,'recovery_data.png'),'Resolution',140);
save(fullfile(resdir,'recovery_data.mat'),'R','nm');
disp('DATA RECOVERY ANALYSIS DONE');

% ================= helpers =================
function [tau,Finf,F0,r2]=fitExpRise(t,F)
  t=t(:)-t(1); F=F(:);
  ok=isfinite(t)&isfinite(F); t=t(ok); F=F(ok);
  if numel(t)<10; tau=NaN;Finf=NaN;F0=NaN;r2=NaN; return; end
  p0=[max(t)/3, max(F), F(1)];
  obj=@(p) sum((F-(p(2)-(p(2)-p(3)).*exp(-t./max(p(1),1e-4)))).^2);
  p=fminsearch(obj,p0,optimset('Display','off','MaxFunEvals',4000,'MaxIter',4000));
  tau=max(p(1),1e-4); Finf=p(2); F0=p(3);
  pred=Finf-(Finf-F0).*exp(-t./tau);
  r2=1-sum((F-pred).^2)/max(sum((F-mean(F)).^2),eps);
end
function o=tern(c,a,b); if c; o=a; else; o=b; end; end
