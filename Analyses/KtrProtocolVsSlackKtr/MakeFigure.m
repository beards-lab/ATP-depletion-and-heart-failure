% MakeFigure.m -- the model-free verdict (trace RMSE, no rate estimator) plus
% the four-panel summary figure.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
Lv=load(fullfile(resdir,'estimator_table.mat'));  M=Lv.M;
Lt=load(fullfile(resdir,'tauscan.mat'));          S=Lt.S;

base=getParams(loadParams('estiff125_opt'),[],true,false);
svt=load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
kvt=load(fullfile(root,'data',base.ktr_velocitytableonfile)).velocitytable;
iS=find(svt(:,2)<-1);
Sd=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
Kd=load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
FisoD=median(Sd.datatable(Sd.datatable(:,1)>svt(iS(1),1)-0.30 & Sd.datatable(:,1)<svt(iS(1),1)-0.01,3));

% ---- model-free: RMSE of the normalised trace over each recovery window ----
tSlk=svt(iS(2)+1,1); tSlkEnd=svt(find(svt(:,2)>1,1)+ (0),1);
iR=find(svt(:,2)>1); tSlkEnd=svt(iR(2),1);
tK0=kvt(end-1,1);    tK1=kvt(end,1);
fprintf('\n=========== MODEL-FREE TRACE RMSE (no rate estimator) ===========\n');
fprintf('normalised by each trace''s own F_iso; window = the recovery itself\n');
fprintf('%-24s %12s %12s\n','','postSlack2','ktr protocol');
rms_ = @(a,b) sqrt(mean((a-b).^2));
% data, normalised
mdS = Sd.datatable(:,1)>=tSlk & Sd.datatable(:,1)<=tSlkEnd;
tdS = Sd.datatable(mdS,1); ydS = Sd.datatable(mdS,3)/FisoD;
mdK = Kd.datatable(:,1)>=tK0 & Kd.datatable(:,1)<=tK1;
tdK = Kd.datatable(mdK,1); ydK = Kd.datatable(mdK,3);
RM=nan(2,2);
for v=1:2
    ymS = interp1(M(v).os.t, M(v).os.Force/M(v).FisoM, tdS,'linear','extrap');
    ymK = interp1(M(v).ok.t, M(v).ok.Force/M(v).FisoK, tdK,'linear','extrap');
    RM(v,:)=[rms_(ymS,ydS), rms_(ymK,ydK)];
    fprintf('%-24s %12.4f %12.4f\n', M(v).name, RM(v,1), RM(v,2));
end
fprintf('%-24s %11.1f%% %11.1f%%\n','change (fast vs incumbent)', ...
    100*(RM(2,1)/RM(1,1)-1), 100*(RM(2,2)/RM(1,2)-1));

% =================== FIGURE ===================
f=figure('Position',[60 60 1350 860],'Color','w');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
C0=[0.85 0.33 0.10]; CF=[0 0.45 0.74]; CD=[0 0 0];

% --- A: ktr protocol ---
nexttile; hold on; box on;
plot((tdK-tK0)*1e3, ydK,'-','Color',CD,'LineWidth',2.2);
for v=1:2
    tt=M(v).ok.t; m=tt>=tK0-0.018 & tt<=tK1;
    plot((tt(m)-tK0)*1e3, M(v).ok.Force(m)/M(v).FisoK,'-', ...
        'Color',(v==1)*C0+(v==2)*CF,'LineWidth',1.8);
end
xlim([-18 59]); ylim([-0.15 1.15]); grid on
xlabel('t since end of the final release (ms)'); ylabel('F / F_{iso}');
title('A   ktr protocol  (HELD OUT of the objective)');
legend({'data 03/27','incumbent \tau_{reg}=26.6 ms','fast gate \tau_{reg}=0.2 ms'}, ...
    'Location','southeast','Box','off');
text(2,1.07,sprintf('k_{63}:  data 49.4   incumbent 83.2   fast 48.6 s^{-1}'),'FontSize',9);

% --- B: post-slack #2 ---
nexttile; hold on; box on;
plot((tdS-tSlk)*1e3, ydS,'-','Color',CD,'LineWidth',2.2);
for v=1:2
    tt=M(v).os.t; m=tt>=tSlk-0.005 & tt<=tSlkEnd;
    plot((tt(m)-tSlk)*1e3, M(v).os.Force(m)/M(v).FisoM,'-', ...
        'Color',(v==1)*C0+(v==2)*CF,'LineWidth',1.8);
end
xlim([-5 150]); grid on
xlabel('t since end of the release ramp (ms)'); ylabel('F / F_{iso}');
title('B   post-slack #2  (SL 2.00 \mum -- same length, IS in the objective)');
text(40,0.25,sprintf('k_{63}:  data 50.0   incumbent 55.2   fast 46.6 s^{-1}'),'FontSize',9);

% --- C: the gate ---
nexttile; hold on; box on;
p=M(1).os.params; iA=p.NumberOfStates*p.ss+8;
for v=1:2
    tt=M(v).ok.t; m=tt>=kvt(2,1)-0.004 & tt<=tK1;
    plot((tt(m)-tK0)*1e3, M(v).ok.PU(m,iA),'-','Color',(v==1)*C0+(v==2)*CF,'LineWidth',1.8);
    tt=M(v).os.t; m=tt>=svt(iS(2),1)-0.004 & tt<=tSlk+0.059;
    plot((tt(m)-tSlk)*1e3, M(v).os.PU(m,iA),'--','Color',(v==1)*C0+(v==2)*CF,'LineWidth',1.8);
end
yline(p.A0,':k'); text(35,p.A0+0.02,'A_0 (isometric availability)','FontSize',9);
xlim([-18 59]); grid on
xlabel('t since the recovery starts (ms)'); ylabel('A_{reg}  (attachment gate)');
title('C   why: the gate remembers how long the protocol shortened');
legend({'incumbent ktr','incumbent slack','fast ktr','fast slack'},'Location','northeast','Box','off');

% --- D: identifiability ---
nexttile; hold on; box on;
tau=[S.tau]*1e3;
yyaxis left
plot(tau,[S.Et],'o-','LineWidth',1.8); ylabel('E_t  (the fitted feature objective)');
ylim([0.85 1.15]);
yyaxis right
plot(tau,[S.kK]./arrayfun(@(s) s.kS(2),S),'s-','LineWidth',1.8);
yline(0.99,':k'); text(0.12,1.03,'data ratio 0.99','FontSize',9);
ylabel('k_{ktr} / k_{postSlack2}');
set(gca,'XScale','log'); xlabel('\tau_{reg}  (ms)'); grid on
title('D   \tau_{reg} is invisible to the objective, decisive for the protocol split');
xline(26.6487,'-','incumbent','Color',C0,'LabelVerticalAlignment','bottom');

exportgraphics(f, fullfile(resdir,'ktr_vs_slack.png'),'Resolution',150);
fprintf('\nwrote %s\n', fullfile(resdir,'ktr_vs_slack.png'));
save(fullfile(resdir,'tracerms.mat'),'RM');
