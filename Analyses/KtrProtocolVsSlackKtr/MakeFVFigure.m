% MakeFVFigure.m -- FV of the two gates on the fine grid, plus the low-velocity
% hump the scored 5-point grid cannot see.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
A=load(fullfile(resdir,'fv_finegrid.mat'));  R=A.R;
B=load(fullfile(resdir,'fv_lowv.mat'));      L=B.L; vlow=-B.VLOW;
C0=[0.85 0.33 0.10]; CF=[0 0.45 0.74]; CG=[0.47 0.67 0.19];

f=figure('Position',[60 60 1350 460],'Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% --- A: full FV ---
nexttile; hold on; box on; grid on
plot(A.AllV, A.D1/A.D1(1),'ks-','MarkerFaceColor','k','LineWidth',1.4,'MarkerSize',6);
plot(A.AllV, A.D2/A.D2(1),'k^--','LineWidth',1.2,'MarkerSize',6);
plot(R(1).v,R(1).fn,'-','Color',C0,'LineWidth',2);
plot(R(2).v,R(2).fn,'-','Color',CF,'LineWidth',2);
plot(A.SCORED, interp1(R(2).v,R(2).fn,A.SCORED),'o','Color',[.3 .3 .3], ...
     'MarkerSize',11,'LineWidth',1.4);
xlabel('shortening velocity (ML/s)'); ylabel('F / F_{iso}'); xlim([0 6.2]); ylim([0 1.3])
title('A   force-velocity, 0.25 ML/s grid');
legend({'data Baker2022','data Isovel2021','incumbent \tau_{reg}=26.6 ms', ...
        'fast gate \tau_{reg}=0.2 ms','the 5 SCORED velocities'},'Box','off','Location','northeast');

% --- B: the unscored low-velocity branch ---
nexttile; hold on; box on; grid on
yline(1,':k');
plot(vlow,L(1).fn,'-','Color',C0,'LineWidth',2.4);
plot(vlow,L(2).fn,'--','Color',CF,'LineWidth',2);
plot(vlow,L(3).fn,'-','Color',CG,'LineWidth',1.6);
plot([0 0.5],[1 A.D1(2)/A.D1(1)],'ks','MarkerFaceColor','k','MarkerSize',8);
plot([0 0.5],[1 A.D2(2)/A.D2(1)],'k^','MarkerSize',8);
xlabel('shortening velocity (ML/s)'); ylabel('F / F_{iso}'); xlim([0 1]); ylim([0.3 1.28])
title('B   between the two lowest scored points — nothing constrains it');
legend({'','incumbent','fast gate','A_0 = 1 control (gate velocity-blind)', ...
        'data Baker2022','data Isovel2021'},'Box','off','Location','southwest');
text(0.12,1.245,sprintf('+21.1%% at 0.10 ML/s'),'FontSize',10,'Color',C0);
text(0.12,1.185,'both gates identical here','FontSize',9,'Color',[.3 .3 .3]);

% --- C: normalised power ---
nexttile; hold on; box on; grid on
ip=A.AllV>0;
plot(A.AllV(ip), A.D1(ip)/A.D1(1).*A.AllV(ip),'ks-','MarkerFaceColor','k','LineWidth',1.4,'MarkerSize',6);
plot(A.AllV(ip), A.D2(ip)/A.D2(1).*A.AllV(ip),'k^--','LineWidth',1.2,'MarkerSize',6);
plot(R(1).v,R(1).pw,'-','Color',C0,'LineWidth',2);
plot(R(2).v,R(2).pw,'-','Color',CF,'LineWidth',2);
xlabel('shortening velocity (ML/s)'); ylabel('normalised power  (F/F_{iso})\cdotv');
xlim([0 6.2]); title('C   normalised power');
text(1.6,0.20,sprintf(['peak  incumbent 0.647 @1.25\n         fast 0.654 @1.25\n' ...
     '         data 0.664 @1.00 / 0.634 @1.00']),'FontSize',9);

exportgraphics(f, fullfile(resdir,'fv_finegrid.png'),'Resolution',150);
fprintf('wrote %s\n', fullfile(resdir,'fv_finegrid.png'));
