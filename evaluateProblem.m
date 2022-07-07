function Etot = evaluateProblem(fcn, g, drawPlots)

LoadData;

%% force x velocity
for k = [1 2 3]
    for j = 1:length(vel)
        F_active(j,k) = evaluateModel(fcn, vel(j), 1, MgATP(k),Pi,MgADP,g);
    end
end



E(1) = sum(abs(F_active(:,1)-Data_ATP(:,2)).^2) + ...
     sum(abs(F_active(:,3)-Data_ATP(:,4)).^2);

%% force x iso MgATP
for k = 1:length(MgATP_iso)
  % Zero velocity:
  F_iso(k) = evaluateModel(fcn, 0, 1, MgATP_iso(k),Pi,MgADP,g);
end


E(2) = sum(abs(F_iso-F_data).^2);
%% time constant of MgATP

Tspan = [0:0.001:0.12];

for k = 1:length(MgATP)

  % Tspan array returns F active array
  F_active_ktr = evaluateModel(fcn,0, Tspan, MgATP(k),Pi,MgADP,g);
  
  Frel(k, :) = F_active_ktr./F_active_ktr(end);
  
  % get the time constant
  Ktr(k) = 1/interp1(Frel(k, :),Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)
end

E(3) = sum(abs(Ktr-Ktr_mean).^2);
%% Return
Etot = sum(E);

%% plot?

if ~drawPlots
    return;
end

%% Force velocity
figure(101); clf; axes('position',[0.1 0.6 0.35 0.35]); hold on;
plot(F_active(:,1), -vel./ML,'b-','linewidth',1.5);
plot(F_active(:,2), -vel./ML,'g-','linewidth',0.5);
plot(F_active(:,3), -vel./ML,'r-','linewidth',1.5);
ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
xlabel('Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
axis([0 65 0 6]);
box on;
plot(Data_ATP(:,2),Data_ATP(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data_ATP(:,3),Data_ATP(:,1),'go','linewidth',1.5,'Markersize',4,'markerfacecolor',[1 1 1]);
plot(Data_ATP(:,4),Data_ATP(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
ldg = legend('8','4','2 mM');
title(ldg,'[MgATP]');

% MgATP
% figure(104); clf; 
axes('position',[0.6 0.6 0.35 0.35]); 
semilogx(MgATP_iso,F_iso,'b-','linewidth',1.5); hold on;
semilogx(MgATP_iso,F_data,'ko','linewidth',2,'markersize',8,'markerfacecolor',[1 1 1]);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',16);
ylabel('Max. Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 100]); 
box on;

% KTR
% figure(106); clf; 
axes('position',[0.1 0.13 0.85 0.33]); hold on;
plot(Tspan,Frel(3, :),'r-','linewidth',1.5);
plot(Tspan,Frel(2, :),'g-','linewidth',1.5);
plot(Tspan,Frel(1, :),'b-','linewidth',1.5);
% plot(Tspan,F_active,'linewidth',1.5);
xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
ylabel('Force (rel.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14,'ylim',[0 1.1]);  box on;
ldg = legend('2','4','8 mM','location','northwest');
title(ldg,'[MgATP]');
axes('position',[0.6 0.2 0.3 0.2]); hold on;
plot(MgATP,Ktr,'k-','linewidth',1.5);
errorbar(MgATP,Ktr_mean,Ktr_err,'ko','linewidth',1.5,'markersize',6);
xlabel('[MgATP] (mM)','interpreter','latex','fontsize',6);
ylabel('$K_{tr}$ (sec.$^{-1}$)','interpreter','latex','fontsize',6);
set(gca,'fontsize',9,'xlim',[0 10],'ylim',[0 45]); box on;
