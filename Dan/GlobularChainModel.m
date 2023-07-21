clear
figure(1);clf;
colors = colormap(lines(5)) ;
% load ..\Data\bakers_passiveStretch_20ms.mat
% datatable5 = datatable;
% load ..\Data\bakers_passiveStretch_100ms.mat
% datatable4 = datatable;
% load ..\Data\bakers_passiveStretch_1000ms.mat
% datatable3 = datatable;
% load ..\Data\bakers_passiveStretch_10000ms.mat
% datatable2 = datatable;
% load ..\Data\bakers_passiveStretch_100000ms.mat
% datatable1 = datatable;
% 
% 
% figure(1); clf; axes('position',[0.15 0.25 0.8 0.25]); hold on; box on;
% plot(datatable1(:,1),datatable1(:,3),'bo-','linewidth',1);
% set(gca,'Xtick',[],'Ytick',[])
% 
% figure(2); clf; axes('position',[0.15 0.25 0.8 0.25]); hold on; box on;
% plot(datatable1(:,1),datatable1(:,3),'bo-','linewidth',1);
% set(gca,'Xtick',[],'Ytick',[])
% 
% figure(3); clf; axes('position',[0.15 0.25 0.8 0.25]); hold on; box on;
% plot(datatable1(:,1),datatable1(:,3),'bo-','linewidth',1);
% set(gca,'Xtick',[],'Ytick',[])
% 
% figure(4); clf; axes('position',[0.15 0.25 0.8 0.25]); hold on; box on;
% plot(datatable1(:,1),datatable1(:,3),'bo-','linewidth',1);
% set(gca,'Xtick',[],'Ytick',[])
% 
% figure(5); clf; axes('position',[0.15 0.25 0.8 0.25]); hold on; box on;
% plot(datatable1(:,1),datatable1(:,3),'bo-','linewidth',1);
% set(gca,'Xtick',[],'Ytick',[])

ks     = 15;
del_U  = 0.02; % (vary del_U and alphaL together)
alphaL = 0.10*3000; % length-dependent unfolding rate

Lmax = 0.4;
N = ceil(Lmax/del_U);
P = zeros(N,1);
x0 = [P; 0];

Vlist = [1 10 100 1000 5000]*Lmax/100; %  half-sarcomere velocity

opts = odeset('MaxStep',0.1);

for i = 1:5
  V = Vlist(i);
  Tend_ramp = Lmax/V; % length of ramp
  [t1,x1] = ode15s(@dPdT,2+[0 Tend_ramp],x0,opts,N,V,del_U,alphaL);
  [t2,x2] = ode15s(@dPdT,2+[Tend_ramp 198],x1(end,:),opts,N,0,del_U,alphaL);

  t = [t1; t2];
  x = [x1; x2];
  L = x(:,end);

  % forces
  PF    = 1 - sum(x(:,1:N)')'; % probability that whole chain is folded
  T     = ks*PF.*L.^1.35; % force due to 100% folded fraction

  for j = 1:length(t)
    T(j) = T(j) + ks*sum(x(j,1:N).* max( 0,L(j) - del_U*(1:N) ).^1.35 );
  end
  T = T + 150*L.^4;

  % viscous force to help get peak in fast ramps
  T(1:length(t1)) = T(1:length(t1)) + 0.3*V;
  t = [t(1); t];
  T = [0; T];

  PeakModel(i,:) = [0.4/V max(T)];

  datatable = load(['../data/bakers_passiveStretch_' num2str(Tend_ramp*1000) 'ms.mat']).datatable;

  % subplot(5, 1, i);cla;
  semilogx(datatable(:,1),datatable(:,3), 'o-','Color', colors(i, :), 'Linewidth', 1);
  hold on;
  semilogx(t,T, '--', 'Color', max(colors(i, :)*0.8, 0), 'LineWidth',3);

  axis([0 202 0 16])
  ylabel('Stress (kPa)')
  set(gca,'Xtick',0:50:200)
  set(gca,'Ytick',0:5:15)
  % set(gca,'Xticklabel',[])
  set(gca,'Fontsize',14)

end

% figure(5); set(gca,'Xticklabel',0:50:200)
xlabel('time (sec.)')


PeakData =[
0.02	14.1969	3.926472727
0.1	    10.6611	3.862672727
1	    7.94194	3.93316
10	    5.9797	3.8093
100	    4.772521951	3.826958537];
PeakModel


figure(6); semilogx(PeakData(:,1),PeakData(:,2),'o',PeakModel(:,1),PeakModel(:,2),'r-','LineWidth',2)
ylabel('Peak stress (kPa)')
xlabel('Ramp time (sec.)')
set(gca,'Fontsize',14)
legend('data','model')
% 
% 
% figure(7); semilogy(t,T - 125*0.4^4)
% title('semilog Y plot tension vs. time')
% 
% figure(8); loglog(t,T - 125*0.4^4)
% title('log log plot tension vs. time')
% 
% figure(2); axis([0 202 0 12])
% ylabel('Stress (kPa)')
% xlabel('time (sec.)')
% set(gca,'Fontsize',16)



