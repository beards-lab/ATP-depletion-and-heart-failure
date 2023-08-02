clear Force
clear Time
clear Length

load ..\Data\bakers_passiveStretch_20ms.mat
datatable5 = datatable;
load ..\Data\bakers_passiveStretch_100ms.mat
datatable4 = datatable;
load ..\Data\bakers_passiveStretch_1000ms.mat
datatable3 = datatable;
load ..\Data\bakers_passiveStretch_10000ms.mat
datatable2 = datatable;
load ..\Data\bakers_passiveStretch_100000ms.mat
datatable1 = datatable;

Lmax = 0.4;
Ls0  = 0.10;
Nx   = 25;          % number of space steps
ds   = (0.36-Ls0)/(Nx-1);      % space step size
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng  = 20;            % number of glubules on globular chain
delU = 0.0125;

kA   = 0;
kD   = 1;
kC   = 103.33;
kS   = 300;         % series element spring constant
alphaU = 2000;       % chain unfolding rate constant
alphaF = 1;
nC = 1.77;
nS = 2.56;
nU = 4;
mu = 2.44; 

% Calculate globular chain force Fc(s,n) for every strain and
% value. 
slack = (0:Ng).*delU;
Fc = kC*(max(0,s-slack)).^nC;

% Calculate the globular chain folding/unfolding probability transition
% rates
RU = alphaU*(max(0,s-slack(1:Ng))).^nU; % unfolding rates from state n to (n+1)
RF = alphaF*(max(0,s-slack(2:(Ng+1)))).^1;  % folding rates from state n+1 to n                 

% Initial state
PU = zeros(1,Ng+1); % initial unfolded probabilities for un-attached rectifier state
PA = zeros(1,Ng+1); % initial unfolded probabilities for attached rectifier state

pu = zeros(Nx,1)*PU;
pa = zeros(Nx,1)*PA;
pu(1,1) = 1/ds; 
x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
x0 = [x0; 0]; 

Vlist = [1 10 100 1000 5000]*Lmax/100; %  half-sarcomere velocity

for j = [1 2 3 4 5]
  V = Vlist(j)

  pu = zeros(Nx,1)*PU;
  pa = zeros(Nx,1)*PA;
  pu(1,1) = 1/ds; 
  x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  x0 = [x0; 0]; 

  opts = odeset('RelTol',1e-1, 'AbsTol',1-3);
  Tend_ramp = Lmax/V; % length of ramp
  [t0,x0] = ode15s(@dXdT,[-100:1:0],x0,[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  [t1,x1] = ode15s(@dXdT,[0 Tend_ramp],x0(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,V);
  [t2,x2] = ode15s(@dXdT,[Tend_ramp 200],x1(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
%   [t3,x3] = ode15s(@dXdT,[2*Tend_ramp:1:200],x2(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  t = [t1; t2];
  x = [x1; x2];

  for i = 1:length(t)
    xi = x(i,:);
    Length{j}(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    Force{j}(i) =  kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pu )) + ...
                   kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pa ));
    Force_pa{j} = kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pa ));
  end
toc
  % add parallel static force
  Force{j} = Force{j} + 130*Length{j}.^4;
  Time{j} = t;

end

figure(1); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable1(:,1)-2,datatable1(:,3),'bo','linewidth',2);
plot(Time{1},Force{1},'linewidth',2); axis([0 200 0 15])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:50:200)
set(gca,'Fontsize',14)

figure(2); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable2(:,1)-2,datatable2(:,3),'bo','linewidth',2);
plot(Time{2},Force{2},'linewidth',2); axis([0 200 0 15])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:50:200)
set(gca,'Fontsize',14)
axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
plot(datatable2(:,1)-2,datatable2(:,3),'bo','linewidth',2);
plot(Time{2},Force{2},'linewidth',2); axis([0 20 0 15])

figure(3); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable3(:,1)-2,datatable3(:,3),'bo','linewidth',2);
plot(Time{3},Force{3},'linewidth',2); axis([0 200 0 15])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:50:200)
set(gca,'Fontsize',14)
axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
plot(datatable3(:,1)-2,datatable3(:,3),'bo','linewidth',2);
plot(Time{3},Force{3},'linewidth',2); axis([0 2 0 15])

figure(4); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable4(:,1)-2,datatable4(:,3),'bo','linewidth',2);
plot(Time{4},Force{4},'linewidth',2); axis([0 200 0 15])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:50:200)
set(gca,'Fontsize',14)
axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
plot(datatable4(:,1)-2,datatable4(:,3),'bo','linewidth',2);
plot(Time{4},Force{4},'linewidth',2); axis([0 0.20 0 15])

figure(5); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable5(:,1)-2,datatable5(:,3),'bo','linewidth',2);
plot(Time{5},Force{5},'linewidth',2); axis([0 200 0 15])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:50:200)
set(gca,'Fontsize',14)
axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
plot(datatable5(:,1)-2,datatable5(:,3),'bo','linewidth',2);
plot(Time{5},Force{5},'linewidth',2); axis([0 0.04 0 15])

PeakData =[
100	    4.772521951	3.826958537
10	    5.9797	3.8093
1	    7.94194	3.93316
0.1	    10.6611	3.862672727
0.02	14.1969	3.926472727];

for j = 1:5
  PeakModel(j) = max(Force{j});
end

figure(6); clf; axes('position',[0.15 0.15 0.8 0.80]);
semilogx(PeakData(:,1),PeakData(:,2),'o',PeakData(:,1),PeakModel,'r-','LineWidth',2)
ylabel('Peak stress (kPa)')
xlabel('Ramp time (sec.)')
set(gca,'Fontsize',14)
legend('data','model')
