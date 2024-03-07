clear Force
clear Time
clear Length

datatable4 = readtable('Data\bakers_passiveStretch_pCa4_100ms.csv');
datatable5 = readtable('Data\bakers_passiveStretch_pCa4_1000ms.csv');
datatable6 = readtable('Data\bakers_passiveStretch_pCa4_10000ms.csv');

g0 = ones(1,11);
kA   = g0(7)*5*200;
kD   = g0(8)*50;
kC   = g0(1)*4*600 ;      % proximal chain force constant
kS   = g0(2)*4*500;        % distal chain force constant
alphaU = g0(6)*(1.0e6);         % chain unfolding rate constant
alphaF = 1;
nC = g0(3)*2;
nS = g0(9)*2;
nU = g0(4)*6;
mu = g0(5)*0.2; 
Ls0  = g0(10)*0.0;

% half-sarcomere ramp height
Lmax = 0.225;
Nx   = 25;          % number of space steps
ds   = 0.90*(Lmax)/(Nx-1);      % space step size
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng = 15; 
delU = 0.010;

% Calculate globular chain force Fc(s,n) for every strain and
% value. 
slack = (0:Ng).*delU;
Fc = kC*(max(0,s-slack-Ls0)).^nC;

% Calculate the globular chain folding/unfolding probability transition
% rates
RU = alphaU*((max(0,s-slack(1:Ng)-Ls0)).^nU).*(ones(Nx,1).*(Ng - (0:Ng-1))); % unfolding rates from state n to (n+1)
RF = alphaF*(max(0,s-slack(2:(Ng+1)))).^1;  % folding rates from state n+1 to n                 

% Initial state
PU = zeros(1,Ng+1); % initial unfolded probabilities for un-attached rectifier state
PA = zeros(1,Ng+1); % initial unfolded probabilities for attached rectifier state

pu = zeros(Nx,1)*PU;
pa = zeros(Nx,1)*PA;
pu(1,1) = 1/ds; 
x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
x0 = [x0; 0]; 

Vlist = 1000*[1/100 1/1000 1/10000]*Lmax; %  half-sarcomere velocity

for j = [1 2 3]
  V = Vlist(j)

  pu = zeros(Nx,1)*PU;
  pa = zeros(Nx,1)*PA;
  pu(1,1) = 1/ds; 
  x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  x0 = [x0; 0]; 

  Tend_ramp = Lmax/V; % length of ramp
  [t0,x0] = ode15s(@dXdT,[-100:1:0],x0,[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  [t1,x1] = ode15s(@dXdT,[0 Tend_ramp],x0(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,V);
  [t2,x2] = ode15s(@dXdT,[Tend_ramp Tend_ramp+30],x1(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
%   [t3,x3] = ode15s(@dXdT,[2*Tend_ramp:1:200],x2(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  t = [t1; t2];
  x = [x1; x2];

  for i = 1:length(t)
    xi = x(i,:);
    Length{j}(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    Force{j}(i) =  kS*ds*sum(sum( (max(0,Length{j}(i) - s - 0*Ls0).^nS).*pu )) + ...
                   kS*ds*sum(sum( (max(0,Length{j}(i) - s - 0*Ls0).^nS).*pa ));
  end

  % add parallel static force
  Force{j} = Force{j} + (0.50e6)*Length{j}.^8;
  Time{j} = t;

end

figure(4); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable4.Time-2,datatable4.F,'b','linewidth',2);
plot(Time{1},Force{1},'linewidth',2); axis([0 30 0 50])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:10:30)
set(gca,'Fontsize',16)
axes('position',[0.525 0.525 0.4 0.4]); hold on; box on;
plot(datatable4.Time-2,datatable4.F,'bo','linewidth',2);
plot(Time{1},Force{1},'linewidth',2); axis([0 0.2 0 50])

figure(5); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable5.Time(1:64100)-2,datatable5.F(1:64100),'b','linewidth',1);
plot(Time{2},Force{2},'linewidth',3); axis([0 30 0 50])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:10:30)
set(gca,'Fontsize',16)
axes('position',[0.525 0.525 0.4 0.4]); hold on; box on;
plot(datatable5.Time-2,datatable5.F,'bo','linewidth',2);
plot(Time{2},Force{2},'linewidth',3); axis([0 2 0 30])

figure(6); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
plot(datatable6.Time-2,datatable6.F,'b','linewidth',2);
plot(Time{3},Force{3},'linewidth',2); axis([0 40 0 50])
ylabel('Stress (kPa)')
xlabel('time (sec.)')
set(gca,'Xtick',0:10:40)
set(gca,'Fontsize',14)


% [max(Force{1}) max(Force{3})]
