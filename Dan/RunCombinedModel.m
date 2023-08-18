clear Force
clear Time
clear Length

rds = fliplr([0.02, 0.1, 1, 10, 100])*1000;
for i_rd = 1:length(rds)
    load(['..\Data\bakers_passiveStretch_' num2str(rds(i_rd)) 'ms.mat']);
    datatables{i_rd} = datatable;
end

Lmax = 0.4;
Ls0  = 0.10*mod(13);
Nx   = 25;          % number of space steps
ds   = (0.36-Ls0)/(Nx-1);      % space step size
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng  = 20;            % number of glubules on globular chain
delU = 0.0125*mod(1);

kA   = 0;
kD   = 1;
kC   = 103.33*mod(2);
kS   = 300*mod(3);         % series element spring constant
alphaU = 2000*mod(4);       % chain unfolding rate constant
alphaF = 1*mod(5);
nC = 1.77*mod(6);
nS = 2.56*mod(7);
nU = 4*mod(8);
mu = 2.44*mod(9); 

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
% PA = zeros(1,Ng+1); % initial unfolded probabilities for attached rectifier state
% PA = []; % initial unfolded probabilities for attached rectifier state

pu = zeros(Nx,1)*PU;
% pa = zeros(Nx,1)*PA;
pu(1,1) = 1/ds; 
% x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
x0 = reshape(pu,[(Ng+1)*Nx,1]);
x0 = [x0; 0]; 

Vlist = [1 10 100 1000 5000]*Lmax/100; %  half-sarcomere velocity
% reducing number of ranges
% Vlist = [10 100 1000]*Lmax/100; %  half-sarcomere velocity
% tic
Force = cell(1, 5); 
Time = cell(1, 5); 
Length = cell(1, 5); 
rampSet = [2 3 4 5];
for j = rampSet
  % tic
  V = Vlist(j);

  pu = zeros(Nx,1)*PU;
  % pa = zeros(Nx,1)*PA;
  pu(1,1) = 1/ds; 
  % x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  x0 = reshape(pu,[(Ng+1)*Nx,1]);
  x0 = [x0; 0]; 
  Tend_ramp = Lmax/V; % length of ramp

  opts = odeset('RelTol',1e-3, 'AbsTol',1e-2);
  % testing tolerances, all with +/- same total cost
  % abstol 1e-6 7.7s,  1e-3 3.6s, 1e-1 2.7s and 1e1 3.1s
  % reltol 1e-3 (normal) 7.7s, 1e-1 6.8s and 8s for 1e-5 
    
  % test solvers - do not touch ode15s
  % ode15s 7.5s, ode45 53s, ode23s 173s, ode23t 9.7s

  % test sparse matrix - takes forever
  % S = [reshape(triu(ones(size(pu))),[(Ng+1)*Nx,1]); 1];
  % Sp = S*S';  % opts = odeset('JPattern',Sp');

  [t0,x0] = ode15s(@dXdT,[-100:1:0],x0,opts,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  [t1,x1] = ode15s(@dXdT,[0 Tend_ramp],x0(end,:),opts,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,V);
  
  % whole decay till the bitter end
  [t2,x2] = ode15s(@dXdT,[Tend_ramp 200],x1(end,:),opts,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  % limited decay
  % [t2,x2] = ode15s(@dXdT,[Tend_ramp min(200, Tend_ramp*4)],x1(end,:),[],Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  % only ramp up, no decay  
  % x2 = [];t2 = []; 

  t = [t1(1:end); t2(2:end)];% prevent overlap at tend_ramp
  x = [x1; x2(2:end, :)];

  for i = 1:length(t)
    xi = x(i,:);
    Length{j}(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    % pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    Force{j}(i) =  kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pu )) ;
    % + ...
                   % kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pa ));
    % Force_pa{j} = kS*ds*sum(sum( (max(0,Length{j}(i) - s - Ls0).^nS).*pa ));
  end

    
    % show state occupation at the end of the ramp
    % clf;
    % i = round(size(x2, 1)*3/4);
    % xi = x2(i,:);    
    % pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    % surf(pu)
    % xlabel('Ng');ylabel('Nx');zlabel('State probability');    
    % hold on;
    % title(['Ramp ' num2str(Tend_ramp) 's, V = ' num2str(V) ' ML/s at t ' num2str(t2(i))]);


  % add parallel static force
  % Force{j} = Force{j} + mod(10)*3;
  % Force{j} = Force{j} + 130*(Length{j}).^4;
  % Force{j} = Force{j} + 0.0045*(1.6 + 2*Length{j}).^7;

  % decay offset is set, now use a nonlinear func to fit the ramp onset
  % a*(-b + Lmax).^c + d = 1.2716*3;
  % a*(-b + Lmax).^c = 1.2716*3 - d;
  b = 0.05*mod(10);
  c = 7*mod(11);
  d = 0.01*mod(12);
  % apply constraints
  if b < 0 || c <= 0 || d < 0 
      cost = inf;
      return;
  end
  % calculate a, so that the max value is the same
  % Fss = mod(10)*3; % optimized previously
  Fss = 1.2716*3*mod(14); % reducing the param space
  a = (Fss - d)/((Lmax -b)^c);
  % calc force
  Force{j} = Force{j} + a*max(Length{j} - b, 0).^c + d; 

  Time{j} = t;
    % ttoc = toc;
  % fprintf("Ramp %0.1fs takes %1.1fs \n", Tend_ramp, ttoc)
    
end

% Get error for the whole ramp-up and decay
t_endFreeware = zeros(1, 5); % time when we start counting the costs
% alternatively, get the error from decay only
% Tend_rampset = Lmax./Vlist + 2;

%% Evaluating all ramps at once
for j = rampSet
    datatable_cur = datatables{j};
    datatable_cur = datatable_cur(datatable_cur(:, 1) >= t_endFreeware(j), :);
    t_int{j} = datatable_cur(:, 1) - 2;
    Ftot_int{j} = interp1(Time{j}, Force{j}, t_int{j}); % total force interpolated
    Es = (Ftot_int{j} - datatable_cur(:,3)).^2; % error set
    Es(isnan(Es)) = 0; % zero outside bounds
    En{j} = 1e3*sum(Es)/length(Es); % normalized error
end
%
PeakData =[
100	    4.772521951	3.826958537
10	    5.9797	3.8093
1	    7.94194	3.93316
0.1	    10.6611	3.862672727
0.02	14.1969	3.926472727];

PeakModel = zeros(1, 5);
for j = 1:5
    m = max(Force{j});
    if ~isempty(m)
        PeakModel(j) = m;
    end
end

Ep = sum((PeakData(2:4, 2) - PeakModel(2:4)').^2);

cost = Ep*100 + sum([En{1:end}], 'all');


if exist('drawPlots', 'var') && ~drawPlots
    return;
end
%%
zoomIns = [0 200 0 10;...
           0 20 0 15;...
           0 2 0 15;...
           0 .2 0 15;...
           0 .04 0 15;...
           ];
clf;
for j = rampSet
    % figure(j); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
    subplot(2, 3, j);hold on;
    % plot(datatables{j}(:,1)-2,datatables{j}(:,3),'bo','linewidth',2);
    % plot(t_int{j},Ftot_int{j},'ro','linewidth',1);
    % plot(Time{j},Force{j},'linewidth',2); axis([0 200 0 15])
    % ylabel('Stress (kPa)')
    % xlabel('time (sec.)')
    % set(gca,'Xtick',0:50:200)
    % set(gca,'Fontsize',14)
    % axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
    plot(datatables{j}(:,1)-2,datatables{j}(:,3),'bo','linewidth',2);
    plot(t_int{j},Ftot_int{j},'-ro','linewidth',1);
    % plot(Time{2},Force{2},'linewidth',2); 
    axis(zoomIns(j, :))
end
%%
% figure(6); clf; axes('position',[0.15 0.15 0.8 0.80]);
subplot(2, 3, 6);
semilogx(PeakData(:,1),PeakData(:,2),'o',PeakData(:,1),PeakModel,'r-','LineWidth',2)
ylabel('Peak stress (kPa)')
xlabel('Ramp time (sec.)')
set(gca,'Fontsize',14)
legend('data','model')
