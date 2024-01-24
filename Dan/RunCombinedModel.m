% Assuming mod = ones(10,1)
% assuming pCa = Inf

% pCa 11 - load Relaxed, do not run the extended, PEVK attachment model
% pCa 10 - load Relaxed, run the PEVK attachment model
% pCa < 10 - load AvgpCa dataset, Ca effect in place
clear Force
clear Time
clear Length
clear outStruct;
if any(mod < 0) 
    cost = inf;
    return;
end
drawAllStates = false;

% rds = fliplr([0.02 0.1, 1, 10 100]);
rds = fliplr([0.1, 1, 10, 100]);
% rds = fliplr([0.1, 10]);
for i_rd = 1:length(rds)
  if isinf(pCa) || pCa >= 10
  %   % hack - the no-Ca noPNB experiments had higher ramps
  %   datatable = readtable(['..\Data\bakers_passiveStretch_' num2str(rds(i_rd)*1000) 'ms.csv']);
  %   datatable.Properties.VariableNames = {'Time'  'ML'  'F'  'SL'};
  %   datatables{i_rd} = datatable;
  % elseif isnan(pCa)
      % newest format of experiments    
    datatables{i_rd} = readtable(['..\Data\AvgRelaxed_' num2str(rds(i_rd)) 's.csv']);
  else
    if exist(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv'], "file")
        datatables{i_rd} = readtable(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv']);
    else
        datatables{i_rd} = [];
    end
  end
  % else
  %   % new format for pCa experiments
  %   datatables{i_rd} = readtable(['..\Data\PassiveCa_2\bakers_passiveStretch_pCa' num2str(pCa) '_' num2str(1000*rds(i_rd)) 'ms.csv']);
  % end
end
% if isinf(pCa)
%     % hack - the no-Ca noPNB experiments had higher ramps
%   Lmax = 0.4;
% else    
    % follow-up experiments had lower ramps
    Lmax = 1.175 - 0.95;
% end

% Ls0  = 0.10*mod(13);
% Nx   = 25;          % number of space steps
% ds   = (0.36-Ls0)/(Nx-1);      % space step size
% s  = (0:1:Nx-1)'.*ds; % strain vector
% Ng  = 20;            % number of glubules on globular chain
% delU = 0.0125*mod(1);

% half-sarcomere ramp height
% Lmax = 0.225;
% Nx   = 25;          % number of space steps
Nx   = 25;          % number of space steps
ds   = 1*(Lmax)/(Nx-1);      % space step size
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng = 14; 
delU = 0.0125*mod(19);
% so all unfolded make Ng*delU slack, i.e. 11*0.137=0.1375um


% propose a function to kA = f(pCa)
% kA   = 1*mod(14);
% kD   = 1*mod(15);
% kC   = 103.33*mod(2);
% kS   = 300*mod(3);         % series element spring constant
% alphaU = 2000*mod(4);       % chain unfolding rate constant
% alphaF = 1*mod(5);
% nC = 1.77*mod(6);
% nS = 2.56*mod(7);
% nU = 4*mod(8);
% mu = 2.44*mod(9); 

% g0 = ones(1,11);
% g0 = mod;
% mod = 0;

kp   = mod(1)*10203*0.7;      % proximal chain force constant
if ~isnan(mod(16))
    kA   = mod(16)*0.1*16.44; % PEVK attachment rate
else
    kA   = mod(7)*16.44;
end

if ~isnan(mod(17))
    kD   = mod(17)*14.977; % PEVK detachment rate
else
    kD   = mod(8)*14.977; % PEVK detachment rate
end

alphaU = mod(6)*(8.4137e5)*0.7;         % chain unfolding rate constant

if pCa < 10
    kp   = mod(9)*10203*4.78*0.7;      % proximal chain force constantkS   = g0(2)*14122;        % distal chain force constant
    kA   = mod(7)*16.44;
    kD   = mod(8)*14.977; % PEVK detachment rate
    if ~isnan(mod(20))
        % cant exceed the no Ca unfolding rate!
        alphaU = min(alphaU, mod(20)*(8.4137e5)*0.7);         % chain unfolding rate constant
    end
end
kd   = mod(2)*14122;        % distal chain force constant
alphaF = 0; % chain folding rate constant - not implemented yet
np = mod(3)*3.27; % proximal chain force exponent
nd = mod(5)*3.25; % distal chain force exponent
nU = mod(4)*6.0; % unfolding rate exponent
nF = 1; % folding rate exponent (not implemented yet)
mu = 1*mod(14); % small enough not to affect the result
Lref  = 0.9*mod(18); % reference sarcomere length (um)
alphaF = 100;
alphaF_0 = 0.05*mod(15);

% Calculate proximal globular chain force Fp(s,n) for every strain and
% value. 
slack = (0:Ng).*delU;
Fp = kp*(max(0,s-slack)/Lref).^(np); 
%% visualizing the Force plot
% clf;
% plot(repmat(s, [1, Ng]), Fp(:, 1:Ng), 'linewidth', 2); 
% legend('\itF_{p,1}', '\itF_{p,2}', '\itF_{p,3}', '\it...', 'Location', 'Northwest');
% xlabel('Strain (um)');ylabel('Tension (kPa)');
% set(gca, 'FontSize', 14)
% set(gcf, 'Position', [500  240  400  300])
% ylim([0 60])
% xlabel('s (\mum)');ylabel(['\itt_p (kPa)']);
% mesh(Fp)
%%
% Calculate the globular chain folding/unfolding probability transition
% rates
% RU = alphaU*(max(0,s-slack(1:Ng))).^nU; % unfolding rates from state n to (n+1)
RU = alphaU*((max(0,s-slack(1:Ng))/Lref).^nU).*(ones(Nx,1).*(Ng - (0:Ng-1))); % unfolding rates from state n to (n+1)

%% visualizing the unfolding rate = fig 1C
% clf;hold on;
% plot(repmat(s, [1 Ng]), RU, '-','linewidth', 2);
% plot(repmat(s(end), [Ng, 1])', RU(end, :)', 'ks', 'linewidth', 2);
% for n = 1:4
%     text(s(end) + 0.01, RU(end, n), ['{\itU}_{' num2str(n)  '\rightarrow' num2str(n+1) '}'], 'Fontsize', 14)
% end
% text(s(end) + 0.01, RU(end, 5), '...', 'FontSize',14, 'FontWeight','bold')
% xlabel('s (\mum)');ylabel(['{\itU_{n\rightarrown+1}  (s^{-1})}']);
% set(gca, 'FontSize', 14)
% 
% axes('position', [0.25 0.5 0.35 0.4], 'YAxisLocation','right')
% plot(1:Ng, RU(end, :), 'ks-', 'linewidth', 2);
% xlabel('\itn')
% text(4, max(RU(end, 1))*0.8, ['\itU_{n\rightarrown+1}' char(10) 'at ' num2str(s(end)) '\mum'], 'FontSize',14, 'FontWeight','bold')
% set(gca, 'FontSize', 14)
% set(gcf, 'Position', [500  240  400  300])
% % set(gca, 'YAxisLocation', 'right');
% set(gca, 'YTickLabel', {})
% clf;mesh(RU)
%% Folding rate design and visualization
% RF = alphaF*(max(0,delU-s))  % folding rates from state n+1 to n            
RF = alphaF_0 + alphaF*(slack(1:Ng) - s).*(slack(1:Ng) > s);  % folding rates from state n+1 to n            
% RF = 0.1 + alphaF*(slack(1:Ng) > s);  % folding rates from state n+1 to n            
% mesh(RF);view(3)

% clf;hold on;
% plot(repmat(s, [1 11]), RF(:, 1:11), 's-', 'linewidth', 2);
% xlabel('Strain (um)');ylabel('Transition rate');
% legend('U_0', 'U_{2\rightarrow1}', 'U_{3\rightarrow2}', '...')

% RF = 0;
%
% Initial state
PU = zeros(1,Ng+1); % initial unfolded probabilities for un-attached rectifier state
PA = zeros(1,Ng+1); % initial unfolded probabilities for attached rectifier state

pu = zeros(Nx,1)*PU;
pa = zeros(Nx,1)*PA;
pu(1,1) = 1/ds; 

if pCa >= 11 
    % no Ca effect assumed
    x0 = reshape(pu,[(Ng+1)*Nx,1]);
else
    % might have some Ca effect
    x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
end
x0 = [x0; 0]; 

% Vlist = [1 10 100 1000 5000]*Lmax/100; %  half-sarcomere velocity
% reducing number of ranges
% Vlist = [10 100 1000]*Lmax/100; %  half-sarcomere velocity (um/s)
Vlist = Lmax./rds;
Force = cell(1, 5); 
Time = cell(1, 5); 
Length = cell(1, 5); 
rampSet = 1:length(rds); %[1 2 3 4 5];
rampSet = [2 4];
for j = rampSet
  if isempty(datatables{j})
      fprintf('Skipping pCa %0.2f %0.0fs dataset\n', pCa, rds(j))
      continue;
  end
  % tic
  V = Vlist(j); % ramp velocity

  pu = zeros(Nx,1)*PU;
  pa = zeros(Nx,1)*PA;
  pu(1,1) = 1/ds; 
  if pCa >= 11
    x0 = reshape(pu,[(Ng+1)*Nx,1]);
  else
    x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  end
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

  [t0,x0] = ode15s(@dXdT,[-100:1:0],x0,opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  [t1,x1] = ode15s(@dXdT,[0 Tend_ramp],x0(end,:),opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,V);
  
  % whole decay till the bitter end
  % [t2,x2] = ode15s(@dXdT,[Tend_ramp 200],x1(end,:),opts,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  % limited decay
  [t2,x2] = ode15s(@dXdT,[Tend_ramp Tend_ramp + 40],x1(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % only ramp up, no decay  
  % x2 = [];t2 = []; 

  t = [t1(1:end); t2(2:end)];% prevent overlap at tend_ramp
  x = [x1; x2(2:end, :)];

  %% ramp downm, wait and up again
  % Vdown = Lmax./0.1;
  % [t3,x3] = ode15s(@dXdT, t2(end)+[0 0.1],x2(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,-Vdown);
  % [t4,x4] = ode15s(@dXdT, t3(end)+[0 100],x3(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % [t5,x5] = ode15s(@dXdT, t4(end)+[0 rds(j)],x4(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,V);
  % [t6,x6] = ode15s(@dXdT, t5(end)+[0 10],x5(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % 
  % 
  % t = [t1(1:end); t2(2:end); t3(2:end); t4(2:end); t5(2:end); t6(2:end)];% prevent overlap at tend_ramp
  % x = [x1; x2(2:end, :); x3(2:end, :); x4(2:end, :); x5(2:end, :); x6(2:end, :)];
  
%%
Time{j} = t;
states{j} = [];states_a{j} = [];    strains{j} = []; i_time_snaps = [];
    
    if drawAllStates
        % save current figure
        g = gcf;
        % open up a new one
        figure(50+j); clf;
        % decide for timepoints
        % fixed time or fraction of ramp durations?
        % time_snaps = [0, 0.1, 1, 10, 30, 40, 100]
        time_snaps = [0, rds(j), rds(j) + 30, 60, 120, 160];
        % i_time_snaps = find(t > time_snaps)
    
        % disable
        % time_snaps = [];i_time_snaps = [];
        for i = 1:length(time_snaps)
            if time_snaps(i) > t
                break;
            end
            i_time_snaps(i) = find(t>=time_snaps(i), 1, 'first');    
        end
    end

  for i = 1:length(t)
    xi = x(i,:);
    Length{j}(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    if pCa >= 11
        pa = 0;
    else
        pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    end
    Fd = kd* max(0,(Length{j}(i) - s)/Lref).^nd; 
    Force_pa{j} = ds*sum(sum(Fd.*pa ));
    Force{j}(i) =  ds*sum(sum(Fd.*pu )) + Force_pa{j};
    states{j}(i, 1:Ng+1) = sum(pu);
    states_a{j}(i, 1:Ng+1) = sum(pa);
    strains{j}(i, 1:Nx) = sum(pu, 2);

    if drawAllStates
        if any(ismember(i_time_snaps, i))
            i_snap = find(i_time_snaps == i);
            % snap{i_snap} = 
            subplot(2, length(i_time_snaps), i_snap);cla;
            surf(s, 0:Ng, pu', 'EdgeColor','none');hold on;
            surf(s, 0:Ng, Fp')
            colormap(1-gray);
            xlim([0, s(end)]);
            shading(gca, 'interp')
            % view(90, -90); 
            xlabel('s'); ylabel('State');
            title(sprintf('U (%fs), S= %0.1f', t(i), sum(pu(:))));
            
    
            % next line
            subplot(2, length(i_time_snaps), length(i_time_snaps) + i_snap);
            plot(repmat(s, [1 Ng+1]), pu(:, 1:Ng+1), 's-')
            legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 'Location', 'best');
            axis([0, s(end), 0, max(max(x))])
            
            % ignore pa for a moment
            % if ~(pa == 0)
            %     subplot(2, length(i_time_snaps), i_snap + length(i_time_snaps));
            %     surf(s, 0:Ng -1, pa, 'EdgeColor','none');
            %     view(90, -90); xlabel('s'); ylabel('State');
            %     title(sprintf('A (%fs)', t(i)))
            % end
        end
    end    

    % if t(i) >= rds(j) && t(max(1, i-1)) < rds(j)
    %     % at the peak
    %     outStruct{j, 1}.pa = pa;
    %     outStruct{j, 1}.pu = pu;
    % elseif t(i) >= rds(j)*2 && t(max(1, i-1)) < rds(j)*2
    %     % after t_rd after the peak
    %     outStruct{j, 2}.pa = pa;
    %     outStruct{j, 2}.pu = pu;
    % elseif t(i) >= rds(j) + 20 && t(max(1, i-1)) < rds(j) + 20
    %     % after 20s after the peak
    %     outStruct{j, 3}.pa = pa;
    %     outStruct{j, 3}.pu = pu;
    % end
    end

%%
if drawAllStates
    figure(60+j); clf;
    set(gcf, 'Name', sprintf('States for pCa %d at %0.1f ramp', pCa, rds(j)) );
    subplot(131);
    h = surf(0:Ng, Time{j}, states{j}, 'EdgeColor','none');
    shading(gca, 'interp')
    view(90, -90)
    ylabel('Time (s)'); xlabel('State occupancy (#)')
    title('States pu');
    colorbar;
    
    subplot(132);
    plot(t, Force{j}, LineWidth=2);
    % subplot(132);
    % surf(0:Ng, Time{j}, states_a{j});
    % shading(gca, 'interp')
    % view(90, -90)
    % ylabel('Time (s)'); xlabel('State occupancy (#)')
    % title('States pa');
    % colorbar;
    
    subplot(133);
    surf((1:Nx)*ds, Time{j}, strains{j});
    shading(gca, 'interp')
    view(90, -90)
    ylabel('Time (s)'); xlabel('Strains (um)')
    title('Strains');
    
    colorbar;
    
    % set to preset figure
    figure(g);
end

%%

    
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
  
  b = 0.05*mod(11);
  c = 7*mod(12);
  d = 1*mod(13);
  % apply constraints
  if b < 0 || c <= 0 || d < 0 
      cost = inf;
      return;
  end
  
  % calculate a, so that the max value is the same  
  Fss = 3.2470*mod(10); % reducing the param space
  a = (Fss - d)/((Lmax -b)^c);
  Force_par{j} = a*max(Length{j} - b, 0).^c + d;
  % calc force
  Force{j} = Force{j} + Force_par{j}; 

  if any(isnan(Force{j}))
      disp('error');
      cost = inf;
      return;
  end

  % Force{j} = Force{j} + (0.55e6)*0.225^8*mod(10); 


  % Fss = mod(10)*3; % optimized previously
  % Force{j} = Force{j} + Fss;

    % ttoc = toc;
  % fprintf("Ramp %0.1fs takes %1.1fs \n", Tend_ramp, ttoc)
    
end

% Get error for the whole ramp-up and decay
t_endFreeware = zeros(1, 5); % time when we start counting the costs
% alternatively, get the error from decay only
% t_endFreeware  = Lmax./Vlist + 2;

%% Evaluating all ramps at once
En = cell(1, length(rampSet));
Es = cell(0);
for j = rampSet
    datatable_cur = datatables{j};
    if isempty(datatable_cur)
        continue;
    end
    inds = find(datatable_cur.Time >= t_endFreeware(j));
    datatable_cur = datatable_cur(inds, :);
    t_int{j} = datatable_cur.Time - 2;
    % t_int{j} = datatable_cur.Time;
    Ftot_int{j} = interp1(Time{j}, Force{j}, t_int{j}); % total force interpolated

    % weighting to fit semilogx
    x = -log10(max(rds(j), t_int{j}));
    % center and scale
    w = (x+abs(min(x)))/mean(x+abs(min(x)));
    
    % semilogx(t_int{j}, x/sum(x));
    % plot(t_int{j}, w);hold on;

    % no weighing, already in the data
    w = 1;
%%
    Es{j} = w.*(Ftot_int{j} - datatable_cur.F).^2; % error set
    Es{j}(isnan(Es{j})) = 0; % zero outside bounds
    En{j} = 1e3*sum(Es{j})/length(Es{j}); % normalized error
    
    PeakData(j, 1) = rds(j);
    PeakData(j, 2) = max(datatable_cur.F);
    PeakModel(j) = max(Ftot_int{j});

end

PeakModel = nan(1, length(PeakData));
for j = 1:length(PeakModel)
    m = max(Force{j});
    if ~isempty(m)
        PeakModel(j) = m;
    end
end

Ep = nansum((PeakData(:, 2) - PeakModel(:)).^2);
% discarding peak fit for high Ca's
if pCa < 10
    Ep = 0;
end

cost = Ep*100 + sum([En{1:end}], 'all');


if exist('drawPlots', 'var') && ~drawPlots
    return;
end
try

%%
% zoomIns = [0 200 0 10;...
%            0 20 0 15;...
%            0 2 0 15;...
%            0 .2 0 15;...
%            0 .04 0 15;...
%            ];

clf;
colors = lines(max(rampSet)+1);

% max out of all
ym = ceil( max(cell2mat(Force)) / 5 ) * 5;

% prepare in advance so that it wont draw over my inset
sp = subplot(212);hold on;
pos = get(sp, 'Position');
ylabel('Tension (kPa)')
xlabel('Time (s)')
set(gca,'Fontsize',14)
title(sprintf('Force response to %.2g ML ramp-up at pCa=%g, costing %1.4e€', Lmax, pCa, cost), 'Parent',sp);
set(sp, 'XLim', [-1 30+max(rds)]);
set(sp, 'YLim', [0 ym*1])

% shift of peaks to have the same tail - just guessed
% shift = [-94, -7.2, -0.35, 0];
% based on pCa 11 shift in data
shift = [5.4, 0.82, 0.22, 0.01] - [100 10 1 0.1];
ymaxScale = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
% figure(j); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
    % subplot(1, 3, j);hold on;
%% primary plot - semilog
    subplot(221)
    semilogx(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    hold on;
    semilogx(t_int{j},Es{j},'--','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    
    semilogx(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    axis([1e-2, 1e2, 0, ym]);  
    set(gca,'Fontsize',14)
    title('Tension response to muscle length ramp-up')
    xlabel('Time (s)')
    ylabel('Tension (kPa)')

%% other view - shifted to see the tail overlap

    subplot(222)
    % Estimating the true offset: Fss = C*(Tss)^-alpha + Fss_true;
    % Fss = Force_par{j}(end); % "steady state" at the end 
    % tss = Time{j}(end) - rds(j);
    % Fss_true = Fss - (4.22*tss^-0.21);
    Fss_true = Force_par{j}(end);

    loglog(datatables{j}.Time-2 + shift(j),datatables{j}.F - Fss_true,'-','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    hold on;
    loglog(Time{j} + shift(j),Force{j} - Fss_true,'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % plot(Time{j} + shift(j),Force{j},styles{j}, 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    % axis([1e-2, 1e2, 0, ym]);  
    xlim([1e-2, 1e2]);
    ymaxScale = max(ymaxScale, max(Force{j} - Fss_true));
    yminScale = (Force{j}(end) - Fss_true); % this should be around the same
    ylim([max(1e-2, 0.8*yminScale), 1.2*ymaxScale])
    if j == 1
        % only after the last one
        legend('Ramp 10s (Data)', 'Ramp 10s (Model)', 'Ramp 1s (Data)', 'Ramp 1s (Model)', 'Ramp 0.1s (Data)', 'Ramp 0.1s (Model)');
    end
    % legend('Ramp 10s, shifted by -8.6s', 'Ramp 1s, shifted by -0.78', 'Ramp 0.1s' );

    set(gca,'Fontsize',14)
    title('Tension response to muscle length ramp-up: shifted peaks')   

%% secondary plots - timebase. Need to cut out

    plot(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',2, 'Color', [colors(j+1, :), 0.15], Parent=sp);
    plot(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8, Parent=sp);
    plot(t_int{j},Es{j},'--|','linewidth',2, 'Color', [colors(j+1, :), 0.3], Parent=sp);

    % plot(t_int{j},Ftot_int{j},'r','linewidth',2.5);
    % max out of all

    
    ylabel('Tension (kPa)')
    xlabel('Time (s)')
    % set(gca,'Xtick',0:50:200)
    
    % zoom-in inset
    
    w = pos(3)*0.18; h = pos(4)*0.5;
    x = pos(1) + pos(3) - (j*1.3 - 0.3)*w; y = pos(2) + pos(4) - h;
    axes('Position',[x, y, w, h]);hold on;
    % axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
    plot(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',3, 'Color', [colors(j+1, :), 0.15]);
    plot(t_int{j},Es{j},'--','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    % plot(Time{j},Force{j}, 'r:', 'linewidth',2); 
    plot(t_int{j},Ftot_int{j},'-','linewidth',2, 'Color', colors(j+1, :)*0.9);
    ym_inset = ceil( max(Force{j}) / 10 ) * 10;
    axis([0, rds(j)*2, 0, ym_inset]);  
    set(gca,'Fontsize',14)    
    
end
% cla;
pos1 = get(subplot(221), 'Position');
w = pos1(3)*0.45; h = pos1(4)*0.4;
x = pos1(1) + pos1(3) - w; y = pos1(2) + pos1(4) - h;
axes('Position',[x, y, w, h]);
semilogx(PeakData(:, 1), PeakData(:, 2), 'ko', LineWidth=2);hold on;
semilogx(PeakData(:, 1), PeakModel, 'x', 'MarkerEdgeColor', [1 1 1]*0.5, LineWidth=2, MarkerSize=8);
axis([1e-1 1e1 0 ym])
semilogx(PeakData(:, 1), PeakModel, '--', Color=[1 1 1]*0.5, LineWidth=1);
legend('Peaks (Data)', 'Peaks (Model)', 'Location', 'southeast')
%%
h = annotation('textbox', [0.07 0.95 0 0], 'String', 'A)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
h = annotation('textbox', [0.5 0.95 0 0], 'String', 'B)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
h = annotation('textbox', [0.07 0.5 0 0], 'String', 'C)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
%%
catch e
    disp(e.message)
end
%%

% fig = gcf;
% set(gcf, 'Position', [50 50 1200 700])
% saveas(fig, ['..\Figures\Fig_' fig.Name], 'png')

return;

%% Overlap plots

%% Draw plots
figure(1001);clf;hold on;legend()
% plot n.1: 
colororder(jet(Ng));
for n = 1:1:Ng
    
    plot(outStruct{1, 1}.pu(:, n), 'x-', LineWidth=2)
    plot(outStruct{1, 2}.pu(:, n), 'x--', LineWidth=2)
    plot(outStruct{1, 3}.pu(:, n), 'x:', LineWidth=2)
end

% plot u to s for different N

% Visualize states in time?
% / Passive resting fit - both semilog and linear?
% / Shift the peaks so the tails overlap?
% / Fit for maximal CA
% Values of the Ca sensitive params - kA, KD?, kd
