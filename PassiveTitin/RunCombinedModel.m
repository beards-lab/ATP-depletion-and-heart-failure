% Assuming mod = ones(10,1)
% assuming pCa = Inf
% pCa 11 - load Relaxed, do not run the extended, PEVK attachment model
% pCa 10 - load Relaxed, run the PEVK attachment model
% pCa < 10 - load AvgpCa dataset, Ca effect in place

% simtype = 'sin';
simtype = 'ramp';
% simtype = 'rampbeat';
% simtype = 'ktr';
% simtype = 'stairs-up'
% simtype = 'slackset'
% simtype = 'atpprotocol'

% rampSet = [1]; % 100s
% rampSet = [2 4]; % 10s and 0.1s
% rampSet = [3]; % nly 1s
% rampSet = [1 2 3 4]; % all
rampSet = [4]; % only 100ms

clear Force
clear Time
clear Length
clear outStruct;
if any(mod < 0) 
    cost = inf;
    return;
end
drawAllStates = false;
drawFig1 = false;
exportRun = false;

figInd = get(groot,'CurrentFigure'); % replace figure(indFig) later without stealing the focus


  
  if strcmp(simtype, 'ramp')
      %% normal
      % rds = fliplr([0.02 0.1, 1, 10 100]);
      rds = fliplr([0.1, 1, 10, 100]);
      for i_rd = 1:length(rds)
          if isinf(pCa) || pCa >= 10
              % newest format of experiments
              datatables{i_rd} = readtable(['..\Data\AvgRelaxed_' num2str(rds(i_rd)) 's.csv']);
          else
              if exist(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv'], "file")
                  datatables{i_rd} = readtable(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv']);
              else
                  datatables{i_rd} = [];
              end
          end
      end
      % tic
      L0 = 0;

      Lmax = 1.175 - 0.95; % Lmax = 0.225; % identified to ramp of this height      
      Vlist = Lmax./rds; %  half-sarcomere velocity (um/s)
      Force = cell(1, 5);
      Time = cell(1, 5);
      Length = cell(1, 5);
      LengthPar = cell(1, 5);
      for j = rampSet
          if isempty(datatables{j})
              fprintf('Skipping pCa %0.2f %0.0fs dataset\n', pCa, rds(j))
              continue;
          end
          V = Vlist(j); % ramp velocity    
          Tend_ramp = Lmax/V; % length of ramp
    
          times = [-100, 0;0 Tend_ramp;Tend_ramp Tend_ramp + 10000];
          velocities = {0 V 0};

          plotOptions = struct();
          if drawAllStates
            plotOptions.time_snaps = [1e-3, rds(j), 0*1 + 2*rds(j), max(1000, rds(j)*10)];
          end

          [Time{j}, Length{j}, Force{j}, Force_par{j}] = evaluateTitinModel(mod, L0, times, velocities, pCa, plotOptions);
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
          Es{j} = ((Ftot_int{j} - datatable_cur.F)/max(Ftot_int{j})).^2; % error set
          Es{j} = Es{j}(~isnan(Es{j})); % zero outside bounds
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

      Ep = nansum(((PeakData(:, 2) - PeakModel(:))/max(Ftot_int{j})).^2);
      % discarding peak fit for high Ca's
      if pCa < 10
          Ep = 0;
      end

      cost = Ep*100 + sum([En{1:end}], 'all');

      if exist('drawPlots', 'var') && ~drawPlots
          return;
      end

      Tarr = t_int;Farr = Ftot_int(1:4);
      if exportRun
          % save data for decay overlay loglog plot
          Tarr = t_int;Farr = Ftot_int(1:4);
          if pCa < 10
              save pcagraphingenv.mat
          else
              save relaxgraphingenv.mat
          end
          save(sprintf('..\\pca%gmodeldata.mat', pCa), 'Tarr', 'Farr')
          % save('..\pca11modeldataDoubleStates.mat', 'Tarr', 'Farr')
      end

  elseif strcmp(simtype, 'sin')
    %% sinusoidal driving
      % cycle time
      Tc = rds(j);
      times = [-100, 0;0 10*Tc];  
      positions = @(t)-Lmax/2*cos(2*pi*t/Tc) + Lmax/2;
      % differentiating positions
      V = @(t)2/2*Lmax*pi/Tc *sin(2*pi*t/Tc);
    
      % syms fpos(t);   % fpos(t) = @(t)Lmax*sin(2*pi*t/Tc);
    
      velocities = {0, V};
      L0 = 0;%Lmax/2;
      x_ = 0:Tc/100:Tc;
      plot(x_, positions(x_), '-', x_, V(x_), x_, cumsum(V(x_)).*[1 diff(x_)], '--');

    if drawAllStates
        time_snaps = [Tc/2, Tc, 3*Tc/2, 2*Tc];
    end

    %% Plot sinusoidal outcome
    aspect = 2;
    figure(900 + j*10 + round(pCa));clf;    tiledlayout(2,2, TileSpacing="compact");
    set(gcf, 'Position', [500  300  7.2*96 7.2*96/aspect])
    rng = Time{j} < 20;
    t = Time{j}(rng)';
    x = Length{j}(rng) + 0.95;
    y = Force{j}(rng);
    y2 = Force{j} - Force_par{j};
    z = zeros(size(t));
    col = linspace(0, t(end), length(t));   
    x_ = @(t) t - floor(t/Tc)*Tc;
    lw = 1

    nexttile;
    plot(t, x, 'k', linewidth = lw);
    % plot(t, x, t, y);legend('Length', 'Force');
    xlabel('$t$ (s)', Interpreter='latex');    ylabel('$L$ ($\mu$m)', Interpreter='latex');
    ylim([0.95, 1.2])

    nexttile(2, [2 1]);
    surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw);
    clim([0 50 100])
    colormap(turbo);    
    ylabel('$\Theta$ (kPa)', Interpreter='latex');     xlabel('$L$ ($\mu$m)', Interpreter='latex');
    cb = colorbar;  cb. title(cb, 't (s)');

    nexttile;
    % plot(x_(t), x, x_(t), y);legend('Length', 'Force');
    plot(t, y, 'k', linewidth = lw);
    xlabel('$t$ (s)', Interpreter='latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex');
    fontsize(12, 'points');
    exportgraphics(gcf,sprintf('../Figures/FigHysteresis%g_%gs.png', pCa, rds(j)),'Resolution',150)
      
  elseif strcmp(simtype, 'ktr')
      %%
      dt1 = 1e-3; % drop down
      dt2 = 100e-3; % hold down
      dt3 = 10e-3; % ramp up
      dt4 = 5e-3; % hold up
      dt5 = 1e-3; % drop to base
      t = [cumsum([0 dt1 dt2 dt3 dt4 dt5]) 1]; % ending times
      times = [[-100 t(1:end-1)]',t']  % generate starting times and merge two columns
      level0 = 1.1;
      level1 = 0.95;
      level2 = 1.11;      
      velocities = {0 (level1-level0)/dt1 0*dt2 (level2-level1)/dt3 0*dt4 (level0-level2)/dt5 0};
      L0 = 0.05; % starting at L0 = 1 ML   
  elseif strcmp(simtype, 'stairs-up')
      %%
      levels = linspace(0, 0.1, 10) + 0.05;
      levels = levels([1:8 10]);
      t_hold = 20e-3;
      t_ramp = 20e-3;
      times = [-100, 0];
      times = [times;0 t_hold];
      velocities = {0 0};
      L0 = levels(1); % starting at L0 = 1 ML   
      for i_lev = 2:length(levels)
          times = [times;times(end) times(end) + t_ramp];
          velocities{end + 1} = (levels(i_lev) - levels(i_lev - 1))/t_ramp;
          times = [times;times(end) times(end) + t_hold];
          velocities{end + 1} = 0;
      end      
  elseif strcmp(simtype, 'slackset')
      %%
      level0 = 1.1;
      % 8 10 12 14 16 % depresion
      depLevels = level0*(1 - [8 10 12 14 16]/100); % depresion
      topLevels = [repmat(level0, [1 length(depLevels) - 1]) 1]; % rebound
      dt1 = 1e-3; % drop time
      dt2 = 100e-3; % hold time
      dt3 = 20e-3; % rise time
      dt4 = 0.4; % wait time
      t = cumsum([dt1 dt2 dt3 dt4]);
      
      times = [-100, 0;0 0.5;0.5 1];
      velocities = {0 0.1/0.5 0};
      for i = 1:length(depLevels)
        tims_ = [[0 t(1:end-1)]',t'] + times(end);
        times = [times;tims_];
        vels = {(depLevels(i) - level0)/dt1, 0, (topLevels(i) - depLevels(i))/dt3, 0};
        velocities = [velocities,vels];
      end
      L0 = 0.05;
  elseif strcmp(simtype, 'atpprotocol')
      %% copypaste all three above

      % ktr 
      dt1 = 1e-3; % drop down
      dt2 = 100e-3; % hold down
      dt3 = 10e-3; % ramp up
      dt4 = 5e-3; % hold up
      dt5 = 1e-3; % drop to base
      t = [cumsum([0 dt1 dt2 dt3 dt4 dt5]) 1]; % ending times
      times = [[-100 t(1:end-1)]',t']  % generate starting times and merge two columns
      level0 = 1.1;
      level1 = 0.95;
      level2 = 1.11;      
      velocities = {0 (level1-level0)/dt1 0*dt2 (level2-level1)/dt3 0*dt4 (level0-level2)/dt5 0};
      L0 = 0.05; % starting at L0 = 1 ML   

      % stairs
      levels = linspace(0, 0.1, 10) + 0.05;
      levels = levels([1:8 10]);
      t_hold = 20e-3;
      t_ramp = 20e-3;
      for i_lev = 2:length(levels)
          times = [times;times(end) times(end) + t_ramp];
          velocities{end + 1} = (levels(i_lev) - levels(i_lev - 1))/t_ramp;
          times = [times;times(end) times(end) + t_hold];
          velocities{end + 1} = 0;
      end      

      % slack
      level0 = 1.1;
      % 8 10 12 14 16 % depresion
      depLevels = level0*(1 - [8 10 12 14 16]/100); % depresion
      topLevels = [repmat(level0, [1 length(depLevels) - 1]) 1]; % rebound
      dt1 = 1e-3; % drop time
      dt2 = 100e-3; % hold time
      dt3 = 20e-3; % rise time
      dt4 = 0.4; % wait time
      t = cumsum([dt1 dt2 dt3 dt4]);      
      for i = 1:length(depLevels)
        tims_ = [[0 t(1:end-1)]',t'] + times(end);
        times = [times;tims_];
        vels = {(depLevels(i) - level0)/dt1, 0, (topLevels(i) - depLevels(i))/dt3, 0};
        velocities = [velocities,vels];
      end
     
  end
  %% repeated - refolding
  % times = [-100, 0;0 Tend_ramp;Tend_ramp Tend_ramp + 40;... % normal ramp-up
  %     Tend_ramp + 40 Tend_ramp + 41;... % rampdown in 1s
  %     Tend_ramp + 41 Tend_ramp + 41 + 60;... % hold down - wait for refold
  %     Tend_ramp + 41 + 60 Tend_ramp + 41 + 60 + Tend_ramp;... ramp-up
  %     Tend_ramp + 41 + 60 + Tend_ramp Tend_ramp + 41 + 60 + Tend_ramp + 40]; % final hold
  % velocities = [0 V 0 ...
  %     -Lmax/1 0 V 0];
  %% 10 heartbeats
  % times = [-100, 0];velocities = [0];
  % HR = 600;Tc = 60/HR;
  % sf = 0.2; % systole fraction
  % df = 0.4; % diastole fraction
  % Lmax = (2.2 - 1.9)/2;
  % Vc = Lmax/(Tc*sf);Vr = Lmax/(Tc*df);
  % for b = 1:10
  %   % relaxation - blowing up
  %   times = [times; times(end) times(end)+Lmax/Vr];
  %   velocities = [velocities Vr];
  %   % hold
  %   times = [times; times(end) times(end)+Tc*(1-sf - df)];
  %   velocities = [velocities 0];
  %   % contraction
  %   times = [times; times(end) times(end) + Lmax/Vc];
  %   velocities = [velocities -Vc];
  %   drawAllStates = true;
  % end
%%



  if any(isnan(Force{j}))
      disp('error');
      cost = inf;
      return;
  end

  

try
% tic
set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
cf = clf;

aspect = 2;
% normal size of 2-col figure on page is 7.2 inches
% matlab's pixel is 1/96 of an inch
% f.Position = [300 200 7.2*96 7.2*96/aspect];

% f.Position = [300 200 7.2*96 7.2*96/aspect];
% colors = lines(max(rampSet)+1); colors(1:end-1, :) = colors(2:end, :);
colors = gray(5);
% colors(1:end-1, :) = colors(2:end, :);
fs = 12;
% tl = tiledlayout(3, 4, 'TileSpacing', 'compact');
rws_s = 4;rws_l = rws_s*4;rws_loglog = 6;
tl = tiledlayout(3, rws_loglog+rws_l, "Padding","compact", "TileSpacing","compact");
tile_semilogx = nexttile(1, [2, rws_l]);
tile_loglog = nexttile(rws_l + 1, [3, rws_loglog]);
tile_r(1) = nexttile((rws_loglog+rws_l)*2 + 1, [1 rws_s]);
tile_r(2) = nexttile((rws_loglog+rws_l)*2 + 1 + rws_s, [1 rws_s]);
tile_r(3) = nexttile((rws_loglog+rws_l)*2 + 1 + 2*rws_s, [1 rws_s]);
tile_r(4) = nexttile((rws_loglog+rws_l)*2 + 1 + 3*rws_s, [1 rws_s]);
tile_positions = [tile_semilogx.Position;tile_loglog.Position;...
    tile_r(1).Position;tile_r(2).Position;tile_r(3).Position;tile_r(4).Position];
clf;

reportCosts = false;
if reportCosts 
    title(tl, sprintf('pCa %g costs %g', pCa, round(cost, 3)));    
end
% tile_semilogx = nexttile(1, [2, rws_l]);hold on;
% tile_semilogx = subplot(3, rws_l + rws_loglog, );hold on;
tile_semilogx = axes('Position', tile_positions(1, :));hold on;

ym = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    hd{j} = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([hd{j}.Bar, hd{j}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [hd{j}.Line.ColorData(1:3); 255*0.4])
    hm{j} = semilogx(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    hph{j} = plot(nan, nan, 'x',Color=[1 1 1]); % just a placeholder
    % set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
    ym = max([ym Force{j}, datatables{j}.F']);
    set(tile_semilogx, 'TickLength', [0.0125 0.05]);
    set(tile_semilogx, 'TickLabelInterpreter', 'latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex')
    
    if rds(j) > 1
        tp = [rds(j)*1.4, max(Force{j})*1.2];
    else
        tp = [rds(j)*0.95, max(Force{j})*0.95];
    end
    text(tp(1), tp(2), ...
        sprintf('$t_r$ = %g s', rds(j)),...
        'horizontalAlignment', 'right', VerticalAlignment='bottom', Interpreter='latex');    
end
xlim([1e-2, 160])
ylim([0 ceil((ym)/5)*5])
yl = ylim();
tile_semilogx.XScale='log';
box on;
% hl = legend([hph{1} hph{2} hph{3} hph{4} hd{4} hd{3} hd{2} hd{1} hm{4} hm{3} hm{2} hm{1}], ...
%     't_r 0.1 s:', '  t_r 1 s:',' t_r 10 s:','t_r 100 s:',...
%     'Data', 'Data', 'Data', 'Data', ...    
%     'Model', 'Model', 'Model', 'Model',...
%     NumColumns=3);
hl = legend([hd{1} hm{1}], ...
    'Data',...    
    'Model',...
    NumColumns=3);

hl.ItemTokenSize = [30, 20];
hl.Box = 'off';

% title('Model fit', Interpreter='latex')

    % 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s'...
    % )
    %%

% tile_semilogx.Position = tile_positions(1, :).*[1 1 1 1]
%% nexttile;
% semilogx(PeakData(:, 1), PeakData(:, 2), 'ko', LineWidth=2);hold on;
% semilogx(PeakData(:, 1), PeakModel, 'x', 'MarkerEdgeColor', [1 1 1]*0.5, LineWidth=2, MarkerSize=8);
% axis([1e-1 1e2 0 ym])
% semilogx(PeakData(:, 1), PeakModel, '--', Color=[1 1 1]*0.5, LineWidth=1);
% hl = legend('Data', 'Model', 'Location', 'best')
% title(hl, 'Peaks')
% ylabel('Tension (kPa)', Interpreter='latex')
%% loglog plot
% tile_loglog = nexttile(rws_l + 1, [3, rws_loglog]);
tile_loglog = axes('Position', tile_positions(2, :).*[1.05 1.8 1 1] + [0 0 0 -0.0440] );box on;
% tile_loglog.Position = tile_positions(2, :).*[1.05 1.8 1 1] + [0 0 0 -0.0440] 
return;
%% best fit from FigFitDecayOverlay
if pCa > 10
    x = [4.7976    0.2392    4.8212];
    [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], false);
else
    x = [17.9381    0.2339    4.3359];
    [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], true);
end


%%
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    % sp = nexttile((rws_loglog+rws_l)*2 + 1 + (4-j)*rws_s, [1 rws_s]);
    sp = axes('Position', tile_positions(7-j, :).*[1 1.8 1.2 0.8] + (4-j).*[-0.008 0 0 0]);box on;
    hold on;
    if reportCosts 
        title(sp, sprintf('Ramp %g costs %g', rds(j), round(En{j}, 3)));    
    end
    h = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255])
    plot(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    xlim([0 min(rds(j)*3, 160)]);
    xl = xlim;
    ylim(yl)
    yticks([0 20 40]);
    xticks([0 rds(j), rds(j)*2])
    if j == 4
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
    else
        yticklabels([]);
    end
    xlabel('$t$ (s)', Interpreter='latex', HorizontalAlignment='left');
    set(sp, 'TickLength', [0.05 0.05]);
    set(sp, 'TickLabelInterpreter', 'latex')
    tit = text(xl(2), min(yl(2)-5, max(Force{j})), sprintf('$t_r$ = %g', rds(j)), 'Interpreter', 'latex', ...
        HorizontalAlignment='right', VerticalAlignment='bottom');
    % tit.Position = tit.Position + [0 max(Force{j}) 0];
    % pos = sp.Position;
    % sp.Position = pos + [0 0 0.1 0];
    fontsize(12, 'points');
end    


aspect = 1.5;
set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);
if exportRun
    exportgraphics(cf,sprintf('../Figures/ModelFitpCa%g.png', pCa),'Resolution',150)
    exportgraphics(cf,sprintf('../Figures/ModelFitpCa%g.eps', pCa))
end
% toc
return
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
tile_semilogx = subplot(212);hold on;
pos = get(tile_semilogx, 'Position');
ylabel('Tension (kPa)')
xlabel('Time (s)')
set(gca,'Fontsize',14)
title(sprintf('Force response to %.2g ML ramp-up at pCa=%g, costing %1.4e€', Lmax, pCa, cost), 'Parent',tile_semilogx);
set(tile_semilogx, 'XLim', [-1 300]);
set(tile_semilogx, 'YLim', [0 ym*1])

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
    
    yyaxis right;
    semilogx(t_int{j},cumsum(Es{j}),'--','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    yyaxis left;
    
    semilogx(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    axis([1e-2, 3e2, 0, ym]);  
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
    Fss_true = Force{j}(end) - Force_par{j}(end);
% shift(j) = 0;
    loglog(datatables{j}.Time-2 + shift(j),datatables{j}.F - Fss_true,'-','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    hold on;
    loglog(Time{j} + shift(j),Force{j} - Fss_true,'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % plot(Time{j} + shift(j),Force{j},styles{j}, 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    % axis([1e-2, 1e2, 0, ym]);  
    xlim([1e-2, 3e2]);
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

    plot(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',2, 'Color', [colors(j+1, :), 0.15], Parent=tile_semilogx);
    plot(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8, Parent=tile_semilogx);
    plot(t_int{j},Es{j},'--|','linewidth',2, 'Color', [colors(j+1, :), 0.3], Parent=tile_semilogx);

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
axis([1e-1 1e2 0 max([ym PeakData(:, 2)'])])
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
return;
Tss =  (0.50e6)*0.225^8;
figure(11); 
loglog(Time{4}-0.1,Force{4}-Tss, ...
       Time{3}-1.0,Force{3}-Tss, ...
       Time{2}-10,Force{2}-Tss, Time{1}-10,Force{1}-Tss); 
grid
save('../modeltesting2.mat', 'Time', 'Force', 'Tss')
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
