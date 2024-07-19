function [Time, Length, Force_INST, Force_par, x0, Force_UNF] = evaluateTitinModel(mod, L0X0, times, velocities, pCa, plotOptions)

if isfield(plotOptions, 'time_snaps') && ~isempty(plotOptions.time_snaps)
    time_snaps = plotOptions.time_snaps;
    drawAllStates = true;
else
    drawAllStates = false;
end

if isfield(plotOptions, 'drawFig1')
    drawFig1 = plotOptions.drawFig1;
else
    drawFig1 = false;
end


Lmax = 1.175 - 0.95; % Lmax = 0.225; % identified to ramp of this height
Nx   = 25;          % number of space steps
ds   = 1*(Lmax)/(Nx-1);      % space step size
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng = 15; 
% so all unfolded make mod(19) = Ng*delU slack, i.e. ~0.15 um. individual delU is about 0.01
delU = mod(19)/Ng;

kp   = mod(1); % proximal chain force constant
kd   = mod(2); % distal chain force constant

if ~isnan(mod(16))
    kA   = mod(16); % PEVK attachment rate
else
    kA   = mod(7);
end

if ~isnan(mod(17))
    kD   = mod(17); % PEVK detachment rate
else
    kD   = mod(8); % PEVK detachment rate
end

alphaU = mod(6);         % chain unfolding rate constant
Fss = mod(10);
if pCa < 10
    kp   = mod(9);      % proximal chain force constantkS   = g0(2)*14122;        % distal chain force constant
    kA   = mod(7);
    kD   = mod(8); % PEVK detachment rate
    if ~isnan(mod(20))
        % cant exceed the no Ca unfolding rate!
        alphaU = min(alphaU, mod(20));         % chain unfolding rate constant
    end
    if ~isnan(mod(21))
        Fss = mod(21);
    end
    if ~isnan(mod(22))
        kd   = mod(22);      % proximal chain force constant high Cabist
    else
        kd = mod(2)*mod(9)/mod(1); % scaled to kp_ca in the same ratio as kd_relaxed to kp_relaxed
    end
        
end
kDf = mod(23);
alphaF = 0; % chain folding rate constant - not implemented yet
np = mod(3); % proximal chain force exponent
nd = mod(5); % distal chain force exponent
nU = mod(4); % unfolding rate exponent
nF = 1; % folding rate exponent (not implemented yet)
mu = mod(14); % small enough not to affect the result
L_0  = mod(18); % reference sarcomere length (um)
alphaF = 0;
alphaF_0 = mod(15);

% Calculate proximal globular chain force Fp(s,n) for every strain and
% value. 
slack = (0:Ng).*delU;
% Fp = kp*(max(0,s-slack)/Lref).^(np); 
Fp = kp*(max(0,s-slack)/L_0).^(np); 
% Calculate the globular chain folding/unfolding probability transition
% rates
% RU = alphaU*(max(0,s-slack(1:Ng))).^nU; % unfolding rates from state n to (n+1)
RU = alphaU*((max(0,s-slack(1:Ng))/L_0).^nU).*(ones(Nx,1).*(Ng - (0:Ng-1))); % unfolding rates from state n to (n+1)
%% visualizing the Force plot - fig 1A&B
if drawFig1
    % Fp = kp*(max(0,s-slack)/Lref).^(np); 

    % do not plot zeros
    % Fp(Fp < 1e-1) = NaN;
    % sx = s(1):0.1:s(end);
    sx = s;
    % Fp = interp1(s, Fp(:, :), sx);
    f = figure(3); clf;
    set(gcf, 'Position', [500  300  7.2*96 3.5*96])
    gc = axes('Position', [0.1, 0.2, 0.35, 0.7]);
    pl1 = plot(sx, Fp(:, 1), 'k-', 'linewidth', 2); hold on;
    pl2 = plot(sx, Fp(:, 2), 'k-.', 'linewidth', 2); hold on;
    pl3 = plot(sx, Fp(:, 3), 'k--', 'linewidth', 2); hold on;
    plN = plot(repmat(sx, [1, Ng-3]), Fp(:, 4:Ng), ':', 'linewidth', 2, 'Color', 'k'); 
    plSS = plot(sx, kd*max(0,(max(sx)-sx)/L_0).^nd, 'r--');
    
    % arrow
    % t0 = 5;
    % plot(s([14, 24]), [t0 t0], 'k');scatter(s(14), [t0], 'k>', 'filled');
    % sx_20 = zeros(1, Ng);
    % for n = 1:Ng
    %     i_s = find(Fp(:, n) > 0, 1);
    %     sx_20(n) = interp1(Fp(i_s:end, n), s(i_s:end), t0);
    % end
    % scatter(sx_20, repmat(t0, [1 Ng]), 'kx', linewidth=1.5)   
    legend([pl1, pl2, pl3, plN(1)], '$\sigma_p(0,s)$', '$\sigma_p(1,s)$', '$\sigma_p(2,s)$', '$\sigma_p(3..14,s)$', 'Location', 'Northwest', Interpreter='latex');
    
    set(gca, 'FontSize', 12)
    aspect = 1.5;
    
    ylim([0 max(Fp(:))]);xlim([0 inf])
    xlabel('$s$ ($\mu$m)', Interpreter='latex');ylabel('$\sigma_p(n,s)$ (kPa)', Interpreter='latex');
    % visualizing the unfolding rate = fig 1C
    % clf;hold on;
    axes('position', [0.55 0.2 0.35, 0.7]);hold on;
    
    plot(s, RU(:, 1), 'k-','linewidth', 2);
    plot(s, RU(:, 2), 'k-.','linewidth', 2);
    plot(s, RU(:, 3), 'k--','linewidth', 2);
    plot(repmat(s, [1 Ng-3]), RU(:, 4:end), 'k:','linewidth', 2);
    
    plot(repmat(s(end), [Ng, 1])', RU(end, :)', 'ks', 'linewidth', 2);
    for n = 1:3
        text(s(end) + 0.01, RU(end, n), ['$U_{' num2str(n)  '\rightarrow ' num2str(n+1) '}$'], 'Fontsize', 12, Interpreter='latex')
    end
    text(s(end) + 0.01, RU(end, 4), '...', 'FontSize',14, 'FontWeight','bold')
    xlabel('$s$ ($\mu$m)', Interpreter='latex');ylabel(['$U_{n\rightarrow n+1}(n, s)$  (s$^{-1}$)'], Interpreter='latex');
    gc = gca;
    set(gca, 'FontSize', 12);
    aspect = 1.5;
    % set(gcf, 'Position', [500  300  3.5*96 3.5*96/aspect])
    % gc.Position = gc.Position + [0 0 -0.1 0];
    set(gca, 'TickLength',[0.025 0.025]);
    
    % inset 
    axes('position', [0.62 0.55 0.21 0.35], 'YAxisLocation','right')
    plot(1:Ng, RU(end, :), 'ks-', 'linewidth', 2);
    xlabel('$n$', Interpreter='latex')
    text(4, max(RU(end, 1))*0.7, ["$U_{n\rightarrow n+1}$ at" sprintf("%0.3f $\\mu$m", max(s))], 'FontSize',14, 'FontWeight','normal', Interpreter='latex')
    set(gca, 'FontSize', 12);
    % set(gcf, 'Position', [500  240  400  300])
    % set(gca, 'YAxisLocation', 'right');
    set(gca, 'YTick',[0 200 400], 'XTick', [0 5 10], 'TickLength',[0.05 0.025]);
    exportgraphics(f,'../Figures/ModelRates2.png','Resolution',150)
    exportgraphics(f,'../Figures/ModelRates.eps','Resolution',150)

    % reconstruct Fp again without the NaN's
    % Fp = kp*(max(0,s-slack)/Lref).^(np); 
    % Fp(isnan(Fp() == NaN) = 0;

end
%% Folding rate design and visualization
% RF = alphaF_0 + alphaF*(slack(1:Ng) - s).*(slack(1:Ng) > s);  % folding rates from state n+1 to n            
RF = alphaF_0;

% Initial state
if size(L0X0) == 1
    % only L0 provided
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
    x0 = [x0; L0X0];
else
    % full initial vector provided
    x0 = L0X0;
end

  opts = odeset('RelTol',1e-3, 'AbsTol',1e-2);          
  t = []; x = [];
  F_alt = [];
  for i_section = 1:size(times, 1)
      [t1,x1] = ode15s(@dXdT,times(i_section, :),x0,opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf,velocities{i_section});
      t = [t; t1(2:end)];
      x = [x; x1(2:end, :)];
      % prep the init vector again
      x0 = x1(end,:);
      F_alt1 = zeros(size(t1));
      for i_t = 1:size(t1)
          [~, outpits] = dXdT(t1(i_t),x1(i_t, :)',Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf, velocities{i_section});
          F_alt1(i_t) = outpits(1);
      end
      F_alt = [F_alt;F_alt1(2:end)];
  end
  validRng = t >= 0;
  t = t(validRng);x = x(validRng, :); F_alt  = F_alt(validRng)';
  if ~any(validRng)
      % important only for the init x0
      Time = []; Length = []; Force_INST = []; Force_par = [];
      return;
  end

%%
Time = t;
states = [];states_a = []; strains = []; i_time_snaps = [];
maxPu = 0; maxPa = 0;
% save pca4statesenv
% load pca4statesenv
% save relaxstatesenv
% load relaxstatesenv
% drawAllStates = true;
    if drawAllStates
        % save current figure
        g = gcf;
        % open up a new one
        layout_x = 2 + (pCa < 10);
        f = figure(50+pCa*10+j); clf; tiledlayout(layout_x, 4, 'TileSpacing','compact', Padding='loose');
        fontsize(12, 'points')
        aspect = 2;
        % normal size of 2-col figure on page is 7.2 inches
        % matlab's pixel is 1/96 of an inch
        f.Position = [300 200 7.2*96 7.2*96/aspect];

        % decide for timepoints
        % fixed time or fraction of ramp durations?
        % time_snaps = [0, 0.1, 1, 10, 30, 40, 100]
        % time_snaps = [0, rds(j), rds(j) + 30, 60, 120, 160];
    
        % disable
        i_time_snaps = [];
        for i = 1:length(time_snaps)
            if time_snaps(i) > t
                break;
            end
            i_time_snaps(i) = find(t>=time_snaps(i), 1, 'first');    
        end
        % custom colormap
        load SoHot.mat;
    end

  for i = 1:length(t)
    xi = x(i,:);
    Length(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    if pCa >= 11
        pa = 0;
    else
        pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    end
    Fd = kd*max(0,(Length(i) - s)/L_0).^nd; 
    Force_pa = ds*sum(sum(Fd.*pa ));
    Force_INST(i) =  0*ds*sum(sum(Fd.*pu )) + Force_pa;
    states(i, 1:Ng+1) = sum(pu);
    states_a(i, 1:Ng+1) = sum(pa);
    strains(i, 1:Nx) = sum(pu, 2);

    if drawAllStates
        if any(ismember(i_time_snaps, i))
            % sum of all states, converting to %
            ss = (sum(pu(:)) + sum(pa(:)))*ds;
            maxPu = max(maxPu, sum(pu(:))*ds*100);
            maxPa = max(maxPa, sum(pa(:))*ds*100);
            i_snap = find(i_time_snaps == i, 1, 'first');


            %% 2D shaded plot - unattached
            % next plot
            tl = nexttile(i_snap);            
            imagesc(s, 0:Ng, pu'*ds*100);hold on; box on;
            set(gca,'YDir','normal');
            mesh([s  - ds/2; s(end) + ds/2], (0:Ng+1) - 0.5, zeros(size(pu') + [1 1]), LineStyle='-', EdgeColor=[1 1 1]*0.9, FaceColor='none')
            contour(s, 0:Ng, pu'*ds*100, [0.01 0.01], 'k-', LineWidth=2)
            colormap(SoHot);
            
            clim([0 max(1, round(maxPu))]);
            if i_snap == 1
                ylabel('$n$', Interpreter='latex');
            elseif i_snap == length(i_time_snaps)
                cb = colorbar;
                cb.Ticks = sort(unique([0 50 round(maxPu)]));
                title(cb, 'Unatt %')
                set(tl, "YTick", [])
            else
                set(tl, "YTick", [])
            end
            if pCa >= 10
                % xlabel is done by attached panels
                xlabel('$s$ ($\mu$m)', Interpreter='latex');
            else
                ga = gca;
                ga.XTick = [];
            end 
            fontsize(12, 'points')
            %% 2D shaded plot - attached
            if pCa < 10
                aspect = 1.4; % make the plots rectangular for pca as well
                f.Position = [300 200 7.2*96 7.2*96/aspect];
                tl = nexttile(i_snap + 4);
                imagesc(s, 0:Ng, pa'*ds*100);hold on; box on;
                set(gca,'YDir','normal');
                mesh([s  - ds/2; s(end) + ds/2], (0:Ng+1) - 0.5, zeros(size(pa') + [1 1]), LineStyle='-', EdgeColor=[1 1 1]*0.9, FaceColor='none')
                contour(s, 0:Ng, pa'*ds*100, [0.01 0.01], 'k-', LineWidth=2)

                colormap(SoHot);
                clim([0 round(maxPa)]);
                if i_snap == 1
                    ylabel('$n$', Interpreter='latex')
                elseif i_snap == length(i_time_snaps)
                    cb = colorbar;
                    cb.Ticks = [0 round(maxPa)];
                    title(cb, 'Att %')
                    set(tl, "YTick", [])
                else
                    set(tl, "YTick", [])
                end
                xlabel('$s$ ($\mu$m)', Interpreter='latex');
                box on;
                fontsize(12, 'points')
            end
        end
    end    
  end
  % add parallel static force
  b = mod(11);
  c = mod(12);
  d = mod(13);
  
  % calculate a, so that the max value is the same    
    a = (Fss - d)/((Lmax -b)^c);
  Force_par = a*max(Length - b, 0).^c + d;
  
  % calc instantenous force
  Force_INST = Force_INST + Force_par; 

  % unfolding force
  Force_UNF = F_alt + Force_par;

  if drawAllStates
        nexttile([1 4]);
        semilogx(Time, Force_INST, 'k-');hold on;
        semilogx(Time, Force_UNF, 'r--');hold on;
        scatter(Time(i_time_snaps), Force_INST(i_time_snaps), 'ko', 'filled');
        % xlim([1e-3 max(1300, rds(j)*10) + 30]);
        % nexttile([1 1]);
        % loglog(Time, Force);hold on;
        % scatter(Time(i_time_snaps), Force(i_time_snaps), 'o', 'filled');
        % xlim([1e-3 inf])
        xlabel('$t$ (s)', Interpreter='latex');
        legend('Tension', 'Insets', 'Location','northwest');
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
        fontsize(12, 'points')
        vmax = max([velocities{:}]);
        % exportgraphics(f,sprintf('../Figures/States%g_%gs.png', pCa, vmax),'Resolution',150)
    end