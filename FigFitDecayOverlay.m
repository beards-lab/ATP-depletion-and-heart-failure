
% Data daved by save('pca4dataAdj60s.mat', 'Tarr', 'Farr') from
% AverageRamps[Ca].m
% load pca4data.mat
% load pca4dataAdj.mat
% load pca4dataAdj60s.mat
load pca11data.mat


fitfun = @(init) evalPowerFit(init, Farr, Tarr, true);
figure(2);clf;
% best fit for pca 11
x = [4.4271    0.2121    4.8964]

% best fit for pca 4.4 - 12.7 t^-1.16, using only first 10s
% x = [12.7405    1.1640   15.9898];

% best fit for pca4.4: limit the Fss to that of pca 11, using only first 10s.
% Surprisingly, the time constant is nearly same!
% note: run the AverageRampsCa to get the Ca Farr
% x = [18.8677    0.2162   4.8964];

% experiment: free Fss, using only 0.1 and 1s ramps for 10s to avoid
% messing with the remaining force buildup. The optimizer pushes the Fss
% high
% x = [16.9548    0.2781    7.0702];

% experiment: fixed Fss, using only 0.1 and 1s ramps for 10s. No big
% difference from the all ramps fit.
% x = [19.0119    0.2407    4.8964];

% Best fit for pca4.4 with Frem correction, no assumptions
% x = [10.4953    0.5869   9];
% x = [10.4953    0.21   11];

% pCa 4.4 with Frem correction, for only 0.1 and 1s ramps for 30s
% x = [11.2095    0.6120   12.2428]

% pCa 4.4 with Frem correction, assuming Fss
% x= [17.0511    0.2176    4.8964];

% % pCa 4.4, optimizing just the tail
% x = [4.3209    0.2100   13.6246];

% limit pCa time constant to the relaxed time constant
% x = [8.6226    0.2100   11.3698];
% no Ca
aspect = 1.5;
f = figure(2);
% normal size of 2-col figure on page is 7.2 inches
% matlab's pixel is 1/96 of an inch
f.Position = [300 200 7.2*96 7.2*96/aspect];
x = [4.4271    0.2121    4.8964];
[c rampShift] = fitfun(x)
f = gcf();
% saveas(f, 'Figures/FigDecayOverlay.png')
% exportgraphics(f,'Figures/FigDecayOverlay7.2.png','Resolution',300)
% saveas(f, 'Figures/FigDecayOverlaypCa4.png')
% exportgraphics(f,'Figures/FigDecayOverlaypCa4.4_7.2.png','Resolution',300)
% saveas(f, 'Figures/FigDecayOverlaypCa4Corr2.png')
%% pCa
aspect = 1.5;
% rampShift = [5.3980    0.8234    0.2223   0.0100];
load('pCa4dataNoAdj60sFremCorr.mat')
x = [4.1240    0.2121   12.0286];
% load('pCa4dataNoAdj60sFremCorrShifted.mat')
% x = [4.0648    0.2121   12.0505];
f = figure(4);
f.Position = [300 200 7.2*96 7.2*96/aspect];
% keep the b, fit a and Tss
% x = [4.1237    0.2121   12.0289];
% x = [4.1233    0.2121   12.0290];
[c rspca] = evalPowerFit(x, Farr, Tarr, true, rampShift, true)

%%
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 150);
% init(1:3) = [0 0 0]
%% fit the no Ca
init = x;
fitfunOpt = @(init) evalPowerFit(init, Farr, Tarr, false);
x = fminsearch(fitfunOpt, init, options)

[c rampShift] = fitfun(x);
%% fit pCa - fixed power exponent b
% init = [x(1) x(3)]
init = x;
pcaFitFun = @(x)evalPowerFit([x(1), init(2), x(2)], Farr, Tarr, false, rampShift, true);
x = fminsearch(pcaFitFun, [init(1) init(3)], options)
init([1 3]) = x;
x = init;
c = evalPowerFit(x, Farr, Tarr, true, rampShift, true)


% pcaFitFun = @(x)evalPowerFit(x, Farr, Tarr, false, rampShift, true);
% x = fminsearch(pcaFitFun, init, options)
% x = init;


function [cost_c0 rampShift] = evalPowerFit(params, Farr, Tarr, plotResults, rampShift, pCa)
%%
if any(params < 0)
    cost_c0 = Inf;
    rampShift = [];
    return;
end
% disp('ev..')
if nargin < 4
    plotResults = false;
    rampShift = [];
    pCa = false;
elseif nargin < 5 
    rampShift = [];
    pCa = false;
end

rampShiftFixed = true;
if isempty(rampShift)
    rampShift = [0 0 0 0];
    rampShiftFixed = false;
end

a = params(1);
b = params(2);
c = params(3);

%% biuld axes
if plotResults
    % set(groot,'CurrentFigure',2); % replace figure(indFig) without stealing the focus
    fug = gcf();
    clf;
    % set width adn height and offsets
    w = 0.255; h = 0.85; x0 = 0.08; y0 = 0.11; gap = 0.003; centerShift = 0.06;
    % data ranges - pca11
    xrng1 = [-10 19]; xrng2 = [-2 27];xrng_add = [55 65]; yrng = [0.8 21];
    % ylogtick = [1e0 1e1];
    xlintick = 0:10:20;xlinticklab  = [string(floor((xlintick(1):10:xlintick(end))/5)*5) ""];
    ylogtick = [1 ceil((yrng(1):5:yrng(end))/5)*5];
    % data ranges - pca4
    if pCa
        xrng1 = [-10 60]; xrng2 = [-2 60]; xrng_add = [0 0]; yrng = [0.8 35];
        xlintick = 0:20:60;xlinticklab = string(floor((xlintick(1):20:xlintick(end))/10)*10);
        ylogtick = [1 ceil((yrng(1):10:yrng(end))/10)*10];
    end
    
    % Font size
    fs = 12;
    
    [w1, w2] = getW1W2(w, gap, xrng1, xrng_add);
    
    % left
    % axes left main
    a_lm = axes(fug, 'Position', [x0 y0 w1 h], 'XLim',xrng1,'XTick',floor((xrng1(1):20:xrng1(end))/5)*5, ... 
        'YLim',[0 yrng(2)], 'FontSize',fs, 'YScale','linear', 'TickLabelInterpreter','latex'); 
    ylabel("$T$ (kPa)", 'Interpreter','latex'); xlabel("$t - t_r$ (s)", 'Interpreter', 'latex')
    hold on; box on; % keep the settings from overwritting
    % axes left adjacent
    if w2 > 0
        a_la = axes(fug, 'Position', [x0+w1+gap y0 w2 h], 'Xlim', xrng_add,'XTick',floor((xrng_add(1):10:xrng_add(end))/5)*5, 'YTick', [],... 
            'YLim',[0 yrng(2)], 'FontSize',fs, 'YScale','linear', 'TickLabelInterpreter','latex'); 
        hold on;box on;
    else
        % invalid, lets make it outside valid range
        a_ca = axes('Position', [-1, 0, 0, 0]);
    end
    
    % center
    [w1, w2] = getW1W2(w, gap, xrng2, xrng_add);
    a_cm = axes('Position', [0.5-w/2 + centerShift y0 w1 h], 'XLim',xrng2,'XTick', xlintick, 'XTickLabel',xlinticklab, ... 
        'YLim',yrng, 'FontSize',fs, 'YScale','log', 'TickLabelInterpreter','latex',...
        YTick=ylogtick); 
    ylabel("$T-T_{ss}$ (kPa)", 'Interpreter','latex');xlabel("$t - t_r + \tau_i$ (s)", 'Interpreter', 'latex');
    hold on; box on;
    
    if w2 > 0
        a_ca = axes('Position', [0.5-w/2 + w1 + gap + centerShift y0 w2 h], 'Xlim', xrng_add,'XTick',floor((xrng_add(1):10:xrng_add(end))/5)*5, 'YTick', [],... 
        'YLim',yrng, 'FontSize',fs, 'YScale','log', 'TickLabelInterpreter','latex');
        hold on; box on;
    else
        % invalid, lets make it outside valid range
        a_la = axes('Position', [-1, 0, 0, 0]);
    end
    
    
    % right
    a_r = axes('Position', [1-x0-w+x0/2  y0 w h],'XLim',[1e-2 1e4], ... 
        'YLim',yrng, 'FontSize',fs, 'XScale', 'log','YScale','log', 'TickLabelInterpreter','latex',...
        YTick=ylogtick, YTickLabel=[], XTick=[1e-1 1e1 1e3]); 
    % ylabel("$T-T_{ss}$ (kPa)", 'Interpreter','latex');
    xlabel("$t - t_r + \tau_i$ (s)", 'Interpreter', 'latex');
    hold on; box on;
end
%
cost_c0 = 0;
cg = gray(5);
% clin = lines(5);
rds = [100.0000   10.0000    1.0000    0.1000];
% set(groot,'CurrentFigure',2); % replace figure(indFig) without stealing the focus
% figure(2);
param_a = string();
for i_rds = [4 3 2 1]
        %%
        Favg = movmean(Farr{i_rds}, [0 0]) - c;
        % loglog(Tarr{i_rds} + rampShift(i_rds), Favg -c, Color=[cg(5-i_rds, :)], LineWidth=5-i_rds);hold on;
    
        % resample log equally from the peak - till the end
        t_s = logspace(log10(rds(i_rds)), log10(Tarr{i_rds}(end)), (Tarr{i_rds}(end) - rds(i_rds))*10);

        if pCa
            % only for 10s to avoid remaining tension build-up
            % t_s = logspace(log10(rds(i_rds)), log10(rds(i_rds) + 30), 100);
            % limit 10s up to fit just the tail for pca
            t_s = logspace(log10(rds(i_rds) + 15-rampShift(i_rds)), log10(min(Tarr{i_rds}(end), rds(i_rds) + 60 - rampShift(i_rds))), 100);
        end

        % force interpolation
        Fint = interp1(Tarr{i_rds}, Favg, t_s, "pchip", 'extrap');
        % t_s_ru = [logspace(log10(1e-2), log10(rds(i_rds)), 10) t_s(2:end)];% including ramp-up
        % t_s_ru = [linspace(0, rds(i_rds), 20) t_s(2:end)];% including ramp-up
        t_s_ru = [0:0.01:(t_s(end))];% including ramp-up
        Fint_ru = interp1(Tarr{i_rds}, Favg, t_s_ru, "pchip", 'extrap');
    
        try
            if rampShiftFixed
                % prescribed rampShift - using the fit only to
                % calculate the costs
                tau0 = rampShift(i_rds);
                pf = @(tau, x) a*max(1e-9, (x-rds(i_rds) + tau0)).^(-b) + tau*0;
                [ae, goodness] = fit(t_s(2:end)', Fint(2:end)', pf, 'StartPoint', [rds(i_rds)/10], 'Lower',[0], 'Upper',[rds(i_rds)]);
            else
                % PowerFunction
                pf = @(tau, x) a*max(1e-9, (x-rds(i_rds) + tau)).^(-b);
                [ae, goodness] = fit(t_s(2:end)', Fint(2:end)', pf, 'StartPoint', [rds(i_rds)/10], 'Lower',[0], 'Upper',[rds(i_rds)]);
                rampShift(i_rds) = ae.tau;
            end
        catch e
            disp(e)
        end        
        cost_c0 = cost_c0 + goodness.rmse;%/length(t_s);
        if plotResults
            cc = [cg(5-i_rds, :)]; % current color 
            clw = 4.5-i_rds; % current line width
            axes(a_lm);
            l_lm(i_rds) = plot(t_s_ru + rampShift(i_rds)*0 - rds(i_rds), Fint_ru + c, '-', Color=cc, LineWidth=clw);
            axes(a_la); plot(t_s_ru + rampShift(i_rds)*0 - rds(i_rds), Fint_ru + c, '-', Color=cc, LineWidth=clw);
            % plot(t_ext+rampShift(i_rds)-rds(i_rds), c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
            axes(a_r);
            loglog(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru +0*c, '-', Color=[cc], LineWidth=clw);
            loglog(t_s + rampShift(i_rds) - rds(i_rds), Fint +0*c, '-', Color=cc, LineWidth=clw);
            xp = t_s + rampShift(i_rds) - rds(i_rds);            
            axes(a_cm);
            l_cm(i_rds) = semilogy(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru + c*0, '-', Color=cc, LineWidth=clw);
            % semilogy(t_s + rampShift(i_rds) - rds(i_rds), Fint +0*c, 'b-', LineWidth=clw);

            axes(a_ca);
            semilogy(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru + c*0, '-', Color=cc, LineWidth=clw);
            % semilogy(t_ext+rampShift(i_rds)-rds(i_rds), 0*c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
        end
        % param_a(i_rds) = sprintf('%0.1f: %0.3fs', rds(i_rds), ae.tau);
    end
    if ~plotResults
        return;
    end
    % plot the fit
    axes(a_r);
    
    % fit area
    if pCa
        l_fitarr = fill([xp(1) xp(end) xp(end) xp(1)], [0.1 0.1 60 60], [0.1 0.1 0.1], 'FaceAlpha',0.14, EdgeColor='none');
    end

    pf = @(x) a*(x).^(-b);
    if true || ~pCa
        % fit area right from the beginning
        t_ext = logspace(log10(1e-2), log10(1*60*60), 100);    
    else
        % fit area from middle
        t_ext = logspace(log10(t_s(2) + rampShift(i_rds) - rds(i_rds)), log10(1*60*60), 100);    
    end

    pf_v = pf(t_ext);
    cl = lines(2);
    l_f = loglog(t_ext, pf_v, 'r:', ... %Color=cl(2, :)
        LineWidth=3);
    
    % already in master legend
    % if ~pCa
    %     legend(l_f, sprintf('$%0.2f(t - t_r + \\tau_i)^{-%0.2f}$',a,b), ...
    %     'Interpreter','latex', 'FontSize',fs, 'Location', 'northeast');
    % else
    %     legend([l_fitarr l_f], 'Fit area', sprintf('$%0.2f(t - t_r + \\tau_i)^{-%0.2f}$',a,b), ...
    %     'Interpreter','latex', 'FontSize',fs, Location='northeast');
    % end

    %
    axes(a_la);
    plot([-10 100], [c, c], 'k:', LineWidth=3);    
    axes(a_lm);
    l_tss = plot([-10 100], [c, c], 'k:', LineWidth=3);
    text(2, c-1.5, sprintf('$T_{ss} = %0.2f$', c), 'Interpreter', 'Latex', 'FontSize',fs+2, 'FontWeight','bold')

    % title('A: Linear axis plot of a_i(x-r_d - \tau_{0, i})^{\tau_i} + Fss')
    valids = isgraphics(l_lm);
    legends = {'$t_r = 100$ s','$t_r = 10$ s','$t_r = 1$ s','$t_r = 0.1$ s','$T_{ss}$'};
    leg = legend([l_lm(valids) l_tss], legends([valids true]), 'Interpreter','latex', 'FontSize',fs);
    %adjust position
    if w1 > 0
        % get top-right corner
        ca_pos = get(a_la, 'Position');
    else
        ca_pos = get(a_lm, 'Position');
    end    
    l_pos = leg.Position;
    leg.Position = [ca_pos(1) + ca_pos(3) - l_pos(3) - 0.01 l_pos(2:4)];
    

    % title('C: Log-log plot a_i(x-r_d - \tau_{0, i})^{\tau_i}, extrapolated to 24hrs')

    
    % fit area - only pCa
    if pCa
        axes(a_cm)
        l_fitarr = fill([xp(1) xp(end) xp(end) xp(1)], [0.1 0.1 60 60], [0.1 0.1 0.1], 'FaceAlpha',0.14, EdgeColor='none');
    end
    
    % show the fit in central too
    axes(a_cm);
    l_f = semilogy(t_ext, pf_v, 'r:', ... %Color=cl(2, :)
        LineWidth=3);

    valids = isgraphics(l_cm);
    legends = { sprintf('$t_r = 100$ s,$\\tau_i=%0.2f$ s', rampShift(1)),...
                sprintf('$t_r = 10$ s, $\\tau_i=%0.2f$ s', rampShift(2)),...
                sprintf('$t_r = 1$ s,  $\\tau_i=%0.2f$ s', rampShift(3)),...
                sprintf('$t_r = 0.1$ s,$\\tau_i=%0.2f$ s', rampShift(4))};
    
    if ~pCa
        leg_gr = [l_cm(valids) l_f];
        [legends(valids) sprintf('$%0.2f(t - t_r + \\tau_i)^{-%0.2f}$',a,b)];
        reorder = [1 3 5 2 4];
    else
        leg_gr = [l_cm(valids) l_fitarr l_f];
        leg_txt = [legends(valids) "Power-law tail" sprintf('$%0.2f(t - t_r + \\tau_i)^{-%0.2f}$',a,b)];
        reorder = [1 3 5 2 4 6];
    end
    leg = legend(leg_gr(reorder), leg_txt(reorder), 'Interpreter','latex', 'FontSize',fs, Location='northwest', NumColumns=2);
    leg.ItemTokenSize=[15; 18];
    %adjust position
    if w2 > 0
        % get top-right corner
        ca_pos = get(a_ca, 'Position');
    else
        ca_pos = get(a_cm, 'Position');
    end    
    l_pos = leg.Position;
    leg.Position = [l_pos(1) l_pos(2) + 0.03 l_pos(3:4)];

    axes(a_ca);
    semilogy(t_ext, pf_v, 'r:', LineWidth=3);

    uistack(a_cm, "top"); % make sure the central is not overlapped y

    % loglog(xl, [c, c], 'r--');
    % xlim(xl);
    % yl = ylim();
    % ylim([yl(1)*0.9, yl(2)])
    
    % legend([l_f], sprintf('%0.1ft^{-%0.2f}', a, b))
    
    % legend([l_d l_f], 'Ramp-up 100s', 'Ramp-up 10s', 'Ramp-up 1s', 'Ramp-up 0.1s', num2str(param_a(1)), num2str(param_a(2)),num2str(param_a(3)),num2str(param_a(4)), 'True F_{ss}')
    % legend([l_d(3:4) l_f(3:4)], 'Ramp-up 1s', 'Ramp-up 0.1s', num2str(param_a(3)), num2str(param_a(4)), 'True F_{ss}')
    % title(num2str(cost_c0))

    % rampShift
end

function [w1, w2] = getW1W2(w_total, gap, xrng_1, xrng_2)
% calculates width of two plots with the same x-scale
    dpi = (w_total-gap)/(diff(xrng_1) + diff(xrng_2)); % 
    w1 = diff(xrng_1)*dpi; w2 = diff(xrng_2)*dpi;
end
