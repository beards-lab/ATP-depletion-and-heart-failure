
% load pca4data.mat
load pca4dataAdj.mat
% load pca11data.mat

fitfun = @(init) evalPowerFit(init, Farr, Tarr);
figure(2);clf;
% best fit for pca 11
% x = [4.4271    0.2121    4.8964]

% best fit for pca 4.4 - 12.7 t^-1.16, using only first 10s
% x = [12.7405    1.1640   15.9898];

% best fit for pca4.4: limit the Fss to that of pca 11, using only first 10s.
% Surprisingly, the time constant is nearly same!
% note: run the AverageRampsCa to get the Ca Farr
x = [18.8677    0.2162   4.8964];

% experiment: free Fss, using only 0.1 and 1s ramps for 10s to avoid
% messing with the remaining force buildup. The optimizer pushes the Fss
% high
% x = [16.9548    0.2781    7.0702];

% experiment: fixed Fss, using only 0.1 and 1s ramps for 10s. No big
% difference from the all ramps fit.
% x = [19.0119    0.2407    4.8964];

% Best fit for pca4.4 with Frem correction, no assumptions
x = [10.4953    0.5869   12.1112];

% pCa 4.4 with Frem correction, for only 0.1 and 1s ramps for 30s
% x = [11.2095    0.6120   12.2428]

% pCa 4.4 with Frem correction, assuming Fss
% x= [17.0511    0.2176    4.8964];


c = fitfun(x)
f = gcf();
% saveas(f, 'Figures/FigDecayOverlay.png')
% saveas(f, 'Figures/FigDecayOverlaypCa4.png')
saveas(f, 'Figures/FigDecayOverlaypCa4Corr2.png')

%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.001, 'PlotFcns', @optimplotfval, 'MaxIter', 50);
% init(1:3) = [0 0 0]
init = x;
x = fminsearch(fitfun, init, options)


function cost_c0 = evalPowerFit(params, Farr, Tarr)
%%
plotResults = true;
% rampShift = [params(1:4) ]; 
a = params(1);
b = params(2);
c = params(3);
% c = 4.8964;
% rampShift = [0 0 0 0];

%% biuld axes
if plotResults
    set(groot,'CurrentFigure',2); % replace figure(indFig) without stealing the focus
    clf;
    % set width adn height and offsets
    w = 0.26; h = 0.85; x0 = 0.05; y0 = 0.1; gap = 0.003; centerShift = 0.04;
    % data ranges - pca11
    xrng1 = [-10 19]; xrng2 = [-2 27];xrng_add = [55 60]; yrng = [0.8 18];
    % ylogtick = [1e0 1e1];
    ylogtick = [1 ceil((yrng(1):5:yrng(end))/5)*5];
    % data ranges - pca4
    xrng = [-10 17]; xrng_add = [0 0]; yrng = [0.5 50];
    ylogtick = [1 ceil((yrng(1):10:yrng(end))/10)*10];
    
    % Font size
    fs = 18;
    
    [w1, w2] = getW1W2(w, gap, xrng1, xrng_add);
    
    % left
    % axes left main
    a_lm = axes('Position', [x0 y0 w1 h], 'XLim',xrng1,'XTick',floor((xrng1(1):5:xrng1(end))/5)*5, ... 
        'YLim',[0 yrng(2)], 'FontSize',fs, 'YScale','linear', 'TickLabelInterpreter','latex'); 
    ylabel("$T$ (kPa)", 'Interpreter','latex'); xlabel("$t - t_r$ (s)", 'Interpreter', 'latex')
    hold on; box on; % keep the settings from overwritting
    % axes left adjacent
    if w2 > 0
        a_la = axes('Position', [x0+w1+gap y0 w2 h], 'Xlim', xrng_add,'XTick',floor((xrng_add(1):5:xrng_add(end))/5)*5, 'YTick', [],... 
            'YLim',[0 yrng(2)], 'FontSize',fs, 'YScale','linear', 'TickLabelInterpreter','latex'); 
        hold on;box on;
    else
        % invalid, lets make it outside valid range
        a_ca = axes('Position', [-1, 0, 0, 0]);
    end
    
    % center
    [w1, w2] = getW1W2(w, gap, xrng2, xrng_add);
    a_cm = axes('Position', [0.5-w/2 + centerShift y0 w1 h], 'XLim',xrng2,'XTick', -5:5:25, 'XTickLabel',[string(floor((xrng2(1):5:xrng2(end))/5)*5) ""], ... 
        'YLim',yrng, 'FontSize',fs, 'YScale','log', 'TickLabelInterpreter','latex',...
        YTick=ylogtick); 
    ylabel("$T-T_{ss}$ (kPa)", 'Interpreter','latex');xlabel("$t - t_r + \tau_i$ (s)", 'Interpreter', 'latex');
    hold on; box on;
    
    if w2 > 0
        a_ca = axes('Position', [0.5-w/2 + w1 + gap + centerShift y0 w2 h], 'Xlim', xrng_add,'XTick',floor((xrng_add(1):5:xrng_add(end))/5)*5, 'YTick', [],... 
        'YLim',yrng, 'FontSize',fs, 'YScale','log', 'TickLabelInterpreter','latex');
        hold on; box on;
    else
        % invalid, lets make it outside valid range
        a_la = axes('Position', [-1, 0, 0, 0]);
    end
    
    
    % right
    a_r = axes('Position', [1-x0-w  y0 w h],'XLim',[1e-2 1e4], ... 
        'YLim',yrng, 'FontSize',fs, 'XScale', 'log','YScale','log', 'TickLabelInterpreter','latex',...
        YTick=ylogtick, YTickLabel=[], XTick=[1e-1 1 1e1 1e2 1e3 1e4]); 
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
        % only for 10s to avoid remaining tension build-up
        % t_s = logspace(log10(rds(i_rds)), log10(rds(i_rds) + 30), 100);
        
        % force interpolation
        Fint = interp1(Tarr{i_rds}, Favg, t_s, "pchip", 'extrap');
        % t_s_ru = [logspace(log10(1e-2), log10(rds(i_rds)), 10) t_s(2:end)];% including ramp-up
        % t_s_ru = [linspace(0, rds(i_rds), 20) t_s(2:end)];% including ramp-up
        t_s_ru = [0:0.01:rds(i_rds) t_s(2:end)];% including ramp-up
        Fint_ru = interp1(Tarr{i_rds}, Favg, t_s_ru, "pchip", 'extrap');
    
        % tau = rampShift(i_rds);
        % PowerFunction
        pf = @(tau, x) a*max(1e-9, (x-rds(i_rds) + tau)).^(-b);
        % ef = @(a, b, x) a*b^(x-rds(i_rds) + rampShift(i_rds));
        try
            [ae goodness] = fit(t_s(2:end)', Fint(2:end)', pf, 'StartPoint', [rds(i_rds)/10], 'Lower',[0], 'Upper',[Inf]);
            % [ae goodness] = fit(t_s(2:end)', Fint(2:end)', ef, 'StartPoint', [1, 2]);%, 'Lower',[0 0 0 -Inf], 'Upper',[Inf 1 10, Inf]);
        catch e
            disp(e)
        end
        rampShift(i_rds) = ae.tau;
        cost_c0 = cost_c0 + goodness.rmse;%/length(t_s);
        % l_f(i_rds) = loglog(t_ext+rampShift(i_rds)-rds(i_rds), 0*c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
        if plotResults
            cc = [cg(5-i_rds, :)]; % current color 
            clw = 4.5-i_rds; % current line width
            axes(a_lm);
            l_lm(i_rds) = plot(t_s_ru + rampShift(i_rds)*0 - rds(i_rds), Fint_ru + c, '-', Color=cc, LineWidth=clw);
            axes(a_la); plot(t_s_ru + rampShift(i_rds)*0 - rds(i_rds), Fint_ru + c, '-', Color=cc, LineWidth=clw);
            % plot(t_ext+rampShift(i_rds)-rds(i_rds), c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
            axes(a_r);
            loglog(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru +0*c, '-', Color=cc, LineWidth=clw);
            axes(a_cm);
            l_cm(i_rds) = semilogy(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru + c*0, '-', Color=cc, LineWidth=clw);
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
    pf = @(x) a*(x).^(-b);
    t_ext = logspace(log10(1e-2), log10(1*60*60), 100);    
    pf_v = pf(t_ext);
    cl = lines(2);
    l_f = loglog(t_ext, pf_v, ':', Color=cl(2, :), LineWidth=4);
    
    legend(l_f, sprintf('$%0.2f(t - t_r + \\tau_i)^{-%0.2f}$',a,b), ...
        'Interpreter','latex', 'FontSize',fs);

    %
    axes(a_la);
    plot([-10 100], [c, c], 'k--', LineWidth=4);    
    axes(a_lm);
    l_tss = plot([-10 100], [c, c], 'k--', LineWidth=4);
    text(2, c-1.5, sprintf('$T_{ss} = %0.2f kPa$', c), 'Interpreter', 'Latex', 'FontSize',fs+5, 'FontWeight','bold')

    % title('A: Linear axis plot of a_i(x-r_d - \tau_{0, i})^{\tau_i} + Fss')
    valids = isgraphics(l_lm);
    legends = {'$t_r = 100s$','$t_r = 10s$','$t_r = 1s$','$t_r = 0.1s$','$T_{ss}$'};
    legend([l_lm(valids) l_tss], legends([valids true]), 'Interpreter','latex', 'FontSize',fs);
    

    % title('C: Log-log plot a_i(x-r_d - \tau_{0, i})^{\tau_i}, extrapolated to 24hrs')

    axes(a_cm);
    valids = isgraphics(l_cm);
    legends = {sprintf('$t_r = 100s, \\tau_i=%0.2fs$', rampShift(1)),...
                sprintf('$t_r = 10s, \\tau_i=%0.2fs$', rampShift(2)),...
                sprintf('$t_r = 1s, \\tau_i=%0.2fs$', rampShift(3)),...
                sprintf('$t_r = 0.1s, \\tau_i=%0.2fs$', rampShift(4))};
    legend([l_cm(valids)], legends(valids), 'Interpreter','latex', 'FontSize',fs);

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
