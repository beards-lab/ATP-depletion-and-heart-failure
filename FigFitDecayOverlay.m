
load pca4data.mat

fitfun = @(init) evalPowerFit(init, Farr, Tarr);
figure(2);clf;
% best fit for pca 11
x = [4.4271    0.2121    4.8964]
% best fit for pca 4.4 - 12.7 t^-1.16, using only first 10s
x = [12.7405    1.1640   15.9898];
% best fit for pca4.4: limit the Fss to that of pca 11, using only first 10s.
% Surprisingly, the time constant is nearly same!
% note: run the AverageRampsCa to get the Ca Farr
x = [18.8677    0.2162   4.8964];
% experiment: free Fss, using only 0.1 and 1s ramps for 10s to avoid
% messing with the remaining force buildup. The optimizer pushes the Fss
% high
x = [16.9548    0.2781    7.0702];
% experiment: fixed Fss, using only 0.1 and 1s ramps for 10s. No big
% difference from the all ramps fit.
x = [19.0119    0.2407    4.8964];
c = fitfun(x)
%% biuld axes
clf;
w = 0.25;
xrng = [-3 17];
gap = 0.005;
xrng_add = [25 30];
yrng = [1e-2 20];
h = 0.85;
dpi = (w-gap)/(diff(xrng) + diff(xrng_add)); % 
w1 = diff(xrng)*dpi
w2 = diff(xrng_add)*dpi
% 20 + 5 = dpi*(w1+w2)
% w1 + w2 + g = w
% 20*dpi = w1
% 5*dpi = w2
x0 = 0.05;
y0 = 0.08;
fs = 14;

% left
a_lm = axes('Position', [x0 y0 w1 h], 'XLim',xrng, 'YLim',[0 yrng(2)]); % axes left main
set(gca, 'FontSize', fs);hold on;
a_la = axes('Position', [x0+w1+gap y0 w2 h], 'Xlim', xrng_add, 'YLim',[0 yrng(2)]); % axes left adjacent
set(gca, 'FontSize', fs);hold on;
% center
axes('Position', [0.5-w/2 y0 w1 h])
set(gca, 'FontSize', fs)
axes('Position', [0.5-w/2 + w1 + gap y0 w2 h])
set(gca, 'FontSize', fs)
% right
axes('Position', [1-x0-w y0 w h])
set(gca, 'FontSize', fs)

% 
axes(a_lm);
plot(outT, outFr);
axes(a_la);
plot(outT, outFr)
%% delete the properties of broken axis
xlim([-2 18])
xticks([0 5 10 15])
xlim([55 60])
ylabel('')
xlabel('')
yticks([]);
xticks([55 60])

%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.001, 'PlotFcns', @optimplotfval, 'MaxIter', 50);
% init(1:3) = [0 0 0]
init = x;
x = fminsearch(fitfun, init, options)

function cost_c0 = evalPowerFit(params, Farr, Tarr)
%%
sx = 1; sy =3;
plotResults = true;
% rampShift = [params(1:4) ]; 
a = params(1);
b = params(2);
c = params(3);
% c = 4.8964;
% b = params(4); 
% c = params(5);
% rampShift = [0 0 0 0];
    cost_c0 = 0;
    cg = gray(6);
    clin = lines(5);
    rds = [100.0000   10.0000    1.0000    0.1000];
    set(groot,'CurrentFigure',2); % replace figure(indFig) without stealing the focus
    % figure(2);
    clf;
    param_a = string();
    for i_rds = 4:-1:3
        %%
        Favg = movmean(Farr{i_rds}, [0 0]) - c;
        % loglog(Tarr{i_rds} + rampShift(i_rds), Favg -c, Color=[cg(5-i_rds, :)], LineWidth=5-i_rds);hold on;
    
        % resample log equally from the peak
        % t_s = logspace(log10(rds(i_rds)), log10(Tarr{i_rds}(end)), (Tarr{i_rds}(end) - rds(i_rds))*10);
        t_s = logspace(log10(rds(i_rds)), log10(rds(i_rds) + 10), 100);
        % force interpolation
        Fint = interp1(Tarr{i_rds}, Favg, t_s, "pchip", 'extrap');
        % t_s_ru = [logspace(log10(1e-2), log10(rds(i_rds)), 10) t_s(2:end)];% including ramp-up
        t_s_ru = [linspace(0, rds(i_rds), 20) t_s(2:end)];% including ramp-up
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
            set(groot,'CurrentFigure',2); % replace figure(indFig) without stealing the focus
            subplot(sx, sy,1);
            l_d(i_rds) = plot(t_s_ru + rampShift(i_rds)*0 - rds(i_rds), Fint_ru + c, '-', Color=[cg(6-i_rds, :)], LineWidth=5-i_rds);hold on;
            % plot(t_ext+rampShift(i_rds)-rds(i_rds), c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
            subplot(sx, sy,3);
            loglog(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru +0*c, '-', Color=[cg(6-i_rds, :)], LineWidth=5-i_rds);hold on;
            subplot(sx, sy,2);
            semilogy(t_s_ru + rampShift(i_rds) - rds(i_rds), Fint_ru + c*0, '-', Color=[cg(6-i_rds, :)], LineWidth=5-i_rds);hold on;
            % semilogy(t_ext+rampShift(i_rds)-rds(i_rds), 0*c+pf(ae.a, ae.b, t_ext), '--', Color=[clin(5-i_rds, :)], LineWidth=4);
        end
        % param_a(i_rds) = sprintf('%0.1f: %0.3fs', rds(i_rds), ae.tau);
    end
    subplot(sx, sy,3);
    pf = @(x) a*(x).^(-b);
    t_ext = logspace(log10(1e-2), log10(24*60*60), 100);    
    pf_v = pf(t_ext);
    l_f(i_rds) = loglog(t_ext, pf_v, ':', Color=[clin(1, :)], LineWidth=4);

    %%just for the legend
    if ~plotResults
        return;
    end
    rampShift
    subplot(sx, sy,1);
    xlim([-5 30])
    ylim([0 18])
    xl = xlim();
    l_g = loglog(xl, [c, c], 'r--');    
    % title('A: Linear axis plot of a_i(x-r_d - \tau_{0, i})^{\tau_i} + Fss')

    % legend([l_d l_g], 'Ramp-up 100s', 'Ramp-up 10s', 'Ramp-up 1s', 'Ramp-up 0.1s', 'Fss');
    
    xlabel("Time (s)");
    ylabel("Tension (kPa)")
    set(gca, 'FontSize', 14);    
    
    subplot(sx, sy,3);
    xlim([1e-2 1e5])
    ylim([min(pf_v) max(pf_v)])
    yl = ylim();    
    % title('C: Log-log plot a_i(x-r_d - \tau_{0, i})^{\tau_i}, extrapolated to 24hrs')
    % loglog(xl, [c, c], 'r--');
    % xlim(xl);
    % yl = ylim();
    % ylim([yl(1)*0.9, yl(2)])
    
    % legend([l_f], sprintf('%0.1ft^{-%0.2f}', a, b))
    
    % legend([l_d l_f], 'Ramp-up 100s', 'Ramp-up 10s', 'Ramp-up 1s', 'Ramp-up 0.1s', num2str(param_a(1)), num2str(param_a(2)),num2str(param_a(3)),num2str(param_a(4)), 'True F_{ss}')
    % legend([l_d(3:4) l_f(3:4)], 'Ramp-up 1s', 'Ramp-up 0.1s', num2str(param_a(3)), num2str(param_a(4)), 'True F_{ss}')
    % title(num2str(cost_c0))
    xlabel('Time + \tau_i (s)');
    ylabel('Tension - F_{ss} (kPa)') 
    set(gca, 'FontSize', 14);

    subplot(sx, sy, 2);
    % l_f(length(l_f)+1) = loglog(NaN, NaN, 'r--');
    % xl = xlim();
    % l_f(length(l_f)+1) = loglog(xl, [c, c], 'r--');
    xlim(xl)
    ylim(yl)
    % title('B: Linear-log a_i(x-r_d - \tau_{0, i})^{\tau_i}')
    xlabel('Time + \tau_i (s)');
    ylabel('Tension - F_{ss} (kPa)')    
    set(gca, 'FontSize', 14);
end