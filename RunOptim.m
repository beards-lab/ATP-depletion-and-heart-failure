%% set g0
load g0;
g = g0;
% original Dan's g0
g = [0.7009    1.2243    1.0965    1.8390    0.4718 2.3357    0.3960    0.2372    0.1465    0.9817 0.8737    1.8333    0.2916    0.9513    1.0085]

% optimized g0
% g = [0.8713    3.5481    1.0423   -1.3882    0.0349 1.3336    0.2065    0.6945    0.1177    2.1158    1.2210    2.0750   -0.9861    1.1191    3.2453]
g = [0.8713    3.5481    1.0423   -1.3882    0.0349 1.3336    0.2065    0.6945    0.1177    2.1158    1.2210    2.0750   0.9861    1.1191    3.2453]

g0 = g;
save g0 g0;
%% set up the problem
fcn = @dPUdT;
E0 = evaluateProblem(fcn, g, true)
%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 1000, 'OutputFcn', @myoutput);
optimfun = @(g)evaluateProblem(fcn, g, false)

x = fminsearch(optimfun, g, options)

%% plot g
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma0", "kmsr", "alpha3", "k3", "K_T1", "s3", "kstiff1", "kstiff2", "K_T3"};
figure(2);clf;
subplot(211);cla;hold on;

bar(g);
plot([0, 16], [1 1], '--r')

title('Optimized values of G');
ylabel('g modifier value');
xticks(1:15)
xticklabels(g_names);
xtickangle(45)

%% sensitivity
delta = 0.1; % 10% difference
% clear E;
for k = 1:length(g)
    disp(num2str(k) + ": " + g_names{k} + '...')
    g_s = g;
    g_s(k) = g(k)*(1 - delta);
    E(k, 1) = evaluateProblem(fcn, g_s, false);
    E(k, 2) = E0;
    g_s(k) = g(k)*(1 + delta);
    E(k, 3) = evaluateProblem(fcn, g_s, false);
    % does this even matter? 
    % i.e. does using zero produces error lower than double?
    g_s(k) = 0;
    try
        E_g0(k) = evaluateProblem(fcn, g_s, false);
    catch
        E_g0(k) = NaN;
    end

end
disp("Done!")

%%
figure(2);subplot(212);cla;hold on;
title('Relative sensitivity to one-at-a-time perturbation of G by \delta');
bar(E/E0);
plot([0, 16], [1 1], '--r')

ylabel('E / E0');
ylim([0, 10])
xticks(1:15)
xticklabels(g_names);
xtickangle(45);

% is the param zeroable?
for k = 1:length(g)
    if ~isnan(g(k)) && E_g0(k) < 2*E0
        text(k, E_g0(k)/E0,'*', 'FontSize', 18);
    end
end
plot(nan, nan, '*k');
legend('g(x) - \delta', 'baseline', 'g(x) + \delta', 'baseline', 'g(x) = 0');
%%

function stop = myoutput(x,optimvalues,state);
    stop = false;    
    if ~isequal(state,'iter') || mod(optimvalues.iteration, 10) > 0
        return;
    end
    fcn = @dPUdT;
    evaluateProblem(fcn, x, true);

end
