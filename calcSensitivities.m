function E = calcSensitivities(fcn, g, g0, g_names, drawPlots, tryZeros, evalParts, sens)
E0 = evaluateProblem(fcn, g, false, evalParts);
if nargin < 8
delta = 0.1; % 10% difference
% clear E;
for k = 1:length(g)
    disp(num2str(k) + ": " + g_names{k} + '...')
    g_s = g;
    g_s(k) = g(k)*(1 - delta);
    E(k, 1) = evaluateProblem(fcn, g_s, false, evalParts);
    g_s(k) = g0(k)
    E(k, 2) = evaluateProblem(fcn, g_s, false, evalParts);;
    g_s(k) = g(k)*(1 + delta);
    E(k, 3) = evaluateProblem(fcn, g_s, false, evalParts);
    % does this even matter? 
    % i.e. does using zero produces error lower than double?
    if tryZeros
        g_s(k) = 0;
        try
            E_g0(k) = evaluateProblem(fcn, g_s, false);
        catch
            E_g0(k) = NaN;
        end
    end

end
disp("Sensitivities Done!")
else
    E = sens;
    disp('sensitivities provided, plotting..')
end

if ~drawPlots
    return;
end
%%
figure();
subplot(211);cla;hold on;

bar([g0;g]');
plot([0, 16], [1 1], '--r')
legend('Original', 'Single-optimized')

title(['Optimized values of G for '  num2str(evalParts)]);
ylabel('g modifier value');
xticks(1:15)
xticklabels(g_names);
xtickangle(45);


subplot(212);cla;hold on;
title('Relative sensitivity to one-at-a-time perturbation of G by \delta');
bar(E/E0);
plot([0, 16], [1 1], '--r')

ylabel('E / E0');
ylim([0, 10])
xticks(1:15)
xticklabels(g_names);
xtickangle(45);

% is the param zeroable?
if tryZeros
for k = 1:length(g)
    if ~isnan(g(k)) && E_g0(k) < 2*E0
        text(k, E_g0(k)/E0,'*', 'FontSize', 18);
    end
end
end
plot(nan, nan, '*k');
legend('g(x) - \delta', 'baseline', 'g(x) + \delta', 'baseline', 'g(x) = 0');
% ylim([0.95 1.05])