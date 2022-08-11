%% set g0//
% clear
load g0;
g = g0;
% original Dan's g0
g = [0.7009    1.2243    1.0965    1.8390    0.4718 2.3357    0.3960    0.2372    0.1465    0.9817 0.8737    1.8333    0.2916    0.9513    1.0085]

% optimized g0
% negative kstiff1
% g = [0.8713    3.5481    1.0423   -1.3882    0.0349 1.3336    0.2065    0.6945    0.1177    2.1158    1.2210    2.0750   -0.9861    1.1191    3.2453]
% negative kstiff1 compensated
g = [0.8713    3.5481    1.0423   1.0    0.0349 1.3336    0.2065    0.6945    0.1177    2.1158    1.2210    2.0750   0.1    1.1191    3.2453]
% % rerun
% g = [1.0695    3.4444    1.1394    0.8159    0.0215 1.0595    0.1061    1.3011    0.0920    2.6525 1.1631    2.2622    0.0997    1.0826    3.7106];
% g0 = g;
% save g0 g0;

% dPUdT2 
g = [1.6841, 1.1324, 1.5813, 1.0150, 1.4815, 4.1527, 0.3237, 1.3563, 0.9649, 1.5560, 0.7972, 1.0657, 2.0000, 0.5240, 0.9312, 0.9737, 1.3799];
% dPUdt2, option 2;
g = [1.8682, 1.1595, 2.1354, 1.0000, 2.0445, 2.8134, 0.4285, 1.6496, 0.9649, 1.5970, 0.7838, 1.0657, 2.0000, 0.5132, 0.9800, 0.8845, 1.4225];
% dPUdt3
g = [3.8413, 1.4807, 3.3446, 1.1211, 3.1461, 0.91776, 0.51122, 1.3703, 0.76326,  2.741, 0.85679, 1.3433, 2.3047, 0.51121, 1.9806, 0.75292, 1.0673];
%% set up the problem
fcn = @dPUdT;
% fcn = @dPUdT_D;
[Etot, E1] = evaluateProblem(fcn, g, true, [1 0 0 0])

[Etot, E1] = evaluateProblem(fcn, g, true, [0 0 1 0])

%% test different setups

g = [1.6841, 1.1324, 1.5813, 1.0150, 1.4815, 4.1527, 0.3237, 1.3563, 0.9649, 1.5560, 0.7972, 1.0657, 2.0000, 0.5240, 0.9312, 0.9737, 1.3799];
[~, E1] = evaluateProblem(@dPUdT2, g, true)
g = [1.8682, 1.1595, 2.1354, 1.0000, 2.0445, 2.8134, 0.4285, 1.6496, 0.9649, 1.5970, 0.7838, 1.0657, 2.0000, 0.5132, 0.9800, 0.8845, 1.4225];
[~, E1_2] = evaluateProblem(@dPUdT2, g, true)
g = [3.8413, 1.4807, 3.3446, 1.1211, 3.1461, 0.91776, 0.51122, 1.3703, 0.76326,  2.741, 0.85679, 1.3433, 2.3047, 0.51121, 1.9806, 0.75292, 1.0673];
[~, E2] = evaluateProblem(@dPUdT3, g, true)

sum([E1;E1_2;E2]')

%% optim for each of plots separately, but together with km
g = [1.6841, 1.1324, 1.5813, 1.0150, 1.4815, 4.1527, 0.3237, 1.3563, 0.9649, 1.5560, 0.7972, 1.0657, 2.0000, 0.5240, 0.9312, 0.9737, 1.3799];
fcn = @dPUdT2;
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
x0101 = fminsearch(@(g)evaluateProblem(fcn, g, false, [0 1 0 1]), g, options);
x0011 = fminsearch(@(g)evaluateProblem(fcn, g, false, [0 0 1 1]), g, options);
x1001 = fminsearch(@(g)evaluateProblem(fcn, g, false, [1 0 0 1]), g, options);

evaluateProblem(@dPUdT2, x0101, true);
h1 = figure(101); h2 = figure(40101);copyobj(allchild(h1),h2);
evaluateProblem(@dPUdT2, x0011, true);
h1 = figure(101); h2 = figure(40011);copyobj(allchild(h1),h2);
evaluateProblem(@dPUdT2, x1001, true);
h1 = figure(101); h2 = figure(41001);copyobj(allchild(h1),h2);

g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma0", "kmsr", "alpha3", "k3", "K_T1", "s3", "kstiff1", "kstiff2", "K_T3", "16", "17"};
s1001 = calcSensitivities(fcn, x1001, g_names, true, false);title('Optim for 1001');
s0101 = calcSensitivities(fcn, x0101, g_names, true, false);title('Optim for 0101');
s0011 = calcSensitivities(fcn, x0011, g_names, true, false);title('Optim for 0011');

% 0101: not significant kd, k16, extreme significant k1, ksr, s3, kstiff2
% 1001: not sig: kd, kt3, k16
% 0011: not sig: kd, k1, k16

% only evaluated sensitivities
se1001 = calcSensitivities(fcn, x1001, g_names, true, false, [1 0 0 1]);title('Eval Optim for 1001');
se0101 = calcSensitivities(fcn, x0101, g_names, true, false, [0 1 0 1]);title('Eval Optim for 0101');
se0011 = calcSensitivities(fcn, x0011, g_names, true, false, [0 0 1 1]);title('Eval Optim for 0011');

evaluateProblem(@dPUdT2, x0001, true);
h1 = figure(101); h2 = figure(40001);copyobj(allchild(h1),h2);
evaluateProblem(@dPUdT2, x0010, true);
h1 = figure(101); h2 = figure(40010);copyobj(allchild(h1),h2);
evaluateProblem(@dPUdT2, x0100, true);
h1 = figure(101); h2 = figure(40100);copyobj(allchild(h1),h2);
evaluateProblem(@dPUdT2, x1000, true);
h1 = figure(101); h2 = figure(41000);copyobj(allchild(h1),h2);
%%
% single ident
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 800);
% g = [1.6841, 1.1324, 1.5813, 1.0150, 1.4815, 4.1527, 0.3237, 1.3563, 0.9649, 1.5560, 0.7972, 1.0657, 2.0000, 0.5240, 0.9312, 0.9737, 1.3799];
x0001 = fminsearch(@(g)evaluateProblem(fcn, g, false, [0 0 0 1]), g, options);
x0010 = fminsearch(@(g)evaluateProblem(fcn, g, false, [0 0 1 0]), g, options);
x0100 = fminsearch(@(g)evaluateProblem(fcn, g, false, [0 1 0 0]), g, options);
x1000 = fminsearch(@(g)evaluateProblem(fcn, g, false, [1 0 0 0]), g, options);

x0001 = [1.6822    1.1376    1.5914    1.0087    1.5093    4.1669    0.3229    1.3498    0.9787    1.5737    0.7892    1.0799    1.9628    0.5371    0.9223    0.9751    1.3855]
x0010 = [0.1274    1.2493    1.3222    1.3170    1.0162    4.6460    0.4426    1.9234    1.0820    1.2623    0.5493    0.7691    2.1239    0.5530    1.0370    1.0055    1.8911]
x0100 = [2.3993    0.9786    0.9433    0.5044    1.6071    1.8642    0.3983    0.9183    1.8131    1.9488    1.4771    0.6349    1.0298    0.9803    1.2062    0.7588    0.7980]
x1000 = [1.0422   -1.4125    3.7460    1.9151    1.7821    2.1777    0.3901    1.7396   -0.2467    3.1627    0.1900   -1.1967    2.5191    1.4679    1.5302    0.8636    2.2706]


% only evaluated sensitivities
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma0", "kmsr", "alpha3", "k3", "K_T1", "s3", "kstiff1", "kstiff2", "K_T3", "16", "17" ,"18"};

se0001 = calcSensitivities(fcn, x0001, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
se0010 = calcSensitivities(fcn, x0010, g, g_names, true, false, [0 0 1 0]);title('Eval Optim for 0010');
se0100 = calcSensitivities(fcn, x0100, g, g_names, true, false, [0 1 0 0]);title('Eval Optim for 0100');
se1000 = calcSensitivities(fcn, x1000, g, g_names, true, false, [1 0 0 0]);title('Eval Optim for 1000');

% save singleOptRes
%%
% only evaluated sensitivities
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma0", "kmsr", "alpha3", "k3", "K_T1", "s3", "kstiff1", "kstiff2", "K_T3", "16", "17" ,"18"};

se0001 = calcSensitivities(fcn, x0001, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
se0010 = calcSensitivities(fcn, x0010, g, g_names, true, false, [0 0 1 0]);title('Eval Optim for 0010');
se0100 = calcSensitivities(fcn, x0100, g, g_names, true, false, [0 1 0 0]);title('Eval Optim for 0100');
se1000 = calcSensitivities(fcn, x1000, g, g_names, true, false, [1 0 0 0]);title('Eval Optim for 1000');


%%
% Lets Run sensitivities for the result per partes
calcSensitivities(fcn, g, g, g_names, true, false, [1 1 1 1]);title('Eval Optim for 0001');

calcSensitivities(fcn, g, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
calcSensitivities(fcn, g, g, g_names, true, false, [0 0 1 0]);title('Eval Optim for 0010');
calcSensitivities(fcn, g, g, g_names, true, false, [0 1 0 0]);title('Eval Optim for 0100');
calcSensitivities(fcn, g, g, g_names, true, false, [1 0 0 0]);title('Eval Optim for 1000');


%% find a parameters that give a neagtive force
fcn = @dPUdT;
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 5000, 'OutputFcn', @myoutput);
x = fminsearch(@EvaluateNegativeForce, g0, options)

%% parameter search
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 5000, 'OutputFcn', @myoutput);
optimfun = @(g)evaluateProblem(fcn, g, false)
x = fminsearch(optimfun, g, options)
% E0 = evaluateProblem(fcn, x, true)
%% reduced g
g_selection = [1 3 5:15];
% leftovers = setdiff(1:length(g),g_selection);

gr = g(g_selection);
g_all = [gr(1) 3 gr(2) 0.8 gr(3:end)];

optimfun = @(gr)evaluateProblem(fcn, g_all, false);


x = fminsearch(optimfun, gr, options)


% x = fmincon(optimfun,g,[],[],[],[],ones(1, 15)*1e-3,[], [],options) 
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

%
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

function stop = myoutput(gr,optimvalues,state);
    stop = false;    
    if ~isequal(state,'iter') || mod(optimvalues.iteration, 20) > 0
        return;
    end
%     g_all = [gr(1) 3 gr(2) 0.8 gr(3:end)];
%     evaluateProblem(fcn, gr_all, true);
    evaluateProblem(@dPUdT, gr, true);
end
