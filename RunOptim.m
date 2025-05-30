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
% update the length
g = [g 1 1]
% 19 pcs from 26-08 incl the mu viscosity optim
g = [0.4156    0.0600    5.3243    2.7286 0.5186    3.9465    0.5116   10.6276 1.3384    0.3296    0.6091    1.1200 1.7357    1.5442    0.8284    1.4614 2.3482    0.6206    1.2782]
% 10/08
g = [1.6469, 1.0000, 1.1611, 0.0053, 1.5031, 1.5605, 0.9276, 0.0220, 0.3455, 1.0114, 0.2624, 1.9307, 1.5611, 1.2518, 1.6784, 0.6298, 1.1186, 0.4442, 1.3752, 1.0959]
g = load('gopt.csv');
%% Eval the model at vel 0
g = ones(1, 22);
% params0 = struct('Pi', 0, 'MgATP', 8, 'MgADP', 0, 'Ca', 1000,'Velocity', 0,'UseCa', false,'UseOverlap', false);
% opts0 = struct('N', 20, 'Slim', 0.04, 'PlotProbsOnFig', 0, 'ValuesInTime', true, ...
%     'SL0', 2.2*1.1, ...
%     'MatchTimeSegments', 1, 'ML', 2.2, 'PlotProbsOnStep', false, 'ReduceSpace', false,...
%     'OutputAtSL', Inf, 'LSE0', 0.0);
% A = [8 4 2];
% params0.Velocity = -6;
g(19) = 10;g(20) = 10;
params0 = getParams(struct('SL0', 1.1*2.0), g);
t_change = 1e-5;
changePercent = -[8, 10, 12, 14, 16];
changeSL = params0.SL0*(changePercent/100);
velocity = changeSL/t_change;


params0.ValuesInTime = true;
params0.UseSlack = true;
params0.UseSerialStiffness = true;
[F, out] = evaluateModel(@dPUdTCa, [0 3], params0);
% out.Force(end)
PU0 = out.PU(end, :);
params0.PU0 = PU0;
%
k = 1;
figure(10);clf;
subplot(211);hold on;
subplot(212);hold on;
for k = 1:length(velocity)
%     params0.MgATP = A(k);
    params0.Velocity = [velocity(k), 0, 0]/params0.ML;
    [F(k), out] = evaluateModel(@dPUdTCa, [0 t_change, 0.5], params0);
    subplot(211);plot(out.t, out.SL, out.t, out.LXB, ':', 'Linewidth', 2);
    subplot(212);plot(out.t, out.Force, ':', 'Linewidth', 2);
end


%% set up the problem
% original fit from tuesday
% g = [0.542876541225473,0.914586975331451,0.147826886651247,7.16847703087333,2.59787534877577,0.0322817577373749,0.794984662510472,0.00676842719352362,0.050540767958257,0.617238554706245,1,2.08884046119593,2.0054310807672,0.921193092113734,1,1.00122403975111e-06,2.27995446259569,0.864010537605939,1,1.0463623046875,0.347826086956522,1.8];
% refit for SL length based sim
g = load('gopt.csv');
% Test the increased dr
% g = [0.5429    0.9146    0.1478    7.1685    2.5979    0.0323    0.7950    0.0068 0.0505    0.6172    1.0000    6.0000    2.0054    0.3067    1.0000    0.0000    2.2800    0.8640    1.0000    1.0464    0.3478    1.8000]

fcn = @dPUdTCa;
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma_0", "kmsr", "\alpha_3", "k3", "K_{T1}", "s3", "k_{stiff1}", "k_{stiff2}", "K_{T3}", "\alpha_1", "\alpha_2" ,"A_{max0}", "\mu_{v0}", "\mu_h0", "k_pas"};
% g = ones(21, 1);
% fcn = @dPUdT_D;
tic
g(20) = 0.5;
g(19) = 500;

% [Etot, E1] = evaluateProblem(fcn, g, true, [0 0 0 0 0 1 1 1 1])
% [Etot, E1] = evaluateProblem(fcn, g, true, [1 0 0 0 0 1 1 1 1])
[Etot, E1] = evaluateProblem(fcn, g, false, [0 0 0 0 0 0 0 0 0 1])
toc
E1
% writematrix(g, 'gopt.csv')
%%
params0 = getParams([], g);
params0.ValuesInTime  =true;
params0.SL0 = 1.4;
params0.Velocity = [0 0.01];
params0 = getParams(params0, g); % need to update the init which is based on the SL0 !!!
[~, out] = evaluateModel(fcn, [0 1 100], params0);

figure(88);clf; hold on;
plot(out.t, out.SL, out.t, out.FXB, out.t, out.OV)
figure(89);clf; hold on;
plot(outBase.OV, outBase.FXB,out.OV, out.FXB)

%% Run through params to eval their importance
close all;
evaluateProblem(fcn, g, true, [1 1 1 1 0 1])
set(gcf, 'Name', 'Baseline' , 'NumberTitle', 'off')

for i = 1:length(g)
    disp(i)
    gt = g(i);
    g(i) = g(i)*2;
    evaluateProblem(fcn, g, true, [1 1 1 1 0 1])
    g(i) = gt;
    set(gcf, 'Name', ['Param ' num2str(i) ' to 150%'], 'NumberTitle', 'off')
end

% [Etot, E1] = evaluateProblem(fcn, g, true, [0 0 0 0 1])
%%
% [Etot, E1] = evaluateProblem(fcn, g, true, [0 0 1 0])

% Run basic sensitivity with updated cost func
% only evaluated sensitivities
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma0", "kmsr", "alpha3", "k3", "K_T1", "s3", "kstiff1", "kstiff2", "K_T3", "16", "17" ,"18", "19", "u", "kpas"};
se1111 = calcSensitivities(fcn, g, g0, g_names, true, false, [0 0 0 0 0 1]);title('Eval Optim for 1111');
saveas(gcf, 'sensitivities.png')
%% extension
se0001 = calcSensitivities(fcn, g, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
se0010 = calcSensitivities(fcn, g, g, g_names, true, false, [0 0 1 0]);title('Eval Optim for 0010');
se0100 = calcSensitivities(fcn, g, g, g_names, true, false, [0 1 0 0]);title('Eval Optim for 0100');
se1000 = calcSensitivities(fcn, g, g, g_names, true, false, [1 0 0 0]);title('Eval Optim for 1000');


%%
fcn = @dPUdTExt;
se0001 = calcSensitivities(fcn, g, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
se1111 = calcSensitivities(fcn, g, g, g_names, true, false, [1 1 1 1]);title('Eval Optim for 1111');
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

se0001 = calcSensitivities(fcn, x0001, g, g_names, true, false, [0 0 0 1]);title('Eval Optim for 0001');
se0010 = calcSensitivities(fcn, x0010, g, g_names, true, false, [0 0 1 0]);title('Eval Optim for 0010');
se0100 = calcSensitivities(fcn, x0100, g, g_names, true, false, [0 1 0 0]);title('Eval Optim for 0100');
se1000 = calcSensitivities(fcn, x1000, g, g_names, true, false, [1 0 0 0]);title('Eval Optim for 1000');

% save singleOptRes
%%
% only evaluated sensitivities

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
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 1500, 'OutputFcn', @myoutput);
% tic
% estart = evaluateProblem(fcn, g, true, [0 0 0 1])
% toc

optimfun = @(g)evaluateProblem(fcn, g, false, [1 1 1 1 1 1]);
x = fminsearch(optimfun, g, options)
g = x;
save gopt1010 g;
% E0 = evaluateProblem(fcn, x, true)
%% Attempt on GA
% parpool
ga_Opts = optimoptions('ga', ...
    'PopulationSize',64, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',4, ...  % 10
    'UseParallel',true);

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(optimfun,size(gr0, 2), ...
    [],[],[],[],ones(size(gr0))*0.01,ones(size(gr0))*20,[],ga_Opts);

save env;

optimfun(p_OptimGA)
x = fminsearch(optimfun, p_OptimGA, options)
save x
optimfun(x)
%% reduced g
fcn = @dPUdTCa;
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'TolX', 1e-3, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% g = load('gopt.csv');
% g = ones(20, 1);
% g(21) = 80/230; % fix the passive force
% g([3,4 5]) = 10;
% , 'OutputFcn', @myoutput);
% g = load('gopt.csv')';

% g_selection = [1 3:10 12:14 16:19 21];

exclude = [2 3 4 11 15 17 18 19 20 21 22];
exclude = [2:5 11 12 15:22];
exclude = [11 15 21:22];
% exclude = [1  2   6     7     8     9    10    12    13    14    16];
g_selection = setdiff(1:length(g), exclude);
% leftovers = setdiff(1:length(g),g_selection);
% g_selection = [1 3 4 10 16]

gr0 = g(g_selection);
% g_all = [gr(1) 3 gr(2) 0.8 gr(3:end)];
% optimfun = @(gr)evaluateProblem(fcn, g_all, false);
% optimfun = @(gr)evaluateProblem(fcn, [gr(1) 1 gr(2:end)], false, [1 1 0 0 0 1]);
% optimfun = @(gr)evaluateProblem(fcn, insertAt(g, gr, g_selection), false, [0 0 1 0 0 1 1 1]);
optimfun = @(gr)evaluateProblem(fcn, insertAt(g, gr, g_selection), false, [0 0 0 0 0 0 0 0 0 1]);
% optimfun = @(gr)evaluateProblem(fcn, gr, false, [1 0 0 0 0 1]);
% tic
% optimfun(gr0)
% toc
%% optimfun(g)
% x = fminsearch(optimfun, gr0, options)
% x = fminsearch(optimfun, p_OptimGA, options)
x = fminsearch(optimfun, gr0, options)
save x


% x = fmincon(optimfun,g,[],[],[],[],ones(1, 15)*1e-3,[], [],options) 
g = insertAt(g, x, g_selection)
% save gopt1_eval6 g;

% to commit as plaintext
writematrix(g, 'gopt.csv_slack')
%%
tic
% [Etot, E1] = evaluateProblem(fcn, g, true, [1 0 1 1 0 0 0 0])
[Etot, E1] = evaluateProblem(fcn, g, true, [1 0 0 0 0 0 0 0 1 1])
toc
saveas(gcf, 'ProblemEval.png')

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

