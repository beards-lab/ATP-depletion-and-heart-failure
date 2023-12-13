x0 = ones(1, 10);

% % param modifiers
% 1  delU
% 2  kC
% 3  kS
% 4  alphaU
% 5  alphaF
% 6  nC
% 7  nS
% 8  nU
% 9 mu
% 10 Fss
% // 10 11 12 - params for ramp-up
% 13 Ls0
% 14 kA
% 15 kD
modNames = {'delU', 'kC', 'kS', 'alphaU', 'alphaF', 'nC', 'nS', 'nU', 'mu', 'b', 'c', 'd/Fps', 'Ls0', 'kA', 'kD', 'N/A'};
x0 = [ 1.0057    1.5981    1.3605    0.0917    1.4696    1.1553    1.0838    0.6758    1.1417]; % first shot

x0 = [ 0.9973    2.9945    2.8851    0.1308    -0.8659    1.3337    1.2446    0.7553     0.9829]; % negative param

% now, fixing the offset to 2.0 at 0.4 offset
x0 = [    1.0394    1.1020    0.8118    0.9860   0.7941    0.9136    0.9657    1.0193    1.4571]; % not really wokring

% optimizing with variable offset, not evaluating ramp-up, cost = 300.7
x0 = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345    1.2716];

% optimized for all params, excl. Fss
x0 = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345   0.6833    1.4778    0.8111 1 1];

% optimizing everything at once as a last step
x0 = [0.9736, 1.1519, 1.0295, 0.9306, 0.4439, 1.0758, 1.0011, 0.8704, 1.0829, 0.6774, 1.4252, 0.8157, 1.0302, 1.0118];

x0 = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];

% optim for -log10 weighting
mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1];
%% set data set
% pCa = Inf; % First round of passive ramp-ups data, 200s long, dML 0.4
% pCa = 11; % Ca&PNB experiments, resting ramp-ups, 30s decay, dML = 0.225

mod = [0.9736, 1.1519, 1.0295, 0.9306, 0.4439, 1.0758, 1.0011, 0.8704, 1.0829, 0.6774, 1.4252, 0.8157, 1.0302, 1.0118 1 1 1 1];
mod = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];


% test for pCa 11
mod(1) = 2.2;
mod(2) = 4.4;
mod(3) = 2.2;
% mod(5) = 2;
mod(6) = 1.03;
% mod(8) = 1.2;
mod(13) = -2.8;
%%
% mod = [2.0398    1.3113    3.8942    1.3500    0.4784 0.7398    0.8176    0.7869    0.8629    0.7200 1.3634    1   -2.8000    1.0150    0.6382 -0.5199];
% mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988, 1, 1, 1];
% mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
% mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];

tic
saSet = 1:15;
[c0_11, cost11_sap, cost11_sam] = runSa(11, saSet);
[c0_4, cost4_sap, cost4_sam] = runSa(4.4, mod, saSet);
% normalize
cost11_sam = cost11_sam/c0_11;
cost11_sap = cost11_sap/c0_11;
cost4_sam = cost4_sam/c0_4;
cost4_sap = cost4_sap/c0_4;
toc
%% plot the result
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0'};

c0_11 = 1;
c0_4 = 1;

% identify tradeoffs
eps = (c0_4 + c0_11)*1e-3;
% better result
bett = (cost11_sap + cost4_sap + eps < c0_11 + c0_4) ...
    | (cost11_sam + cost4_sam + eps < c0_11 + c0_4);
betts = strings(1, length(cost4_sap));betts(bett) = "$$$";
% tradeoffs between low and high Ca - suggesting some dependency?
tdoCa = (cost11_sap + eps < c0_11 & cost4_sap > c0_4 + eps) | (cost11_sam + eps < c0_11 & cost4_sam > c0_4 + eps) ...
    | (cost11_sap > c0_11 + eps & cost4_sap + eps < c0_4) | (cost11_sam > c0_11 + eps & cost4_sam + eps < c0_4);
tdoCas = strings(1, length(cost4_sap));tdoCas(tdoCa) = "*";

figure(101);clf; 
b = bar([cost11_sap; cost4_sap; zeros(size(cost4_sap)); cost11_sam;cost4_sam]'); hold on;
title('Sensitivity Anal')
plot([0.5, length(cost4_sap)], [c0_11 c0_11], 'b--')
plot([0.5, length(cost4_sap)], [c0_4 c0_4], 'm--')
title('Grouped sensitivity analysis (* indicates Ca trade-off, $$$ potential of improvement)')
legend('10xpCa11+','pCa4+','', '10xpCa11-','pCa4-','pCa11_0','pCa4_0')
set(gca,'XTick',saSet);
set(gca,'XTickLabels',strcat(string(1:length(cost4_sap)), ':', modNames(1:length(cost4_sap)),tdoCas, betts));
set(gca, "FontSize", 14)
%%
% close all
tic
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0'};

% optim for -log10 weighting
% mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
% reoptim on avg data
% mod = [0.0165    0.7035    0.4403    1.0234    1.0077    0.5754    0.9379    1.1950    0.9099    0.8988    1.0000    1.1421    1.4792    1.1156    2.9834];
% reoptim on avg data incl. 20231102
% mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
% modSel = [2 3 4 6 7 8 9 12 14 15];
% mod(16) = 0;
% evalCombined(mod(modSel))
% mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% modSel = [7, 9];
% mod_pca6 = [0.1657, 0.3895];
% mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];
% for absolute average - relaxed
% mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    0.9379 1.1950    0.9099    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
% for absolute average - pca 4.4
% mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    1.5266 0.8291    0.0356    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
% for absolute average excluding the ramp-up
% mod = [0.0142    0.4204    0.4267    1.0234    0.7567 0.3025    0.1535    1.0318    0.0133    0.8988 0.8679    4.5429    0.7510    1.2811    1.6208];
%% optim for filtering out remaining force (attempt 01)
mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    1.2811    1.6663];
mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    5.2811    1.6663];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% optim for pCa 11 only, from original estimates, without ramp up
% modSel = [1 2 3 5 6 10 15];
% mod = 0.0278    0.7953    0.4768    1.0000    1.0770    1.7512    1.0000    1.0000    1.0000    1.3607    1.0000    1.0000    1.0000    1.0000   -0.1149
% fixing the refolding to non-negative
% mod =[0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    1.0000    0.50000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% retuned for reduced mu
% mod =  [0.0287    0.7984    0.4766    1.0000    1.0715    1.8054    1.0000    1.0000    1.0000    1.3254    1.0000    1.0000    1.0000    1.0000    0.5000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% Retuned for ramp-ups, strange bump in 10s
% mod =  [0.0293    0.8574    0.4871    1.0000    0.9528    1.7817    1.0000    1.0000    1.0000    1.3369    3.6372    0.2425    0.0030    0.1000    0.5000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 11 and 4.4, from scratch except for ramp-up and fixed mu and alphaR0 
% Candidate 1
tic
modSel = [1 2 3 4 5 6 7 8 9 10]; mod = [0.4852    0.2070    1.0403    1.1617    0.7393    1.2668    1.3151    1.5592    1.1949    1.6778    3.6372    0.2425    0.0030    0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, [mod 1 1], 1:15)
toc
%% yet another set
mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    1.0000    0.50000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for reduced mu
modSel = [1 2 3 4 5 6 7 8 9 10]; mod =  [0.0287    0.7984    0.4766    1.0000    1.0715    1.8054    1.0000    1.0000    1.0000    1.3254    1.0000    1.0000    1.0000    1.0000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for both pCas in log space
% candidate 2
mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 4 in log space only, incl. ramp-up
mod = [0.4804    0.0661    0.9735    2.0109    0.5917 1.6109e+05 0.1694    0.2900    0.7105    1.4348    0.3534    0.1366    0.0002 0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 4 in log space only, w/o. the ramp-up
mod = [0.4804    0.0667    0.9711    1.9971    0.5672 2.0654e+05 0.0919    0.3252    0.6417    2.0434    0.3652    0.0866    0.0002 0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, [mod 1 1], 1:15);
%% Tuned for data readjusted for remaining force and scaled by Fmax - with 1.4e4
mod = [0.6108    0.0578    1.0426    1.1110    0.5432    1.3586    0.9776    1.3503 1.2387    1.8919    3.6372    0.2425    0.0030    0.7421    0.4327    0.9964 1.0091];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% FineTuning for b, c and d
mod = [0.6108    0.0578    1.0426    1.1110    0.5432    1.3586    0.9776    1.3503    1.2387    1.8919    0.001    1.0877    0.001    0.7421    0.4327    0.9964    1.0091];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% Retuned for Fmax normalized data without Frem adjustment, incl. mu
mod = [0.4847    0.0514    1.0389    1.1202    0.5064    1.3293    0.9757    1.3062    1.3521    1.9492    0.0003    1.0841    0.0007    0.7551    0.6101    0.9964    1.0090];
% [0.4777    0.0503    1.0404    1.1199    0.5045    1.3213     0.9756    1.2949    1.3656    1.9521    0.0003    1.0843    0.0008    0.7595    0.6170    0.9964    1.0089];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% unfinished finetuning with free PEVK_A_Relaxed, costs 3e4
mod = [0.2651    0.0171    1.0557    1.1643    0.3515    1.4502    0.9884    1.1832    1.0980    2.0724    0.0000    0.8893    0.0027    0.6528    0.6272    0.9964    1.0068];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%%
modSel = 1:15;
i = 14;
mod(i) = mod(i)*1;
disp(modNames(i))
%%
tic
% modSel = [1 2 3 5 6 10];
% default is 403
% modSel = [7, 9];
mod = [mod ones(1, 17 - length(mod))]
modSel = 1:length(mod)
evalCombined(mod(modSel), mod, modSel)
% evalCombined(mod_pca6)
% evalCombined(mod)
% evalCombined([1 1 1 1 1 1])
toc
%%
clf;
modNames{16} = 'k_{PEVK,A} (low Ca)';
modNames{17} = 'k_{PEVK,D} (low Ca)';
mm = min(floor(log10(mod)));

n = length(mod);
polarplot(linspace(2*pi/n, 2*pi, n), log10(mod) - mm, 'x-', LineWidth=2);
% evalCombined(mod, mod, 1:n);
% set(gca, 'ThetaAxisUnits', 'radians', 'thetatick', linspace(2*pi/n, 2*pi, n));
set(gca, 'ThetaAxisUnits', 'radians', ...
    'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
    'Rlim', [-3 1]-mm, 'RTick', [-2 -1 0 1]-mm, 'RTickLabel', [0.01 0.1 1 10]);
hold on;

%%
% mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    1.2811    1.6663];
% mod = ones(1, 15);
% modSel = [1 2 3 5 6 11 12 15];
% modSel = [7, 8, 9];
% modSel = [1 2 3 5 6 7 8 9 11 12 15];
% 
% modSel = [1 2 3 5 6 10 11 12 13 15];
% modSel = [1 2 3 5 6 8 10 11 12 13 15];
% mod = [    0.0140    0.3630    0.4241    1.0234    0.7953    0.2941    0.2514     0.8819    0.0135    0.8988    0.7571    4.3872    0.7510    1.2811    1.6652];
% cost = 301.7
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 5000);
% x0 = x0([1:4 6:10 13]);
% x0 = [0.5795    1.6029    0.9047    1.0569    0.9236    0.6627    1.0497    0.9465    1.0641    0.7124];

% best in the worst so far optMods =     1.3705    1.4770    1.1905    1.1205    0.8670    0.3074    0.6812    0.9407    1.3393    0.8915    1.1628    1.0166
% modSel = [2 3 4 6 7 8 9 12 14 15];

% init = mod(modSel);
% init = ones(1, 10);
% modSel = [7, 9];
% mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    0.1    0.50000];
% reducing the mu to nonsensitive value, 
% modSel = [1 2 3 5 6 10]; mod =  [0.0293    0.8574    0.4871    1.0000    0.9528    1.7817    1.0000    1.0000    1.0000    1.3369    1.0000    1.0000    1.0000    0.1000    0.5000];
% optimizing for the ramp-ups in logspace
% modSel = [1 2 3 4 5 6 7 8 9 10];mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];
% optimizing for high pCa
modSel = [2 3 4 5 6 7 8 9 10 11 12 13]; mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];


%%
modSel = [1:10 14:17];
init = mod(modSel);
% evalFunc = @(optMods) evalCombined(optMods, mod, modSel);
x = fminsearch(@evalCombined, init, options, mod, modSel);

mod(modSel) = x;
% mod = x;

save mod;
%% optim in log param space
% mod = ones(1, 17);
mod = [mod ones(1, 17 - length(mod))];
modSel = [1:10 14:17];
modSel = [11 12 13];
modSel = 1:17;
init = log10(mod(modSel));
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
x = fminsearch(evalLogCombined, init, options);
mod(modSel) = 10.^x;

%% optim with bounded params
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

modSel = [11 12 13];
init = log10(mod(modSel));
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
[x, FV, EF, OUT] = fmincon(evalLogCombined, init, [], [], [], [], log10(lb(modSel)), log(ub(modSel)), [], options);
%% test in GA
% parpool
ga_Opts = optimoptions('ga', ...
    'PopulationSize',64, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',8, ...  % 10
    'UseParallel',true);
modSel = 1:14;
% Ng = length(modSel);
ub = 100*ones(1, length(mod));
ub([3, 4, 5, 10, 14]) = [10 10 10, 10, 10];
lb = 0.001*ones(1, length(mod));
lb([3, 4, 5, 10, 14]) = [.1 .1 .1 .1 .1];
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
init = log10(mod(modSel));

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(evalLogCombined,length(modSel), ...
    [],[],[],[],...
    log10(lb(modSel)),log10(ub(modSel)),[],ga_Opts);

mod(modSel) = 10.^p_OptimGA;
% use fminserach afterwards
%% result
mod = [2.0398    0.9359    4.3424    2.3068    0.4784 0.6410    0.8054    0.8220    0.2923    0.7200 1.3634    1.0000   -2.8000    1.3022    0.8551];
mod = [2.0398    0.9359    4.3397    2.2756    0.4784 0.6412    0.8053    0.8333    0.2923    0.7200 1.3634    1.0468   -2.8000    1.3200    0.8759];

% proximal chain Ca dependency
mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% Fix for nS
mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];

% mods for pCa 6: modSel = [7, 9], x = [0.1657, 0.3895];
mod([7,9]) = [0.1657, 0.3895];
% thus, the profile is:
pCa = [11, 6, 4];
kC   = [10203,10203*4.78*0.3895,10203*4.78*0.9889];      % proximal chain force constant
kA   = [0, 16.44*0.1657,16.44*0.9403];
figure(44);clf;
% subplot(211);
plot(pCa, kA, 'x:', LineWidth=2);
ylabel('kA Value');
yyaxis right;
plot(pCa, kC, '+:', LineWidth=2);
ylabel('kC Value');

set ( gca, 'xdir', 'reverse' );
xlabel('pCa');
legend('kA (PEVK attachment)', 'kC (stiffening proximal chain)')
%%
function totalCost = evalCombined(optMods, mod, modSel)

    
    %normal - optimizing for all
    % modSel = 1:15;

    % modSel = [1 2 3 5 6 10];
    % optimizing only subset of mods
    % mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];
    % mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
    % mod = [0.0165    0.7035    0.4403    1.0234    1.0077    0.5754    0.9379    1.1950    0.9099    0.8988    1.0000    1.1421    1.4792    1.1156    2.9834];
    % mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
    % for absolute average
    % mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    0.9379 1.1950    0.9099    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
    % optim for -log10 weighting
    % mod = [0.928351480405723	0.928351480405723	1.01367539052550	1.02336567490158	1.01367539052550	1.10319213342611	0.937882365838957	1.19500150970587	0.909890571615859 1 1 1 1];
    % reset the search
    % mod = ones(1, 15);
    % pCa 11 decay with reduced mu
    % mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    0.1    0.50000];
    % modSel = [1 2 3 5 6 10];

    % % mod([1:4 6:10 13]) = optMods;
    % modSel = [11, 12, 13];
    % modSel = [1 2 3 5 6 11 12 15];
    % modSel = [7, 8, 9];
    % modSel = [1 2 3 5 6 10];
    
    % store the init
    baseMods = mod;
    mod(modSel) = optMods;

    drawPlots = true;
    % drawPlots = false;
    totalCost = 0;

    if drawPlots && 1
        figInd = 100;
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = 'Parameter modifiers';
        end      
        clf;
        plotParams(baseMods);hold on;
        plotParams(mod);hold off;
    end

    %% pCa 6
    % pCa = 6;
    % mod([7,9]) = [0.1657, 0.3895];
    % if drawPlots
    %     f = figure(35);f.Name = 'pCa 6';
    % end
    % 
    % cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    % totalCost = totalCost + cost;
    % % return;


    %% no Ca
    pCa = 10;
    figInd = 110;
    if drawPlots
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = 'pCa 11';
        end        
    end    
    % RunCombinedModel;
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    totalCost = totalCost + cost*10;
% return
    %% pCa 4
    pCa = 4.4;
    figInd = 104;
    % mod([7,9]) = [0.9403    0.9889];
    if drawPlots
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = 'pCa 4.4';
        end
    end
    % RunCombinedModel;
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    totalCost = totalCost + cost;

end

function cost = isolateRunCombinedModel(mod, pCa, drawPlots)
% just to isolate the script, so the variables can't intervene
    % drawPlots = true;
    RunCombinedModel;
end

function plotParams(mod, resetGca)

    modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)'};
    mm = min(floor(log10(mod))); n = length(mod);
    
    polarplot(linspace(2*pi/n, 2*pi, n), log10(mod) - mm, 'x-', LineWidth=2);
    
    if strcmp('degrees', get(gca, 'ThetaAxisUnits')) || nargin <= 2 || ~resetGca
        % degrees means it is not adjusted yet 
        % or the param has not been provided 
        % or we want to reset
        return
    end
    set(gca, 'ThetaAxisUnits', 'radians', ...
        'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
        'Rlim', [-3 1]-mm, 'RTick', [-2 -1 0 1]-mm, 'RTickLabel', [1 10 100 1000]*10^mm);
end

function [c0, cost_sap, cost_sam] = runSa(pCa, mod, saSet)
    drawPlots = false;
    cost_sap = []; % SA plus
    cost_sam = []; % SA minus

    % pCa = 4;
    if drawPlots
        figure(100);
    end
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    c0 = cost;
    %
    % saSet = [1:9 12:15];
    % saSet = [1:8 13];
    % saSet = [14 15];
    
    SAFact = 1.05;
    for i_m = saSet
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('Mod %g is up to %g..', i_m, mod(i_m));
        % figure(i_m)
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sap(i_m) = cost;
        
        mod(i_m) = mod(i_m)/SAFact/SAFact;
        fprintf('costing %1.4e€ and down to %g...', cost, mod(i_m));
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sam(i_m) = cost;
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('costing %1.4e€. \n', cost);
    end
end