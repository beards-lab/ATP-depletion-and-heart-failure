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
drawPlots = true;
cost_sap = []; % SA plus
cost_sam = []; % SA minus
pCa = 11;
figure(100);
RunCombinedModel;
cost
c0 = cost;
%
saSet = [1:9 12:15];
% saSet = [1:8 13];
% saSet = [14 15];
for i_m = saSet
    mod(i_m) = mod(i_m)*1.1;
    fprintf('Mod %g is up to %g..', i_m, mod(i_m));
    % figure(i_m)
    RunCombinedModel;    
    cost_sap(i_m) = cost;
    
    mod(i_m) = mod(i_m)/1.1/1.1;
    fprintf('costing %1.4e€ and down to %g...', cost, mod(i_m));
    RunCombinedModel;    
    cost_sam(i_m) = cost;
    mod(i_m) = mod(i_m)*1.1;
    fprintf('costing %1.4e€. \n', cost);
end
%% plot the result
% %%
% cost11_sap = cost_sap;
% cost11_sam = cost_sam;
% c0_11 = c0;
% %%
% cost4_sap = cost_sap;
% cost4_sam = cost_sam;
% c0_4 = c0;

% identify tradeoffs
eps = (c0_4 + c0_11)*1e-3;
% better result
bett = (cost11_sap*10 + cost4_sap + eps < c0_11*10 + c0_4) ...
    | (cost11_sam*10 + cost4_sam + eps < c0_11*10 + c0_4);
betts = strings(1, length(cost4_sap));betts(bett) = "+++";
% tradeoffs
tdo = (cost11_sap + eps < c0_11 & cost4_sap > c0_4 + eps) | (cost11_sam + eps < c0_11 & cost4_sam > c0_4 + eps) ...
    | (cost11_sap > c0_11 + eps & cost4_sap + eps < c0_4) | (cost11_sam > c0_11 + eps & cost4_sam + eps < c0_4);
tdos = strings(1, length(cost4_sap));tdos(tdo) = "*";
figure(101);clf; 
b = bar([cost11_sap*10; cost4_sap; zeros(size(cost4_sap)); cost11_sam*10;cost4_sam]'); hold on;
title('Sensitivity Anal')
plot([0, 16], [c0_11 c0_11]*10, 'b--')
plot([0, 16], [c0_4 c0_4], 'm--')
title('Grouped sensitivity analysis (* indicates trade-off, +++ potential of improvement)')
legend('10xpCa11+','pCa4+','', '10xpCa11-','pCa4-','pCa11_0','pCa4_0')
set(gca,'XTick',saSet);
set(gca,'XTickLabels',strcat(string(1:length(cost4_sap)), ':', modNames(1:15),tdos, betts   ));
set(gca, "FontSize", 14)
%%
tic
% modSel = [2 3 4 6 7 8 9 12 14 15];
% mod(16) = 0;
% evalCombined(mod(modSel))
evalCombined(mod)
toc
%%
% cost = 301.7
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% x0 = x0([1:4 6:10 13]);
% x0 = [0.5795    1.6029    0.9047    1.0569    0.9236    0.6627    1.0497    0.9465    1.0641    0.7124];

% best in the worst so far optMods =     1.3705    1.4770    1.1905    1.1205    0.8670    0.3074    0.6812    0.9407    1.3393    0.8915    1.1628    1.0166
modSel = [2 3 4 6 7 8 9 12 14 15];

% init = mod(modSel);
init = ones(1, 10);
x = fminsearch(@evalCombined, init, options);

% mod(modSel) = x;
mod = x;

save mod;
%% result
mod = [2.0398    0.9359    4.3424    2.3068    0.4784 0.6410    0.8054    0.8220    0.2923    0.7200 1.3634    1.0000   -2.8000    1.3022    0.8551];
mod = [2.0398    0.9359    4.3397    2.2756    0.4784 0.6412    0.8053    0.8333    0.2923    0.7200 1.3634    1.0468   -2.8000    1.3200    0.8759];

% proximal chain Ca dependency
mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
%%
function totalCost = evalCombined(optMods)
    %normal - optimizing for all
    % mod = optMods;
    
    % optimizing only subset of mods
    % mod = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345    1.2716];
    % mod([11, 12, 13]) = optMods(:);
    % x0 = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];    
    % mod = [2.2000    4.4000    2.2000    1.0226    0.4784    1.0300    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411   -2.8000    1.0150    1.0000    1.0000];

    % mod = [2.0398    0.9359    4.3424    2.3068    0.4784 0.6410    0.8054    0.8220    0.2923    0.7200 1.3634    1.0000   -2.8000    1.3022    0.8551];
    % 
    % modSel = [2 3 4 6 7 8 9 12 14 15];
    % % mod([1:4 6:10 13]) = optMods;
    % mod(modSel) = optMods;

    mod = optMods;

    drawPlots = true;
    %% no Ca
    pCa = 11;
    if drawPlots
        figure(33);title('pCa 11');
    end    
    % RunCombinedModel;
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    totalCost = cost*10;

    %% pCa 4
    pCa = 4;
    if drawPlots
        figure(34);title('pCa 4');
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
