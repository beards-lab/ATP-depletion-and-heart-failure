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

%%

mod = [0.9736, 1.1519, 1.0295, 0.9306, 0.4439, 1.0758, 1.0011, 0.8704, 1.0829, 0.6774, 1.4252, 0.8157, 1.0302, 1.0118 1 1 1 1];
mod = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];

pCa = Inf;
figure(3)
RunCombinedModel;

%%
tic
modSel = [1:4 6:9 11 13 15 16]
evalCombined(modSel)
toc
%%
% cost = 301.7
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);
% x0 = x0([1:4 6:10 13]);
% x0 = [0.5795    1.6029    0.9047    1.0569    0.9236    0.6627    1.0497    0.9465    1.0641    0.7124];

% best in the worst so far optMods =     1.3705    1.4770    1.1905    1.1205    0.8670    0.3074    0.6812    0.9407    1.3393    0.8915    1.1628    1.0166
init = x0(modSel);
x = fminsearch(@evalCombined, init, options);

save x;
%%
function cost = evalCombined(optMods)
    %normal - optimizing for all
    % mod = optMods;
    
    % optimizing only subset of mods
    % mod = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345    1.2716];
    % mod([11, 12, 13]) = optMods(:);
    x0 = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];    
    mod = x0;
    % mod([1:4 6:10 13]) = optMods;
    mod([1:4 6:9 11 13 15 16]) = optMods;

    drawPlots = false;

    pCa = 11;
    
    if drawPlots
        figure(33);title('pCa 11');
    end
    
    RunCombinedModel;
    totalCost = cost;

    pCa = 4;
    if drawPlots
        figure(34);title('pCa 4');
    end

    drawPlots = true;
    RunCombinedModel;
    totalCost = totalCost + cost;
end
