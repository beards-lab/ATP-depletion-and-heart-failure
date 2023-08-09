x0 = ones(1, 10);
x0 = [ 1.0057    1.5981    1.3605    0.0917    1.4696    1.1553    1.0838    0.6758    1.1417]; % first shot

x0 = [ 0.9973    2.9945    2.8851    0.1308    -0.8659    1.3337    1.2446    0.7553     0.9829]; % negative param

% now, fixing the offset to 2.0 at 0.4 offset
x0 = [    1.0394    1.1020    0.8118    0.9860   0.7941    0.9136    0.9657    1.0193    1.4571]; % not really wokring

% optimizing with variable offset, not evaluating ramp-up, cost = 300.7
mod = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345    1.2716];


%%
tic
drawPlots =false;
evalCombined(x0)
toc
%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);

x = fminsearch(@evalCombined, x0, options);

save x;
%%
function cost = evalCombined(mod)
    drawPlots = false;
    RunCombinedModel;
end
