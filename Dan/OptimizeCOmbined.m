x0 = ones(1, 9);
x0 = [ 1.0057    1.5981    1.3605    0.0917    1.4696    1.1553    1.0838    0.6758    1.1417]; % first shot

x0 = [ 0.9973    2.9945    2.8851    0.1308    -0.8659    1.3337    1.2446    0.7553     0.9829]; % negative param

% now, fixing the offset to 2.0 at 0.4 offset
x0 = [    1.0394    1.1020    0.8118    0.9860

  Columns 5 through 8

    0.7941    0.9136    0.9657    1.0193

  Column 9

    1.4571];
%%
evalCombined(x0)
%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 1500);

x = fminsearch(@evalCombined, x0, options);

save x;
%%
function cost = evalCombined(mod)
    drawPlots = false;
    RunCombinedModel;
end
