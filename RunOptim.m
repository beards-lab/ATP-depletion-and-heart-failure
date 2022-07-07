%% set g0
load g0;
g = g0;
% original Dan's g0
g = [0.7009    1.2243    1.0965    1.8390    0.4718 2.3357    0.3960    0.2372    0.1465    0.9817 0.8737    1.8333    0.2916    0.9513    1.0085]

% optimized g0
g = [0.8713    3.5481    1.0423   -1.3882    0.0349 1.3336    0.2065    0.6945    0.1177    2.1158    1.2210    2.0750   -0.9861    1.1191    3.2453]

g0 = g;
save g0 g0;
%% set up the problem
fcn = @dPUdT;
E0 = evaluateProblem(fcn, g, true)
%%
options = optimset('Display','iter', 'TolFun', 1e-3, 'TolX', 1, 'PlotFcns', @optimplotfval, 'MaxIter', 1000, 'OutputFcn', @myoutput);
optimfun = @(g)evaluateProblem(fcn, g, false)

x = fminsearch(optimfun, g, options)

function stop = myoutput(x,optimvalues,state);
    stop = false;    
    if ~isequal(state,'iter') && mod(optimvalues.iteration, 10) == 1
        return;
    end
    evaluateProblem(fcn, g, true);

end