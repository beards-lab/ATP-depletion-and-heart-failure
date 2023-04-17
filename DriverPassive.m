% driver passive

% set the params

% set the modifiers
% let it optim


mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma'}
opt_mods = ones(7, 1);
plotEach = true;

tic
datastruct = load('data/bakers_passiveStretch_1000ms.mat');
datatable = datastruct.datatable;
toc 
tic
evaluatePassive;
toc

figure(1010);clf;
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.1, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];
optimfun = @(g)evalPassiveCost(opt_mods, datatable);
x = fminsearch(optimfun, g, options);
x


%%
function Es = evalPassiveCost(opt_mods, datatable)

mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma'};
opt_mods = ones(7, 1);
plotEach = false;

% datastruct = load('data/bakers_passiveStretch_1000ms.mat');
% datatable = datastruct.datatable;
evaluatePassive;


end


