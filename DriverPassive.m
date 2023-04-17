% driver passive

% set the params

% set the modifiers
% let it optim
clear;

mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma'};
% optimized for 1s
opt_mods = [  1.1797,     0.6314,    1.1477,    0.5181,    0.5833,    1.9550,    1.6055, 1];
% optimized for 100s
opt_mods =   [  0.0298    0.0290    0.3837    1.2774    1.0610    2.0034     1.7758, 1];
% hand tuned for 100ms
opt_mods =   [  0.0298    0.10290    0.006837    8.2774    1.0610    2.0034     1.7758, 1];
% opt_mods =   [  0.06298    0.060    0.3837    1.2774    4.0610    2.0034     1.7758, 1];

% optimized for 0.1s
% opt_mods = [0.0298    0.0290    200   0.1   1.0610    2.0034    1.7758, 1];
% optimized for 1
plotEach = true;

% ramp duration
rd = 0.100; 

tic
datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
time_end = datatable(end, 1);
toc 
tic

evaluatePassive;
toc
%%
figure(1010);clf;
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% g = [1, 1, 1, 1, 1, 1, 1, 1];
% g = [1.2539    0.4422];

% ident all params
% optimfun = @(g)evalPassiveCost(g, datatable, rd);

% ident just a subset
ps = 3:4;
disp('Optimizing for ') 
disp(mods(ps))
x0 = opt_mods(ps);
optimfun = @(g)evalPassiveCost([opt_mods(1:2),g(1, :),opt_mods(5:end)], datatable, rd);
%
x = fminsearch(optimfun, x0, options);
x


%%
function Es = evalPassiveCost(opt_mods, datatable, rd)
if any(mods) < 0 
    Es = Inf;
    return;
end
plotEach = false;

% datastruct = load('data/bakers_passiveStretch_1000ms.mat');
% datatable = datastruct.datatable;
evaluatePassive;


end


