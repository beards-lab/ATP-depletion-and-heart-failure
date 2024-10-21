dsc = load('DataStruct20241010.mat').dsc;

%% all input combinations

% relaxed
% ids = 1; arr = [2 3 4 5];
% x = [5.1073    0.1500    4.2836];
% ids = 1; arr = [9 8 7 6];
% x = [3.0652    0.2952    5.4493];

% PNB and MAVA, no Ca
% ids = 3; arr = [2 3 4 5];
% x = [5.2542    0.1328    3.1200];

% max Ca, PNB and Mava
ids = 3; arr = [7 8 9 11];
% last 10s
% x = [4.5449    0.2347    4.5626  -10.0684];
% first 10s
x = [12.4807    0.4731    5.0766    3.4441];
x = [12.4807    0.6731    5.0766    10.4441];

% 10C - no Ca, just PNB mava
% ids = 4; arr = [2 3 4 5];
% x = [3.3884    0.2061    4.6010];

% 10C - max Ca, PNB and Mava
% ids = 4; arr = [7 8 9 11];
% x = [4.2542    0.27    4.1200];

% extracted relaxed
% ids = 5; arr = [2 3 4 5];
% x = [2.7368    0.1600    2.4991];

% extracted activated
% ids = 6; arr = [2 3 4 5];
% x = [3.5907    0.1182    2.1046];

[Tarr Farr] = getArrs(dsc, ids, arr);

c = evalPowerFit(x, Farr, Tarr, true, [], false)

%%
fitfun = @(init) evalPowerFit(init, Farr, Tarr, true);
%%
%%
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 50);

% Tarr{4} = [];Farr{4} = [];Tarr{3} = [];Farr{3} = [];
% Tarr{2} = [];Farr{2} = [];Tarr{1} = [];Farr{1} = [];
init = x;
% init = [x 10];
% init = [10 0.2, 4];
fitfunOpt = @(init) evalPowerFit(init, Farr, Tarr, false);
x = fminsearch(fitfunOpt, init, options)

%%

% Tarr{3} = [];Tarr{2} = [];Tarr{1} = [];
% x = [15.7900    0.3980    2.0187   16.2781];
% x = [15.7900    0.73980    5.0187];

[c] = evalPowerFit(x, Farr, Tarr, true, [], false)


function [Tarr Farr] = getArrs(dsc, ids, arr)

    for iarr = 1:length(arr)
        Farr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.F;
        Tarr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.t - 10;
    
        i_cutoff = find(Farr{iarr} > 4 & Tarr{iarr} > 50, 1, 'last');
        Farr{iarr} = Farr{iarr}(1:i_cutoff);
        Tarr{iarr} = Tarr{iarr}(1:i_cutoff);
    end
end