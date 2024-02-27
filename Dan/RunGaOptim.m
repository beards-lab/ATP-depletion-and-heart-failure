mod = [0.1730    0.0718    0.8944    1.2103    0.4733    1.1266    0.1235    0.8796    0.2051    1.9012         0    1.4225    0.0027    0.5112    2.5843       NaN       NaN    1.0000    1.0000       NaN       NaN    0.9521];
%%  better mods

logmod = [-0.7508   -1.1525   -0.0499    0.0828   -0.3219    0.0526   -0.8812   -0.0560   -0.7104    0.2722   -0.2797   -0.0222];
modSel = [1     2     3     4     5     6     7     8     9    10    14    22];
mod(modSel) = 10.^logmod;    

%%
ga_Opts = gaoptimset('InitialPopulation', mod(modSel), ...
    'PopulationSize',64, ...            % 250
    'Display','iter', ...
    'UseParallel',true, 'TimeLimit', 60*60*3);
%%
tic
% mod = [2.6973    0.0053    1.5531    0.1331    0.2561    0.2523    3.4075   0.0220   31.3090    2.7587    3.8009    0.1083    0.0110    2.7344    2.5843       NaN       NaN    1.0000    1.0000       NaN       NaN    0.9502];
% mod(1) = 1;
evalLogCombined(log10(mod(modSel)))
toc
%%
modSel = [1:10 14 22];
% Ng = length(modSel);
ub = 100*ones(1, length(mod));
ub([3, 4, 5, 10, 14]) = [10 10 10, 10, 10];
lb = 0.001*ones(1, length(mod));
lb([3, 4, 5, 10, 14]) = [.1 .1 .1 .1 .1];
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel, [4.4 11]);
init = log10(mod(modSel));
%%
[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(evalLogCombined,length(modSel), ...
    [],[],[],[],...
    log10(lb(modSel)),log10(ub(modSel)),[],ga_Opts);

mod(modSel) = 10.^p_OptimGA;

%%
function totalCost = evalCombined(optMods, mod, modSel, pCas)
    if nargin < 4
        % evaluate pCa 4.4 and 10 by default
        pCas = [4.4 10];
    end
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
        % clf;
        % plotParams(baseMods, [-2, 1]);hold on;
        plotParams(mod, [-2, 1]);hold on;
    end


    %% pCa 4
    if ismember(4.4, pCas)
        pCa = 4.4;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end
    %% pCa 5.5
    if ismember(5.5, pCas)
        pCa = 5.5;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end    
    %% pCa 5.75
    if ismember(5.75, pCas)
        pCa = 5.75;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% pCa 6
    if ismember(6, pCas)
        pCa = 6;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% no Ca, but evaluate PEVK binding
    if ismember(10, pCas)    
        pCa = 10;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end
    %% no Ca, PEVK binding not evaluated
    if ismember(11, pCas)    
        pCa = 11;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end    
% return
end 

function plotOnBackground(drawPlots, pCa)
    figInd = 100 + round(pCa*10);
    if drawPlots
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = ['pCa ' num2str(pCa)];
        end
    end
end

function cost = isolateRunCombinedModel(mod, pCa, drawPlots)
% just to isolate the script, so the variables can't intervene
    % drawPlots = true;
    RunCombinedModel;
end

function plotParams(mod, mm, resetGca)

    % modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU'};
    modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa'};

    indxs = strcat(string(1:length(modNames)), ':');
    modNames = strcat(indxs, modNames);

    if nargin < 2 || isempty(mm)
        % try
        %     mima = get(gca, 'RLim');
        %     mi = mima(1);ma=mima(2);
        % catch
            mi = min(floor(log10(mod)));
            ma = mi + 3;
        % end
    else
        mi = mm(1);
        ma = mm(2);
    end
    n = length(mod);
    
    polarplot(linspace(2*pi/n, 2*pi, n), max(1e-3, log10(mod) - mi), 'x-', LineWidth=2);
    
    if ~strcmp('degrees', get(gca, 'ThetaAxisUnits')) || (nargin >= 3 && resetGca)
        % degrees means it is not adjusted yet 
        % or the param has not been provided 
        % or we want to reset
        return
    end
    rng = mi:ma;
    set(gca, 'ThetaAxisUnits', 'radians', ...
        'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
        'Rlim', [min(rng) max(rng)]-mi, 'RTick', rng - mi, 'RTickLabel', 10.^rng);
end
