function cost = captureOptimIter(params0, iter_num, fval, fval_pre, namePrefix)
% Run a visualization pass and save the plotFeatures figure + cost to ../params/.
% namePrefix e.g. 'ModelOptParams_TL3' — iter number is appended automatically.

    p = params0;
    p.PlotEachSeparately = 1;

    prevVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanup = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));

    outBase = fullfile('..', 'params', sprintf('%s_iter_%g', namePrefix, iter_num));

    params0 = p;
    figure(1);clf;
    RunBakersExp;
    print(figure(1), [outBase '_timecourse.png'], '-dpng', '-r150');

    cost = plotFeatures(features_data, features_model, [], p.fn);
    % plotFeatures always renders into figure 80085

    print(figure(80085), [outBase '_feats.png'], '-dpng', '-r150');
    save([outBase '_cost.mat'], 'cost', 'fval', 'fval_pre');
    fprintf('[captureOptimIter] iter %d  cost=%.3f  -> %s\n', iter_num, sum(cost), outBase);
end
