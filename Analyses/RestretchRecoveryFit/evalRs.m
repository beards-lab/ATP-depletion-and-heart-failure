function R = evalRs(params0, opts)
%EVALRS Evaluate the model and summarise the post-restretch recovery fit.
%
%   R = evalRs(PARAMS0)          full battery (FV + slack + passive)
%   R = evalRs(PARAMS0, opts)    opts.slackOnly = true -> slack protocol only
%
%   Returns a struct with the feature-cost breakdown and the model/data
%   post-restretch rates, so a lever sweep can see BOTH what it buys on rsK
%   and what it costs everywhere else.

    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'slackOnly'); opts.slackOnly = false; end

    if opts.slackOnly
        params0.RunForceVelocity = 0;
        params0.RunSlackPassive  = 0;
    end
    params0.RunKtr    = 0;
    params0.RunStairs = 0;
    params0.RunForceVelocityTime = 0;
    params0.EvalFeatures       = 1;
    params0.PlotEachSeparately = 0;
    params0.PlotFeatureFitting = 0;

    ATP_c         = params0.MgATP;          %#ok<NASGU>  (read by RunBakersExp)
    FV_velocities = params0.FV_velocities;  %#ok<NASGU>

    R = struct('ok', false, 'err', '', 'E', NaN, 'featTotal', NaN);
    tt = tic;
    try
        RunBakersExp;   % script: reads params0/ATP_c/FV_velocities from THIS workspace
    catch ME
        R.err = ME.message; R.time = toc(tt); return;
    end
    R.time = toc(tt);

    % Only score fn entries whose features actually exist in this run
    % (slackOnly drops the FV_* and PS_* features).
    fn = params0.fn;
    keep = true(1, numel(fn));
    if opts.slackOnly
        for i = 1:numel(fn)
            if contains(fn{i}, 'FV_') || contains(fn{i}, 'PS_'); keep(i) = false; end
        end
        fn = fn(keep);
    end
    [ct, w, c] = evalFeatureCost(features_data, features_model, fn, 1);

    R.ok        = true;
    R.E         = E;
    R.fn        = fn;
    R.cost      = c;
    R.weight    = w;
    R.costw     = ct;
    R.featTotal = sum(ct);
    R.features_model = features_model;
    R.features_data  = features_data;

    % headline observables
    R.rsK_m   = getf(features_model, 'rsK');
    R.rsK_d   = getf(features_data,  'rsK');
    R.ktr_m   = getf(features_model, 'ktr');
    R.ktr_d   = getf(features_data,  'ktr');
    R.steady_m= getf(features_model, 'steady');
    R.steady_d= getf(features_data,  'steady');
    R.rsA_m   = getf(features_model, 'rsA');
    R.rsA_d   = getf(features_data,  'rsA');
    R.rsK_ratio  = mean(R.rsK_m, 'omitnan') / mean(R.rsK_d, 'omitnan');
    R.ktr_ratio  = mean(R.ktr_m, 'omitnan') / mean(R.ktr_d, 'omitnan');
    R.Fiso_m     = mean(R.steady_m(1:min(2,end)), 'omitnan');

    % named per-feature costs, for the sweep table
    R.byName = struct();
    for i = 1:numel(fn)
        nm = regexprep(strtok(fn{i}, '|'), '[^A-Za-z0-9_]', '');
        if isempty(nm); nm = sprintf('f%d', i); end
        if isfield(R.byName, nm); nm = sprintf('%s_%d', nm, i); end
        R.byName.(nm) = ct(i);
    end
end

function v = getf(s, f)
    if isfield(s, f); v = s.(f); else; v = NaN; end
end
