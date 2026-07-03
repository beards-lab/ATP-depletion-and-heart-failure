function [tc, cv, fm, fd] = costOfSnap(snapfile, fn, mrt, extra)
% COSTOFSNAP  Evaluate a param snapshot against a feature list, return total cost.
%   [tc, cv, fm, fd] = costOfSnap(SNAPFILE, FN, MRT, EXTRA)
%   Loads SNAPFILE (a writeParamsToMFile snapshot) into params0, runs FV+slack,
%   and returns total feature cost TC (=sum of evalFeatureCost with FN), the per-
%   feature cost vector CV, model features FM and data features FD.
%   MRT = MaxRunTime (s). EXTRA (optional struct) overrides fields on params0
%   after load (e.g. struct('UseVelGaussAttachment',1,'v_att_center',-1)).
    if nargin < 3 || isempty(mrt); mrt = 120; end
    if nargin < 4; extra = struct(); end
    params0 = getParams();
    run(snapfile);
    if ~isfield(params0,'RegAvailShorteningOnly'); params0.RegAvailShorteningOnly = true; end
    ef = fieldnames(extra);
    for i = 1:numel(ef); params0.(ef{i}) = extra.(ef{i}); end
    params0 = getParams(params0, [], true, false);
    params0.RunForceVelocity=true; params0.RunKtr=false; params0.RunStairs=false;
    params0.RunForceVelocityTime=false; params0.RunSlack=true; params0.RunSlackSegments='All';
    params0.FV_velocities=-[0 0.5 1 2 4]; params0.EvalFeatures=true;
    params0.PlotEachSeparately=0; params0.PlotFeatureFitting=0; params0.BreakOnODEUnstable=false;
    params0.MaxRunTime=mrt; params0.fn=fn;
    RunBakersExp;
    [cv,~,~] = evalFeatureCost(features_data, features_model, fn, 2);
    tc = sum(cv,'omitnan');
    fm = features_model; fd = features_data;
end
