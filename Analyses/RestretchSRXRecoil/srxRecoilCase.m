function C = srxRecoilCase(pBase, name, over, amp, rsK_d)
%SRXRECOILCASE  One variant of the SRX-recoil test, scored on the full battery.
%
%   C = srxRecoilCase(PBASE, NAME, OVER, AMP, RSK_D)
%
%   OVER is a struct of parameter overrides applied to PBASE before getParams is
%   re-run (so expression fields and PU0 are rebuilt consistently). AMP is the
%   per-cycle re-stretch amplitude (ML) and RSK_D the data rsK, used for the
%   LEVEL / SLOPE decomposition that is the point of this analysis.
%
%   Returns the feature costs in BOTH norms (the optimiser minimises L2; the
%   reports elsewhere quote L1), the post-restretch rate per cycle, and the
%   amplitude slope.
%
%   See also: evalRs, AssessSnapshot (Analyses/RestretchRecoveryFit)

    p = pBase;
    f = fieldnames(over);
    for i = 1:numel(f)
        p.(f{i}) = over.(f{i});
    end
    p = getParams(p, [], true, false);

    C = struct('name', name, 'over', over, 'ok', false, 'err', '');

    r = evalRs(p, struct('slackOnly', false));
    C.ok = r.ok;
    if ~r.ok
        C.err = r.err;
        return;
    end

    [ct2, ~, ~] = evalFeatureCost(r.features_data, r.features_model, p.fn, 2);
    C.L1     = r.featTotal;
    C.L2     = sum(ct2);
    C.ct2    = ct2;
    C.ct1    = r.costw;
    C.fn     = p.fn;
    C.rsK_m  = r.rsK_m;
    C.rsK_x  = r.rsK_ratio;
    C.ktr_x  = r.ktr_ratio;
    C.steady = mean(r.steady_m);
    C.rsA_x  = mean(r.rsA_m, 'omitnan') / mean(r.rsA_d, 'omitnan');
    C.time   = r.time;
    C.features_model = r.features_model;

    % LEVEL vs SLOPE: the data's mild NEGATIVE amplitude dependence is the
    % structural residual this mechanism is meant to supply.
    ok = isfinite(r.rsK_m(:)) & isfinite(amp(:));
    if sum(ok) >= 2
        pm = polyfit(amp(ok), r.rsK_m(ok), 1);
        C.slope = pm(1);
        C.icept = pm(2);
    else
        C.slope = NaN; C.icept = NaN;
    end
    C.rsK_d = rsK_d;
end
