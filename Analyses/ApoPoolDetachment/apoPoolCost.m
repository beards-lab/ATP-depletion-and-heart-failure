function [C, D] = apoPoolCost(g, base, cfg)
%APOPOOLCOST Objective for the ATP-ceiling / apo-pool refit.
%
%   [C, D] = apoPoolCost(G, BASE, CFG)
%
%   G    - multiplicative modifiers, one per CFG.mods entry (see RunOptimApoPool).
%   BASE - seed params struct (already through getParams).
%   CFG  - struct from RunOptimApoPool: .mods, .targets, .w, .bounds.
%
%   C    - scalar cost. D is a breakdown struct (for logging / plotting).
%
%   WHY THIS EXISTS RATHER THAN evaluateBakersExp:
%   the existing cost scores force TRACES. The defect being fixed is a RATE
%   (post-restretch recovery, 73-100 /s in the model vs 44-59 /s in the data),
%   and a trace MSE is almost blind to it — a fit can be 2x wrong on the
%   recovery rate for a few percent of MSE. Handed only the trace cost, the
%   optimiser's cheapest move is to switch the new mechanisms straight back off
%   (k2rip -> 0, k2max -> inf) and recover the old fit. So the recovery rates
%   measured by recoveryWindows are scored EXPLICITLY here, alongside the trace
%   costs and the restretch peak/valley that the ceiling is known to distort.
%
%   PHYSIOLOGY BOUNDS ARE PART OF THE OBJECTIVE, DELIBERATELY.
%   evalPhysiologyCost is included with the alpha-MHC bounds active. The
%   mechanism is only known to work at k2max = 150-300 /s (below the alpha-MHC
%   ADP-release range) and k2rip = 3000-10000 (6-20x over bounds.k2rip.ub).
%   If it can only be fitted by leaving those ranges, this objective makes that
%   show up as cost rather than as an unnoticed extrapolation.
%   See conclusions.md, "CRITICAL REVIEW OF THE PHYSIOLOGY".
%
%   Every experiment is run ONCE and both its trace cost and its window rates
%   are taken from the same output, so an evaluation costs one slack + one ktr
%   (+ one FV), not two of each.

    C = NaN; D = struct();
    p = base;

    % ---- apply modifiers ----
    for i = 1:numel(cfg.mods)
        p.(cfg.mods{i}) = cfg.seed.(cfg.mods{i}) * g(i);
    end
    p.PlotEachSeparately = 0; p.ghostSave = ''; p.ghostLoad = '';
    p.EvalFeatures = 0; p.PlotFeatureFitting = 0;

    try
        p = getParams(p, p.g, true);
    catch
        C = 1e9; return;
    end

    T = cfg.targets;
    try
        % ================= SLACK =================
        ps = p; ps.RunSlackSegments = 'All';
        [E_slack, os] = runSlackExperiment(ps);
        Fiso = interp1(os.t, os.Force, T.svt(T.iS(1),1) - 0.02, 'linear', 'extrap');
        if ~isfinite(Fiso) || Fiso <= 0; C = 1e9; return; end
        ws = recoveryWindows(os.t, os.Force, T.svt, 'slack', Fiso);

        kS = arrayfun(@(c) pick(ws,'postSlack',c),     1:4);
        kR = arrayfun(@(c) pick(ws,'postRestretch',c), 1:4);

        % restretch peak and valley — the ceiling is known to inflate both
        gm = [true; diff(os.t(:)) > 0]; tm = os.t(gm); Fm = os.Force(gm)/Fiso;
        pk1 = mean(arrayfun(@(c) max(Fm(tm >= T.svt(T.iR(c),1) & ...
                                        tm <= T.svt(T.iR(c)+1,1)+0.005)), 1:4));
        vy  = mean(arrayfun(@(c) getf(ws,'postRestretch',c,'B'), 1:4));

        % ================= KTR =================
        [E_ktr, ok] = runKtrExperiment(p, []);
        FisoK = interp1(ok.t, ok.Force, T.kvt(2,1), 'linear', 'extrap');
        if ~isfinite(FisoK) || FisoK <= 0; C = 1e9; return; end
        wk = recoveryWindows(ok.t, ok.Force, T.kvt, 'ktr', FisoK);

        % ================= FORCE-VELOCITY =================
        E_fv = 0;
        if cfg.w.fv > 0
            E_fv = runFVExperiment(p, T.ATP_c, T.Data_ATP);
        end
    catch e
        % a diverged / stalled solve is a bad point, not a crash
        D.err = e.message; C = 1e8; return;
    end

    % ================= assemble =================
    rel = @(x, t) ((x - t) ./ t).^2;              % relative squared error

    D.E_slack = E_slack(1);
    D.E_ktr   = E_ktr;
    D.E_fv    = E_fv;
    D.kS      = mean(kS);   D.kR = mean(kR);   D.kK = wk.k63;
    D.peak1   = pk1;        D.valley = vy;

    D.c_rate  = cfg.w.kS * mean(rel(kS, T.kS)) ...
              + cfg.w.kR * mean(rel(kR, T.kR)) ...
              + cfg.w.kK * rel(wk.k63, T.kK);
    D.c_shape = cfg.w.peak * rel(pk1, T.peak1) + cfg.w.vall * rel(vy, T.valley);
    D.c_trace = cfg.w.slack * (E_slack(1)/T.E_slack_ref) ...
              + cfg.w.ktr   * (E_ktr    /T.E_ktr_ref) ...
              + cfg.w.fv    * (E_fv     /T.E_fv_ref);

    D.c_phys = 0;
    if cfg.w.phys > 0
        D.c_phys = cfg.w.phys * evalPhysiologyCost(p, cfg.bounds);
    end

    C = D.c_rate + D.c_shape + D.c_trace + D.c_phys;
    if ~isfinite(C); C = 1e9; end
    D.C = C;
end

% ---------------------------------------------------------------
function k = pick(ws, ty, c)
    w = ws(strcmp({ws.type}, ty) & [ws.cyc] == c);
    if isempty(w) || ~isfinite(w.k63); k = 1e3; else; k = w.k63; end
end

function v = getf(ws, ty, c, f)
    w = ws(strcmp({ws.type}, ty) & [ws.cyc] == c);
    if isempty(w); v = NaN; else; v = w.(f); end
end
