function a = extractSawtoothAttributes(out, meta, srsStrainWin)
%EXTRACTSAWTOOTHATTRIBUTES Outcome measures for one sawtooth simulation.
%
%   a = extractSawtoothAttributes(out, meta, srsStrainWin)
%
%   out          evaluateModel output struct (needs t, Force, FXB, FXBPassive,
%                p1_0, p2_0, SR).
%   meta         struct from buildSawtoothProtocol.
%   srsStrainWin strain window (fraction of L0) over which the STRAIN-MATCHED
%                short-range stiffness is regressed. The skeletal 10 ms window
%                at 0.45 L0/s spans 0.0045 L0; passing that value makes the two
%                preparations comparable even when ramp speed differs.
%
%   Short-range stiffness follows the paper: slope of the NORMALISED
%   force-length curve, by linear regression, over a window starting at stretch
%   onset. Two versions are returned because they diverge badly for cardiac:
%     SRS_t10  - the literal recipe: first 0.010 s of the ramp.
%     SRS_eps  - the same slope over the strain window instead of the time one.
%   Units of both: F0 per L0 (normalised force per muscle length).

    t = out.t(:); F = out.Force(:);
    F0 = interp1(t, F, meta.t0, 'linear', 'extrap');   % isometric at onset
    if F0 == 0; F0 = 1; end

    a.F0_kPa = F0;

    % Commanded length (ML). The triangle is exact, so no numeric integration.
    Lof = @(tt) meta.L0 + meta.v * max(0, min(tt, meta.t1) - meta.t0) ...
                        - meta.v * max(0, min(tt, meta.t2) - meta.t1);

    %% Short-range stiffness -------------------------------------------------
    a.SRS_t10 = iSlope(t, F/F0, Lof, meta.t0, min(meta.t0 + 0.010, meta.t1));
    a.SRS_eps = iSlope(t, F/F0, Lof, meta.t0, ...
                       min(meta.t0 + srsStrainWin*meta.L0/meta.v, meta.t1));
    a.srs_t10_dur    = min(0.010, meta.dt);                     % s actually used
    a.srs_t10_strain = a.srs_t10_dur * meta.v / meta.L0;        % L0 covered
    a.srs_t10_nmHS   = a.srs_t10_strain * meta.L0 * meta.ML * 1000 / 2;
    a.srs_eps_strain = min(srsStrainWin, meta.fracAmp);
    a.srs_eps_dur    = min(srsStrainWin*meta.L0/meta.v, meta.dt);

    % How many RAW solver points fall inside each window. ode15s picks its own
    % steps, and a 15 ms ramp can be resolved by only a handful of them; if
    % these counts are small the slopes above are interpolation, not data.
    a.n_srs_t10 = nnz(t >= meta.t0 & t <= meta.t0 + a.srs_t10_dur);
    a.n_srs_eps = nnz(t >= meta.t0 & t <= meta.t0 + a.srs_eps_dur);
    a.n_ramp    = nnz(t >= meta.t0 & t <= meta.t1);

    %% Force landmarks (all normalised by isometric F0) ----------------------
    inRamp = t >= meta.t0 & t <= meta.t1;
    inRel  = t >= meta.t1 & t <= meta.t2;
    inPost = t >= meta.t2;

    % Peak on a uniform 20 kHz resample, not on the raw solver points. The ramp
    % force peak is BROAD (dF/dt ~ 0 for several ms), so a raw argmax hops
    % between neighbouring solver steps: perturbing the commanded velocity by
    % one ulp moved the reported peak time by 0.5 ms while leaving the peak
    % amplitude identical to 4 digits. Resampling makes it a property of the
    % trace rather than of the solver's step placement. Residual uncertainty on
    % peak_t is ~0.4 ms (see RunVerify.m, sampling convergence); peak_Fn is
    % determined to < 0.02%.
    tq = linspace(meta.t0, meta.t1, max(2001, ceil((meta.t1 - meta.t0)*2e4)))';
    Fq = interp1(t, F/F0, tq, 'linear', 'extrap');
    [a.peak_Fn, ip] = max(Fq);
    a.peak_t      = tq(ip) - meta.t0;
    a.peak_strain = a.peak_t * meta.v / meta.L0;              % L0
    a.peak_nmHS   = a.peak_strain * meta.L0 * meta.ML * 1000 / 2;
    a.yields      = a.peak_t < 0.995 * meta.dt;   % peak before apex = yield

    a.Fn_apex     = interp1(t, F/F0, meta.t1, 'linear', 'extrap');
    a.Fn_minRel   = min(F(inRel) / F0);
    a.Fn_return   = interp1(t, F/F0, meta.t2, 'linear', 'extrap');
    a.Fn_end      = F(end) / F0;
    a.RFE         = a.Fn_end - 1;                 % residual force enhancement
    a.slackedRel  = a.Fn_minRel <= 1e-6;          % went slack during release

    % Post-release transient. Force does not simply relax back to F0 after the
    % muscle is returned to L0 - it can OVERSHOOT, because the release leaves
    % more bridges attached than isometric. Broad peak, so resampled like the
    % ramp peak above.
    tp = linspace(meta.t2, min(meta.t2 + 0.30, t(end)), 3001)';
    Fp = interp1(t, F/F0, tp, 'linear', 'extrap');
    [a.Fn_postMax, ipm] = max(Fp);
    a.Fn_postMax_t = tp(ipm) - meta.t2;

    %% Work loop -------------------------------------------------------------
    % W = int F dL over each half (kPa*ML). Stretch does work ON the muscle;
    % release returns part of it. The difference is the hysteresis area, i.e.
    % the energy dissipated in one cycle.
    a.W_stretch = iWork(t, F, Lof, meta.t0, meta.t1);
    a.W_release = -iWork(t, F, Lof, meta.t1, meta.t2);
    a.W_hyst    = a.W_stretch - a.W_release;
    a.elasticReturn = a.W_release / max(a.W_stretch, eps);

    %% State-variable snapshots ---------------------------------------------
    att = out.p1_0(:) + out.p2_0(:);
    keys = {'t0', 't1', 't2'};
    for k = 1:numel(keys)
        tk = meta.(keys{k});
        a.(['att_' keys{k}]) = interp1(t, att, tk, 'linear', 'extrap');
        a.(['SR_'  keys{k}]) = interp1(t, out.SR(:), tk, 'linear', 'extrap');
    end
    a.att_min      = min(att(inRamp | inRel));
    a.att_end      = att(end);
    a.Fpas_apex    = interp1(t, out.FXBPassive(:), meta.t1, 'linear', 'extrap');
    a.Fact_apex    = interp1(t, out.FXB(:), meta.t1, 'linear', 'extrap');
    a.pasFrac_apex = a.Fpas_apex / max(a.Fpas_apex + a.Fact_apex, eps);
    a.recovered    = (F(end) - min(F(inPost))) / max(F0 - min(F(inPost)), eps);
end

% -------------------------------------------------------------------------
function b = iSlope(t, Fn, Lof, ta, tb)
%ISLOPE Regression slope dFn/dL over [ta, tb], on a uniform resample.
%   Uniform resampling matters: ode15s clusters its output where the solution
%   is stiff, which would otherwise weight the window unevenly.
    if tb <= ta; b = NaN; return; end
    tq = linspace(ta, tb, 201)';
    Fq = interp1(t, Fn, tq, 'linear', 'extrap');
    Lq = Lof(tq);
    p  = polyfit(Lq, Fq, 1);
    b  = p(1);
end

function W = iWork(t, F, Lof, ta, tb)
%IWORK int F dL over [ta, tb] (kPa*ML), trapezoidal on a uniform resample.
    if tb <= ta; W = 0; return; end
    tq = linspace(ta, tb, 2001)';
    Fq = interp1(t, F, tq, 'linear', 'extrap');
    Lq = Lof(tq);
    W  = trapz(Lq, Fq);
end
