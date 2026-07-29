function L = windowFluxLedger(out, w, params)
%WINDOWFLUXLEDGER Integrate every recorded flux over one recovery window.
%
%   L = windowFluxLedger(OUT, W, PARAMS)
%
%   OUT    output struct from evaluateModel (needs the out.R* flux fields)
%   W      one element of the struct array from captureWindows (.i0 .i63 .i1)
%   PARAMS the params struct actually used
%
%   Both recovery windows are pure holds (vel = 0), so the entire k63 difference
%   between them is a property of the INITIAL STATE, not of the protocol. This
%   function attributes the force rise to its underlying fluxes.
%
%   The recorded rates are already strain-integrated scalars: the ODE assembles
%   `rates = [RTD, RDT, RD1, sum([R1D,R12,R21,R2,XB_Ripped],1)*dS, ...]`, so
%   every out.R* below is a total flux in units of (fraction of heads)/s.
%
%   Key fields of L:
%     .nu     new-head share at t63.  ->1 recovery is fresh attachment
%                                     ->0 recovery is redistribution of survivors
%     .I      internal-traffic index (p1<->p2 churn / attach+detach traffic)
%     .T      turnover ratio, new heads delivered per head present at window start
%     .N_*    cumulative fluxes over [i0,i63]
%     .srx    SRX sub-ledger (entry, exit, peak excursion)
%     .closure  mass-balance residuals; must be ~0 or the ledger is invalid

    i0 = w.i0; i63 = w.i63; i1 = w.i1;
    seg = i0:i63;  t = out.t(:);

    %% ---- guards: if any fails the labels below are wrong ----
    % The protocol's "hold" rows carry a ~0.002 um/s drift, not exact zero, and
    % kA2hop is gated on vel>0 (not on magnitude), so a positive drift can fire an
    % UNRECORDED flux. Report the range; the mass-closure check below is the real
    % arbiter of whether anything unrecorded actually moved.
    vw = out.v(i0:i1);
    L.guard.vmax = max(vw); L.guard.vmin = min(vw);
    L.guard.holdOK  = max(abs(vw)) < 0.01;      % vs 72-214 um/s during the ramps
    L.guard.hopRisk = any(vw > 0) && params.kA2hop > 0;
    L.guard.satOff   = ~params.UseTargetZoneSaturation && ~params.UseVernierVelocity;
    L.guard.pathways = params.k2rip == 0 && params.kA2re == 0 && ...
                       params.ksrd == 0 && params.kmsr == 0;
    if ~L.guard.holdOK
        warning('windowFluxLedger:vel','%s: |vel| up to %.3g um/s — not a hold.', w.name, max(abs(vw)));
    end
    if ~L.guard.satOff
        warning('windowFluxLedger:sat','%s: per-bin f_sat is active — out.RD1 is NOT the delivered attachment flux.', w.name);
    end

    %% ---- cumulative fluxes over [i0, i63] ----
    % out.R* are stored as ROW vectors while t is a column; force columns or the
    % elementwise ops below implicitly expand into matrices.
    col = @(fld) reshape(out.(fld)(seg), [], 1);
    ig  = @(fld) trapz(t(seg), col(fld));
    L.N_att  = ig('RD1');      % PD -> p1 : the ONLY source of new bound heads
    L.N_det1 = ig('R1D');      % p1 -> PD
    L.N_12   = ig('R12');      % p1 -> p2
    L.N_21   = ig('R21');      % p2 -> p1
    L.N_det2 = ig('R2T');      % p2 -> PT (stroke completes)
    L.N_hyd  = ig('RTD');      % PT -> PD (hydrolysis)
    L.N_hydr = ig('RDT');      % PD -> PT

    %% ---- mass-balance closure (validity check) ----
    d = @(fld) out.(fld)(i63) - out.(fld)(i0);
    L.closure.p1 = d('p1_0') - (L.N_att - L.N_det1 - L.N_12 + L.N_21);
    L.closure.p2 = d('p2_0') - (L.N_12 - L.N_21 - L.N_det2);
    L.closure.ok = max(abs([L.closure.p1 L.closure.p2])) < 1e-3;

    %% ---- the fractions ----
    Pb  = out.p1_0(:) + out.p2_0(:);
    % survivor fraction: mean-field first-order loss of the pool present at i0
    loss = (col('R1D') + col('R2T')) ./ max(Pb(seg), 1e-12);
    L.S   = exp(-trapz(t(seg), loss));
    L.nu  = 1 - L.S * Pb(i0) / max(Pb(i63), 1e-12);
    L.T   = L.N_att / max(Pb(i0), 1e-12);
    L.I   = (L.N_12 + L.N_21) / max(L.N_att + L.N_det1 + L.N_det2, 1e-12);

    %% ---- SRX sub-ledger: the H_SRX test ----
    L.srx.N_PT2SR  = ig('RPT2SR');    % PT -> SR   (entry, force-SUPPRESSED)
    L.srx.N_SRD2PD = ig('RSRD2PD');   % SRD -> PD  (exit,  force-ACTIVATED)
    L.srx.net      = L.srx.N_PT2SR - L.srx.N_SRD2PD;
    SRX = out.SR(:) + out.SRD(:);
    L.srx.at_i0    = SRX(i0);
    L.srx.peak     = max(SRX(i0:i1));
    L.srx.excursion= L.srx.peak - SRX(i0);
    % instantaneous gate rates at the window start, for the narrative
    F0 = out.Force(i0);
    L.srx.entryRate0 = params.ksr0  * exp(-F0/params.sigma2);
    L.srx.exitRate0  = params.kmsrd * exp( F0/params.sigma_srd1);

    %% ---- pool state at the window start ----
    L.state.F      = out.Force(i0);
    L.state.PD     = out.PuR(i0);
    L.state.PT     = out.PuATP(i0);
    L.state.p1     = out.p1_0(i0);
    L.state.p2     = out.p2_0(i0);
    L.state.Pb0    = Pb(i0);
    L.state.Pb63   = Pb(i63);
    L.dur_ms       = 1000*(t(i63) - t(i0));
    L.name         = w.name;
    L.k63          = w.k63;
end
