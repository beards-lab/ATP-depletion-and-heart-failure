function [features_model, features_data, out_pas] = runPassiveExperiment(params0)
%RUNPASSIVEEXPERIMENT Run the passive (blebbistatin/mavacamten) slack protocol.
%
%   [features_model, features_data, out_pas] = runPassiveExperiment(params0)
%
%   Simulates the passive slack-release protocol with active cross-bridge
%   attachment switched OFF (ka = 0), reproducing the PNB/Mava condition in
%   which only titin + serial-elastic mechanics bear load. It exists so the
%   PASSIVE observations can be added to the JOINT feature cost (RunBakersExp
%   -> evalFeatureCost on params0.fn): the passive fit alone is non-unique (a
%   sloppy k_pas/gamma/Lsc0/kSE/kSE_M/eta_M manifold), so identifying it in
%   isolation lands on an arbitrary point that may be pessimal for the active
%   fit. Scoring passive INSIDE the joint cost lets the active data select,
%   along the passive soft directions, the point also consistent with active.
%
%   Protocol / data are loaded from a .mat named by
%   params0.passive_velocitytableonfile (default the ActivePNBMava slack file),
%   holding {datatable, velocitytable}. The velocitytable has one 4-row block
%   per slack cycle: [shorten; hold-short; restretch; hold-long].
%
%   Returns ONLY feature structs (no direct cost): the cost is taken by
%   evalFeatureCost from params0.fn exactly like every other feature, so the
%   passive terms are weighted/selected the same way. All feature names are
%   prefixed PS_ (passive-slack). Features (per cycle, from BOTH model & data):
%     PS_restretchPeak  (1xNc) - peak force over the restretch+hold window
%     PS_rampupRMSE     (1xNc) - RMSE(model,data) over the restretch ramp;   data=0
%     PS_holdDecayRMSE  (1xNc) - RMSE(model,data) over the post-restretch hold; data=0
%     PS_steady22       (1x1)  - mean force in the long pre-slack hold  (SL 2.2)
%     PS_steady20       (1x1)  - mean force in the final hold           (SL 2.0)
%   The RMSE features are model-vs-data residuals, so their data value is 0 by
%   definition; in the feature cost (costExp=2) RMSE^2 recovers the window MSE.

    modelFcn = str2func(params0.modelFcn);

    features_model = struct();
    features_data  = struct();

    %% --- Load protocol + data ---
    if isfield(params0, 'passive_velocitytableonfile') && ~isempty(params0.passive_velocitytableonfile)
        fname = params0.passive_velocitytableonfile;
    else
        fname = 'protocol_03_27_2026_ActivePNBMava_slack.mat';
    end
    datastruct    = load(['data/' fname]);
    datatable     = datastruct.datatable;         % [t, ML, force]
    velocitytable = datastruct.velocitytable;     % [t, vel, 2*vel, ML]
    dt_t = datatable(:, 1);
    dt_F = datatable(:, 3);

    %% --- Passive condition: no active attachment ---
    params = params0;
    ka_passive = 0;
    if isfield(params0, 'PassiveKa'); ka_passive = params0.PassiveKa; end
    params.ka = ka_passive;
    params.SL0 = 2 * velocitytable(1, 4);         % ML(1) -> SL0 (col4 is ML)

    seg = find(velocitytable(:, 2) < -1)';        % one row per slack cycle
    Nc  = numel(seg);
    if isempty(seg)
        error('runPassiveExperiment:noSlacks', 'No slack cycles (vel<-1) in %s', fname);
    end
    velocitytable(1, 1) = velocitytable(seg(1), 1) - 2;   % short warm-start

    params = getParams(params, params.g, true);
    params.Velocity = velocitytable(:, 2);
    [~, out_pas] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    out_pas.datatable = datatable;

    % Monotone model time for interpolation
    mono = makeMonotonous(out_pas.t);
    m_t  = out_pas.t(mono);
    m_F  = out_pas.Force(mono);
    Fmod = @(tq) interp1(m_t, m_F, tq, 'linear', 'extrap');

    %% --- Per-cycle feature extraction ---
    peak_m = nan(1, Nc); peak_d = nan(1, Nc);
    ramp_m = nan(1, Nc); hold_m = nan(1, Nc);
    for i = 1:Nc
        s = seg(i);
        t_rs   = velocitytable(s + 2, 1);                          % restretch start
        t_hold = velocitytable(s + 3, 1);                          % hold-long start
        if s + 4 > size(velocitytable, 1)
            t_end = velocitytable(end, 1);
        else
            t_end = velocitytable(s + 4, 1);                       % hold-long end
        end

        % restretch peak: max force over restretch ramp + hold (model & data)
        wp_d = dt_t >= t_rs & dt_t <= t_end;
        if any(wp_d); peak_d(i) = max(dt_F(wp_d)); end
        wp_m = m_t >= t_rs & m_t <= t_end;
        if any(wp_m); peak_m(i) = max(m_F(wp_m)); end

        % rampup / hold-decay RMSE (model vs data)
        ramp_m(i) = iWindowRMSE(dt_t, dt_F, Fmod, t_rs, t_hold);
        hold_m(i) = iWindowRMSE(dt_t, dt_F, Fmod, t_hold, t_end);
    end

    %% --- Steady states ---
    t_sh1 = velocitytable(seg(1), 1);
    wb_d  = dt_t >= t_sh1 - 0.40 & dt_t < t_sh1 - 0.005;
    wb_m  = m_t  >= t_sh1 - 0.40 & m_t  < t_sh1 - 0.005;
    s22_d = iMeanOr(dt_F(wb_d)); s22_m = iMeanOr(m_F(wb_m));
    t_last = velocitytable(end, 1);
    we_d   = dt_t >= t_last - 0.50 & dt_t <= t_last;
    we_m   = m_t  >= t_last - 0.50 & m_t  <= t_last;
    s20_d  = iMeanOr(dt_F(we_d)); s20_m = iMeanOr(m_F(we_m));

    %% --- Assemble feature structs (PS_ = passive-slack; RMSE data = 0 by def) ---
    features_data.PS_restretchPeak = peak_d;
    features_data.PS_rampupRMSE    = zeros(1, Nc);
    features_data.PS_holdDecayRMSE = zeros(1, Nc);
    features_data.PS_steady22      = s22_d;
    features_data.PS_steady20      = s20_d;

    features_model.PS_restretchPeak = peak_m;
    features_model.PS_rampupRMSE    = ramp_m;
    features_model.PS_holdDecayRMSE = hold_m;
    features_model.PS_steady22      = s22_m;
    features_model.PS_steady20      = s20_m;

    %% --- Plot ---
    if params.PlotEachSeparately
        if params.PlotFullscreen; clf; else; nexttile; end
        hold on; box on;
        plot(dt_t, dt_F, 'k-', 'LineWidth', 0.75);
        plot(m_t, m_F, 'b-', 'LineWidth', 1.5);
        for i = 1:Nc
            s = seg(i); t_rs = velocitytable(s + 2, 1);
            hp = min(s + 3, size(velocitytable, 1));
            plot(velocitytable(hp, 1), peak_m(i), 'b^', 'MarkerSize', 8, 'LineWidth', 1.5);
            plot(velocitytable(hp, 1), peak_d(i), 'kv', 'MarkerSize', 8, 'LineWidth', 1.5);
            xline(t_rs, ':', 'Color', [0.6 0.6 0.6]);
        end
        yline(s22_d, 'k--'); yline(s22_m, 'b--');
        xlabel('$t$ (sec.)', 'interpreter', 'latex', 'fontsize', 14);
        ylabel('Passive force (kPa)', 'interpreter', 'latex', 'fontsize', 14);
        title('Passive (ka=0): titin + SE mechanics');
        legend('data', 'model', 'Location', 'best');
        xlim([t_sh1 - 0.3, t_last]);
    end
end

% ----------------------------------------------------------------------------
function r = iWindowRMSE(dt_t, dt_F, Fmod, ta, tb)
%IWINDOWRMSE RMSE of model vs data force over [ta, tb], on the data time grid
%   (data is the reference; its self-error is 0 by definition). RMSE (not MSE)
%   so that the feature-cost square, RMSE^2, recovers the window MSE.
    w = dt_t >= ta & dt_t <= tb;
    if ~any(w)
        r = 0; return;
    end
    resid = dt_F(w) - Fmod(dt_t(w));
    r     = sqrt(mean(resid.^2, 'omitnan'));
end

function m = iMeanOr(v)
%IMEANOR Mean of v, or NaN if empty.
    if isempty(v); m = NaN; else; m = mean(v, 'omitnan'); end
end
