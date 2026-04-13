%% CompareProtocols.m
% Universal cross-dataset comparison of slack, ktr, step perturbation, and
% staircase features.
%
% HOW TO ADD A DATASET
%   Append one struct entry to DS below with the relevant fields set.
%   Feature sections that have no corresponding data_* range are skipped
%   automatically — no other code changes needed.
%
% UNITS
%   All lengths are converted to µm inside loadDataset() via Lo_ref_um.
%   Baker data is already in µm → set Lo_ref_um = 1.0 (no-op multiply).

close all; clc;
addpath(genpath('..'));

%% ══════════════════════════════════════════════════════════════════════════
%% Section 1 — Dataset descriptor array
%% ══════════════════════════════════════════════════════════════════════════

DSC = {};  % cell array — comment out / reorder any block freely

% % ── Baker 8mM  (slack only; datatable+velocitytable in .mat)
% d = struct();
% d.label      = 'Baker 8mM';
% d.matFile    = '../data/bakers_slack8mM_all.mat';
% d.Lo_ref_um  = 1.0;        % SL already in µm
% d.ATP        = 8;
% d.color      = [0.00 0.45 0.70];
% d.marker     = 'o';   % 8 mM
% d.lineStyle  = '-';
% d.fillMarker = true;
% d.data_slack  = [1.0 3.1];   % Baker time axis: slack window is ~1–3 s
% d.ktr_matFile = '../data/bakers_ktr_8.mat';
% d.ktr_tstart  = 0.0;   % s — recovery begins here
% DSC{end+1} = d;
% 
% % ── Baker 2mM  (slack from .mat; ktr from raw txt — no processed .mat exists)
% % Raw ktr txt: t in ms, L in Lo (×2.0 → µm), 4 header lines.
% % Restretch ends at ~96 ms raw → with −100 ms offset: t_ktr_start ≈ 0 s.
% d = struct();
% d.label          = 'Baker 2mM';
% d.matFile        = '../data/bakers_slack2mM.mat';
% d.Lo_ref_um      = 1.0;        % SL already in µm in the .mat
% d.ATP            = 2;
% d.color          = [0.00 0.45 0.70];
% d.marker         = 's';   % 2 mM
% d.lineStyle      = '--';
% d.fillMarker     = false;
% d.data_slack     = [1.0 3.1];
% d.ktr_txtFile    = '../data/2 mM ATP ktr.txt';
% d.ktr_header_lines = 4;   % ASI format: 4 text header lines
% d.ktr_t_scale    = 1e-3;  % ms → s
% d.ktr_t_offset   = -0.1;  % −100 ms offset applied in LoadBakersExp
% d.ktr_Lo_ref_um  = 2.0;   % raw txt L is in Lo; ML = 2.0 µm
% d.ktr_tstart     = 0.0;   % s — restretch ends near t = 0 after offset
% DSC{end+1} = d;

% ── 03/27 8mM
d = struct();
d.label      = '03/27 8mM';
d.dataFile   = '../data/03 27 2026 M/02_Merged_8mM_Active.txt';
d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
d.Lo_ref_um  = 2.0;
d.ATP        = 8;
d.color      = [0.84 0.10 0.11];
d.marker     = 'o';   % 8 mM
d.lineStyle  = '-';
d.fillMarker = true;
d.data_slack = [74.0 77.5];
d.vt_slack   = [74.4 Inf];
d.data_ktr   = [70.8 71.5];
d.data_step  = [69.2 70.5];
d.data_stair = [72.4 73.2];
DSC{end+1} = d;

% % ── 03/27 2mM
% d = struct();
% d.label      = '03/27 2mM';
% d.dataFile   = '../data/03 27 2026 M/03_Merged_2mM_Active.txt';
% d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
% d.Lo_ref_um  = 2.0;
% d.ATP        = 2;
% d.color      = [0.84 0.10 0.11];
% d.marker     = 's';   % 2 mM
% d.lineStyle  = '--';
% d.fillMarker = false;
% d.data_slack = [74.0 77.5];
% d.vt_slack   = [74.4 Inf];
% d.data_ktr   = [70.8 71.5];
% d.data_step  = [69.2 70.5];
% d.data_stair = [72.4 73.2];
% DSC{end+1} = d;

% ── 03/27 8mM repeat
d = struct();
d.label      = '03/27 8mM repeat';
d.dataFile   = '../data/03 27 2026 M/04_Merged_8mM_Active_repeat.txt';
d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
d.Lo_ref_um  = 2.0;
d.ATP        = 2;
d.color      = [0.10 0.84 0.11];
d.marker     = 'p';
d.lineStyle  = '-';
d.fillMarker = false;
d.data_slack = [74.0 77.5];
d.vt_slack   = [74.4 Inf];
d.data_ktr   = [70.8 71.5];
d.data_step  = [69.2 70.5];
d.data_stair = [72.4 73.2];
DSC{end+1} = d;


% ── 03/27 8mM repeat - rundown compensation Model 1
d = struct();
d.label      = '03/27 8mM repeat model 1';
d.dataFile   = '../data/03 27 2026 M/Merged_RunDownCorrection_Model1.txt';
d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
d.Lo_ref_um  = 2.0;
d.ATP        = 2;
d.color      = [0.00 0.45 0.70];
d.marker     = '<';
d.lineStyle  = '-';
d.fillMarker = false;
d.data_slack = [74.0 77.5];
d.vt_slack   = [74.4 Inf];
d.data_ktr   = [70.8 71.5];
d.data_step  = [69.2 70.5];
d.data_stair = [72.4 73.2];
DSC{end+1} = d;

% ── 03/27 8mM repeat - rundown compensation Model 2
d = struct();
d.label      = '03/27 8mM repeat model 2';
d.dataFile   = '../data/03 27 2026 M/Merged_RunDownCorrection_Model2.txt';
d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
d.Lo_ref_um  = 2.0;
d.ATP        = 2;
d.color      = [0.47 0.67 0.19];
d.marker     = '^';
d.lineStyle  = '-';
d.fillMarker = false;
d.data_slack = [74.0 77.5];
d.vt_slack   = [74.4 Inf];
d.data_ktr   = [70.8 71.5];
d.data_step  = [69.2 70.5];
d.data_stair = [72.4 73.2];
DSC{end+1} = d;

% ── 03/27 8mM repeat - rundown compensation Model 3
d = struct();
d.label      = '03/27 8mM repeat model 3';
d.dataFile   = '../data/03 27 2026 M/Merged_RunDownCorrection_Model3.txt';
d.vtFile     = '../data/protocol_03_27_2026_velocitytable.mat';
d.Lo_ref_um  = 2.0;
d.ATP        = 2;
d.color      = [0.19 0.47 0.67];
d.marker     = '>';
d.lineStyle  = '-';
d.fillMarker = false;
d.data_slack = [74.0 77.5];
d.vt_slack   = [74.4 Inf];
d.data_ktr   = [70.8 71.5];
d.data_step  = [69.2 70.5];
d.data_stair = [72.4 73.2];
DSC{end+1} = d;

% % ── 04/03 8mM  (preamble slack before 74.8 → vt_slack offset)
% d = struct();
% d.label      = '04/03 8mM';
% d.dataFile   = '../data/04 03 2026 F/02_Merged_8mM_Active.txt';
% d.vtFile     = '../data/protocol_04_03_2026_velocitytable.mat';
% d.Lo_ref_um  = 2.0;
% d.ATP        = 8;
% d.color      = [0.47 0.67 0.19];
% d.marker     = 'o';   % 8 mM
% d.lineStyle  = '-';
% d.fillMarker = true;
% d.data_slack = [74.0 78.0];
% d.vt_slack   = [74.8 Inf];
% d.data_ktr   = [70.8 71.5];
% d.data_step  = [69.2 70.5];
% d.data_stair = [72.4 73.2];
% DSC{end+1} = d;
% 
% % ── 04/03 2mM
% d = struct();
% d.label      = '04/03 2mM';
% d.dataFile   = '../data/04 03 2026 F/03_Merged_2mM_Active.txt';
% d.vtFile     = '../data/protocol_04_03_2026_velocitytable.mat';
% d.Lo_ref_um  = 2.0;
% d.ATP        = 2;
% d.color      = [0.47 0.67 0.19];
% d.marker     = 's';   % 2 mM
% d.lineStyle  = '--';
% d.fillMarker = false;
% d.data_slack = [74.0 78.0];
% d.vt_slack   = [74.8 Inf];
% d.data_ktr   = [70.8 71.5];
% d.data_step  = [69.2 70.5];
% d.data_stair = [72.4 73.2];
% DSC{end+1} = d;

% Unify fields across all structs so we can build a struct array.
% Each DSC entry may have different fields — collect the superset,
% fill missing ones with [], then concatenate.
allFields = {};
for k = 1:numel(DSC)
    allFields = union(allFields, fieldnames(DSC{k})', 'stable');
end
allFields = allFields(:);  % column vector — matches fieldnames() orientation
for k = 1:numel(DSC)
    for fi = 1:numel(allFields)
        if ~isfield(DSC{k}, allFields{fi}); DSC{k}.(allFields{fi}) = []; end
    end
    DSC{k} = orderfields(DSC{k}, allFields);
end
DS = [DSC{:}];

%% ══════════════════════════════════════════════════════════════════════════
%% Section 3 — Slack feature extraction
%% ══════════════════════════════════════════════════════════════════════════

for di = 1:numel(DS)
    if ~isfield(DS(di), 'data_slack') || isempty(DS(di).data_slack); continue; end

    fprintf('=== Slack: %s ===\n', DS(di).label);
    [t, F, L_um, vt_full] = loadDataset(DS(di));

    % Subset data to slack window
    win = t >= DS(di).data_slack(1) & t <= DS(di).data_slack(2);

    % Subset vt to the slack velocity-table range
    if isfield(DS(di), 'vt_slack') && ~isempty(DS(di).vt_slack)
        vt = vt_full(vt_full(:,1) >= DS(di).vt_slack(1) & ...
                     vt_full(:,1) <= DS(di).vt_slack(2), :);
    else
        vt = vt_full;
    end

    figure(400 + di); clf;
    set(gcf, 'Name', [DS(di).label ' – slack fits'], ...
        'Units', 'centimeters', 'Position', [1 1 20 12]);
    DS(di).slackFeats = extractSlackAttributes( ...
        t(win), F(win), L_um(win), vt, struct(), [], true, true);  % useSmoothing=true for experimental data
    rmsed = DS(di).slackFeats.ktr_rmse - DS(di).slackFeats.ktr2_rmse;
    if any(rmsed < 0)
        warning('Single exponential RMSE is better!', di);
    end
    title([DS(di).label ' – slack fits ']);
    fprintf('RMSE is better by %f\n' , sum(rmsed));
end

%% ══════════════════════════════════════════════════════════════════════════
%% Section 4 — KTR extraction
%% ══════════════════════════════════════════════════════════════════════════
figure(420); clf;hold on;
for di = 1:numel(DS)
    has_ktr_mat = isfield(DS(di), 'ktr_matFile') && ~isempty(DS(di).ktr_matFile);
    has_ktr_txt = isfield(DS(di), 'ktr_txtFile') && ~isempty(DS(di).ktr_txtFile);
    has_ktr_win = isfield(DS(di), 'data_ktr')    && ~isempty(DS(di).data_ktr);
    if ~has_ktr_mat && ~has_ktr_txt && ~has_ktr_win; continue; end

    fprintf('=== KTR: %s ===\n', DS(di).label);

    if has_ktr_mat
        % Baker 8mM-style: processed .mat, L already in µm, t already in s.
        s = load(DS(di).ktr_matFile, 'datatable');
        M = s.datatable;
        datatable_fr = [M(:,1), M(:,2) * DS(di).Lo_ref_um, M(:,3)];
        t_ktr_start  = DS(di).ktr_tstart;
    elseif has_ktr_txt
        % Baker 2mM-style: raw ASI txt, t in ms, L in Lo. Convert and offset.
        M = readmatrix(DS(di).ktr_txtFile, 'NumHeaderLines', DS(di).ktr_header_lines);
        t_s  = M(:,1) * DS(di).ktr_t_scale + DS(di).ktr_t_offset;
        L_um = M(:,2) * DS(di).ktr_Lo_ref_um;
        datatable_fr = [t_s, L_um, M(:,3)];
        t_ktr_start  = DS(di).ktr_tstart;
    else
        % Protocol file: find restretch end via velocity table.
        [t, F, L_um, vt_full] = loadDataset(DS(di));
        datatable_fr = [t, L_um, F];
        rng = DS(di).data_ktr;
        % The restretch end is the last big positive velocity segment before rng(2)
        ktr_segs = findVtSegments(vt_full, rng(1), rng(2), +1, 10, 0);
        if isempty(ktr_segs)
            fprintf('  [WARN] no ktr restretch found in vt for %s\n', DS(di).label);
            continue;
        end
        t_ktr_start = ktr_segs(end, 2);  % end of restretch = start of recovery
    end

    zone_ms = [(t_ktr_start * 1000) + 10, (t_ktr_start * 1000) + 500];

    
    nexttile();
    % subplot(1, numel(DS), di); hold on;
    title(DS(di).label, 'Interpreter', 'none');
    try
        [~, ktr_val, df_val, ~, E_val, SL_val, ~, ktr2] = fitRecovery(datatable_fr, zone_ms, 0, [], true);
        xlim([t_ktr_start - 0.05 t_ktr_start + 0.3]);
        DS(di).ktrFeats.ktr  = ktr_val;
        DS(di).ktrFeats.SL  = SL_val;

        DS(di).ktrFeats.df   = df_val;
        DS(di).ktrFeats.rmse = E_val;

DS(di).ktrFeats.k     = ktr2.k    ;
DS(di).ktrFeats.df    = ktr2.df   ;
DS(di).ktrFeats.B     = ktr2.B    ;
DS(di).ktrFeats.gam   = ktr2.gam  ;
DS(di).ktrFeats.omega = ktr2.omega;
DS(di).ktrFeats.rmse  = ktr2.rmse ;

        fprintf('  ktr = %.1f s⁻¹, df = %.1f kPa, rmse = %.2f\n', ktr_val, df_val, E_val);
    catch ME
        fprintf('  [WARN] fit failed: %s\n', ME.message);
        DS(di).ktrFeats.ktr  = NaN;
        DS(di).ktrFeats.SL   = NaN;
        DS(di).ktrFeats.df   = NaN;
        DS(di).ktrFeats.rmse = NaN;
    end
end
if any(arrayfun(@(d) isfield(d,'ktrFeats'), DS))
    figure(420); sgtitle('KTR Fits');
end

%% ══════════════════════════════════════════════════════════════════════════
%% Section 5 — Step perturbation extraction
%% ══════════════════════════════════════════════════════════════════════════

% ── Step perturbation model: decoupled 2-exponential ─────────────────────
%   F(t) = Fss + A1*exp(-k1*t) + A2*exp(-k2*t)
%
%   Rationale: the classic coupled form  Fss + A*(r*exp(-k1*t) - (1-r)*exp(-k2*t))
%   suffers a strong k1-k2 tradeoff — the optimiser can trade rate against
%   amplitude freely.  Decoupling amplitudes and bounding each from the data
%   fixes this:
%
%     A1  = "decay"  amplitude  → bounded by  (max(F) – Fss_est)  for step-up
%                                             (min(F) – Fss_est)  for step-down
%     A2  = "undershoot" amplitude (opposite sign to A1 if any overshoot present)
%                                → bounded by  (min(F) – Fss_est)  for step-up
%                                             (max(F) – Fss_est)  for step-down
%     k1  = fast rate  (k1 > k2 enforced by lower/upper bounds)
%     k2  = slow rate
%
%   Both rates are constrained independently of the amplitudes, eliminating
%   the k1/k2 swap ambiguity present in the r-parameterisation.
mdl_step = fittype('Fss + A1*exp(-k1*x) + A2*exp(-k2*x)', ...
    'coeff', {'Fss','A1','A2','k1','k2'}, 'independent', 'x');

for di = 1:numel(DS)
    if ~isfield(DS(di), 'data_step') || isempty(DS(di).data_step); continue; end

    fprintf('=== Step perturbation: %s ===\n', DS(di).label);
    [t, F, L_um, vt_full] = loadDataset(DS(di));

    rng      = DS(di).data_step;
    step_ups = findVtSegments(vt_full, rng(1), rng(2), +1, 0.5, 1);
    step_dns = findVtSegments(vt_full, rng(1), rng(2), -1, 0.5, 1);
    nPairs   = min(size(step_ups,1), size(step_dns,1));
    fprintf('  Found %d step pairs\n', nPairs);

    sf.label   = DS(di).label;
    sf.k1_up   = nan(1, nPairs);  sf.k2_up   = nan(1, nPairs);
    sf.k1_dwn  = nan(1, nPairs);  sf.k2_dwn  = nan(1, nPairs);
    sf.A1_up   = nan(1, nPairs);  sf.A2_up   = nan(1, nPairs);
    sf.A1_dwn  = nan(1, nPairs);  sf.A2_dwn  = nan(1, nPairs);
    sf.Fss_up  = nan(1, nPairs);  sf.Fss_dwn = nan(1, nPairs);
    sf.dL      = nan(1, nPairs);

    clr = DS(di).color;
    fig_step = figure(430 + di); clf;
    set(fig_step, 'Name', [DS(di).label ' – step fits'], ...
        'Units', 'centimeters', 'Position', [1 1 24 max(nPairs*5, 8)]);
    tl_step = tiledlayout(nPairs, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_step, [DS(di).label ' — step perturbation fits'], 'Interpreter', 'none');

    t_maxFrame = 0.1;  % 100 ms fit window
    t_fit = linspace(0, t_maxFrame*2, 300)';
    for s = 1:nPairs
        msk_up  = t > step_ups(s,1) & t < step_ups(s,2) & t < step_ups(s,1) + t_maxFrame;
        msk_dwn = t > step_dns(s,1) & t < step_dns(s,2) & t < step_dns(s,1) + t_maxFrame;

        L_before = L_um(find(t <= step_ups(s,1), 1, 'last'));
        L_after  = L_um(find(t >= step_ups(s,2), 1));
        sf.dL(s) = (L_after - L_before) / DS(di).Lo_ref_um;  % store in Lo

        % ── Step-up fit ───────────────────────────────────────────────────
        % A1 > 0: fast decay from peak down to Fss.  Bounded by (max-Fss_est).
        % A2 ≤ 0: slow undershoot below Fss (if present). Bounded by (min-Fss_est).
        t0_up   = t(find(msk_up, 1));
        tf_up   = t(msk_up) - t0_up;  Ff_up = F(msk_up);
        Fss_est = median(Ff_up(max(1, end-round(0.2*numel(Ff_up))):end));
        A1_est  = max(Ff_up) - Fss_est;          % decay amplitude (> 0)
        A2_est  = min(Ff_up) - Fss_est;          % undershoot amplitude (≤ 0)
        fu = [];
        try
            fu = fit(tf_up, Ff_up, mdl_step, ...
                'StartPoint', [Fss_est,  A1_est,          A2_est,         300, 30], ...
                'Lower',      [Fss_est*0.7, 0,            min(A2_est*1.5, -1e-3), 50,  1], ...
                'Upper',      [Fss_est*1.3, A1_est*2,     0,              3000, 200]);
            sf.Fss_up(s) = fu.Fss; sf.k1_up(s) = fu.k1; sf.k2_up(s) = fu.k2;
            sf.A1_up(s) = fu.A1;   sf.A2_up(s) = fu.A2;
        catch ME
            fprintf('  [WARN] step-up s=%d: %s\n', s, ME.message);
        end
        ax = nexttile(tl_step); hold(ax, 'on');
        plot(ax, tf_up, Ff_up, 'k.', 'MarkerSize', 4);
        if ~isempty(fu)
            plot(ax, t_fit, fu(t_fit), '-', 'Color', clr, 'LineWidth', 1.5);
            title(ax, sprintf('Up %d  k1=%.0f k2=%.0f', s, fu.k1, fu.k2), 'Interpreter', 'none');
        else
            title(ax, sprintf('Up %d  (failed)', s), 'Interpreter', 'none');
        end
        xlabel(ax, 't (s)'); ylabel(ax, 'F (kPa)'); box(ax, 'on');

        % ── Step-down fit ─────────────────────────────────────────────────
        % A1 < 0: fast drop below Fss.  Bounded by (min-Fss_est).
        % A2 ≥ 0: slow recovery/overshoot above Fss (if present). Bounded by (max-Fss_est).
        t0_dwn  = t(find(msk_dwn, 1));
        tf_dwn  = t(msk_dwn) - t0_dwn;  Ff_dwn = F(msk_dwn);
        Fss_est = median(Ff_dwn(max(1, end-round(0.2*numel(Ff_dwn))):end));
        A1_est  = min(Ff_dwn) - Fss_est;         % drop amplitude (< 0)
        A2_est  = max(Ff_dwn) - Fss_est;         % overshoot amplitude (≥ 0)
        fd = [];
        try
            fd = fit(tf_dwn, Ff_dwn, mdl_step, ...
                'StartPoint', [Fss_est,  A1_est,          A2_est,         300, 30], ...
                'Lower',      [Fss_est*0.7, A1_est*1.5,   0,              50,  1], ...
                'Upper',      [Fss_est*1.3, 0,             max(A2_est*1.5, 1e-3), 3000, 200]);
            sf.Fss_dwn(s) = fd.Fss; sf.k1_dwn(s) = fd.k1; sf.k2_dwn(s) = fd.k2;
            sf.A1_dwn(s) = fd.A1;   sf.A2_dwn(s) = fd.A2;
        catch ME
            fprintf('  [WARN] step-down s=%d: %s\n', s, ME.message);
        end
        ax = nexttile(tl_step); hold(ax, 'on');
        plot(ax, tf_dwn, Ff_dwn, 'k.', 'MarkerSize', 4);
        if ~isempty(fd)
            plot(ax, t_fit, fd(t_fit), '-', 'Color', clr, 'LineWidth', 1.5);
            title(ax, sprintf('Down %d  k1=%.0f k2=%.0f', s, fd.k1, fd.k2), 'Interpreter', 'none');
        else
            title(ax, sprintf('Down %d  (failed)', s), 'Interpreter', 'none');
        end
        xlabel(ax, 't (s)'); ylabel(ax, 'F (kPa)'); box(ax, 'on');
    end
    DS(di).stepFeats = sf;
end

%% ══════════════════════════════════════════════════════════════════════════
%% Section 6 — Staircase extraction
%% ══════════════════════════════════════════════════════════════════════════

% ── Staircase model: decoupled 2-exponential (same rationale as step) ────
%   F(t) = Fss + A1*exp(-k1*t) + A2*exp(-k2*t)
%
%   After each shortening step force is at a minimum and recovers to Fss.
%   A1 < 0: fast recovery component (k1 large, dominant).
%   A2 ≥ 0: slow overshoot component (k2 small), may be ~0 for monotone data.
%   Bounds are data-derived (min/max of the recovery window) to prevent
%   the k1/k2 tradeoff seen in the r-parameterisation.
%
%   TODO: staircase 2-exp fits are unreliable when the recovery window is
%   short (< 50 ms) or when the signal has significant noise.  Possible
%   improvements:
%     - Use single-exp fallback when A2 is near zero
%     - Jointly fit across all staircase steps (shared k1, k2; free Fss, A)
%     - Increase window duration in vt_stair descriptor
mdl_stair = fittype('Fss + A1*exp(-k1*x) + A2*exp(-k2*x)', ...
    'coeff', {'Fss','A1','A2','k1','k2'}, 'independent', 'x');

for di = 1:numel(DS)
    if ~isfield(DS(di), 'data_stair') || isempty(DS(di).data_stair); continue; end

    fprintf('=== Staircase: %s ===\n', DS(di).label);
    [t, F, L_um, vt_full] = loadDataset(DS(di));

    rng       = DS(di).data_stair;
    stair_ups = findVtSegments(vt_full, rng(1), rng(2), +1, 0.3, 1);
    nStairs   = size(stair_ups, 1);
    fprintf('  Found %d staircase steps\n', nStairs);

    dt = median(diff(t));
    fs = 1 / dt;
    isHighRes = fs > 1000;

    sf.label   = DS(di).label;
    sf.Fss     = nan(1, nStairs);
    sf.k1      = nan(1, nStairs);
    sf.k2      = nan(1, nStairs);
    sf.A1      = nan(1, nStairs);
    sf.A2      = nan(1, nStairs);
    sf.dL      = nan(1, nStairs);
    sf.L_after = nan(1, nStairs);

    clr = DS(di).color;
    ncols_st = ceil(sqrt(nStairs));
    nrows_st = ceil(nStairs / ncols_st);
    fig_stair = figure(440 + di); clf;
    set(fig_stair, 'Name', [DS(di).label ' – staircase fits'], ...
        'Units', 'centimeters', 'Position', [1 1 ncols_st*8 nrows_st*6]);
    tl_stair = tiledlayout(nrows_st, ncols_st, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_stair, [DS(di).label ' — staircase fits'], 'Interpreter', 'none');

    for s = 1:nStairs
        t_start = stair_ups(s, 1);
        t_end   = stair_ups(s, 2);
        if s < nStairs
            t_next = stair_ups(s+1, 1);
        else
            t_next = t_start + 0.15;
        end

        L_before   = L_um(find(t <= t_start, 1, 'last'));
        L_after_v  = mean(L_um(t >= t_start+0.025 & t <= t_start+0.035));
        sf.dL(s)      = (L_after_v - L_before) / DS(di).Lo_ref_um;
        sf.L_after(s) = L_after_v / DS(di).Lo_ref_um;

        mask_raw = t >= t_end & t <= t_next - 0.005;

        ax = nexttile(tl_stair); hold(ax, 'on');

        if sum(mask_raw) < 10
            title(ax, sprintf('Stair %d  (no data)', s), 'Interpreter', 'none');
            continue;
        end

        if ~isHighRes
            sf.Fss(s) = mean(F(mask_raw));
            plot(ax, t(mask_raw) - t_end, F(mask_raw), 'k.', 'MarkerSize', 4);
            title(ax, sprintf('Stair %d  Fss=%.1f', s, sf.Fss(s)), 'Interpreter', 'none');
            xlabel(ax, 't (s)'); ylabel(ax, 'F (kPa)'); box(ax, 'on');
            continue;
        end

        % Downsample to ~1000 Hz
        ds_factor = max(1, round(fs / 1000));
        t_raw = t(mask_raw);  F_raw = F(mask_raw);
        n_bins = floor(numel(t_raw) / ds_factor);
        t_ds = arrayfun(@(k) mean(t_raw((k-1)*ds_factor+1:k*ds_factor)), 1:n_bins)';
        F_ds = arrayfun(@(k) mean(F_raw((k-1)*ds_factor+1:k*ds_factor)), 1:n_bins)';
        tf = t_ds - t_ds(1);
        % Data-derived amplitude bounds:
        % A1 < 0: fast recovery component (force rises from min → Fss, so A1 is negative)
        % A2 ≥ 0: slow overshoot component (may be ~0 for monotone recovery)
        Fss_est = mean(F_ds(max(1,end-round(0.2*n_bins)):end));
        A1_est  = min(F_ds) - Fss_est;   % initial drop below Fss (< 0)
        A2_est  = max(F_ds) - Fss_est;   % overshoot above Fss (≥ 0, often ~0)

        plot(ax, tf, F_ds, 'k.', 'MarkerSize', 4);
        f0 = [];
        try
            f0 = fit(tf, F_ds, mdl_stair, ...
                'StartPoint', [Fss_est,  A1_est,         A2_est,      100, 10], ...
                'Lower',      [Fss_est*0.7, A1_est*1.5,  0,            20,  1], ...
                'Upper',      [Fss_est*1.1, 0,           max(A2_est*2, 1), 800, 50]);
            sf.Fss(s) = f0.Fss; sf.k1(s) = f0.k1; sf.k2(s) = f0.k2;
            sf.A1(s) = f0.A1;   sf.A2(s) = f0.A2;
            fprintf('  Stair %d: Fss=%.1f k1=%.0f k2=%.0f\n', s, f0.Fss, f0.k1, f0.k2);
        catch ME
            fprintf('  [WARN] stair s=%d: %s\n', s, ME.message);
        end
        if ~isempty(f0)
            t_fit_st = linspace(0, tf(end), 300)';
            plot(ax, t_fit_st, f0(t_fit_st), '-', 'Color', clr, 'LineWidth', 1.5);
            title(ax, sprintf('Stair %d  k1=%.0f k2=%.0f', s, f0.k1, f0.k2), 'Interpreter', 'none');
        else
            title(ax, sprintf('Stair %d  (failed)', s), 'Interpreter', 'none');
        end
        xlabel(ax, 't (s)'); ylabel(ax, 'F (kPa)'); box(ax, 'on');
    end
    DS(di).stairFeats = sf;
end

%% ══════════════════════════════════════════════════════════════════════════
%% Section 7 — Universal comparison figures
%% ══════════════════════════════════════════════════════════════════════════

% Slack features
[fc, lbls, clrs, mkrs, lsts, fills] = collectDS(DS, 'slackFeats');
if ~isempty(fc)
    % fn_slack = fieldnames(fc{1});
    % customized
    fn_slack = {'ktr|SLslack', 'ktr2_k|SLslack','ktr2_A|SLslack', 'A|SLslack', 'ktr2_omega', 't0', 'restretchSlopeStart', 'restretchSlopeEnd', 'peak1_y', 'peak1_t', ...
        'peak1_dSL', 'vall_y', 'vall_t', 'peak2', 'steady', 'vall2_dy', 'vall2_t', 'ovrsht_dy', 'ovrsht_t'}
    % Dropped: 'SLslack', 'SLdiff', 'v_restretch', 'peak1_SL', 'XTOR', 'A|Am' 
    figure(403); clf;
    set(gcf, 'Name', 'Slack feature comparison', 'Units', 'centimeters', 'Position', [2 2 32 24]);
    title_str = ['Slack features — ' strjoin(lbls, ' vs ')];
    sgtitle(title_str, 'Interpreter', 'none');
    plotMultipleFeatures(fc, lbls, clrs, mkrs, fn_slack, lsts, fills);
end

%% KTR — overlay ktr-experiment points onto tile 1 of fig 403 (ktr vs SLslack)
    % Each dataset contributes one point: ktr from the restretch experiment
    % at its measured SL.  Use '*' marker (distinct from slack-series markers)
    % with the same per-dataset color so the source is visually clear.
[fc, lbls, clrs, mkrs, lsts, fills] = collectDS(DS, 'slackFeats');
[fc_ktr, lbls_ktr, clrs_ktr, mkrs_ktr, lsts_ktr, fills_ktr] = collectDS(DS, 'ktrFeats');
if ~isempty(fc) && ~isempty(fc_ktr)

    fn_slack = {'ktr|SLslack', };
    figure(420); clf;
    set(gcf, 'Name', 'Slack + ktr feature comparison', 'Units', 'centimeters', 'Position', [2 2 32 24]);
    title_str = ['Slack + ktr features — ' strjoin(lbls, ' vs ')];
    sgtitle(title_str, 'Interpreter', 'none');
    plotMultipleFeatures(fc, lbls, clrs, mkrs, fn_slack, lsts, fills);


    for di = 1:numel(fc_ktr)
        kf = fc_ktr{di};
        if isnan(kf.ktr); continue; end
        nexttile(1); hold on;
        plot(kf.SL, kf.ktr, mkrs{di}, 'Color', clrs_ktr(di, :), 'MarkerFaceColor', clrs_ktr(di, :),'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']);
        
        nexttile(2);hold on;title('k with osc (filled) vs ktr (open)')
        plot(di, kf.k, mkrs{di}, 'Color', clrs_ktr(di, :), 'MarkerFaceColor', clrs_ktr(di, :),'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']); 
        plot(di, kf.ktr, mkrs{di}, 'Color', clrs_ktr(di, :) ,'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']); 
                        
        nexttile(3);hold on;title('B');
        plot(di, kf.B, mkrs{di}, 'Color', clrs_ktr(di, :), 'MarkerFaceColor', clrs_ktr(di, :),'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']);

        nexttile(4);hold on;title('gamma')
        plot(di, kf.gam, mkrs{di}, 'Color', clrs_ktr(di, :), 'MarkerFaceColor', clrs_ktr(di, :),'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']);
        
        nexttile(5);hold on;title('omega')
        plot(di, kf.omega, mkrs{di}, 'Color', clrs_ktr(di, :), 'MarkerFaceColor', clrs_ktr(di, :),'MarkerSize', 8, ...
            'LineWidth', 2, 'DisplayName', [lbls_ktr{di} ' (ktr exp)']);
        
    end
end

%% Step perturbation
[fc, lbls, clrs, mkrs, lsts, fills] = collectDS(DS, 'stepFeats');
if ~isempty(fc)
    fn_step = {'k1_up','k2_up','k1_dwn','k2_dwn','A1_up','A2_up','A1_dwn','A2_dwn','dL'};
    figure(430); clf;
    set(gcf, 'Name', 'Step perturbation comparison', 'Units', 'centimeters', 'Position', [2 2 28 18]);
    sgtitle('Step perturbation features', 'Interpreter', 'none');
    plotMultipleFeatures(fc, lbls, clrs, mkrs, fn_step, lsts, fills);
end

%% Staircase
[fc, lbls, clrs, mkrs, lsts, fills] = collectDS(DS, 'stairFeats');
if ~isempty(fc)
    fn_stair = {'k1|L_after','k2|L_after','Fss|L_after','A1|L_after','A2|L_after'};
    % ,'dL','L_after'
    figure(440); clf;
    set(gcf, 'Name', 'Staircase comparison', 'Units', 'centimeters', 'Position', [2 2 28 18]);
    sgtitle('Staircase features', 'Interpreter', 'none');
    plotMultipleFeatures(fc, lbls, clrs, mkrs, fn_stair, lsts, fills);
end

%% ══════════════════════════════════════════════════════════════════════════
%% Section 8 — Universal print table
%% ══════════════════════════════════════════════════════════════════════════

% Explicit field whitelists per feature type (avoids struct-array
% cross-contamination where sharing field names leaks values across entries)
feat_types  = {'slackFeats', 'ktrFeats', 'stepFeats', 'stairFeats'};
feat_names  = {'Slack', 'KTR', 'Step perturbation', 'Staircase'};
feat_fields = {
    {},   ...  % slack: auto-discover (no contamination issue)
    {'ktr','df','rmse'}, ...
    {'k1_up','k2_up','k1_dwn','k2_dwn','A1_up','A2_up','A1_dwn','A2_dwn','Fss_up','Fss_dwn','dL'}, ...
    {'k1','k2','Fss','A1','A2','dL','L_after'} ...
};

for ft = 1:numel(feat_types)
    fname = feat_types{ft};
    have  = find(arrayfun(@(d) isfield(d, fname) && ~isempty(d.(fname)), DS));
    if isempty(have); continue; end

    if ~isempty(feat_fields{ft})
        all_fn = feat_fields{ft};
    else
        % Auto-discover for slack (struct not shared across feat types)
        all_fn = {};
        for di = have
            fs = DS(di).(fname);
            if ~isstruct(fs); continue; end
            fn_i = fieldnames(fs);
            for kk = 1:numel(fn_i)
                v = fs.(fn_i{kk});
                if isnumeric(v) && ~isempty(v)
                    all_fn = union(all_fn, fn_i(kk), 'stable');
                end
            end
        end
        all_fn = all_fn(~strcmp(all_fn, 'label'));
    end

    COL = 30;
    sep = repmat('-', 1, 22 + (COL+2)*numel(have));
    fprintf('\n── %s ──\n%s\n', feat_names{ft}, sep);
    fprintf('%-22s', 'Field');
    for di = have; fprintf('  %-*s', COL, DS(di).label); end
    fprintf('\n%s\n', sep);

    for k = 1:numel(all_fn)
        fn_k = all_fn{k};
        fprintf('%-22s', fn_k);
        for di = have
            fs = DS(di).(fname);
            if ~isfield(fs, fn_k)
                fprintf('  %-*s', COL, '—');
            else
                v = fs.(fn_k);
                if ischar(v) || isstring(v)
                    fprintf('  %-*s', COL, v);
                elseif isnumeric(v) && ~isempty(v)
                    s = mat2str(round(v, 4, 'significant'));
                    fprintf('  %-*s', COL, s);
                else
                    fprintf('  %-*s', COL, '[]');
                end
            end
        end
        fprintf('\n');
    end
    fprintf('%s\n', sep);
end

%% ══════════════════════════════════════════════════════════════════════════
%% Local functions
%% ══════════════════════════════════════════════════════════════════════════

function [fc, lbls, clrs, mkrs, lsts, fills] = collectDS(DS, fieldname)
%COLLECTDS  Return cell arrays for datasets that have the given field.
    idx   = find(arrayfun(@(d) isfield(d, fieldname) && ~isempty(d.(fieldname)), DS));
    fc    = {DS(idx).(fieldname)};
    lbls  = {DS(idx).label};
    clrs  = vertcat(DS(idx).color);
    mkrs  = {DS(idx).marker};
    lsts  = {DS(idx).lineStyle};
    fills = [DS(idx).fillMarker];
end

function [t, F, L_um, vt] = loadDataset(d)
%LOADDATASET  Load time/force/length and velocity table for one descriptor.
%   Applies Lo_ref_um conversion and removes duplicate near-zero-neg VT rows.

    if isfield(d, 'matFile') && ~isempty(d.matFile)
        % Baker-style .mat with datatable [t, SL, F] and velocitytable
        s = load(d.matFile, 'datatable', 'velocitytable');
        M  = s.datatable;       % [t, SL(already in um or Lo), F]
        vt_raw = s.velocitytable;
    else
        % Text file: cols [t, L_Lo, F]
        M = readmatrix(d.dataFile);
        % Remove non-monotonic timestamps
        mono = [true; diff(M(:,1)) > 0];
        M = M(mono, :);
        s_vt = load(d.vtFile, 'velocitytable');
        vt_raw = s_vt.velocitytable;
    end

    t    = M(:, 1);
    L_um = M(:, 2) * d.Lo_ref_um;   % µm (Baker: ×1, new data: ×2.0)
    F    = M(:, 3);

    % VT: keep only cols [t, vel]; remove second-of-consecutive near-zero-neg
    vt2       = vt_raw(:, 1:2);
    is_nzn    = vt2(:,2) < 0 & vt2(:,2) > -0.1;
    noise_row = is_nzn & [false; is_nzn(1:end-1)];
    vt        = vt2(~noise_row, :);
end

function segs = findVtSegments(vt, tmin, tmax, velSign, velThresh, shift)
%FINDVTSEGMENTS  Return [start, end] time pairs for velocity-table segments.
%   Finds contiguous rows in vt where vt(:,1) ∈ [tmin,tmax] and
%   vt(:,2)*velSign > velThresh. Retuinrns Nx2 matrix.

    mask = vt(:,1) >= tmin & vt(:,1) <= tmax & (vt(:,2) * velSign) > velThresh;
    idx  = find(mask) + shift;
    if isempty(idx)
        segs = zeros(0, 2);
        return;
    end
    breaks = [0; find(diff(idx) > 1); numel(idx)];
    segs   = zeros(numel(breaks)-1, 2);
    for k = 1:numel(breaks)-1
        gi = idx(breaks(k)+1);
        ge = idx(breaks(k+1));
        ge_next = min(ge + 1, size(vt, 1));
        segs(k, :) = [vt(gi, 1), vt(ge_next, 1)];
    end
end
