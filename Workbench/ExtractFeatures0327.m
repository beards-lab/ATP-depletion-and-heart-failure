%% Extract Features from 03/27/2026 Merged Data
close all; clear; clc;
addpath(genpath('..'));

dataDir = '../data/03 27 2026 M/';
SL0 = 2.2; % sarcomere length at Lo (um)

% Load velocity table (from CreateProtocolVelocityTable.m)
load('../data/protocol_03_27_2026_velocitytable.mat', 'velocitytable');
vt = velocitytable; % [time(s), vel(Lo/s), vel_um, L(Lo)]

%% Section 1: Load Data
files = dir([dataDir '*Merged*.txt']);
[~, idx] = sort({files.name});
files = files(idx);

n = length(files);
D = struct('t', {}, 'L', {}, 'F', {}, 'label', {}, 'isActive', {});
for i = 1:n
    M = readmatrix(fullfile(files(i).folder, files(i).name));
    % Remove non-monotonic timestamps (burst/Log boundary artifacts)
    mono = [true; diff(M(:,1)) > 0];
    M = M(mono, :);
    D(i).t = M(:,1);
    D(i).L = M(:,2);
    D(i).F = M(:,3);
    D(i).label = strrep(files(i).name(4:end-4), '_', ' ');
    w = D(i).t >= 69 & D(i).t <= 70;
    D(i).isActive = mean(D(i).F(w)) > 10;
    fprintf('Loaded %s (active=%d, meanF=%.1f kPa)\n', files(i).name, D(i).isActive, mean(D(i).F(w)));
end

colors = lines(n);
activeIdx = find([D.isActive]);

%% fix the steps VT
vt = velocitytable;
% drop lines
dl = [11, 14, 17];
keep = setdiff(1:length(vt)', dl');
vt = vt(keep, :);
% vt = vt()
plot(vt(:, 1), vt(:, 4), 'r-o')

%% Helper: find velocity table segments by time and velocity sign
% Returns Nx2 matrix of [seg_start_time, seg_end_time] for segments
% where velocity exceeds threshold in the given time window
findSegs = @(tmin, tmax, velSign, velThresh) ...
    findVtSegments(vt, tmin, tmax, velSign, velThresh);

%% Section 2: Step Perturbation Decay/Recovery Fitting (69.2 - 70.5s)
fprintf('\n=== Step Perturbation Analysis ===\n');
t_maxFrame = 1; % 100 ms timewindows for the fit
% From velocity table, identify step-up and step-down segments
step_ups = findVtSegments(vt, 69.2, 70.5, +1, 0.5);    % positive velocity > 0.5 Lo/s
step_downs = findVtSegments(vt, 69.2, 70.5, -1, 0.5);   % negative velocity
nStepPairs = min(size(step_ups,1), size(step_downs,1));
fprintf('Found %d step-up and %d step-down segments from velocity table\n', ...
    size(step_ups,1), size(step_downs,1));

mdl = fittype('Fss + A*exp(-k*x) +b*x', 'coeff',{'Fss','A','k', 'b'},'independent', 'x');

stepResults = struct();
figure(101);clf;hold on; 
ax1 = nexttile(1);hold on;xlim([69.4, 70.6])
ax2 = nexttile(2);hold on;

for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;

    stepResults(ci).label = D(ai).label;
    stepResults(ci).k_up = nan(1, nStepPairs);
    stepResults(ci).k_down = nan(1, nStepPairs);
    stepResults(ci).A_up = nan(1, nStepPairs);
    stepResults(ci).A_down = nan(1, nStepPairs);
    stepResults(ci).dL = nan(1, nStepPairs);

        
plot(ax1, t, F);
        
    for s = 1:nStepPairs
        %% Step-up: plateau starts at step_ups(s,2), ends at step_downs(s,1)
        t_up_end = step_ups(s, 2);
        t_dn_start = step_downs(s, 1);
        L_before = L(find(t <= step_ups(s,1), 1, 'last'));
        L_after  = L(find(t >= t_up_end, 1));
        stepResults(ci).dL(s) = L_after - L_before;
        msk_up = t > step_ups(s, 1) & t < step_ups(s, 2) & t < step_ups(s, 1) + t_maxFrame;
        msk_dwn = t > step_downs(s, 1) & t < step_downs(s, 2) & t < step_downs(s, 1)+ t_maxFrame;
        % figure(1);clf;plot(t, L, t(msk), L(msk));

        plot(ax1, t(msk_up), F(msk_up), t(msk_dwn), F(msk_dwn));
        % plot(ax2, t(msk_up), L(msk_up), t(msk_dwn), L(msk_dwn));linkaxes([ax1, ax2], 'x');
        % xlim([69.4, 70.6])

        

        % Coupled 2-exp model: F = Fss + A*(r*exp(-k1*t) - (1-r)*exp(-k2*t))
        % At t=0: F = Fss + A*(2r-1). Fast (k1) drives the spike; slow (k2) drives the undershoot.
        % Amplitude A and split ratio r are coupled → avoids exaggerated individual amplitudes.
        mdl = fittype('Fss + A*(r*exp(-k1*x) - (1-r)*exp(-k2*x))', ...
            'coeff', {'Fss','A','r','k1','k2'}, 'independent', 'x');

        % --- Step-up fit ---
        t0 = t(find(msk_up,1));
        tf = t(msk_up) - t0;
        Ff = F(msk_up);
        Fmed = median(Ff);
        A0_up = Ff(1) - Fmed;   % expected positive for stretch spike

        try
            f0_up = fit(tf, Ff, mdl, ...
                'StartPoint', [Fmed,  abs(A0_up), 0.8, 500, 50], ...
                'Lower',      [Fmed*0.5,  0,   0.5,  10,   1], ...
                'Upper',      [Fmed*2,  500,   1,  3000, 500]);
            stepResults(ci).Fss_up(s) = f0_up.Fss;
            stepResults(ci).k1_up(s)  = f0_up.k1;
            stepResults(ci).k2_up(s)  = f0_up.k2;
            stepResults(ci).A1_up(s)  = f0_up.A;
            stepResults(ci).A2_up(s)  = f0_up.r;
            plot(ax1, t(msk_up), f0_up(t(msk_up)-t0), 'r-', LineWidth=2);
            fprintf('  %s UP   s=%d dL=%.4f: Fss=%.1f A=%.2f r=%.2f k1=%.0f k2=%.0f\n', ...
                D(ai).label, s, stepResults(ci).dL(s), f0_up.Fss, f0_up.A, f0_up.r, f0_up.k1, f0_up.k2);
        catch ME_up
            fprintf('  [WARN] %s UP s=%d fit failed: %s\n', D(ai).label, s, ME_up.message);
        end

        % --- Step-down fit ---
        t0 = t(find(msk_dwn,1));
        tf = t(msk_dwn) - t0;
        Ff = F(msk_dwn);
        Fmed = median(Ff);
        A0_dwn = Ff(1) - Fmed;  % expected negative for release drop

        try
            f0_dwn = fit(tf, Ff, mdl, ...
                'StartPoint', [Fmed, -abs(A0_dwn), 0.8, 500, 50], ...
                'Lower',      [Fmed*0.5, -500,  0.5,  10,   1], ...
                'Upper',      [Fmed*2,      0,   1,  3000, 500]);
            stepResults(ci).Fss_dwn(s) = f0_dwn.Fss;
            stepResults(ci).k1_dwn(s)  = f0_dwn.k1;
            stepResults(ci).k2_dwn(s)  = f0_dwn.k2;
            stepResults(ci).A1_dwn(s)  = f0_dwn.A;
            stepResults(ci).A2_dwn(s)  = f0_dwn.r;
            plot(ax1, t(msk_dwn), f0_dwn(t(msk_dwn)-t0), 'b-', LineWidth=2);
            fprintf('  %s DOWN s=%d dL=%.4f: Fss=%.1f A=%.2f r=%.2f k1=%.0f k2=%.0f\n', ...
                D(ai).label, s, stepResults(ci).dL(s), f0_dwn.Fss, f0_dwn.A, f0_dwn.r, f0_dwn.k1, f0_dwn.k2);
        catch ME_dwn
            fprintf('  [WARN] %s DOWN s=%d fit failed: %s\n', D(ai).label, s, ME_dwn.message);
        end

        % --- Diagnostic plot for this step ---
        figure(30 + ci*10 + s); clf;
        tl_d = tiledlayout(3,1,'TileSpacing','compact');
        title(tl_d, sprintf('%s  step %d  dL=%.4f', D(ai).label, s, stepResults(ci).dL(s)));

        % UP panel
        axU = nexttile(tl_d);
        t0u = t(find(msk_up,1));
        tfu = t(msk_up) - t0u;  Fu = F(msk_up);
        plot(axU, tfu*1e3, Fu, 'k.'); hold(axU,'on');
        if isfield(stepResults(ci),'k1_up') && ~isnan(stepResults(ci).k1_up(s))
            plot(axU, tfu*1e3, f0_up(tfu), 'r-', LineWidth=2);
            text(axU, tfu(end)*1e3*0.6, max(Fu), ...
                sprintf('k1=%.0f  k2=%.0f\nA=%.2f  r=%.2f', f0_up.k1, f0_up.k2, f0_up.A, f0_up.r), ...
                'FontSize',8);
        end
        ylabel(axU,'Force'); title(axU,'Step UP'); grid(axU,'on');

        % DOWN panel
        axD = nexttile(tl_d);
        t0d = t(find(msk_dwn,1));
        tfd = t(msk_dwn) - t0d;  Fd = F(msk_dwn);
        plot(axD, tfd*1e3, Fd, 'k.'); hold(axD,'on');
        if isfield(stepResults(ci),'k1_dwn') && ~isnan(stepResults(ci).k1_dwn(s))
            plot(axD, tfd*1e3, f0_dwn(tfd), 'b-', LineWidth=2);
            text(axD, tfd(end)*1e3*0.6, min(Fd), ...
                sprintf('k1=%.0f  k2=%.0f\nA=%.2f  r=%.2f', f0_dwn.k1, f0_dwn.k2, f0_dwn.A, f0_dwn.r), ...
                'FontSize',8);
        end
        ylabel(axD,'Force'); title(axD,'Step DOWN'); grid(axD,'on');

        % Residuals panel
        axR = nexttile(tl_d);
        hold(axR,'on');
        if exist('f0_up','var'),  plot(axR, tfu*1e3, Fu  - f0_up(tfu),  'r.'); end
        if exist('f0_dwn','var'), plot(axR, tfd*1e3, Fd  - f0_dwn(tfd), 'b.'); end
        yline(axR,0,'k--'); ylabel(axR,'Residual'); xlabel(axR,'Time (ms)'); grid(axR,'on');
        legend(axR,{'UP resid','DOWN resid'},'Location','best','FontSize',7);    
        
    end

   
end

%% Summary figure: fitted parameters per condition
% Row 1 = step-up, Row 2 = step-down; cols = Fss, k1, k2, A, r
quantFields  = {'Fss_up','k1_up','k2_up','A1_up','A2_up', ...
                'Fss_dwn','k1_dwn','k2_dwn','A1_dwn','A2_dwn'};
quantLabels  = {'Fss (mN/mm²)','k1 (s⁻¹)','k2 (s⁻¹)','A (mN/mm²)','r', ...
                'Fss (mN/mm²)','k1 (s⁻¹)','k2 (s⁻¹)','A (mN/mm²)','r'};
nQ = length(quantFields);
nCols = 5; nRows = 2;
figure(20); clf;
tl = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
title(tl, 'Step-up (top row) / Step-down (bottom row)');
for qi = 1:nQ
    ax = nexttile(tl);
    hold(ax,'on');
    for ci = 1:length(activeIdx)
        if isfield(stepResults, quantFields{qi})
            vals = stepResults(ci).(quantFields{qi});   % 1 x nStepPairs
            plot(ax, 1:length(vals), vals, 'o-', 'DisplayName', stepResults(ci).label, linewidth =2);
        end
    end
    ylabel(ax, quantLabels{qi});
    xlabel(ax, 'Step pair #');
    legend(ax, 'Location','best','FontSize',7);
    % grid(ax,'on');
end

%%

%% Section 3: KTR Extraction via fitRecovery (70.8 - 71.5s)
fprintf('\n=== KTR Analysis ===\n');

% From vt: the restretch is the large positive velocity segment near 70.9s
ktr_restretch = findSegs(70.9, 70.95, +1, 10); % big positive vel
t_ktr_start = ktr_restretch(end, 2); % end of restretch

ktrResults = struct();
figure(3); clf;
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    datatable_fr = [D(ai).t, D(ai).L * SL0, D(ai).F];

    zone_ms = [t_ktr_start*1000 + 10, t_ktr_start*1000 + 300];

    subplot(1, length(activeIdx), ci); hold on;
    title(D(ai).label);
    try
        [~, ktr_val, df_val, del_val, E_val] = fitRecovery(datatable_fr, zone_ms, 0, [], true);
        ktrResults(ci).label = D(ai).label;
        ktrResults(ci).ktr = ktr_val;
        ktrResults(ci).df = df_val;
        ktrResults(ci).rmse = E_val;
        fprintf('%s: ktr = %.1f s-1, df = %.1f kPa, rmse = %.2f\n', ...
            D(ai).label, ktr_val, df_val, E_val);
    catch ME
        fprintf('%s: KTR fit failed: %s\n', D(ai).label, ME.message);
        ktrResults(ci).label = D(ai).label;
        ktrResults(ci).ktr = NaN;
        ktrResults(ci).df = NaN;
        ktrResults(ci).rmse = NaN;
    end
end
sgtitle('KTR Fits');

%% Section 4: Staircase Decay Fitting (72.4 - 73.2s)
fprintf('\n=== Staircase Analysis ===\n');

% From vt: staircase step-ups are positive velocity segments in 72.4-73.2
stair_ups = findSegs(72.4, 73.2, +1, 0.3);
nStairs = size(stair_ups, 1);
fprintf('Found %d staircase steps from velocity table\n', nStairs);

% Coupled model — same as section 2 step-up (A > 0, k1 > k2)
% F(t) = Fss + A*(r*exp(-k1*t) - (1-r)*exp(-k2*t))
mdl_stair = fittype('Fss + A*(r*exp(-k1*x) - (1-r)*exp(-k2*x))', ...
    'coeff', {'Fss','A','r','k1','k2'}, 'independent', 'x');

% Diagnostic figure: rows = conditions, cols = stairs
nCI = length(activeIdx);
figure(40); clf;
tl_stair = tiledlayout(nCI, nStairs, 'TileSpacing','compact','Padding','compact');
title(tl_stair, 'Staircase decay fits');

stairResults = struct();
for ci = 1:nCI-1
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;
    dt = median(diff(t));
    fs = 1/dt;
    isHighRes = fs > 1000;

    stairResults(ci).label   = D(ai).label;
    stairResults(ci).Fss     = nan(1, nStairs);
    stairResults(ci).k1      = nan(1, nStairs);
    stairResults(ci).k2      = nan(1, nStairs);
    stairResults(ci).A       = nan(1, nStairs);
    stairResults(ci).r       = nan(1, nStairs);
    stairResults(ci).dL      = nan(1, nStairs);
    stairResults(ci).L_after = nan(1, nStairs);

    for s = 1:nStairs
        t_vt_start = stair_ups(s, 1);
        if s < nStairs
            t_next = stair_ups(s+1, 1);
        else
            t_next = stair_ups(s, 1) + 0.15;  % extended for last (bigger) stair; slack starts at 73.5
        end

        % dL: L_before from just before VT start; L_after from plateau
        L_before = L(find(t <= t_vt_start, 1, 'last'));
        L_after  = mean(L(t >= t_vt_start+0.025 & t <= t_vt_start+0.035));
        stairResults(ci).dL(s)      = L_after - L_before;
        stairResults(ci).L_after(s) = L_after;

        % t0 = force peak after stretch is done (search starts at t_vt_end,
        % not t_vt_start, to exclude the rising ramp where fast-kinetics
        % conditions can peak during the stretch motion itself).
        t_vt_end = stair_ups(s, 2);
        m_post = t >= t_vt_end & t <= t_next - 0.005;
        Fsm = smoothdata(F(m_post), 'gaussian', round(0.003*fs));
        % Cap peak search to first 60% of window so the fit always has data
        max_search = max(1, round(0.6 * length(Fsm)));
        [~, pi_pk] = max(Fsm(1:max_search));
        idx_post = find(m_post);
        t0 = t(idx_post(pi_pk));
        t0 = stair_ups(s, 2);

        ax_s = nexttile(tl_stair, (ci-1)*nStairs + s);
        hold(ax_s, 'on');

        % Fit window from force peak to 5ms before next stretch
        mask_raw = t >= t0 & t <= t_next - 0.005;
        if sum(mask_raw) < 10
            title(ax_s, sprintf('s%d no data', s), FontSize=7, Color='r');
            continue
        end

        if ~isHighRes
            % 100 Hz: only Fss extraction, skip decay fit
            Fss_val = mean(F(mask_raw));
            stairResults(ci).Fss(s) = Fss_val;
            plot(ax_s, (t(mask_raw)-t0)*1e3, F(mask_raw), 'k.-');
            yline(ax_s, Fss_val, 'r-', LineWidth=2);
            title(ax_s, sprintf('s%d 100Hz\nFss=%.1f', s, Fss_val), FontSize=7);
            ylabel(ax_s,'F (kPa)'); xlabel(ax_s,'ms'); grid(ax_s,'on');
            continue
        end

        % Downsample to ~500 Hz (bin-average) to reduce noise for fitting
        ds_factor = max(1, round(fs / 1000));
        t_raw = t(mask_raw);  F_raw = F(mask_raw);
        n_bins = floor(length(t_raw) / ds_factor);
        t_ds = arrayfun(@(k) mean(t_raw((k-1)*ds_factor+1 : k*ds_factor)), 1:n_bins)';
        F_ds = arrayfun(@(k) mean(F_raw((k-1)*ds_factor+1 : k*ds_factor)), 1:n_bins)';
        tf = t_ds - t_ds(1);
        % Fss estimate from tail (last 20% of window); A0 = peak - tail
        Fss_est = mean(F_ds(max(1,end-round(0.2*n_bins)):end));
        A0 = F_ds(1) - Fss_est;
        A0 = (max(F_ds) - min(F_ds))/2;

        % Plot raw (gray) and downsampled (black)
        plot(ax_s, (t_raw-t0)*1e3, F_raw, 'Color',[0.8 0.8 0.8]);
        plot(ax_s, (t_ds-t0)*1e3,  F_ds,  'k.-', MarkerSize=6, LineWidth=1);
        
        % if s == 1

        try
            f0 = fit(tf, F_ds, mdl_stair, ...
                'StartPoint', [Fss_est, abs(A0),      0.8, 100, 10], ...
                'Lower',      [Fss_est*0.7, 0,        0.5,  10,  10], ...
                'Upper',      [Fss_est*1.05, abs(A0)*5, 1, 500, 40]);
            stairResults(ci).Fss(s) = f0.Fss;
            stairResults(ci).k1(s)  = f0.k1;
            stairResults(ci).k2(s)  = f0.k2;
            stairResults(ci).A(s)   = f0.A;
            stairResults(ci).r(s)   = f0.r;
            t_fit = linspace(0, tf(end), 300)';
            plot(ax_s, t_fit*1e3, f0(t_fit), 'r-', LineWidth=2);
            title(ax_s, sprintf('s%d L=%.3f\nk1=%.0f k2=%.0f', s, L_after, f0.k1, f0.k2), FontSize=7);
            fprintf('  %s Stair %d: L=%.4f  Fss=%.1f  k1=%.0f  k2=%.0f  A=%.2f  r=%.2f\n', ...
                D(ai).label, s, L_after, f0.Fss, f0.k1, f0.k2, f0.A, f0.r);
        catch ME
            title(ax_s, sprintf('s%d FAIL\n%s', s, ME.message), FontSize=6, Color='r');
            fprintf('  [WARN] %s Stair %d fit failed: %s\n', D(ai).label, s, ME.message);
        end
        ylabel(ax_s, 'F (kPa)'); xlabel(ax_s, 'ms'); grid(ax_s, 'on');
    end
end

%% Summary figure: k1, k2, Fss vs L per condition (high-res only)
figure(41); clf;
ax_s = nexttile(tl_stair, (ci-1)*nStairs + s);
tl41 = tiledlayout(1,3,'TileSpacing','compact');
% tl41 = tl_stair;
title(tl41,'Staircase summary');
axK1 = nexttile(tl41); hold(axK1,'on'); ylabel(axK1,'k1 (s^{-1})'); xlabel(axK1,'L (Lo)'); title(axK1,'Fast rate k1'); grid(axK1,'on');
axK2 = nexttile(tl41); hold(axK2,'on'); ylabel(axK2,'k2 (s^{-1})'); xlabel(axK2,'L (Lo)'); title(axK2,'Slow rate k2'); grid(axK2,'on');
axFs = nexttile(tl41); hold(axFs,'on'); ylabel(axFs,'Fss (kPa)');    xlabel(axFs,'L (Lo)'); title(axFs,'Fss');            grid(axFs,'on');
for ci = 1:nCI
    dt = median(diff(D(activeIdx(ci)).t));
    if 1/dt < 1000, continue; end   % skip 100Hz
    L_ax = stairResults(ci).L_after;
    plot(axK1, L_ax, stairResults(ci).k1, 'o-', LineWidth=2, DisplayName=stairResults(ci).label);
    plot(axK2, L_ax, stairResults(ci).k2, 'o-', LineWidth=2, DisplayName=stairResults(ci).label);
    plot(axFs, L_ax, stairResults(ci).Fss,'o-', LineWidth=2, DisplayName=stairResults(ci).label);
end
legend(axK1,'Location','best'); legend(axK2,'Location','best'); legend(axFs,'Location','best');

%% Section 5: Slack KTR Extraction (73.5 - 75.5s)
fprintf('\n=== Slack Analysis ===\n');

% From vt: slack events are large negative velocity segments
slack_drops = findSegs(73.5, 78, -1, 10);
% Restretches follow: large positive velocity
slack_restretches = findSegs(73.5, 78, +1, 2);
nSlacks = size(slack_drops, 1);
fprintf('Found %d slack events from velocity table\n', nSlacks);

ai = activeIdx(1);
mask = D(ai).t > 73.5 & D(ai).t < 78;
data_t = D(ai).t(mask); data_L = D(ai).L(mask); data_F = D(ai).F(mask);
vtm = [vt(1,:);vt(vt(:,1) > 73.5, :)];
features = extractSlackAttributes(data_t, data_F, data_L*2, vtm, struct(), [], true);

fn = fieldnames(features)';
plotFeatures(features, features, [], fn);
%%
slackResults = struct();
figure(5); clf;
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    datatable_fr = [D(ai).t, D(ai).L * SL0, D(ai).F];

    slackResults(ci).label = D(ai).label;
    slackResults(ci).ktr = nan(1, nSlacks);
    slackResults(ci).df = nan(1, nSlacks);
    slackResults(ci).SL = nan(1, nSlacks);
    slackResults(ci).rmse = nan(1, nSlacks);

    subplot(1, length(activeIdx), ci); hold on;
    title(sprintf('%s - Slack', D(ai).label));

    for s = 1:nSlacks
        t_drop_end = slack_drops(s, 2);
        % Find the matching restretch (first restretch after this drop)
        rs_idx = find(slack_restretches(:,1) > t_drop_end, 1);
        if ~isempty(rs_idx)
            t_restretch = slack_restretches(rs_idx, 1);
        else
            t_restretch = t_drop_end + 0.3;
        end

        zone_ms = [t_drop_end*1000 + 10, t_restretch*1000 - 10];
        if zone_ms(2) <= zone_ms(1) + 20
            fprintf('  Slack %d: zone too short, skipping\n', s);
            continue;
        end

        try
            [~, ktr_val, df_val, ~, E_val, SL_val] = fitRecovery(datatable_fr, zone_ms, 0, [], true);
            slackResults(ci).ktr(s) = ktr_val;
            slackResults(ci).df(s) = df_val;
            slackResults(ci).SL(s) = SL_val;
            slackResults(ci).rmse(s) = E_val;
            fprintf('  %s Slack %d: ktr=%.1f s-1, df=%.1f kPa, SL=%.3f um\n', ...
                D(ai).label, s, ktr_val, df_val, SL_val);
        catch ME
            fprintf('  Slack %d fit failed: %s\n', s, ME.message);
        end
    end
end
sgtitle('Slack KTR Fits');

%% Section 6: Visualization

% Figure 1: Overview
figure(1); clf;
ax1 = subplot(2,1,1); hold on;
ax2 = subplot(2,1,2); hold on;
for i = 1:n
    plot(ax1, D(i).t, D(i).F, 'Color', colors(i,:));
    plot(ax2, D(i).t, D(i).L, 'Color', colors(i,:));
end
ylabel(ax1, 'Force (kPa)'); ylabel(ax2, 'Length (Lo)');
xlabel(ax2, 'Time (s)');
legend(ax1, {D.label}, 'Interpreter', 'none', 'Location', 'best');
linkaxes([ax1 ax2], 'x');
xlim(ax1, [69, 76]);
for ax = [ax1, ax2]
    xline(ax, 70.5, '--k', 'Alpha', 0.3);
    xline(ax, 72.3, '--k', 'Alpha', 0.3);
    xline(ax, 73.5, '--k', 'Alpha', 0.3);
end
title(ax1, 'Overview: all conditions');

% Figure 2: Step perturbation fits
figure(2); clf;
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;
    win = t >= 69.2 & t <= 70.6;

    subplot(length(activeIdx), 2, (ci-1)*2+1); hold on;
    plot(t(win), F(win), 'Color', colors(ai,:));
    title(sprintf('%s - Force', D(ai).label)); ylabel('F (kPa)');

    subplot(length(activeIdx), 2, (ci-1)*2+2); hold on;
    plot(t(win), L(win), 'Color', colors(ai,:));
    title('Length'); ylabel('L (Lo)');
end
sgtitle('Step Perturbations');

% Figure 4: Staircase
figure(4); clf;
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;
    win = t >= 72.3 & t <= 73.4;

    subplot(length(activeIdx), 2, (ci-1)*2+1); hold on;
    plot(t(win), F(win), 'Color', colors(ai,:));
    title(sprintf('%s - Force', D(ai).label)); ylabel('F (kPa)');

    subplot(length(activeIdx), 2, (ci-1)*2+2); hold on;
    plot(t(win), L(win), 'Color', colors(ai,:));
    title('Length'); ylabel('L (Lo)');
end
sgtitle('Staircase Stretch');

% Figure 6: Comparison bar charts
figure(6); clf;

subplot(1,3,1); hold on;
ktr_vals = [ktrResults.ktr];
bar(ktr_vals);
set(gca, 'XTickLabel', {ktrResults.label}, 'XTickLabelRotation', 30);
ylabel('ktr (s^{-1})'); title('KTR');

subplot(1,3,2); hold on;
for ci = 1:length(activeIdx)
    plot(stepResults(ci).dL, stepResults(ci).k_up, 'o-', 'Color', colors(activeIdx(ci),:), ...
        'DisplayName', [stepResults(ci).label ' up']);
    plot(stepResults(ci).dL, stepResults(ci).k_down, 's--', 'Color', colors(activeIdx(ci),:), ...
        'DisplayName', [stepResults(ci).label ' down']);
end
xlabel('Step dL (Lo)'); ylabel('k (s^{-1})');
title('Step Perturbation Rates'); legend('Location', 'best');

subplot(1,3,3); hold on;
for ci = 1:length(activeIdx)
    if ~isempty(slackResults(ci).SL) && any(~isnan(slackResults(ci).SL))
        plot(slackResults(ci).SL, slackResults(ci).ktr, 'o-', ...
            'Color', colors(activeIdx(ci),:), 'DisplayName', slackResults(ci).label);
    end
end
xlabel('SL (um)'); ylabel('ktr (s^{-1})');
title('Slack ktr vs SL'); legend('Location', 'best');

sgtitle('Feature Comparison');

%% Section 7: Results Table
fprintf('\n\n========== RESULTS SUMMARY ==========\n\n');
fprintf('%-25s', 'Condition');
fprintf('ktr_KTR\t\t');
for s = 1:3; fprintf('k_up_%d\t\t', s); end
for s = 1:3; fprintf('k_dn_%d\t\t', s); end
for s = 1:3; fprintf('ktr_slk_%d\t', s); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 140));

for ci = 1:length(activeIdx)
    fprintf('%-25s', ktrResults(ci).label);
    fprintf('%.1f\t\t', ktrResults(ci).ktr);
    for s = 1:min(3, length(stepResults(ci).k_up))
        fprintf('%.1f\t\t', stepResults(ci).k_up(s));
    end
    for s = 1:min(3, length(stepResults(ci).k_down))
        fprintf('%.1f\t\t', stepResults(ci).k_down(s));
    end
    for s = 1:min(3, length(slackResults(ci).ktr))
        fprintf('%.1f\t\t', slackResults(ci).ktr(s));
    end
    fprintf('\n');
end

save(fullfile(dataDir, 'results_0327.mat'), 'stepResults', 'ktrResults', 'stairResults', 'slackResults', 'D');
fprintf('\nResults saved to %sresults_0327.mat\n', dataDir);

%% Helper function
function segs = findVtSegments(vt, tmin, tmax, velSign, velThresh)
    % Find velocity table segments in [tmin, tmax] where
    % velocity*velSign > velThresh
    % Returns Nx2 matrix [seg_start, seg_end]
    mask = vt(:,1) >= tmin & vt(:,1) <= tmax & (vt(:,2) * velSign) > velThresh;
    idx = find(mask);
    if isempty(idx)
        segs = zeros(0, 2);
        return;
    end
    % Cluster consecutive entries
    breaks = [0; find(diff(idx) > 1); length(idx)];
    segs = zeros(length(breaks)-1, 2);
    for k = 1:length(breaks)-1
        gi = idx(breaks(k)+1);         % first vt row of this high-vel segment
        ge = idx(breaks(k+1));         % last vt row of this high-vel segment
        ge_next = min(ge + 1, size(vt, 1)); % first plateau row = end of motion
        segs(k, :) = [vt(gi, 1), vt(ge_next, 1)];
    end
end
