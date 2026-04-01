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

%% Helper: find velocity table segments by time and velocity sign
% Returns Nx2 matrix of [seg_start_time, seg_end_time] for segments
% where velocity exceeds threshold in the given time window
findSegs = @(tmin, tmax, velSign, velThresh) ...
    findVtSegments(vt, tmin, tmax, velSign, velThresh);

%% Section 2: Step Perturbation Decay/Recovery Fitting (69.2 - 70.5s)
fprintf('\n=== Step Perturbation Analysis ===\n');

% From velocity table, identify step-up and step-down segments
step_ups = findVtSegments(vt, 69.2, 70.5, +1, 0.5);    % positive velocity > 0.5 Lo/s
step_downs = findVtSegments(vt, 69.2, 70.5, -1, 0.5);   % negative velocity
nStepPairs = min(size(step_ups,1), size(step_downs,1));
fprintf('Found %d step-up and %d step-down segments from velocity table\n', ...
    size(step_ups,1), size(step_downs,1));

mdl = fittype('Fss + A*exp(-k*x) +b*x', 'coeff',{'Fss','A','k', 'b'},'independent', 'x');

stepResults = struct();
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;

    stepResults(ci).label = D(ai).label;
    stepResults(ci).k_up = nan(1, nStepPairs);
    stepResults(ci).k_down = nan(1, nStepPairs);
    stepResults(ci).A_up = nan(1, nStepPairs);
    stepResults(ci).A_down = nan(1, nStepPairs);
    stepResults(ci).dL = nan(1, nStepPairs);

    for s = 1:nStepPairs
        % Step-up: plateau starts at step_ups(s,2), ends at step_downs(s,1)
        t_up_end = step_ups(s, 2);
        t_dn_start = step_downs(s, 1);
        L_before = L(find(t <= step_ups(s,1), 1, 'last'));
        L_after  = L(find(t >= t_up_end, 1));
        stepResults(ci).dL(s) = L_after - L_before;
        msk_up = t > step_ups(s, 1) & t < step_ups(s, 2);
        msk_dwn = t > step_downs(s, 1) & t < step_downs(s, 2);
        % figure(1);clf;plot(t, L, t(msk), L(msk));
        figure(101);clf;hold on; ax1 = nexttile(1);hold on;
        plot(t, F, t(msk_up), F(msk_up), t(msk_dwn), F(msk_dwn));
        ax2 = nexttile(2);plot(t, L,t(msk_up), L(msk_up), t(msk_dwn), L(msk_dwn));linkaxes([ax1, ax2], 'x');

        % Fit step-up decay (force spike decaying after stretch)
        buf = 0.003; % 3ms buffer
        % mask = t >= t_up_end + buf & t <= t_dn_start - buf;
        % if sum(msk_up) > 10
        t0 = t(find(msk_up,1));
        tf = t(msk_up) - t0;
        Ff = F(msk_up);
        try
            f0 = fit(tf, Ff, mdl, ...
                'StartPoint', [median(Ff), Ff(1)-median(Ff), 50, 0], ...
                'Lower', [median(Ff)/2, 10, 1, -100], 'Upper', [Inf, Inf, 500, 100]);
            stepResults(ci).k_up(s) = f0.k;
            stepResults(ci).A_up(s) = f0.A;
            plot(ax1, t(msk_up), f0(t(msk_up)-t0), linewidth =1);
            % Ffit = feval(mdl, median(Ff), max(Ff)-median(Ff), 100, 100, tf); plot(tf, Ff, tf, Ffit)
        catch; end

        % figure(101);clf;hold on; ax1 = nexttile(1);plot(t, F, tf + t0, Ff);ax2 = nexttile(2);plot(t, L);linkaxes([ax1, ax2], 'x');
        t0 = t(find(msk_dwn,1));
        tf = t(msk_dwn) - t0;
        Ff = F(msk_dwn);
        try
            f0 = fit(tf, Ff, mdl, ...
                'StartPoint', [median(Ff), Ff(1)-median(Ff), 50, 0], ...
                'Lower', [median(Ff)/2, 10, 1, -100], 'Upper', [Inf, Inf, 100, 100]);
            stepResults(ci).k_up(s) = f0.k;
            stepResults(ci).A_up(s) = f0.A;
            plot(ax1, t(msk_dwn), f0(t(msk_dwn)-t0), linewidth =1);
        catch; end
        % Step-down: plateau starts at step_downs(s,2)
        t_dn_end = step_downs(s, 2);
        if s < nStepPairs
            t_next = step_ups(s+1, 1);
        else
            t_next = 70.6;
        end
        mask = t >= t_dn_end + buf & t <= t_next - buf;
        if sum(mask) > 10
            tf = t(mask) - t(find(mask,1));
            Ff = F(mask);
            try
                f0 = fit(tf, Ff, mdl, ...
                    'StartPoint', [median(Ff), Ff(1)-median(Ff), 50], ...
                    'Lower', [-Inf, -Inf, 0.1], 'Upper', [Inf, Inf, 2000]);
                stepResults(ci).k_down(s) = f0.k;
                stepResults(ci).A_down(s) = f0.A;
            catch; end
        end

        fprintf('  %s Step %d: dL=%.4f, k_up=%.1f, k_down=%.1f s-1\n', ...
            D(ai).label, s, stepResults(ci).dL(s), stepResults(ci).k_up(s), stepResults(ci).k_down(s));
    end
end

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

stairResults = struct();
for ci = 1:length(activeIdx)
    ai = activeIdx(ci);
    t = D(ai).t; L = D(ai).L; F = D(ai).F;

    stairResults(ci).label = D(ai).label;
    stairResults(ci).k_decay = nan(1, nStairs);
    stairResults(ci).A = nan(1, nStairs);
    stairResults(ci).Fss = nan(1, nStairs);
    stairResults(ci).dL = nan(1, nStairs);
    stairResults(ci).L_after = nan(1, nStairs);

    for s = 1:nStairs
        t_step_end = stair_ups(s, 2);
        if s < nStairs
            t_next = stair_ups(s+1, 1);
        else
            t_next = 73.35; % before the big drop
        end

        L_before = L(find(t <= stair_ups(s,1), 1, 'last'));
        L_after  = L(find(t >= t_step_end, 1));
        stairResults(ci).dL(s) = L_after - L_before;
        stairResults(ci).L_after(s) = L_after;

        buf = 0.001; % 1ms buffer to capture fast spike decay
        mask = t >= t_step_end + buf & t <= t_next - buf;
        if sum(mask) > 10
            tf = t(mask) - t(find(mask,1));
            Ff = F(mask);
            try
                f0 = fit(tf, Ff, mdl, ...
                    'StartPoint', [median(Ff), Ff(1)-median(Ff), 30], ...
                    'Lower', [-Inf, -Inf, 0.1], 'Upper', [Inf, Inf, 1000]);
                stairResults(ci).k_decay(s) = f0.k;
                stairResults(ci).A(s) = f0.A;
                stairResults(ci).Fss(s) = f0.Fss;
            catch; end
        end

        fprintf('  %s Stair %d: dL=%.4f, L=%.4f, k_decay=%.1f s-1\n', ...
            D(ai).label, s, stairResults(ci).dL(s), stairResults(ci).L_after(s), stairResults(ci).k_decay(s));
    end
end

%% Section 5: Slack KTR Extraction (73.5 - 75.5s)
fprintf('\n=== Slack Analysis ===\n');

% From vt: slack events are large negative velocity segments
slack_drops = findSegs(73.5, 75.5, -1, 10);
% Restretches follow: large positive velocity
slack_restretches = findSegs(73.5, 75.5, +1, 2);
nSlacks = size(slack_drops, 1);
fprintf('Found %d slack events from velocity table\n', nSlacks);

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
        % gi = idx(breaks(k)+1);
        % ge = idx(breaks(k+1));
        % % Segment end is the NEXT VT row (velocity at ge applies until ge+1)
        % ge_next = min(ge + 1, size(vt, 1));
        % segs(k, :) = [vt(gi, 1), vt(ge_next, 1)];
        gi = idx(breaks(k)+1)+1;
        % Segment end is the NEXT VT row (velocity at ge applies until ge+1)
        ge_next = min(gi+1, size(vt, 1));
        segs(k, :) = [vt(gi, 1), vt(ge_next, 1)];
        
    end
end
