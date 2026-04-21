% BuildUpdatedSlackProtocol.m
% Build an updated Baker slack velocity table and save as
% data/bakers_slack8mM_update.mat.
%
% Changes vs. original (bakers_slack8mM_all.mat):
%   1. All slack hold durations      -> t_ref (uniform)
%   2. All recovery periods          -> t_ref (uniform)
%   3. All slack shortening speeds   -> v = -174.0 (col-2, fastest in original)
%      Shortening duration adjusted per cycle so same SL target is reached.
%   4. All restretch durations       -> 0.020 s (original velocities preserved)
%   5. Two ±5 % ML square perturbations:
%        Perturb 1: ends 200 ms before first slack
%        Perturb 2: starts 500 ms after last restretch end

datastruct = load('../data/bakers_slack8mM_all.mat');
vt_orig   = datastruct.velocitytable;
datatable = datastruct.datatable;

%% Constants
t_ref        = vt_orig(19,1) - vt_orig(18,1);  % 0.3794 s
rsl0_ML      = 2.0;                             % µm
SL0          = vt_orig(1,4);                    % 2.2006 µm
SL_step      = 0.0055 * SL0;                    % ±0.0121 µm  (0.55 % ML)
dt_pfast     = 0.005;                           % 5 ms
dt_phold     = 0.020;                           % 20 ms
v_perturb    = SL_step / dt_pfast / rsl0_ML;   % 11.00
perturb_dur  = 5*dt_pfast + 2*dt_phold;        % 0.065 s
v_slack_fast = -174.0;                          % col-2
t_restretch  = 0.020;                           % s

% Perturbation 1 SL extremes
SL_hi1 = SL0 + SL_step;
SL_lo1 = SL0 - SL_step;

% Timing
t_first_slack = vt_orig(3,1);                          % 1.1594 s
t_p1_start    = t_first_slack - 0.200 - perturb_dur;   % 0.8944 s
t_p1_end      = t_p1_start + perturb_dur;              % 0.9594 s

% Cycle data from original (cols: SL_before, SL_after_slack, v_rest_orig, SL_after_rest)
cycle_info = [
    vt_orig(3,4),  vt_orig(4,4),  vt_orig(5,2),  vt_orig(6,4);
    vt_orig(7,4),  vt_orig(8,4),  vt_orig(9,2),  vt_orig(10,4);
    vt_orig(11,4), vt_orig(12,4), vt_orig(13,2), vt_orig(14,4);
    vt_orig(15,4), vt_orig(16,4), vt_orig(17,2), vt_orig(18,4);
    vt_orig(19,4), vt_orig(20,4), vt_orig(21,2), vt_orig(22,4);
];

%% Build vt_new
vt_new = zeros(0, 4);

% Row 1: initial isometric
vt_new(end+1,:) = [0,          0,          0,          SL0];

% Perturbation 1
t = t_p1_start;
vt_new(end+1,:) = [t,          v_perturb,  2*v_perturb, SL0];    t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_hi1]; t = t + dt_phold;
vt_new(end+1,:) = [t,         -v_perturb, -2*v_perturb, SL_hi1]; t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL0];    t = t + dt_pfast;
vt_new(end+1,:) = [t,         -v_perturb, -2*v_perturb, SL0];    t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_lo1]; t = t + dt_phold;
vt_new(end+1,:) = [t,          v_perturb,  2*v_perturb, SL_lo1]; t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL0];   % t = t_p1_end

% Isometric until first slack (velocity=0 already set above; next row starts slack)

% Slack cycles
t_cursor   = t_first_slack;
t_p2_start = NaN;
SL_base_p2 = NaN;

for c = 1:5
    SL_before = cycle_info(c,1);
    SL_after  = cycle_info(c,2);
    v_rest    = cycle_info(c,3);
    SL_rest   = cycle_info(c,4);

    t_fast = abs(SL_after - SL_before) / (abs(v_slack_fast) * rsl0_ML);

    vt_new(end+1,:) = [t_cursor,          v_slack_fast, 2*v_slack_fast, SL_before];
    t_cursor = t_cursor + t_fast;
    vt_new(end+1,:) = [t_cursor,          0,            0,              SL_after];
    t_cursor = t_cursor + t_ref;

    vt_new(end+1,:) = [t_cursor,          v_rest,       2*v_rest,       SL_after];
    t_cursor = t_cursor + t_restretch;
    vt_new(end+1,:) = [t_cursor,          0,            0,              SL_rest];

    if c < 5
        t_cursor = t_cursor + t_ref;
    else
        t_p2_start = t_cursor + 0.500;
        SL_base_p2 = SL_rest;
    end
end

% Perturbation 2
SL_hi2 = SL_base_p2 + SL_step;
SL_lo2 = SL_base_p2 - SL_step;

t = t_p2_start;
vt_new(end+1,:) = [t,          v_perturb,  2*v_perturb, SL_base_p2]; t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_hi2];     t = t + dt_phold;
vt_new(end+1,:) = [t,         -v_perturb, -2*v_perturb, SL_hi2];     t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_base_p2]; t = t + dt_pfast;
vt_new(end+1,:) = [t,         -v_perturb, -2*v_perturb, SL_base_p2]; t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_lo2];     t = t + dt_phold;
vt_new(end+1,:) = [t,          v_perturb,  2*v_perturb, SL_lo2];     t = t + dt_pfast;
vt_new(end+1,:) = [t,          0,          0,           SL_base_p2];

% Tail
vt_new(end+1,:) = [t + 0.100,  0,          0,           SL_base_p2];

%% Sanity checks
assert(all(diff(vt_new(:,1)) > 0), 'ERROR: non-monotone time entries in vt_new');
fprintf('New table: %d rows, %.3f s  (orig: %d rows, %.3f s)\n', ...
    size(vt_new,1), vt_new(end,1), size(vt_orig,1), vt_orig(end,1));

slack_idx = find(vt_new(:,2) == v_slack_fast);
assert(numel(slack_idx) == 5, 'Expected exactly 5 slack rows');
fprintf('\nCycle verification:\n');
for c = 1:5
    r = slack_idx(c);
    t_fast_d = vt_new(r+1,1) - vt_new(r,1);
    t_hold_d = vt_new(r+2,1) - vt_new(r+1,1);
    t_rest_d = vt_new(r+3,1) - vt_new(r+2,1);
    if c < 5
        t_recov_d = vt_new(r+4,1) - vt_new(r+3,1);
        fprintf('  Cycle %d: fast=%.4fs  hold=%.4fs  restretch=%.4fs  recovery=%.4fs\n', ...
            c, t_fast_d, t_hold_d, t_rest_d, t_recov_d);
    else
        fprintf('  Cycle %d: fast=%.4fs  hold=%.4fs  restretch=%.4fs\n', ...
            c, t_fast_d, t_hold_d, t_rest_d);
    end
end

%% SL traces
to = []; SLo = []; SLc = vt_orig(1,4);
for i = 1:size(vt_orig,1)-1
    ts  = linspace(vt_orig(i,1), vt_orig(i+1,1), max(3, round((vt_orig(i+1,1)-vt_orig(i,1))/5e-4)+1));
    SLs = SLc + vt_orig(i,2)*rsl0_ML*(ts - vt_orig(i,1));
    to  = [to,  ts(1:end-1)];  SLo = [SLo, SLs(1:end-1)];  SLc = SLs(end); %#ok<AGROW>
end
to(end+1) = vt_orig(end,1);  SLo(end+1) = SLc;

tn = []; SLn = []; SLc = vt_new(1,4);
for i = 1:size(vt_new,1)-1
    ts  = linspace(vt_new(i,1), vt_new(i+1,1), max(3, round((vt_new(i+1,1)-vt_new(i,1))/5e-4)+1));
    SLs = SLc + vt_new(i,2)*rsl0_ML*(ts - vt_new(i,1));
    tn  = [tn,  ts(1:end-1)];  SLn = [SLn, SLs(1:end-1)];  SLc = SLs(end); %#ok<AGROW>
end
tn(end+1) = vt_new(end,1);  SLn(end+1) = SLc;

%% Plot comparison
figure(201); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(to, SLo, 'b-', 'LineWidth', 1.5);
ylabel('SL (µm)');
title(sprintf('Original  (%d rows, %.2f s)', size(vt_orig,1), vt_orig(end,1)));
grid on;  xlim([0 vt_orig(end,1)+0.05]);  ylim([1.80 2.35]);

nexttile; hold on;
plot(tn, SLn, 'r-', 'LineWidth', 1.5);
xregion(t_p1_start, t_p1_start+perturb_dur, 'FaceColor',[0.4 0.7 1],  'FaceAlpha',0.4,'EdgeColor','none');
xregion(t_p2_start, t_p2_start+perturb_dur, 'FaceColor',[1.0 0.6 0.2],'FaceAlpha',0.4,'EdgeColor','none');
ylabel('SL (µm)');
title(sprintf('New  (%d rows, %.2f s) — blue=perturb1, orange=perturb2', size(vt_new,1), vt_new(end,1)));
grid on;  xlim([0 vt_new(end,1)+0.05]);  ylim([1.80 2.35]);

nexttile; hold on;
win = [t_p1_start-0.04, t_first_slack+0.04];
plot(to(to>=win(1)&to<=win(2)), SLo(to>=win(1)&to<=win(2)), 'b-', 'LineWidth',1.5, 'DisplayName','Original');
plot(tn(tn>=win(1)&tn<=win(2)), SLn(tn>=win(1)&tn<=win(2)), 'r-', 'LineWidth',2,   'DisplayName','New');
xline(t_first_slack, 'k--', 'LineWidth',1, 'Label','1st slack');
xlabel('Time (s)');  ylabel('SL (µm)');
title(sprintf('Zoom: perturb1  (±5%%ML = ±%.3f µm, ends %.0f ms before first slack)', ...
    SL_step, (t_first_slack - t_p1_end)*1000));
legend('Location','best');  grid on;
sgtitle('Original vs New Baker Slack Protocol');

%% Save
velocitytable = vt_new;
save('../data/bakers_slack8mM_update.mat', 'velocitytable', 'datatable');
fprintf('\nSaved  ../data/bakers_slack8mM_update.mat\n');
fprintf('SL range: [%.4f, %.4f] µm\n', min(SLn), max(SLn));
