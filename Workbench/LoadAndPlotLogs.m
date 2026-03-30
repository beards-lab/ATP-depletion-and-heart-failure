%% Load and plot all Log_*.txt files from a data directory
close all; clear; clc;

dataDir = '../data/03 27 2026 M/';
S = dir([dataDir '*Log_*.txt']);
[~, idx] = sort({S.name});
S = S(idx);

n = length(S);
colors = lines(n);
labels = cell(n, 1);

%% Load all files and detect sync anchor points
% Anchor events (in t-space after initial drop sync):
%   t1: big L drop in 68-72s
%   t2: start of gradual L step-up in 72-73s
%   t3: big L drop in 74-75s (may not exist in all files)

data = struct('t', {}, 'F', {}, 'L', {}, 'anchors', {});
for i = 1:n
    fpath = fullfile(S(i).folder, S(i).name);
    M = readmatrix(fpath, 'NumHeaderLines', 4);

    % Initial sync: t=0 at first big L drop
    idx0 = find(diff(M(:,2)) < -0.01, 1);
    t = (M(:,1) - M(idx0,1)) / 1000;

    % t1: big L drop in 68-72s window
    w1 = find(t >= 68 & t <= 72);
    i1 = find(diff(M(w1,2)) < -0.05, 1);
    t1 = t(w1(i1));

    % t2: first step-up (start of ramp) in 72-73s window
    w2 = find(t >= 72 & t <= 73);
    i2 = find(diff(M(w2,2)) > 0.003, 1);
    t2 = t(w2(i2));

    % t3: big L drop in 74-75s window (optional)
    w3 = find(t >= 74 & t <= 75);
    i3 = find(diff(M(w3,2)) < -0.05, 1);
    if ~isempty(i3)
        t3 = t(w3(i3));
    else
        t3 = NaN;
    end

    data(i).t       = t;
    data(i).F       = M(:,3);
    data(i).L       = M(:,2);
    data(i).anchors = [t1, t2, t3];
    labels{i} = strrep(S(i).name(1:end-4), '_', ' ');
end

%% Use file 06 (last file) as reference time axis
refIdx = n; % 06_Log_8mM_Active_PNB_Mava is last after sorting
t_ref  = data(refIdx).t;
a_ref  = data(refIdx).anchors; % [t1_ref, t2_ref, t3_ref]

figure(1); clf;
ax1 = subplot(2,1,1); hold on;
ax2 = subplot(2,1,2); hold on;

for i = 1:n
    t_i = data(i).t;
    a_i = data(i).anchors;

    % Piecewise shift: shift each segment so its anchor aligns with ref
    t_shifted = t_i; % before 68s: no extra shift (initial sync is fine)

    % Segment 68-72: shift by (t1_ref - t1)
    seg1 = t_i >= 68 & t_i < 72;
    t_shifted(seg1) = t_i(seg1) + (a_ref(1) - a_i(1));

    % Segment 72-74.4: shift by (t2_ref - t2)
    seg2 = t_i >= 72 & t_i < 74.4;
    t_shifted(seg2) = t_i(seg2) + (a_ref(2) - a_i(2));

    % Segment 74.4-end: shift by (t3_ref - t3), only if t3 exists
    if ~isnan(a_i(3)) && ~isnan(a_ref(3))
        seg3 = t_i >= 74.4;
        t_shifted(seg3) = t_i(seg3) + (a_ref(3) - a_i(3));
    end

    % Ensure strict monotonicity (shifts at boundaries can cause overlap/dups)
    [t_shifted, iu] = unique(t_shifted, 'stable');
    F_i = data(i).F(iu);
    L_i = data(i).L(iu);
    ok = [true; diff(t_shifted) > 0];
    t_shifted = t_shifted(ok);
    F_i = F_i(ok);
    L_i = L_i(ok);

    % Resample onto reference time axis
    F_rs = interp1(t_shifted, F_i, t_ref, 'linear', NaN);
    L_rs = interp1(t_shifted, L_i, t_ref, 'linear', NaN);

    plot(ax1, t_ref, F_rs, 'Color', colors(i,:));
    plot(ax2, t_ref, L_rs, 'Color', colors(i,:));
end

ylabel(ax1, 'Force (kPa)');
ylabel(ax2, 'Length (Lo)');
xlabel(ax2, 'Time (s, synced to 06 Log)');
legend(ax1, labels, 'Interpreter', 'none', 'Location', 'best');
title(ax1, strrep(dataDir, '_', ' '));
linkaxes([ax1 ax2], 'x');
