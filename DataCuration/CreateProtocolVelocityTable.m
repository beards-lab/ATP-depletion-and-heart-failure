% TODO:
% - extend to ktr and minor perturbations export
% - fix the velocitytable for those as well

% CreateProtocolVelocityTable.m
% clear;
dataDir    = '../data/03 27 2026 M/';
fpath      = [dataDir '06_Merged_8mM_Active_PNB_Mava.txt'];
fpath      = [dataDir '02_Merged_8mM_Active.txt'];
% exportData = [dataDir '06_Merged_8mM_Active_PNB_Mava_proc.mat'];
% outPath    = '../data/protocol_03_27_2026_velocitytable_slack.mat';
outPath    = '../data/protocol_03_27_2026_8mM_slack.mat';
% outPath    = '../data/protocol_03_27_2026_ActivePNBMava_slack.mat';

% set #02
% dataDir = '../data/04 03 2026 F/';
% fpath = [dataDir '04_Merged_8mM_Active_PNB_Mava.txt'];
% outPath = '../data/protocol_04_03_2026_velocitytable.mat';


M = readmatrix(fpath);

% We only want t >= 0
idx0 = find(M(:,1) >= -0.1, 1);
t = M(idx0:end, 1);
L = M(idx0:end, 2);
% clf;ax1 = nexttile(1);plot(t, L); ax2 = nexttile(2);plot(t(2:end), diff(L)); linkaxes([ax1 ax2], 'x')

F = M(idx0:end, 3);

% Ramer-Douglas-Peucker line simplification to find turning points
epsilon = 0.001; % 0.2% length tolerance
keep = rdp(t, L, epsilon);

t_vt = t(keep);
L_vt = L(keep);

% Post-process RDP output: cluster-based collapse.
% RDP can emit multiple close-in-time points near a sharp transition, where
% raw-L noise causes boundary points to land off the true plateau (so naive
% removal keeps a mid-slope point with the wrong L). Within a "cluster"
% (consecutive points with dt < dt_cluster), we identify the "near-max" and
% "near-min" groups (within tol of the cluster's extrema). The group first in
% time is the before-plateau; keep its LATEST point. The other group is the
% after-plateau; keep its EARLIEST point. All mid-slope points are dropped.
% Flat clusters (range < 4*tol) collapse to a single point (the last one).
% Direction of the step (up or down) is detected from time-order of the groups.
dt_cluster = 0.020;     % 20 ms — cluster window
tol        = epsilon;   % amplitude tolerance for "near plateau"
hard_dt    = 5e-4;      % 0.5 ms absolute minimum (final ODE-stiffness safety)

N = length(t_vt);
keep_mask = true(N, 1);
k = 1;
while k < N
    j = k;
    while j < N && (t_vt(j+1) - t_vt(j)) < dt_cluster
        j = j + 1;
    end
    if j > k
        Lc = L_vt(k:j);
        Lmax = max(Lc); Lmin = min(Lc);
        if (Lmax - Lmin) < 4*tol
            % Flat cluster — keep only the last point
            keep_mask(k:j-1) = false;
        else
            nm = Lc >= Lmax - tol;
            nn = Lc <= Lmin + tol;
            first_nm = find(nm, 1, 'first');
            first_nn = find(nn, 1, 'first');
            keep_mask(k:j) = false;
            if first_nm < first_nn
                % Downward step: latest near-max + earliest near-min
                keep_mask(k - 1 + find(nm, 1, 'last'))  = true;
                keep_mask(k - 1 + find(nn, 1, 'first')) = true;
            else
                % Upward step: latest near-min + earliest near-max
                keep_mask(k - 1 + find(nn, 1, 'last'))  = true;
                keep_mask(k - 1 + find(nm, 1, 'first')) = true;
            end
        end
    end
    k = j + 1;
end
t_vt = t_vt(keep_mask);
L_vt = L_vt(keep_mask);

% Final ODE-stiffness safety: merge any residual pair with dt < hard_dt
k = 1;
while k < length(t_vt)
    if t_vt(k+1) - t_vt(k) < hard_dt
        t_vt(k+1) = [];
        L_vt(k+1) = [];
    else
        k = k + 1;
    end
end

% Compute velocities (Lo/s)
v_vt = diff(L_vt) ./ diff(t_vt);
v_vt = [v_vt; 0]; % last segment velocity is 0

% velocitytable format: [Time(s), Vel(ML/s), Vel_um(um/s), L(ML)]
% (um/s isn't used actively but we will populate it as v_vt * 2 or similar)
vel_um = v_vt * 2; % placeholder scaling

velocitytable = [t_vt, v_vt, vel_um, L_vt];

% Generate datatable (downsampled full trace for cost function evaluation)
dsf = 5; % downsample factor
datatable = [downsample(t, dsf), downsample(L, dsf), downsample(F, dsf)];


% clip it - catch the first slack burst
tWarmStart   = 50;    % warm-start entry time prepended to velocitytable
slackSegment = [74, 78];
i_firstSlack = find(velocitytable(:, 1) >= slackSegment(1) & velocitytable(:, 2) < -10, 1, 'first');
clipping     = [velocitytable(i_firstSlack, 1), slackSegment(2)];  % [t_start, t_end] seconds — adjust as needed

% ── velocitytable ──
velocitytable = velocitytable(velocitytable(:,1) >= clipping(1) & velocitytable(:,1) <= clipping(2), :);
if velocitytable(1,1) > slackSegment(1) || velocitytable(1,1) > tWarmStart
    warmRow       = [tWarmStart, 0, 0, velocitytable(1, 4)];
    velocitytable = [warmRow; velocitytable];
end

if velocitytable(end,1) < slackSegment(2)
    % add an end row to include the ending zero velocity part
    endRow = [slackSegment(2), 0, 0, velocitytable(end, 4)];
    velocitytable = [velocitytable;endRow];
end


% ── datatable: clip only ──
datatable = datatable(datatable(:,1) >= slackSegment(1) & datatable(:,1) <= slackSegment(2), :);

% clf;nexttile;plot(datatable(:,1), datatable(:,2), '|-', velocitytable(:, 1),velocitytable(:, 4), 'x-');

if mean(datatable(:,2)) < 1.5
    % fix ML to SL
    datatable(:,2) = datatable(:,2)*2;
end
%% fit data features in slack
features_data = struct();
features_data = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable, features_data, [], true);
            % % Print extracted values to console so they can be pasted into the else branch below
            % fieldNames = fieldnames(features_data);
            % for k = 1:numel(fieldNames)
            %     fname = fieldNames{k};
            %     val   = features_data.(fname);
            %     if ~isnumeric(val), continue; end
            %     s = mat2str(round(val, 4, 'significant'));
            %     fprintf('features_data.%s = %s;\n', fname, s);
            % end


%%

% Save it
save(outPath, 'velocitytable', 'datatable', 'features_data');
fprintf('Saved new velocitytable (%dx4) and datatable (%dx3) to %s\n', ...
        size(velocitytable,1), size(datatable,1), outPath);

figure(50); clf;
ax1 = nexttile;
plot(t, L, 'Color', [0.7 0.7 0.7]); hold on;
plot(velocitytable(:,1), velocitytable(:,4), 'ro-', 'MarkerSize', 4);
title('Extracted Velocity Table over Original L trace');
ylabel('Length (Lo)'); xlabel('Time (s)');
legend('Raw L', 'Velocity Table Segments');
ax2 = nexttile;
plot(t, F, 'Color', [0.7 0.7 0.7]); hold on;
plot(datatable(:,1), datatable(:,3), 'r-');
linkaxes([ax1 ax2], 'x')

function keep = rdp(x, y, epsilon)
    % Iterative RDP algorithm
    dmax = 0;
    index = 0;
    end_val = length(x);

    % Line from start to end: a*x + b*y + c = 0
    % Line defined by p1=(x(1),y(1)) and p2=(x(end),y(end))
    dx = x(end) - x(1);
    dy = y(end) - y(1);

    % calculate distance of all points to the line
    % Distance formula: |dy*x0 - dx*y0 + x(end)*y(1) - y(end)*x(1)| / sqrt(dx^2 + dy^2)
    c = x(end)*y(1) - y(end)*x(1);
    norm_factor = sqrt(dx^2 + dy^2);

    if norm_factor == 0
        dists = sqrt((x - x(1)).^2 + (y - y(1)).^2);
    else
        dists = abs(dy.*x - dx.*y + c) / norm_factor;
    end

    [dmax, index] = max(dists);

    if dmax > epsilon
        % Recursive call
        keep1 = rdp(x(1:index), y(1:index), epsilon);
        keep2 = rdp(x(index:end), y(index:end), epsilon);

        % Combine results
        keep = [keep1(1:end-1); keep2 + index - 1];
    else
        keep = [1; end_val];
    end
end
