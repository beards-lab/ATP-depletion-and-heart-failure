% RunFVBeard2022Analysis.m -- score the incumbent against the Beard/Baker 2022
% isovelocity dataset: F-V points (8/4/2 mM), F-V shape (Vmax, curvature),
% within-preparation ATP ratios, and the raw traces.
%
% Needs results/sims.mat (RunFVBeard2022.m) and results/steady.mat
% (RunFVBeard2022Steady.m). Seconds to run.
%
% Normalisation (issue #75): ONE force scale per preparation, set so the model's
% 8 mM isometric force at 1.00 Lo equals that preparation's measured 8 mM
% isometric force at 1.00 Lo, and carried unchanged to the same preparation's
% 4 or 2 mM record. Everything is then expressed in units of that reference
% force Fref (data: F/Fref; model: s*F/Fref = F / F_model(8 mM, v=0)), so the
% ATP force ratio stays a prediction.
%
% Both data and model are read by readIsoV on the DATA's motor-length timing
% (the force window opens when the measured length crosses 1.00 Lo), so reader
% bias is common to both.
%
% Output -> results/metrics.mat, figures/*.png, and the tables printed below.

rootP = 'C:\home\git\ATP-depletion-and-heart-failure';
aP = fullfile(rootP, 'Analyses', 'FV_Beard2022'); addpath(aP);
figP = fullfile(aP, 'figures'); if ~isfolder(figP), mkdir(figP); end
D = loadBeard2022FV();
Sm = load(fullfile(aP, 'results', 'sims.mat')); S = Sm.S;
St = load(fullfile(aP, 'results', 'steady.mat'));
assert(all([S.ok]), 'some simulations failed');
F_of = @(k) S([S.k] == k).F;

ATP = [8 4 2]; col = [0 0 0; 0.85 0.45 0; 0.75 0 0];
nS = numel(D);
%% ---- 1. read every record, data and model, on the data's timing --------------
for k = 1:nS
    D(k).Rd  = readIsoV(D(k).t, D(k).L, D(k).F);
    D(k).Rm  = readIsoV(D(k).t, D(k).L, F_of(k));
    D(k).Fm  = F_of(k);
end
% reader check against the lab's own iso-V column
dLab = cell2mat(arrayfun(@(d) d.Rd.Fx - d.summary.Fiso_v, D(:), 'uni', 0));
fprintf('\nReader vs workbook iso-V force (n = %d): median |d| %.2f kPa, 95%% %.2f, max %.2f\n', ...
    numel(dLab), median(abs(dLab)), prctile(abs(dLab), 95), max(abs(dLab)));

%% ---- 2. per-preparation scale (8 mM isometric at 1.00 Lo) ---------------------
preps = unique({D.prep}, 'stable');
for ip = 1:numel(preps)
    ks = find(strcmp({D.prep}, preps{ip}));
    k8 = ks([D(ks).atp] == 8);
    Fref = D(k8).Rd.Fx(D(k8).Rd.vMeas == 0);
    s    = Fref / D(k8).Rm.Fx(D(k8).Rm.vMeas == 0 | isnan(D(k8).Rm.tStart));
    for k = ks, D(k).Fref = Fref; D(k).scale = s; end
end
v0 = D(1).summary.v; [vS, iv] = sort(v0);           % nominal velocities, sorted
nV = numel(vS);
P  = @(R, k) R.Fx(iv) / D(k).Fref;                  % sorted, in Fref units
fprintf('Force scale s = F_data/F_model at 8 mM, 1.00 Lo: median %.3f (range %.3f-%.3f), Fref %.1f +- %.1f kPa\n', ...
    median([D.scale]), min([D.scale]), max([D.scale]), mean(unique([D.Fref])), std(unique([D.Fref])));

%% ---- 3. F-V points per ATP: pooled table --------------------------------------
T = table(); Pd = struct(); Pm = struct();
for a = ATP
    ks = find([D.atp] == a); f = sprintf('a%d', a);
    Pd.(f) = cell2mat(arrayfun(@(k) P(D(k).Rd, k)', ks(:), 'uni', 0));               % data
    Pm.(f) = cell2mat(arrayfun(@(k) (D(k).scale * D(k).Rm.Fx(iv) / D(k).Fref)', ks(:), 'uni', 0));
    for j = 1:nV
        d = Pd.(f)(:, j); m = Pm.(f)(:, j);
        T = [T; table(a, vS(j), numel(d), mean(d), std(d), mean(m), mean(m) - mean(d), ...
            (mean(m) - mean(d)) / std(d), 'VariableNames', ...
            {'ATP', 'v', 'n', 'data', 'sd', 'model', 'err', 'z'})]; %#ok<AGROW>
    end
end
fprintf('\nF-V points in units of the preparation''s 8 mM isometric force (model scaled once per prep):\n');
disp(T);
for a = ATP
    r = T(T.ATP == a & T.v > 0, :);
    fprintf('%d mM: RMS error over v>0 = %.3f Fref; mean |z| = %.2f; %d/%d points within 1 sd, %d within 2 sd\n', ...
        a, rms(r.err), mean(abs(r.z)), sum(abs(r.z) <= 1), height(r), sum(abs(r.z) <= 2));
end

%% ---- 4. F-V shape: Hill fits per record (P/P0 within the record) -------------
hill = @(x, v) max((1 + x(1)) * x(2) ./ (v + x(2)) - x(1), 0);   % x = [a/F0, b]; Vmax = b/(a/F0)
H = table();
for k = 1:nS
    for src = ["data", "model"]
        if src == "data", R = D(k).Rd; else, R = D(k).Rm; end
        p0 = R.Fx(R.vMeas == 0); m = R.vMeas > 0;
        sc = 1; if src == "model", sc = D(k).scale; end   % P0 column in Fref units
        v = R.vMeas(m); p = R.Fx(m) / p0;
        cost = @(lx) sum((hill(exp(lx), v) - p).^2);
        lx = fminsearch(cost, log([0.3 1.5]), optimset('TolX', 1e-8, 'TolFun', 1e-10, 'MaxFunEvals', 4000));
        x = exp(lx);
        pw = max(p .* v);                                  % max power, F0*ML/s units
        H = [H; table(k, D(k).atp, string(src), x(1), x(2) / x(1), pw, sc * p0 / D(k).Fref, ...
            'VariableNames', {'k', 'ATP', 'src', 'curv', 'Vmax', 'Pmax', 'P0'})]; %#ok<AGROW>
    end
end
fprintf('\nHill fits per record (P/P0, v>0): mean +- sd over preparations\n');
HS = table();
for a = ATP
    for src = ["data", "model"]
        h = H(H.ATP == a & H.src == src, :);
        % medians + IQR: single-record Hill fits on 7 noisy points have outliers
        HS = [HS; table(a, src, height(h), median(h.curv), iqr(h.curv), median(h.Vmax), iqr(h.Vmax), ...
            median(h.Pmax), iqr(h.Pmax), mean(h.P0), std(h.P0), 'VariableNames', ...
            {'ATP', 'src', 'n', 'curv_med', 'curv_iqr', 'Vmax_med', 'Vmax_iqr', 'Pmax_med', 'Pmax_iqr', 'P0', 'P0_sd'})]; %#ok<AGROW>
    end
end
disp(HS);

%% ---- 5. within-preparation ATP ratios (paired) --------------------------------
RT = table();
for a = [4 2]
    pr = unique({D([D.atp] == a).prep});
    rd = nan(numel(pr), nV); rm = rd;
    for ip = 1:numel(pr)
        k8 = find(strcmp({D.prep}, pr{ip}) & [D.atp] == 8);
        ka = find(strcmp({D.prep}, pr{ip}) & [D.atp] == a);
        rd(ip, :) = D(ka).Rd.Fx(iv)' ./ D(k8).Rd.Fx(iv)';
        rm(ip, :) = D(ka).Rm.Fx(iv)' ./ D(k8).Rm.Fx(iv)';
        % a ratio involving a slack model record (F < 1 % of its own F0) is undefined
        f0 = @(k) D(k).Rm.Fx(D(k).Rm.vMeas == 0);
        slack = D(ka).Rm.Fx(iv)' < 0.01 * f0(ka) | D(k8).Rm.Fx(iv)' < 0.01 * f0(k8);
        rm(ip, slack) = NaN;
    end
    for j = 1:nV
        RT = [RT; table(a, vS(j), numel(pr), mean(rd(:, j)), std(rd(:, j)), mean(rm(:, j), 'omitnan'), ...
            'VariableNames', {'ATP', 'v', 'n', 'data', 'sd', 'model'})]; %#ok<AGROW>
    end
    RTr.(sprintf('a%d', a)) = struct('rd', rd, 'rm', rm, 'preps', {pr});
end
fprintf('\nWithin-preparation force ratio [ATP]/8 mM at each velocity (model ratio is scale-free):\n');
disp(RT);
% rundown check: the same ratio split by which [ATP] ran first
fprintf('Ratio by running order (mean over 3 preps each; rundown depresses whichever ran second):\n');
for a = [4 2]
    r = RTr.(sprintf('a%d', a));
    first8 = cellfun(@(p) D(strcmp({D.prep}, p) & [D.atp] == 8).slot == 1, r.preps);
    for j = find(ismember(vS, [0 1 2 3]))'
        fprintf('  %d/8 at %g ML/s: 8 first %.3f | %d first %.3f | mean %.3f | model %.3f\n', a, vS(j), ...
            mean(r.rd(first8, j)), a, mean(r.rd(~first8, j)), mean(r.rd(:, j)), mean(r.rm(:, j)));
    end
end

%% ---- 6. traces: ramp-window RMSE and the between-prep envelope ---------------
tq = D(1).t; sm = @(x) movmean(x, 21);               % 1 ms smoothing (white noise ~1.3 kPa)
Rd1 = D(1).Rd;                                        % cycle timing (aligned to <0.7 ms)
inRamp = false(size(tq));
for i = 1:height(Rd1)
    if isnan(Rd1.tStart(i)), continue; end
    inRamp = inRamp | (tq >= Rd1.tStart(i) - 2 & tq <= Rd1.tSlack(i));
end
TR = table(); env = struct();
for a = ATP
    ks = find([D.atp] == a); f = sprintf('a%d', a);
    Xd = cell2mat(arrayfun(@(k) sm(D(k).F) / D(k).Fref, ks, 'uni', 0));
    Xm = cell2mat(arrayfun(@(k) D(k).scale * D(k).Fm / D(k).Fref, ks, 'uni', 0));
    mu = mean(Xd, 2); sd = std(Xd, 0, 2); mm = mean(Xm, 2);
    env.(f) = struct('mu', mu, 'sd', sd, 'mm', mm, 'Xd', Xd, 'Xm', Xm, 'ks', ks);
    rmsePrep = sqrt(mean((Xm(inRamp, :) - Xd(inRamp, :)).^2, 1));
    hold_ = tq > 3800;                                % final isometric hold at 1.00 Lo
    TR = [TR; table(a, numel(ks), mean(rmsePrep), std(rmsePrep), ...
        mean(abs(mm(inRamp) - mu(inRamp)) <= sd(inRamp)), mean(abs(mm(inRamp) - mu(inRamp)) <= 2 * sd(inRamp)), ...
        mean(sd(inRamp)), sqrt(mean((mm(hold_) - mu(hold_)).^2)), ...
        'VariableNames', {'ATP', 'n', 'rampRMSE', 'rampRMSE_sd', 'in1sd', 'in2sd', 'envSD', 'holdRMSE'})]; %#ok<AGROW>
end
fprintf('\nTrace comparison over the ramp windows (units Fref; data smoothed 1 ms):\n');
disp(TR);

% force redevelopment after each slack + re-stretch (the record's other half):
% time from the end of the re-stretch to 50 % and 90 % of the following hold
Lref = D(1).L; Lh = median(Lref(tq < 400));
RD = table();
for a = ATP
    E = env.(sprintf('a%d', a));
    for i = 1:height(Rd1) - 1                         % cycles 1..7 (the last row is v = 0)
        if isnan(Rd1.tSlack(i)), continue; end
        isLast = i == height(Rd1) - 1;                % cycle 7 re-stretches to 1.00 Lo
        Ltarget = Lh; if isLast, Ltarget = median(Lref(tq > 4100)); end
        iR = find(tq > Rd1.tSlack(i) + 2 & abs(Lref - Ltarget) < 0.002, 1);
        if isLast, wP = tq > 4100; else, wP = tq >= Rd1.tStart(i+1) - 60 & tq <= Rd1.tStart(i+1) - 5; end
        % force at fixed delays after the end of the re-stretch, as a fraction of
        % the hold it recovers to (the data dip and recover; the model overshoots
        % and settles, so a threshold-crossing time is not defined for both)
        x = [E.mu, E.mm]; ff = nan(2, 3); dl = [5 20 50];
        for c = 1:2
            pl = mean(x(wP, c)); y = movmean(x(:, c), 41);
            ff(c, :) = interp1(tq, y, tq(iR) + dl) / pl;
        end
        RD = [RD; table(a, v0(i), ff(1, 1), ff(2, 1), ff(1, 2), ff(2, 2), ff(1, 3), ff(2, 3), ...
            'VariableNames', {'ATP', 'v', 'f5_data', 'f5_model', 'f20_data', 'f20_model', 'f50_data', 'f50_model'})]; %#ok<AGROW>
    end
end
fprintf('\nForce 5/20/50 ms after each re-stretch, fraction of the following hold (pooled mean traces, median over 7 cycles):\n');
RDs = groupsummary(RD, 'ATP', 'median', {'f5_data', 'f5_model', 'f20_data', 'f20_model', 'f50_data', 'f50_model'});
disp(RDs);

%% ---- 7. figures -------------------------------------------------------------
set(0, 'DefaultAxesFontSize', 11);
% (a) F-V points + steady model curve, absolute (Fref units)
fig = figure('Visible', 'off', 'Position', [50 50 1200 450]); tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
nexttile(tl); hold on;
s8 = St.Fs(1, 1);                                     % steady model normalised by its own 8 mM F0 at SL 2.0
for ia = 1:3
    a = ATP(ia); r = T(T.ATP == a, :);
    errorbar(r.v, r.data, r.sd, 'o', 'Color', col(ia, :), 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', sprintf('data %d mM (n=%d)', a, r.n(1)));
    plot(r.v, r.model, 'x', 'Color', col(ia, :), 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', sprintf('model, measured protocol %d mM', a));
    plot(St.vFine, St.Fs(ia, :) / s8, '-', 'Color', col(ia, :), 'LineWidth', 1, 'DisplayName', sprintf('model, step-free ramp %d mM', a));
end
xlabel('shortening velocity (ML/s)'); ylabel('F / F_{ref} (8 mM isometric, 1.00 Lo)'); xlim([0 7]); ylim([0 1.5]);
title('Force-velocity, one force scale per preparation'); legend('Location', 'northeast', 'FontSize', 8); box on; grid on;
% (b) shape: P/P0 per ATP
nexttile(tl); hold on;
for ia = 1:3
    a = ATP(ia); f = sprintf('a%d', a);
    pd = Pd.(f) ./ Pd.(f)(:, 1); pm = Pm.(f) ./ Pm.(f)(:, 1);
    errorbar(vS, mean(pd), std(pd), 'o', 'Color', col(ia, :), 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', sprintf('data %d mM', a));
    plot(vS, mean(pm), 'x', 'Color', col(ia, :), 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', sprintf('model (measured protocol) %d mM', a));
    plot(St.vFine, St.Fs(ia, :) / St.Fs(ia, 1), '-', 'Color', col(ia, :), 'LineWidth', 1, 'DisplayName', sprintf('model step-free %d mM', a));
end
xlabel('shortening velocity (ML/s)'); ylabel('P / P_0 (own record)'); xlim([0 7]); ylim([0 1.3]);
title('Force-velocity shape'); legend('Location', 'northeast', 'FontSize', 8); box on; grid on;
exportgraphics(fig, fullfile(figP, 'fv_points.png'), 'Resolution', 150); close(fig);

% (c) within-prep ratios
fig = figure('Visible', 'off', 'Position', [50 50 650 420]); hold on;
for a = [4 2]
    ia = find(ATP == a); r = RT(RT.ATP == a, :);
    errorbar(r.v, r.data, r.sd, 'o', 'Color', col(ia, :), 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'DisplayName', sprintf('data %d/8 mM (n=%d preps)', a, r.n(1)));
    plot(r.v, r.model, 'x-', 'Color', col(ia, :), 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', sprintf('model %d/8 mM, measured protocol (gaps = model slack)', a));
    plot(St.vFine, St.Fs(ia, :) ./ St.Fs(1, :), '--', 'Color', col(ia, :), 'DisplayName', sprintf('model %d/8 mM, step-free', a));
end
yline(1, 'k:', 'HandleVisibility', 'off'); xlabel('shortening velocity (ML/s)'); ylabel('F([ATP]) / F(8 mM), same preparation');
xlim([0 7]); ylim([0 2]); legend('Location', 'northwest', 'FontSize', 8); box on; grid on;
title('ATP effect at each velocity (paired within preparation)');
exportgraphics(fig, fullfile(figP, 'fv_ratio.png'), 'Resolution', 150); close(fig);

% (d) pooled traces per ATP: full record + per-cycle zooms
order = [5 2 7 3 6 4 1];                              % cycles sorted by speed 0.5,1,2,3,4,5,6
for ia = 1:3
    a = ATP(ia); f = sprintf('a%d', a); E = env.(f);
    fig = figure('Visible', 'off', 'Position', [50 50 1400 800]);
    tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, [1 4]); hold on;
    fill([tq; flipud(tq)], [E.mu - E.sd; flipud(E.mu + E.sd)], col(ia, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(tq, E.mu, '-', 'Color', col(ia, :), 'LineWidth', 1);
    plot(tq, E.mm, 'b-', 'LineWidth', 1);
    plot(tq, D(E.ks(1)).L * 1 - 0.1, 'k-', 'LineWidth', 0.5);
    xlabel('time (ms)'); ylabel('F / F_{ref}'); xlim([0 4300]); ylim([-0.15 1.9]);
    title(sprintf('%d mM: data mean \\pm sd (n=%d), model (blue); black = motor length - 0.1 (Lo)', a, numel(E.ks)));
    box on;
    for ic = 1:7
        i = order(ic); nexttile(tl); hold on;
        w = tq >= Rd1.tStart(i) - 5 & tq <= Rd1.tSlack(i) + 3;
        fill([tq(w); flipud(tq(w))], [E.mu(w) - E.sd(w); flipud(E.mu(w) + E.sd(w))], col(ia, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(tq(w), E.mu(w), '-', 'Color', col(ia, :), 'LineWidth', 1.2);
        plot(tq(w), E.mm(w), 'b-', 'LineWidth', 1.2);
        xline(Rd1.tCross(i), 'k:');
        title(sprintf('%g ML/s', v0(i))); xlabel('time (ms)'); ylabel('F / F_{ref}'); box on;
        xlim([Rd1.tStart(i) - 5, Rd1.tSlack(i) + 3]); ylim([-0.1 1.4]);
    end
    nexttile(tl); hold on;                            % the step-to-1.00 isometric hold
    w = tq >= 3630 & tq <= 4300;
    fill([tq(w); flipud(tq(w))], [E.mu(w) - E.sd(w); flipud(E.mu(w) + E.sd(w))], col(ia, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(tq(w), E.mu(w), '-', 'Color', col(ia, :)); plot(tq(w), E.mm(w), 'b-');
    title('v = 0 (hold at 1.00 Lo)'); xlabel('time (ms)'); box on; ylim([-0.1 1.4]);
    exportgraphics(fig, fullfile(figP, sprintf('traces_%dmM.png', a)), 'Resolution', 130); close(fig);
end

% (e) appendix: every record, model vs data
fig = figure('Visible', 'off', 'Position', [50 50 1600 1500]);
tl = tiledlayout(fig, 6, 4, 'TileSpacing', 'tight', 'Padding', 'compact');
for k = 1:nS
    nexttile(tl); hold on;
    plot(tq, sm(D(k).F) / D(k).Fref, 'Color', [0.5 0.5 0.5]);
    plot(tq, D(k).scale * D(k).Fm / D(k).Fref, 'b-');
    title(sprintf('%s  %d mM (slot %d)', D(k).date, D(k).atp, D(k).slot), 'FontSize', 9);
    xlim([0 4300]); ylim([-0.15 1.6]); box on;
end
exportgraphics(fig, fullfile(figP, 'appendix_all_records.png'), 'Resolution', 110); close(fig);

save(fullfile(aP, 'results', 'metrics.mat'), 'T', 'H', 'HS', 'RT', 'TR', 'RD', 'RDs', 'dLab', 'Pd', 'Pm');
fprintf('\nfigures -> %s\n', figP);
