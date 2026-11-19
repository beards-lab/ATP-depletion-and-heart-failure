% RunProtocolByATP  Is the 4 and 2 mM isovelocity protocol really different from
% the 8 mM one, in every record and both running orders, and is it an Lo
% misalignment (a constant length offset) or a different programmed step?  (#75)
%
% Data only, no simulation. Reads the cached records (loadBeard2022FV), measures
% per record and cycle the length before the step, after it, at the slack and
% the timing, prints the tables to results/protocol_by_atp.log and draws
% figures/protocol_by_atp.png.

here = fileparts(mfilename('fullpath'));
addpath(here);
D = loadBeard2022FV();
logf = fullfile(here, 'results', 'protocol_by_atp.log');
if isfile(logf), delete(logf); end
diary(logf);

vNom = [6 1 3 5 0.5 4 2];                 % recorded order of the cycles
nC = numel(vNom); nR = numel(D);
[hold0, Lpre, Lpost, Lmin, tOn, tRead, tRamp, vMeas, Lread] = deal(nan(nR, nC));
[Lfirst, Lend, sdHold] = deal(nan(nR, 1));
for k = 1:nR
    t = D(k).t; L = D(k).L;
    R = readIsoV(t, L, D(k).F);
    assert(height(R) >= nC, 'record %d: %d cycles', k, height(R));
    Lfirst(k) = median(L(t < 400)); sdHold(k) = std(L(t < 400));
    Lend(k) = median(L(t > t(end) - 200));
    for c = 1:nC
        t0 = R.tStart(c);
        Lpre(k, c)  = median(L(t >= t0 - 20 & t < t0 - 1));        % hold before the step
        Lpost(k, c) = mean(L(t >= t0 + 1.25 & t <= t0 + 1.75));    % 1.5 ms after onset
        seg = t >= t0 & t <= R.tSlack(c) + 60;
        Lmin(k, c)  = min(L(seg));
        tOn(k, c) = t0; tRead(k, c) = R.tCross(c) - t0; tRamp(k, c) = R.tSlack(c) - t0;
        vMeas(k, c) = R.vMeas(c); Lread(k, c) = R.Lx(c);
    end
end
[~, ord] = sort(vNom); vS = vNom(ord);
atp = [D.atp]'; slot = [D.slot]'; lab = strcat({D.date}', " ", {D.order}', " ", string(atp), " mM, slot ", string(slot));

fprintf('\n== 1. Hold lengths per record (Lo): first hold, sd within it, hold before each cycle (min..max), final v=0 hold, slack minimum (min..max)\n');
T1 = table(lab, Lfirst, sdHold, min(Lpre, [], 2), max(Lpre, [], 2), Lend, min(Lmin, [], 2), max(Lmin, [], 2), ...
    'VariableNames', {'record', 'hold_first', 'sd', 'hold_min', 'hold_max', 'hold_v0', 'slack_min', 'slack_max'});
disp(T1);

fprintf('\n== 2. Step size per record (hold before the step minus length 1.5 ms after onset, Lo), columns v = %s ML/s\n', mat2str(vS));
step = Lpre - Lpost;
disp(array2table(round(step(:, ord), 4), 'RowNames', cellstr(lab), 'VariableNames', compose('v%g', vS)));

fprintf('\n== 3. By [ATP] and slot: step size mean and full range over records (Lo)\n');
for a = [8 4 2]
    for s = 1:2
        m = atp == a & slot == s; if ~any(m), continue; end
        x = step(m, ord);
        fprintf('%d mM slot %d (n=%d): mean %s\n                    range %s\n', a, s, sum(m), ...
            mat2str(round(mean(x, 1), 4)), mat2str(round(max(x, [], 1) - min(x, [], 1), 4)));
    end
end

fprintf('\n== 4. Within one preparation: low-ATP record minus its 8 mM record\n');
preps = unique({D.prep}, 'stable');
dStep = nan(numel(preps), nC); dHold = dStep; dRead = dStep; dRamp = dStep; dLread = dStep; dMin = dStep; dOn = dStep;
plab = strings(numel(preps), 1);
for p = 1:numel(preps)
    i8 = find(strcmp({D.prep}, preps{p}) & [D.atp] == 8); iL = find(strcmp({D.prep}, preps{p}) & [D.atp] ~= 8);
    plab(p) = sprintf('%s %s', D(iL).date, D(iL).order);
    dStep(p, :) = step(iL, ord) - step(i8, ord);   dHold(p, :) = Lpre(iL, ord) - Lpre(i8, ord);
    dRead(p, :) = tRead(iL, ord) - tRead(i8, ord); dRamp(p, :) = tRamp(iL, ord) - tRamp(i8, ord);
    dLread(p, :) = Lread(iL, ord) - Lread(i8, ord); dMin(p, :) = Lmin(iL, ord) - Lmin(i8, ord);
    dOn(p, :) = tOn(iL, ord) - tOn(i8, ord);
end
vn = compose('v%g', vS);
fprintf('step size difference (Lo):\n');            disp(array2table(round(dStep, 4), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('hold-before-step difference (Lo):\n');     disp(array2table(round(dHold, 4), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('slack-minimum difference (Lo):\n');        disp(array2table(round(dMin, 4), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('length-at-read difference (Lo):\n');       disp(array2table(round(dLread, 4), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('cycle-onset time difference (ms):\n');     disp(array2table(round(dOn, 2), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('onset-to-read difference (ms):\n');        disp(array2table(round(dRead, 2), 'RowNames', cellstr(plab), 'VariableNames', vn));
fprintf('onset-to-slack (ramp) difference (ms):\n'); disp(array2table(round(dRamp, 2), 'RowNames', cellstr(plab), 'VariableNames', vn));

fprintf('\n== 5. Ramp design by [ATP] (medians): step, ramp duration, speed, and length at the slack release = hold - step - v*duration\n');
for a = [8 4 2]
    m = atp == a;
    st = median(step(m, ord), 1); du = median(tRamp(m, ord), 1); vv = median(vMeas(m, ord), 1);
    fprintf('%d mM: step   %s\n      dur ms %s\n      v ML/s %s\n      L end  %s\n', a, mat2str(round(st, 3)), ...
        mat2str(round(du, 1)), mat2str(round(vv, 3)), mat2str(round(median(Lpre(m, ord), 1) - st - vv .* du / 1000, 3)));
end

fprintf('\n== 6. The lab''s own summary columns: iso-V Length per v, by record (Lo)\n');
for k = 1:nR
    S = sortrows(D(k).summary, 'v');
    fprintf('%-34s v %s  Liso %s\n', lab(k), mat2str(S.v'), mat2str(round(S.Liso', 4)));
end
diary off;

%% Figure: commanded length around each cycle onset, all 24 records
col = containers.Map({8, 4, 2}, {[0 0 0], [0.85 0.33 0.10], [0 0.45 0.74]});
fig = figure('Visible', 'off', 'Position', [50 50 1500 900], 'Color', 'w');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:nC
    c = ord(j); nexttile; hold on;
    span = median(tRamp(:, c), 'omitnan');
    for k = 1:nR
        t = D(k).t - tOn(k, c); w = t > -5 & t < span + 8;
        ls = '-'; if D(k).slot == 2, ls = '--'; end
        plot(t(w), D(k).L(w), ls, 'Color', [col(D(k).atp) 0.6], 'LineWidth', 0.8);
    end
    yline(1.00, ':'); ylim([0.97 1.105]); xlim([-5 span + 8]);
    title(sprintf('%g ML/s', vS(j))); xlabel('time from cycle onset (ms)'); ylabel('length (L_o)');
end
% zoom on the step, 1 ML/s
nexttile; hold on; c = find(vNom == 1);
for k = 1:nR
    t = D(k).t - tOn(k, c); w = t > -3 & t < 6;
    ls = '-'; if D(k).slot == 2, ls = '--'; end
    plot(t(w), D(k).L(w), ls, 'Color', [col(D(k).atp) 0.6], 'LineWidth', 0.8);
end
title('1 ML/s, the step'); xlabel('time from cycle onset (ms)'); ylabel('length (L_o)'); ylim([1.06 1.105]);
% step size by v
nexttile([1 2]); hold on;
for a = [8 4 2]
    for s = 1:2
        m = find(atp == a & slot == s); mk = 'o'; if s == 2, mk = 's'; end
        for k = m'
            plot(vS + (a == 4) * 0.06 - (a == 2) * 0.06, step(k, ord), mk, 'Color', col(a), 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
        plot(nan, nan, mk, 'Color', col(a));   % legend entry
    end
end
legend({'8 mM run first', '8 mM run second', '4 mM run first', '4 mM run second', '2 mM run first', '2 mM run second'}, 'Location', 'northwest', 'Box', 'off');
xlabel('ramp speed (ML/s)'); ylabel('pre-ramp step (L_o)'); title('Step size, every record'); grid on;
% hold before the step: is there an Lo offset?
nexttile([1 2]); hold on;
for k = 1:nR
    mk = 'o'; if D(k).slot == 2, mk = 's'; end
    plot(k + zeros(1, nC), Lpre(k, :), mk, 'Color', col(D(k).atp), 'MarkerSize', 4);
end
yline(1.10, ':'); ylim([1.09 1.11]); xlim([0 nR + 1]);
xlabel('record (two per preparation)'); ylabel('hold before each step (L_o)'); title('Hold length before the 7 steps, every record'); grid on;
title(tl, 'Beard2022 isovelocity: length around each cycle. Black 8 mM, orange 4 mM, blue 2 mM; solid = run first, dashed = run second');
exportgraphics(fig, fullfile(here, 'figures', 'protocol_by_atp.png'), 'Resolution', 150);
close(fig);
fprintf('done\n');
