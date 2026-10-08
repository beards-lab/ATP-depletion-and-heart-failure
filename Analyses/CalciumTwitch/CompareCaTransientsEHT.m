%% CompareCaTransientsEHT.m -- (1) audit: which parameters the EHT fits changed; (2) Ca transients
%
% (1) For fits B, C and M, the params struct actually used by the simulator (simEHTTwitch's S.P)
%     is compared field by field with getParams of the unmodified snapshot params/ExampleHighATP.m.
%     Every differing field is listed. Expected: only the Ca layer (CaModel, CaInput, ..., tf_*),
%     the derived Ca struct, SL0 and what SL0 sets (initial state, strain pivot), and the
%     simulator's isometric settings (Vums, v). No cross-bridge parameter.
% (2) Ca transients of fits B, C and M against the two untuned transients:
%       'twoState' as configured before tuning: caTransient 'rice2008' defaults
%       (rat 24 C: 0.09 -> 1.45 uM, tau 20/110 ms, onset 5 ms)
%       MyoSim: the pCa trace of MATMyoSim's simple_twitch demo (generate_twitch_protocol
%       two-compartment SR pulse; stimulus at 0.1 s), loaded from ../MATMyoSim if present
% Output: results/eht_param_audit.txt, results/eht_ca_transients.png
% Report: EHT_identification_report.md, sections 2.2 and 3.6.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
res  = fullfile(here, 'results');
tfine = (0:0.0025:1)';
snap = loadParams(fullfile(root, 'params', 'ExampleHighATP.m'));
Pref = getParams(snap, [], true, false);                       % the snapshot as calibrated (SL 2.2)

%% (1) parameter audit ---------------------------------------------------------------------------
G = load(fullfile(res, 'eht_fit.mat'));  FC = load(fullfile(res, 'eht_fitC.mat'));  FM = load(fullfile(res, 'eht_myosim.mat'));
base = snap;  base.SL0 = 2.0;
fitTh = {'fit B', G.fits(2).th;  'fit C', FC.thC;  'fit M', FM.fits{2}.th};
fid = fopen(fullfile(res, 'eht_param_audit.txt'), 'w');
xbFields = {'kah', 'kamh', 'ka', 'kd', 'k1', 'k_1', 'k2', 'kstiff1', 'kstiff2', 'ksr0', 'kmsrd', 'ksrd', 'k_pas', 'gamma', 'Lsc0', 'kSE', 'A0', 'xrate'};
for k = 1:size(fitTh, 1)
    [~, S] = simEHTTwitch(base, fitTh{k, 2}, tfine(1:10));
    Pf = S.P;
    fn = union(fieldnames(Pref), fieldnames(Pf));
    diffs = {};
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(Pref, f), diffs(end+1, :) = {f, '(absent)', short(Pf.(f))}; continue; end %#ok<SAGROW>
        if ~isfield(Pf, f), diffs(end+1, :) = {f, short(Pref.(f)), '(absent)'}; continue; end %#ok<SAGROW>
        if ~isequaln(Pref.(f), Pf.(f)), diffs(end+1, :) = {f, short(Pref.(f)), short(Pf.(f))}; end %#ok<SAGROW>
    end
    nXB = nnz(ismember(diffs(:, 1), xbFields));
    fprintf(fid, '== %s: %d of %d fields differ from the snapshot; cross-bridge/SRX/passive/SE fields among them: %d\n', ...
        fitTh{k, 1}, size(diffs, 1), numel(fn), nXB);
    fprintf('A %s: %d of %d fields differ; listed XB/SRX/passive/SE fields among them: %d\n', fitTh{k, 1}, size(diffs, 1), numel(fn), nXB);
    for i = 1:size(diffs, 1)
        fprintf(fid, '   %-26s %-40s -> %s\n', diffs{i, :});
        fprintf('     %-24s %-38s -> %s\n', diffs{i, :});
    end
end
fclose(fid);

%% (2) Ca transients -----------------------------------------------------------------------------
tr = struct('name', {}, 't', {}, 'Ca', {}, 'note', {});
tr(end+1) = struct('name', 'fit B (twoState)', 't', tfine, 'Ca', G.fits(2).S.Ca(:), 'note', 'fitted, onset +18.7 ms');
tr(end+1) = struct('name', 'fit C (twoState)', 't', tfine, 'Ca', FC.S.Ca(:), 'note', 'fitted, onset -2.7 ms');
tr(end+1) = struct('name', 'fit M (myosim)', 't', tfine, 'Ca', FM.fits{2}.S.Ca(:), 'note', 'fitted, onset -3.9 ms');
[~, ro] = caTransient('rice2008', 0, struct('T', 1));
tr(end+1) = struct('name', 'original twoState input: Rice 2008 defaults', 't', tfine, 'Ca', caTransient('rice2008', tfine, ro), 'note', 'untuned, onset +5 ms');
pf = fullfile(root, '..', 'MATMyoSim', 'code', 'demos', 'getting_started', 'simple_twitch', 'protocol.txt');
if exist(pf, 'file')
    Tm = readtable(pf, 'FileType', 'text', 'Delimiter', '\t');
    tm = [0; cumsum(Tm.dt(1:end-1))] - 0.1;                         % stimulus at 0.1 s -> t = 0
    tr(end+1) = struct('name', 'original MyoSim demo (simple\_twitch)', 't', tm, 'Ca', 10.^(6 - Tm.pCa), 'note', 'untuned, stimulus at 0');
end
for k = 1:numel(tr)
    c0 = min(tr(k).Ca);                                         % diastolic level (fits C, M rise just before t = 0)
    t = tr(k).t;  c = tr(k).Ca;  w = t >= 0 & t <= 1;  t = t(w);  c = c(w);
    [cm, ip] = max(c);  half = c0 + 0.5*(cm - c0);
    j = find(c(ip:end) <= half, 1) + ip - 1;
    tr(k).dia = c0;  tr(k).pk = cm;  tr(k).ttp = t(ip);
    if isempty(j), tr(k).t50 = NaN; else, tr(k).t50 = t(j) - t(ip); end
    fprintf('C %-44s diastolic %.3f uM, peak %.3f uM at %.0f ms, peak -> 50 %% decay %.0f ms (%s)\n', ...
        strrep(tr(k).name, '\', ''), c0, cm, 1e3*t(ip), 1e3*tr(k).t50, tr(k).note);
end

cm_ = [42 120 214; 0 0 0; 235 104 52; 150 150 150; 27 175 122]/255;  ls = {'-', '-', '--', ':', '-.'};
fig = figure('Color', 'w', 'Position', [40 40 1400 480], 'Visible', 'off');
tl_ = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_, 'Ca transients: fitted (B, C, M) and the untuned inputs (Rice 2008 default; MyoSim demo)', 'FontWeight', 'normal');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(tr), plot(ax, 1e3*tr(k).t, tr(k).Ca, ls{k}, 'Color', cm_(k,:), 'LineWidth', 1.8, 'DisplayName', tr(k).name); end
xlim(ax, [-20 800]);  style(ax, 'A  Absolute [Ca^{2+}]_i (fitted levels are a fixed gauge, 0.10 \rightarrow 1.0 \muM)', '[Ca^{2+}]_i (\muM)');
legend(ax, 'Box', 'off', 'FontSize', 8);  xlabel(ax, 'time from stimulus (ms)');
ax = nexttile; hold(ax, 'on');
for k = 1:numel(tr)
    t = tr(k).t;  c = tr(k).Ca;  w = t >= -0.02;
    plot(ax, 1e3*t(w), (c(w) - tr(k).dia)/(tr(k).pk - tr(k).dia), ls{k}, 'Color', cm_(k,:), 'LineWidth', 1.8);
end
xlim(ax, [-20 800]);  ylim(ax, [-0.05 1.05]);
style(ax, 'B  Normalised shape: what the twitch fit determines', '(Ca - Ca_{dia}) / (Ca_{peak} - Ca_{dia})');  xlabel(ax, 'time from stimulus (ms)');
exportgraphics(fig, fullfile(res, 'eht_ca_transients.png'), 'Resolution', 150);
save(fullfile(res, 'eht_ca_transients.mat'), 'tr');
fprintf('figure written\n');

function s = short(v)
if isnumeric(v) || islogical(v)
    if isscalar(v), s = sprintf('%.6g', v); else, s = sprintf('[%s, %d values]', class(v), numel(v)); end
elseif ischar(v) || isstring(v), s = ['''' char(v) ''''];
elseif isstruct(v), s = '(struct)';
else, s = ['(' class(v) ')'];
end
end

function style(ax, ttl, yl)
title(ax, ttl, 'FontWeight', 'normal');  ylabel(ax, yl);  box(ax, 'off');  grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'XColor', [0.32 0.32 0.31], 'YColor', [0.32 0.32 0.31], 'TickDir', 'out');
end
