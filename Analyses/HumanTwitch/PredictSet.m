%% PredictSet.m -- case B prediction validity: predictions of every acceptable tB* fit (profile points, Ca-error refits)
% Acceptable = SSR within the profile-likelihood threshold of the best tB fit (chi2 with 1 dof at 95 %, scaled by
% the best fit's residual variance, i.e. SSR <= SSR_min*(1 + 3.84/(n - p))). For each acceptable fit:
%   muscle-level (humanProbe, saturating Ca, SL 2.2): F0, attached, SRX, ATP use R2T, ktr, a/F0
%   twitch-level (1 Hz, high ATP): F_dev, ATP/beat, and the responses to Ca amplitude x1.5 and /1.5,
%   2 Hz pacing (same transient shape), and a mavacamten-like ka x0.5 -- each as % change of peak force.
% Output: results/predict_caseB.mat, results/predict_caseB.csv, results/predict_caseB.png
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
L = dir(fullfile(res, 'tB*_fit.mat'));  names = erase({L.name}, '_fit.mat');
F = cellfun(@(n) load(fullfile(res, [n '_fit.mat'])), names);
ssr = arrayfun(@(f) f.R.ssr, F);  n = numel(F(1).R.r);  pB = numel(tierConfigs('tB').par);
[ssrMin, ib] = min(ssr);  thr = ssrMin*(1 + 3.84/(n - pB));
acc = ssr <= thr;
fprintf('%d tB fits, best %s SSR %.5f (rms %.3f %%), threshold %.5f -> %d acceptable\n', numel(F), names{ib}, ssrMin, 100*sqrt(ssrMin/n), thr, nnz(acc));
tb = (0:0.0025:1)';
if isempty(gcp('nocreate')), parpool('local', 8); end
P = cell(numel(F), 1);
parfor k = 1:numel(F)
    cfg = F(k).cfg;  x = F(k).R.x;  [base, th] = buildModel(cfg, x);
    q = struct('name', names{k}, 'ssr', ssr(k), 'rms', 100*F(k).R.rms, 'acc', acc(k));
    % which parameter is profiled (if any), its value
    q.prof = '';  q.profVal = NaN;  if isfield(cfg, 'profFixed'), q.prof = cfg.profFixed.par;  q.profVal = 10^cfg.profFixed.value; end
    nm = {cfg.par.name};
    for f = ["r_net" "r_Pi" "r_SRX" "r_kSE" "tf_k_on" "tf_k_off" "tf_K_coop"]
        j = find(strcmp(nm, f), 1);
        if ~isempty(j), v = x(j); if cfg.par(j).log, v = 10^v; end, q.(f) = v;
        elseif strcmp(q.prof, f), q.(f) = q.profVal;
        else, q.(f) = NaN; end
    end
    % muscle level
    hp = humanProbe(base, struct('Tiso', 4, 'Tktr', 2, 'vels', [0.5 2]));   % 2 velocities: FV curvature not used, keeps the probe short
    q.F0 = hp.F0;  q.attached = hp.attached;  q.SRX = hp.SRX;  q.R2T = hp.R2T;  q.ktr = hp.ktr;  q.aF0 = hp.a_over_F0;
    % twitch level
    [fm, S] = simHumanTwitch(base, th, tb);  q.Fdev = max(fm);  q.sigma = F(k).R.sigma;  q.ATPbeat = trapz(tb, S.R2T);
    o = cfg.ca.opt;  thc = th;
    thc.CaOpt = o;  thc.CaOpt.Camax = o.Cadia + 1.5*(o.Camax - o.Cadia);  q.dF_Ca15 = 100*(max(simHumanTwitch(base, thc, tb))/q.Fdev - 1);
    thc.CaOpt.Camax = o.Cadia + (o.Camax - o.Cadia)/1.5;                    q.dF_Ca067 = 100*(max(simHumanTwitch(base, thc, tb))/q.Fdev - 1);
    th2 = th;  th2.T = 0.5;  th2.nb = 4;  f2 = simHumanTwitch(base, th2, tb(tb < 0.5));  q.dF_2Hz = 100*(max(f2)/q.Fdev - 1);
    bm = base;  bm.ka = 0.5*bm.ka;  q.dF_mava = 100*(max(simHumanTwitch(bm, th, tb))/q.Fdev - 1);
    P{k} = q;
end
P = [P{:}];
T = struct2table(P);  disp(T);  writetable(T, fullfile(res, 'predict_caseB.csv'));
save(fullfile(res, 'predict_caseB.mat'), 'P', 'thr', 'ssrMin', 'names');
% prediction intervals over the acceptable set
pv = {'r_net', 'r_SRX', 'r_kSE', 'tf_k_off', 'tf_K_coop', 'F0', 'ktr', 'R2T', 'SRX', 'aF0', 'Fdev', 'sigma', 'ATPbeat', 'dF_Ca15', 'dF_Ca067', 'dF_2Hz', 'dF_mava'};
fprintf('\n%-10s %10s %10s %10s   (acceptable set, n = %d)\n', 'quantity', 'best', 'min', 'max', nnz(acc));
for i = 1:numel(pv)
    v = [P.(pv{i})];  va = v(acc);
    fprintf('%-10s %10.4g %10.4g %10.4g\n', pv{i}, v(ib), min(va), max(va));
end
%% figure: profiles and prediction spread
fig = figure('Color', 'w', 'Position', [20 20 1500 900], 'Visible', 'off');
tl_ = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for pr = {'r_net', 'r_SRX', 'r_kSE'}
    ax = nexttile; hold(ax, 'on');
    m = strcmp({P.prof}, pr{1});  v = [P(m).profVal];  s = 100*sqrt([P(m).ssr]/n);  [v, o] = sort(v);  s = s(o);
    plot(ax, v, s, 'o-', 'LineWidth', 1.4);  plot(ax, P(ib).(pr{1}), 100*sqrt(ssrMin/n), 'rs', 'MarkerFaceColor', 'r');
    yline(ax, 100*sqrt(thr/n), '--', '95 % threshold');  set(ax, 'XScale', 'log');  grid(ax, 'on');
    title(ax, sprintf('profile: %s', strrep(pr{1}, '_', '\_')), 'FontWeight', 'normal');  xlabel(ax, 'fixed value (multiplier)');  ylabel(ax, 'rms (% of peak)');
end
for pr = {'ktr', 'R2T', 'dF_Ca15'}
    ax = nexttile; hold(ax, 'on');
    m = strcmp({P.prof}, 'r_net');  v = [P(m).profVal];  y = [P(m).(pr{1})];  a = [P(m).acc];  [v, o] = sort(v);  y = y(o);  a = a(o);
    plot(ax, v, y, 'o-', 'Color', [.6 .6 .6]);  plot(ax, v(a), y(a), 'o', 'MarkerFaceColor', [0 .45 .75], 'Color', [0 .45 .75]);
    e = strncmp({P.name}, 'tB_caE', 6);  plot(ax, [P(e).r_net], [P(e).(pr{1})], 's', 'MarkerFaceColor', [.85 .33 .1], 'Color', [.85 .33 .1]);
    set(ax, 'XScale', 'log');  grid(ax, 'on');  title(ax, sprintf('prediction: %s along the r\\_net profile (filled: acceptable; squares: Ca error)', strrep(pr{1}, '_', '\_')), 'FontWeight', 'normal');  xlabel(ax, 'r\_net');
end
exportgraphics(fig, fullfile(res, 'predict_caseB.png'), 'Resolution', 130);  close(fig);
fprintf('prediction figure written\n');
