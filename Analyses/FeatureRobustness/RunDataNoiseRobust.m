%% RunDataNoiseRobust.m -- RunDataNoise's surrogates, every definition in robustSlackFeatures.
%
% Regenerates the identical surrogate recordings (same smooth part, same
% residual pool, same seeds) and reads them with robustSlackFeatures only --
% its current-definition replicas are bit-identical to extractSlackAttributes,
% so they double as a check against RunDataNoise. Adds, per feature:
%   s_phase : SD over the 4 phases of the record decimated to 1 kHz
%   b_noise : (mean over surrogates - value on the smooth trace) / |mean| --
%             the BIAS the transducer noise puts into the data value, which a
%             noise-free model trace can never reproduce
% t0 (Fit-1 onset, the data side of t0_crossing) comes from RunDataNoise.
if ~exist('NB', 'var'), NB = 40; end
if ~exist('NW', 'var'), NW = 4; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root)); res = fullfile(here, 'results');
DAYS = {'03_27', '04_03', '04_10'}; COND = {'8mM', '2mM'};
SEG14 = {'rsK80','rsKfit','rsT20','rsTfit','rsT20e','rsTfitE','rsKls','rsT0ls'};
if isempty(gcp('nocreate')), parpool('Threads', NW); end
RUN = struct();
for d = 1:numel(DAYS)
    for c = 1:numel(COND)
        S = load(fullfile(root, 'data', sprintf('protocol_%s_2026_%s_slack.mat', DAYS{d}, COND{c})));
        t = S.datatable(:,1); SL = S.datatable(:,2); y = S.datatable(:,3); vt = S.velocitytable;
        dt = median(diff(t));
        segs = find(vt(:,2) < -1)'; rq = {};
        for s = segs
            if s + 4 > size(vt,1), break; end
            m = t >= vt(s+4,1) - 0.15 & t <= vt(s+4,1) - 0.005;
            tq = t(m); yq = y(m); rq{end+1} = yq - polyval(polyfit(tq - tq(1), yq, 2), tq - tq(1)); %#ok<SAGROW>
        end
        rpool = vertcat(rq{:});
        ysm = sgolayfilt(y, 3, 9);
        g = var(sgolayfilt(rpool, 3, 9)) / var(rpool); cfac = sqrt(max(0, 1 - g));
        BL = round(0.010 / dt); N = numel(y);
        R = struct('day', DAYS{d}, 'cond', COND{c});
        R.base  = robustSlackFeatures(t, y, SL, vt, true);
        R.clean = robustSlackFeatures(t, ysm, SL, vt, true);
        Xr = cell(NB, 1);
        parfor b = 1:NB
            st = RandStream('mt19937ar', 'Seed', 1000*d + 100*c + b);
            nz = zeros(N,1); k = 0;
            while k < N
                j0 = randi(st, numel(rpool) - BL); L = min(BL, N - k);
                nz(k+1:k+L) = rpool(j0:j0+L-1); k = k + L;
            end
            Xr{b} = robustSlackFeatures(t, ysm + cfac * nz, SL, vt, true);
        end
        R.Xr = Xr;
        q = round(1e-3 / dt); R.Xp = cell(q, 1);
        for ph = 1:q, ii = ph:q:N; R.Xp{ph} = robustSlackFeatures(t(ii), y(ii), SL(ii), vt, true); end
        RUN(d,c).R = R;
        fprintf('%s %s done\n', DAYS{d}, COND{c});
    end
end
save(fullfile(res, 'data_noise_robust.mat'), 'RUN', 'DAYS', 'COND', 'NB', 'SEG14');

%% ---- summary ------------------------------------------------------------
ALL = fieldnames(RUN(1,1).R.base)';
ALL = ALL(~strcmp(ALL, 't0x'));                     % model-only
fprintf('\n%-20s %8s %8s %8s %8s %8s %8s %8s %8s\n', 'feature', 's_betw', 's_noise', 'noise/b', 's_phase', 'b_noise', 'r_ATP', 'sr', 'z_ATP');
T = table();
for f = ALL
    sg = 1:5; if any(strcmp(f{1}, SEG14)), sg = 1:4; end
    V8 = nan(3,5); V2 = nan(3,5);
    for d = 1:3, V8(d,:) = row(RUN(d,1).R.base, f{1}); V2(d,:) = row(RUN(d,2).R.base, f{1}); end
    mu8 = mean(V8, 1, 'omitnan'); sb = relSD(V8(:, sg));
    sn = nan(3,1); sp = nan(3,1); bn = nan(3,1);
    for d = 1:3
        X = cell2mat(cellfun(@(z) row(z, f{1}), RUN(d,1).R.Xr, 'uni', 0));
        sn(d) = rmsv(std(X(:, sg), 0, 1, 'omitnan') ./ abs(mu8(sg)));
        bn(d) = mean((mean(X(:, sg), 1, 'omitnan') - row2(RUN(d,1).R.clean, f{1}, sg)) ./ abs(mu8(sg)), 'omitnan');
        Xq = cell2mat(cellfun(@(z) row(z, f{1}), RUN(d,1).R.Xp, 'uni', 0));
        sp(d) = rmsv(std(Xq(:, sg), 0, 1, 'omitnan') ./ abs(mu8(sg)));
    end
    s_noise = rmsv(sn); s_phase = rmsv(sp); b_noise = mean(bn, 'omitnan');
    Rp = nan(3,1);
    for d = 1:3
        qq = V2(d,:) ./ V8(d,:); qq = qq(isfinite(qq) & qq > 0);
        if ~isempty(qq), Rp(d) = exp(mean(log(qq))); end
    end
    r = mean(Rp, 'omitnan'); sr = std(Rp, 'omitnan') / abs(r); z = abs(r - 1) / sqrt(sb^2 + sr^2);
    fprintf('%-20s %7.1f%% %7.1f%% %8.2f %7.1f%% %+7.1f%% %8.3f %8.3f %8.2f\n', f{1}, 100*sb, 100*s_noise, s_noise/sb, 100*s_phase, 100*b_noise, r, sr, z);
    T = [T; table(string(f{1}), sb, s_noise, s_noise/sb, s_phase, b_noise, r, sr, z, mu8, ...
        'VariableNames', {'feature','s_between','s_noise','noise_over_between','s_phase','b_noise','r_ATP','sr','z_ATP','mean8'})]; %#ok<AGROW>
end
save(fullfile(res, 'data_noise_robust.mat'), 'T', '-append');
writetable(T(:, 1:9), fullfile(res, 'data_noise_robust_summary.csv'));

function v = row(s, f)
    if isfield(s, f), v = double(s.(f)(:))'; else, v = nan(1,5); end
    if numel(v) < 5, v(end+1:5) = NaN; end; v = v(1:5);
end
function v = row2(s, f, sg), v = row(s, f); v = v(sg); end
function r = rmsv(x), r = sqrt(mean(x(isfinite(x)).^2)); end
function s = relSD(V)
    mu = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
    g = isfinite(mu) & abs(mu) > 0 & isfinite(sd);
    if ~any(g); s = NaN; return; end
    s = sqrt(mean((sd(g)./abs(mu(g))).^2));
end
