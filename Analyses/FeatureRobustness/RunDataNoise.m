%% RunDataNoise.m -- how much of each slack feature is the transducer's noise?
%
% For each 2026 preparation x ATP condition: extract the features (current
% definitions via extractSlackAttributes, exactly as poolLoadRuns does, and the
% robust replacements via robustSlackFeatures), then re-extract them from NB
% surrogate recordings = the trace's smooth part + fresh noise with the
% recording's own level and short-range correlation (block bootstrap of the
% residuals in the flat end of each post-restretch hold). The SD over the
% surrogates is the feature's MEASUREMENT-noise SD -- what it would scatter by
% if the same muscle were recorded again. Compared against
%   s_between : the between-preparation relative SD (buildPooledTarget's relSD)
%   z_ATP     : the 2 mM / 8 mM discriminability (buildPooledTarget's r, sr, z)
% A feature whose noise SD is a large part of s_between is mostly measuring the
% transducer; one whose z_ATP is small does not carry the ATP effect either way.
%
% Also a SAMPLING-PHASE test: the 4 kHz record decimated to 1 kHz (the rate the
% model trace is read at) at its 4 possible phases.
if ~exist('NB', 'var'), NB = 60; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root)); res = fullfile(here, 'results');
DAYS = {'03_27', '04_03', '04_10'}; COND = {'8mM', '2mM'};
CUR = {'rsK80','rsKfit','rsT20','rsTfit','peak1_t','peak1_dSL','vall_t','ktr80','ktrfit','t0','restretchSlopeStart','steady','peak1_y','vall_y'};
ROB = {'rsT20e','rsTfitE','peak1_tq','peak1_dSLq','vall_tq','tFall50','tOn20','rsSlope1ms'};
SEG14 = {'rsK80','rsKfit','rsT20','rsTfit','rsT20e','rsTfitE'};   % scored on cycles 1-4
if ~exist('NW', 'var'), NW = 5; end
if isempty(gcp('nocreate')), parpool(NW); end

rng(5);
RUN = struct();
for d = 1:numel(DAYS)
    for c = 1:numel(COND)
        S = load(fullfile(root, 'data', sprintf('protocol_%s_2026_%s_slack.mat', DAYS{d}, COND{c})));
        t = S.datatable(:,1); SL = S.datatable(:,2); y = S.datatable(:,3); vt = S.velocitytable;
        dt = median(diff(t));

        % ---- noise model: residuals in the flat end of each post-restretch hold
        segs = find(vt(:,2) < -1)'; rq = {};
        for s = segs
            if s + 4 > size(vt,1), break; end
            m = t >= vt(s+4,1) - 0.15 & t <= vt(s+4,1) - 0.005;
            tq = t(m); yq = y(m); rq{end+1} = yq - polyval(polyfit(tq - tq(1), yq, 2), tq - tq(1)); %#ok<SAGROW>
        end
        rpool = vertcat(rq{:});
        sig = std(rpool);
        ac = arrayfun(@(L) corr(rpool(1:end-L), rpool(1+L:end)), 1:4);
        ysm = sgolayfilt(y, 3, 9);
        g = var(sgolayfilt(rpool, 3, 9)) / var(rpool);   % noise the smoother lets through
        cfac = sqrt(max(0, 1 - g));
        BL = round(0.010 / dt);                            % 10 ms blocks keep the correlation

        base = extractSlackAttributes(t, y, SL, vt, [], [], false, true);
        rb   = robustSlackFeatures(t, y, SL, vt, true);

        Xc = cell(NB,1); Xr = cell(NB,1);
        N = numel(y);
        parfor b = 1:NB
            st = RandStream('mt19937ar', 'Seed', 1000*d + 100*c + b);
            nz = zeros(N,1); k = 0;
            while k < N
                j0 = randi(st, numel(rpool) - BL); L = min(BL, N - k);
                nz(k+1:k+L) = rpool(j0:j0+L-1); k = k + L;
            end
            ys = ysm + cfac * nz;
            fc = extractSlackAttributes(t, ys, SL, vt, [], [], false, true);
            Xc{b} = fc; Xr{b} = robustSlackFeatures(t, ys, SL, vt, true);
        end

        % ---- sampling phase: decimate to 1 kHz at each phase ----------------
        q = round(1e-3 / dt); Xp = cell(q,1);
        for ph = 1:q
            ii = ph:q:N;
            Xp{ph} = robustSlackFeatures(t(ii), y(ii), SL(ii), vt, true);
        end

        R = struct('day', DAYS{d}, 'cond', COND{c}, 'sig', sig, 'ac', ac, 'g', g, ...
                   'base', base, 'rb', rb);
        R.Xc = Xc; R.Xr = Xr; R.Xp = Xp;
        RUN(d,c).R = R;
        fprintf('%s %s: noise sigma %.2f kPa, autocorr lag1-4 %s, smoother lets through %.2f\n', ...
                DAYS{d}, COND{c}, sig, mat2str(round(ac,2)), g);
    end
end
save(fullfile(res, 'data_noise.mat'), 'RUN', 'CUR', 'ROB', 'SEG14', 'DAYS', 'COND', 'NB');

%% ---- summary ------------------------------------------------------------
% value of feature f for run (d,c) from source: 'base' | surrogate b | phase
getv = @(R, f, src, b) pick(R, f, src, b);
ALL = [CUR ROB];
fprintf('\n%-20s %8s %8s %8s %8s %8s %8s %8s\n', 'feature', 's_betw', 's_noise', 'noise/b', 's_phase', 'r_ATP', 'sr', 'z_ATP');
T = table();
for f = ALL
    sg = 1:5; if any(strcmp(f{1}, SEG14)), sg = 1:4; end
    V8 = nan(3,5); V2 = nan(3,5);
    for d = 1:3
        V8(d,:) = getv(RUN(d,1).R, f{1}, 'base', 0);
        V2(d,:) = getv(RUN(d,2).R, f{1}, 'base', 0);
    end
    mu8 = mean(V8, 1, 'omitnan');
    sb = relSD(V8(:, sg));
    % measurement noise: per run, per cycle SD over surrogates / |pooled 8 mM mean|
    sn = nan(3,1); sp = nan(3,1);
    for d = 1:3
        X = cell2mat(arrayfun(@(b) getv(RUN(d,1).R, f{1}, 'sur', b), (1:NB)', 'uni', 0));
        sn(d) = sqrt(mean((std(X(:, sg), 0, 1, 'omitnan') ./ abs(mu8(sg))).^2, 'omitnan'));
        if isfield(RUN(d,1).R.Xp{1}, f{1})
            Xq = cell2mat(cellfun(@(z) z.(f{1})(:)', RUN(d,1).R.Xp, 'uni', 0));
            sp(d) = sqrt(mean((std(Xq(:, sg), 0, 1, 'omitnan') ./ abs(mu8(sg))).^2, 'omitnan'));
        end
    end
    s_noise = sqrt(mean(sn.^2, 'omitnan')); s_phase = sqrt(mean(sp.^2, 'omitnan'));
    Rp = nan(3,1);
    for d = 1:3
        qq = V2(d,:) ./ V8(d,:); qq = qq(isfinite(qq) & qq > 0);
        if ~isempty(qq), Rp(d) = exp(mean(log(qq))); end
    end
    r = mean(Rp, 'omitnan'); sr = std(Rp, 'omitnan') / abs(r); z = abs(r - 1) / sqrt(sb^2 + sr^2);
    fprintf('%-20s %7.1f%% %7.1f%% %8.2f %7.1f%% %8.3f %8.3f %8.2f\n', f{1}, 100*sb, 100*s_noise, s_noise/sb, 100*s_phase, r, sr, z);
    T = [T; table(string(f{1}), sb, s_noise, s_noise/sb, s_phase, r, sr, z, ...
        'VariableNames', {'feature','s_between','s_noise','noise_over_between','s_phase','r_ATP','sr','z_ATP'})]; %#ok<AGROW>
end
save(fullfile(res, 'data_noise.mat'), 'T', '-append');
writetable(T, fullfile(res, 'data_noise_summary.csv'));

function v = pick(R, f, src, b)
    switch src
        case 'base', if isfield(R.base, f), v = R.base.(f)(:)'; else, v = R.rb.(f)(:)'; end
        case 'sur',  if isfield(R.Xc{b}, f), v = R.Xc{b}.(f)(:)'; else, v = R.Xr{b}.(f)(:)'; end
    end
    v = double(v); if numel(v) < 5, v(end+1:5) = NaN; end; v = v(1:5);
end

function s = relSD(V)
    mu = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
    g = isfinite(mu) & abs(mu) > 0 & isfinite(sd);
    if ~any(g); s = NaN; return; end
    s = sqrt(mean((sd(g)./abs(mu(g))).^2));
end
