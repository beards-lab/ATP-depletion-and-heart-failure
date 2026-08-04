% RunR2CeilingSlope.m — does capping the strain-dependent detachment remove
% the model's spurious amplitude dependence?
%
% RunAmplitudeVsData established the signature to explain:
%   MODEL  d(rsK)/d(amplitude) = +714 s^-1 per ML   (r = +0.97)
%   DATA   d(rsK)/d(amplitude) = -126 s^-1 per ML   (r = -0.31, 3 preps)
% i.e. the model speeds up the harder you stretch it, and the real muscle
% does not. The leading suspect is the unbounded strain acceleration of the
% strong-bridge detachment R2, which reaches 2100-4250 s^-1 at re-stretch
% strains (RestretchVsKtrRecovery Part 3): a bigger stretch pushes bridges to
% higher strain, they detach faster, and the pool turns over faster.
%
% `UseR2Ceiling` (from ApoPoolDetachment) imposes
%     1/R2_eff = 1/R2_strain(s) + 1/R2max
% so R2 saturates instead of running away. That analysis FALSIFIED the ceiling
% on the ktr overstretch peak — but it never tested it against the amplitude
% SLOPE, which is a different and more specific observable. Doing that here.
%
% Writes results/r2ceiling_slope.mat + .png

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if isempty(gcp('nocreate')); parpool('local', 3); end

SNAP   = 'params/optfull7_opt_mov.m';
DEPTHS = [0.080 0.060 0.040 0.030 0.020 0.015 0.010];
K2MAX  = [Inf 1500 800 400 200];    % Inf = ceiling OFF (the current model)

p0 = getParams(loadParams(SNAP), [], true, false);
p0.MaxRunTime = 300;
fcn = str2func(p0.modelFcn);

V_REL = -36; V_STR = 3.9; T_ACT = 2.0; T_HOLD = 0.298; T_REC = 0.400;

R = struct('k2max',{},'depth',{},'k',{},'A',{},'Fss',{});
for ik = 1:numel(K2MAX)
    km = K2MAX(ik);
    fprintf('\n--- k2max = %g %s---\n', km, repmat('(ceiling OFF) ', 1, isinf(km)));
    for id = 1:numel(DEPTHS)
        d = DEPTHS(id);
        t_rel = T_ACT; dt_rel = d/abs(V_REL); t_h0 = t_rel+dt_rel;
        t_str = t_h0+T_HOLD; dt_str = d/V_STR; t_h1 = t_str+dt_str; t_end = t_h1+T_REC;
        vt = [-2.0,0,0,1.10; t_rel,V_REL,0,1.10; t_h0,0,0,1.10-d; ...
              t_str,V_STR,0,1.10-d; t_h1,0,0,1.10; t_end,0,0,1.10];

        p = p0;
        p.UseR2Ceiling = ~isinf(km);
        if ~isinf(km); p.k2max = km; end
        p.Velocity = vt(:,2);
        p = getParams(p, p.g, true);

        try
            [~, out] = evaluateModel(fcn, vt(:,1), p);
        catch ME
            fprintf('  depth %.3f FAILED: %s\n', d, ME.message); continue;
        end
        nr = makeMonotonous(out.t); T = out.t(nr); F = out.Force(nr);
        s = fitWin(T, F, t_h1, t_end, 0.080);

        n = numel(R)+1;
        R(n).k2max = km; R(n).depth = d; R(n).k = s.k; R(n).A = s.A;
        R(n).Fss = interp1(T, F, t_rel-0.005);
        fprintf('  depth %.3f  k %6.1f  A %6.2f  F_iso %6.2f\n', d, s.k, s.A, R(n).Fss);
    end
end

save(fullfile(here,'results','r2ceiling_slope.mat'), 'R', 'DEPTHS', 'K2MAX', 'SNAP');

%% ---- slopes -------------------------------------------------------------
fprintf('\n%-12s %10s %10s %10s %12s %10s\n', 'k2max', 'slope', 'r', 'k@0.08', 'k@0.01', 'F_iso');
rows = [];
for ik = 1:numel(K2MAX)
    m = [R.k2max] == K2MAX(ik) | (isinf(K2MAX(ik)) & isinf([R.k2max]));
    a = [R(m).depth]'; k = [R(m).k]';
    ok = isfinite(a) & isfinite(k);
    if nnz(ok) < 3; continue; end
    pf = polyfit(a(ok), k(ok), 1); rr = corr(a(ok), k(ok));
    kb = k(find(a==0.080,1)); ks = k(find(a==0.010,1));
    fs = mean([R(m).Fss], 'omitnan');
    if isempty(kb); kb = NaN; end; if isempty(ks); ks = NaN; end
    fprintf('%-12g %10.0f %10.2f %10.1f %12.1f %10.1f\n', K2MAX(ik), pf(1), rr, kb, ks, fs);
    rows(end+1,:) = [K2MAX(ik), pf(1), rr, kb, ks, fs]; %#ok<SAGROW>
end
fprintf('%-12s %10.0f %10.2f %10s %12s %10s\n', 'DATA (3 preps)', -126, -0.31, '~52', 'n/a', '~80');

fprintf(['\nREADING: the ceiling works on this signature only if the slope drops\n' ...
         'toward (or below) zero WITHOUT collapsing F_iso. A ceiling that flattens\n' ...
         'the slope by simply killing force is not a fix.\n']);

fig = figure(704); clf; set(fig,'Position',[70 70 950 400],'Color','w');
subplot(1,2,1); hold on; box on; co = lines(numel(K2MAX));
for ik = 1:numel(K2MAX)
    m = ([R.k2max] == K2MAX(ik)) | (isinf(K2MAX(ik)) & isinf([R.k2max]));
    lab = sprintf('k2max=%g', K2MAX(ik)); if isinf(K2MAX(ik)); lab = 'ceiling OFF'; end
    plot([R(m).depth], [R(m).k], 'o-', 'Color', co(ik,:), 'LineWidth',1.5, 'DisplayName',lab);
end
yline(43.7,'k--','data ~44 s^{-1}');
xlabel('re-stretch amplitude (ML)'); ylabel('rsK (s^{-1})');
legend('Location','northwest','FontSize',8); title('rate vs amplitude, by R2 ceiling');

subplot(1,2,2); hold on; box on;
if ~isempty(rows)
    x = 1:size(rows,1);
    yyaxis left;  bar(x, rows(:,2)); ylabel('slope (s^{-1} per ML)');
    yline(-126,'r--','data slope');
    yyaxis right; plot(x, rows(:,6), 'ks-','LineWidth',1.6); ylabel('F_{iso} (kPa)');
    set(gca,'XTick',x,'XTickLabel',arrayfun(@(v) sprintf('%g',v), rows(:,1),'UniformOutput',false));
    xlabel('k2max'); title('slope vs force cost');
end
exportgraphics(fig, fullfile(here,'results','r2ceiling_slope.png'), 'Resolution',140);
fprintf('wrote results/r2ceiling_slope.png\n');

% ---------------------------------------------------------------- helper
function s = fitWin(T, F, ta, tb, vallWin)
    m = T >= ta & T <= tb; t = T(m); y = F(m);
    s = struct('k',NaN,'A',NaN,'t0',NaN,'r2',NaN);
    if numel(t) < 20; return; end
    nW = max(3, nnz(t <= t(1)+vallWin)); [B, iB] = min(y(1:nW));
    tt = t(iB:end) - t(iB); yy = y(iB:end) - B;
    if numel(tt) < 20; return; end
    A = median(yy(tt >= tt(end) - 0.15*(tt(end)-tt(1))));
    if ~isfinite(A) || A <= 0; return; end
    obj = @(q) sum((yy - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
    q = fminsearch(obj, [45, 0.004], optimset('Display','off'));
    s.k = abs(q(1)); s.t0 = max(q(2),0); s.A = A;
    pred = A.*(1-exp(-s.k.*max(tt-s.t0,0)));
    s.r2 = 1 - sum((yy-pred).^2)/max(sum((yy-mean(yy)).^2),eps);
end
