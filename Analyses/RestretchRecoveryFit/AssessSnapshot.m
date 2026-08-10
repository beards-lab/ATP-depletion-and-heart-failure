% AssessSnapshot.m — evaluate one parameter snapshot on the full battery and
% report the post-restretch diagnosis: where rsK sits, whether the residual is
% LEVEL (rate too fast everywhere) or SLOPE (wrong amplitude dependence), and
% what the fit paid elsewhere.
%
% Set SNAPS (cell array of .m paths) before running, or it uses the default
% three-point comparison.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

if ~exist('SNAPS', 'var')
    SNAPS = {'params/optfull7_opt_mov.m', 'params/rsk_w025_opt.m', 'params/rskR2_w025_opt.m'};
end

ds  = load('data/protocol_03_27_2026_8mM_slack.mat');
vt  = ds.velocitytable;
iR  = find(vt(:,2) > 1);
amp = (vt(iR+1,4) - vt(iR,4))';          % re-stretch amplitude per cycle (ML)
rsK_d = ds.features_data.rsK;

R = struct();
for is = 1:numel(SNAPS)
    snap = SNAPS{is};
    fprintf('\n================ %s ================\n', snap);
    p = getParams(loadParams(snap), [], true, false);
    p.MaxRunTime = 600;
    p.RunForceVelocity = 1; p.RunSlack = 1; p.RunSlackPassive = 1;
    p.RunKtr = 0; p.RunStairs = 0; p.RunForceVelocityTime = 0;
    p.EvalFeatures = 1; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;

    r = evalRs(p, struct('slackOnly', false));
    if ~r.ok; fprintf(2, 'FAILED: %s\n', r.err); continue; end

    key = matlab.lang.makeValidName(erase(erase(snap,'params/'),'.m'));
    R.(key) = r;
    R.(key).snap = snap;
    R.(key).params = p;

    % L2 breakdown (what the optimiser actually minimises)
    [ct2, w2, c2] = evalFeatureCost(r.features_data, r.features_model, p.fn, 2);
    R.(key).ct2 = ct2; R.(key).fnList = p.fn;

    fprintf('feature total  L1 %.3f | L2 %.3f (optimiser''s norm)\n', r.featTotal, sum(ct2));
    fprintf('rsK model %s\n', mat2str(round(r.rsK_m,1)));
    fprintf('rsK data  %s   ratio %.2f\n', mat2str(round(rsK_d,1)), r.rsK_ratio);
    fprintf('ktr  x%.2f | steady %.1f kPa (data %.1f) | rsA x%.2f\n', ...
        r.ktr_ratio, mean(r.steady_m), mean(r.steady_d), ...
        mean(r.rsA_m,'omitnan')/mean(r.rsA_d,'omitnan'));

    % LEVEL vs SLOPE decomposition of the rsK residual
    pm = polyfit(amp(:), r.rsK_m(:), 1);
    pd = polyfit(amp(:), rsK_d(:),  1);
    fprintf('rsK vs re-stretch amplitude:  model slope %+7.0f  intercept %6.1f\n', pm(1), pm(2));
    fprintf('                              data  slope %+7.0f  intercept %6.1f\n', pd(1), pd(2));
    R.(key).slope_m = pm(1); R.(key).slope_d = pd(1);

    % top L2 contributors
    [~, o] = sort(ct2, 'descend');
    fprintf('\n  top L2 terms:\n');
    for i = o(1:min(8,numel(o)))
        nm = regexprep(strtok(p.fn{i},'|'), '\[.*$','');
        fprintf('    %-20s %8.3f\n', nm, ct2(i));
    end
end

save(fullfile(here,'results','assess_snapshot.mat'), 'R', 'SNAPS', 'amp', 'rsK_d');

%% ---- side-by-side --------------------------------------------------------
k = fieldnames(R);
fprintf('\n\n%-26s %8s %8s %8s %8s %8s %9s\n', 'snapshot', 'L2tot', 'rsK x', ...
        'ktr x', 'steady', 'rsA x', 'slope');
for i = 1:numel(k)
    r = R.(k{i});
    fprintf('%-26s %8.3f %8.2f %8.2f %8.1f %8.2f %9.0f\n', k{i}, sum(r.ct2), ...
        r.rsK_ratio, r.ktr_ratio, mean(r.steady_m), ...
        mean(r.rsA_m,'omitnan')/mean(r.rsA_d,'omitnan'), r.slope_m);
end
fprintf('%-26s %8s %8.2f %8.2f %8.1f %8.2f %9.0f\n', 'DATA', '-', 1.00, 1.00, ...
        mean(R.(k{1}).steady_d), 1.00, R.(k{1}).slope_d);

%% ---- parameter drift -----------------------------------------------------
if numel(k) >= 2
    pA = R.(k{1}).params; pB = R.(k{end}).params;
    f = fieldnames(pB); rows = {};
    for i = 1:numel(f)
        a = pA.(f{i}); b = pB.(f{i});
        if ~isnumeric(a) || ~isnumeric(b) || numel(a)~=numel(b) || isempty(a); continue; end
        if any(a==0); continue; end
        rr = b(:)./a(:);
        if any(abs(log(abs(rr))) > 0.02)
            rows(end+1,:) = {f{i}, a(1), b(1), max(rr)}; %#ok<SAGROW>
        end
    end
    fprintf('\nparameters that moved %s -> %s (|ln ratio| > 0.02):\n', k{1}, k{end});
    fprintf('%-34s %13s %13s %8s\n','param','first','last','ratio');
    for i = 1:size(rows,1)
        fprintf('%-34s %13.5g %13.5g %8.3f\n', rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4});
    end
end
