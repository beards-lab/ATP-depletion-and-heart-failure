%% ScoreStrategies.m -- the slack objective's own noise under each strategy.
% Feature noise matters through the COST, which amplifies it: a term is a
% squared residual minus a deadband, so near the deadband a 1 % feature wobble
% is a large fraction of what is left. For every configuration x point saved by
% RunModelNoise the slack terms of the objective are scored three ways:
%   S0 pipe   : as the objective does now
%   S1 reader : same definitions, same targets and floors, but every
%               trace-read feature read on the DATA's own 4 kHz grid anchored
%               on absolute time (robustSlackFeatures replicas); t0_crossing
%               from the interpolated LSE crossing
%   S2 replace: S1 plus the replacements (REPL below); a replacement's target
%               is the 3-preparation 8 mM mean and its floor nEl*s^2/3 with s
%               its between-preparation spread -- the construction the replaced
%               terms were given (RunDataNoiseRobust); its weight is the
%               replaced term's
% Reported per strategy x configuration: slack cost of the unperturbed point
% and its relative SD over the 1e-6 points.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');
L = load(fullfile(res, ['model_noise_' PARAMSET '.mat'])); M = L.M; vt = L.vt;
DR = load(fullfile(res, 'data_noise_robust.mat'), 'T');
DT = DR.T;
P = getParams(loadParams(PARAMSET), [], true, false);
% trace-read features the reader fix re-reads (current definitions)
REREAD = {'rsK80','rsKfit','rsT20','rsTfit','peak1_t','peak1_dSL','vall_t','ktr80','ktrfit', ...
          'restretchSlopeStart','steady','peak1_y','vall_y'};
% replacements: {scored name in fn, replacement feature}
if ~exist('REPL', 'var')
REPL = { 'rsK80',     'rsKls'
         'rsT20',     'rsT20e'
         'peak1_dSL', 'peak1_dSLq'
         'vall_t',    'tFall50'
         'restretchSlopeStart', 'rsSlope1ms' };
end
base = @(e) regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '');
cf = fieldnames(M)'; OUT = struct();
for c = cf
    O = M.(c{1}); n = numel(O.TR);
    fn = P.fn; FD = O.FD;
    keep = cellfun(@(e) ~contains(e, ',') && isfield(O.FM{1}, base(e)), fn); fnS = fn(keep);
    % S2 objective: swap the replaced entries, new targets from the data
    fnR = fnS; FDR = FD;
    for k = 1:size(REPL, 1)
        j = find(strcmp(cellfun(base, fnS, 'uni', 0), REPL{k,1}), 1);
        if isempty(j), continue; end
        tok = split(fnS{j}, '|'); w = str2double(tok{end});
        sel = regexp(tok{1}, '\[.*\]', 'match'); sel = [sel{:}];
        row = DT(DT.feature == REPL{k,2}, :);
        nEl = 5; if ~isempty(sel), nEl = numel(eval(sel)); end
        fl = nEl * row.s_between^2 / 3;
        fnR{j} = sprintf('%s%s|fl:%g|%g', REPL{k,2}, sel, fl, w);
        FDR.(REPL{k,2}) = row.mean8;
    end
    C = nan(n, 3); CT = cell(n, 3);
    for r = 1:n
        fm = O.FM{r};
        CT{r,1} = evalFeatureCost(FD, fm, fnS); C(r,1) = sum(CT{r,1});
        tr = O.TR{r}; [t, iu] = unique(tr.t);
        g = (ceil(t(1)/2.5e-4) : floor(t(end)/2.5e-4))' * 2.5e-4;
        rb = robustSlackFeatures(g, interp1(t, tr.F(iu), g), interp1(t, tr.SL(iu), g), vt, false, t, tr.LSE(iu));
        f1 = fm;
        for q = REREAD, f1.(q{1}) = rb.(q{1}); end
        f1.t0_crossing = rb.t0x;
        f1.ktrfit_shape = rb.ktrfit / mean(rb.ktrfit, 'omitnan');
        f1.ktr80_shape  = rb.ktr80 / mean(rb.ktr80, 'omitnan');
        CT{r,2} = evalFeatureCost(FD, f1, fnS); C(r,2) = sum(CT{r,2});
        f2 = f1; for k = 1:size(REPL,1), f2.(REPL{k,2}) = rb.(REPL{k,2}); end
        CT{r,3} = evalFeatureCost(FDR, f2, fnR); C(r,3) = sum(CT{r,3});
    end
    OUT.(c{1}) = struct('name', O.name, 'C', C, 'CT', {CT}, 'fnS', {fnS}, 'fnR', {fnR});
    fprintf('%-44s', O.name);
    for s = 1:3, fprintf('  S%d %.4f +- %5.2f %%', s-1, C(1,s), 100*std(C(:,s))/mean(C(:,s))); end
    fprintf('\n');
end
save(fullfile(res, ['strategies_' PARAMSET '.mat']), 'OUT', 'REPL');

%% per-term breakdown: unperturbed cost and SD over the points, per strategy
lab = @(e) regexprep(e, '\|.*$', '');
for c = cf
    O = OUT.(c{1}); fprintf('\n--- %s ---\n%-26s %18s %18s %-26s %18s\n', O.name, 'term (S0/S1)', 'S0 cost (SD)', 'S1 cost (SD)', 'term (S2)', 'S2 cost (SD)');
    K = cell(1,3); for s = 1:3, K{s} = vertcat(O.CT{:,s}); end
    for j = 1:numel(O.fnS)
        if max(abs([K{1}(:,j); K{2}(:,j); K{3}(:,j)])) < 1e-6, continue; end
        fprintf('%-26s %9.4f (%6.4f) %9.4f (%6.4f) %-26s %9.4f (%6.4f)\n', lab(O.fnS{j}), K{1}(1,j), std(K{1}(:,j)), ...
            K{2}(1,j), std(K{2}(:,j)), lab(O.fnR{j}), K{3}(1,j), std(K{3}(:,j)));
    end
end
save(fullfile(res, ['strategies_' PARAMSET '.mat']), 'OUT', 'REPL');
