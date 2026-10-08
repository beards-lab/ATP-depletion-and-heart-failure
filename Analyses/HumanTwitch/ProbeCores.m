%% ProbeCores.m -- humanProbe (isometric F0, SRX, ATP use, Brenner ktr, FV) of the fitted cores at saturating Ca
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
names = {'incumbent', 't1_isoT', 't2c_SE', 't2_all', 's3_free', 's3_free_b'};
if exist(fullfile(res, 't2_allc_fit.mat'), 'file'), names{end+1} = 't2_allc'; end
if exist(fullfile(res, 's3_free_c_fit.mat'), 'file'), names{end+1} = 's3_free_c'; end
cores = cell(size(names));
for k = 1:numel(names)
    if strcmp(names{k}, 'incumbent')
        c0 = tierConfigs('t0_v0');  cores{k} = c0.base;
    else
        F = load(fullfile(res, [names{k} '_fit.mat']));  cores{k} = buildModel(F.cfg, F.R.x);
    end
end
if isempty(gcp('nocreate')), parpool('local', min(numel(names), 6)); end
Rp = cell(size(names));
parfor k = 1:numel(names), Rp{k} = humanProbe(cores{k}, struct('Tiso', 4, 'Tktr', 2)); end
q = {'F0', 'FXB0', 'attached', 'SRX', 'PT', 'PD', 'R2T', 'XTOR', 'dwell_ms', 'ktr', 'ktr_t50', 'Vmax_hill', 'Vmax_lin', 'a_over_F0', 'Pmax'};
fprintf('%-10s', 'quantity');  fprintf('%12s', names{:});  fprintf('\n');
T = zeros(numel(q), numel(names));
for i = 1:numel(q)
    fprintf('%-10s', q{i});
    for k = 1:numel(names), v = NaN; if isfield(Rp{k}, q{i}), v = Rp{k}.(q{i}); end, T(i, k) = v; fprintf('%12.4g', v); end
    fprintf('\n');
end
save(fullfile(res, 'probe_cores.mat'), 'Rp', 'names', 'q', 'T');
