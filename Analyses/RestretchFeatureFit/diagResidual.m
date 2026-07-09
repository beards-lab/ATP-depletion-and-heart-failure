% diagResidual.m
% =========================================================================
% Diagnostic: decompose the feature cost of the current best 2-state snapshot
% (params/opt2state_opt.m, reported total ~4.48) into per-feature
% contributions, and dump model-vs-data values for the named features so we
% can see WHERE the residual lives (ktr, restretch peaks, physiology bounds...).
%
% Writes a plaintext report to Analyses/RestretchFeatureFit/diagResidual_report.txt
% and saves the raw feature structs to diagResidual_feats.mat
%
% Run:  cd(root); addpath(genpath('.')); diagResidual
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_opt.m';
REPORT = fullfile('Analyses','RestretchFeatureFit','diagResidual_report.txt');

% --- read the fn (cost list) straight off the snapshot for exact fidelity ---
params0 = getParams(); run(SNAP); FN = params0.fn;

% --- evaluate cost + collect features ---
[tc, cv, fm, fd] = costOfSnap(SNAP, FN, 120);

fid = fopen(REPORT, 'w');
fprintf(fid, '==== diagResidual: %s ====\n', SNAP);
fprintf(fid, 'TOTAL feature cost = %.4f\n\n', tc);

% --- per-feature cost, sorted descending ---
[~, order] = sort(cv, 'descend', 'MissingPlacement', 'last');
fprintf(fid, '---- per-feature cost (sorted) ----\n');
fprintf(fid, '%-72s %10s\n', 'feature spec', 'cost');
for i = 1:numel(order)
    j = order(i);
    fprintf(fid, '%-72s %10.4f\n', FN{j}, cv(j));
end
fprintf(fid, '\nsum(cv) = %.4f\n\n', sum(cv, 'omitnan'));

% --- model vs data for named features (arrays over 5 slack segments) ---
fprintf(fid, '---- model vs data, key features (per slack segment) ----\n');
names = {'ktr','ktr_rmse','A','peak1_y','peak1_dSL','peak2','vall_y','vall_t', ...
         'vall2_dy','ovrsht_dy','steady','restretchSlopeStart','t0_crossing', ...
         'XTOR','XTOR_vmax','SRX_ss','attached_ss','PT_ss'};
for i = 1:numel(names)
    dumpfeat_local(fid, names{i}, fm, fd);
end

% --- FV block ---
fprintf(fid, '\n---- force-velocity ----\n');
if isfield(fm, 'FV_fnorm'); fprintf(fid, '  FV_fnorm model = [%s]\n', num2str(fm.FV_fnorm(:)', '%7.3f')); end
if isfield(fd, 'FV_fnorm'); fprintf(fid, '  FV_fnorm data  = [%s]\n', num2str(fd.FV_fnorm(:)', '%7.3f')); end
if isfield(fm, 'FV_v');     fprintf(fid, '  FV_v           = [%s]\n', num2str(fm.FV_v(:)',     '%7.3f')); end

fclose(fid);
save(fullfile('Analyses','RestretchFeatureFit','diagResidual_feats.mat'), 'fm','fd','cv','FN','tc');

type(REPORT);
disp('DONE diagResidual');

% ---- local function (must be at end of script) ----
function dumpfeat_local(fid, name, fm, fd)
    m = NaN; d = NaN;
    if isfield(fm, name) && ~isempty(fm.(name)); m = fm.(name); end
    if isfield(fd, name) && ~isempty(fd.(name)); d = fd.(name); end
    fprintf(fid, '  %-18s model=[%s]\n', name, num2str(m(:)', '%9.4g'));
    fprintf(fid, '  %-18s data =[%s]\n', '',   num2str(d(:)', '%9.4g'));
end
