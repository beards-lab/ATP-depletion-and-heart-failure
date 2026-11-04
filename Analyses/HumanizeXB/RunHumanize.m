% RunHumanize.m -- score candidate "human" scalings of the incumbent
% (robust1_seed, mouse/rat alpha-MHC, ~20-25 C) with humanProbe, and write
% one table: results/<TAG>.txt (+ .mat).
%
%   SETS = {'name', struct('ka',0.2,...); ...}; TAG = 'explore'; run('.../RunHumanize.m')
%
% Every set is applied to the incumbent's *base* rates (PCHIP strain shapes are
% multipliers with central value ~1, so scaling the base rate scales the whole
% strain-dependent curve). Observables are read out with the same probe for the
% incumbent and every set.
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
if ~exist('TAG','var'), TAG = 'explore'; end
if ~exist('CORE','var'), CORE = 'robust1_seed'; end
if ~exist('OPTS','var'), OPTS = struct('Tiso', 3, 'Tktr', 2); end
if ~exist('SETS','var'), error('define SETS'); end
resH = fullfile(root, 'Analyses', 'HumanizeXB', 'results');
diary(fullfile(resH, [TAG '.log'])); diary on;
if isempty(gcp('nocreate'))
    cl = parcluster('local'); jsl = fullfile(tempdir, sprintf('mljob_hum_%d', feature('getpid')));
    if ~isfolder(jsl), mkdir(jsl); end
    cl.JobStorageLocation = jsl; parpool(cl, 5);
end
base = loadParams(CORE);
names = [{'incumbent'}; SETS(:,1)];
H = struct('name', {}, 'S', {}, 'R', {}, 'sec', {});
% resume: keep sets already recorded under this TAG (a dropped session must not cost the whole sweep)
if isfile(fullfile(resH, [TAG '.mat']))
    Hold = load(fullfile(resH, [TAG '.mat']), 'H'); Hold = Hold.H;
    Hold = Hold(arrayfun(@(h) isfield(h.R, 'ktr'), Hold));
    H = Hold(:)';
end
for i = 1:numel(names)
    t0 = tic;
    if any(strcmp({H.name}, names{i})), fprintf('%-14s kept from %s.mat\n', names{i}, TAG); continue; end
    if i == 1, S = struct(); else, S = SETS{i-1, 2}; end
    q = applyScaling(base, S);
    q = getParams(q, [], true, false);
    try
        R = humanProbe(q, OPTS);
    catch ME
        warning('%s failed: %s', names{i}, ME.message); R = struct();
    end
    H(end+1) = struct('name', names{i}, 'S', S, 'R', R, 'sec', toc(t0)); %#ok<SAGROW>
    fprintf('%-14s done in %.0f s\n', names{i}, H(end).sec);
    save(fullfile(resH, [TAG '.mat']), 'H', 'SETS', 'CORE', 'OPTS');
end
% ---- table (in SETS order)
[~, io] = ismember(names, {H.name}); H = H(io(io > 0));
cols = {'F0','kPa'; 'FXB0','kPa'; 'attached','-'; 'p2','-'; 'SRX','-'; 'PT','-'; 'PD','-'; ...
        'XTOR','1/s'; 'XTOR_net','1/s'; 'R2T','1/s'; 'R2D','1/s'; 'k2eff','1/s'; 'dwell_ms','ms'; 'ktr','1/s'; 'ktr_t50','ms'; ...
        'Vmax_hill','ML/s'; 'Vmax_lin','ML/s'; 'a_over_F0','-'; 'Pmax','F0*ML/s'; 'v_Pmax','ML/s'; 'iso_settled','-'};
fid = fopen(fullfile(resH, [TAG '.txt']), 'w');
fprintf(fid, '%-12s %-8s', 'quantity', 'unit'); fprintf(fid, ' %12s', names{:}); fprintf(fid, '\n');
for c = 1:size(cols,1)
    fprintf(fid, '%-12s %-8s', cols{c,1}, cols{c,2});
    for i = 1:numel(H)
        if isfield(H(i).R, cols{c,1}), fprintf(fid, ' %12.4g', H(i).R.(cols{c,1})); else, fprintf(fid, ' %12s', 'NaN'); end
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\nF/F0 at v (ML/s):\n');
vv = H(1).R.FV_v;
for k = 1:numel(vv)
    fprintf(fid, '%-12s %-8s', sprintf('v=%.2f', vv(k)), '-');
    for i = 1:numel(H)
        if isfield(H(i).R, 'FV_fnorm') && numel(H(i).R.FV_fnorm) >= k, fprintf(fid, ' %12.4f', H(i).R.FV_fnorm(k)); else, fprintf(fid, ' %12s', 'NaN'); end
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\nscalings:\n');
for i = 2:numel(H)
    f = fieldnames(H(i).S); s = '';
    for j = 1:numel(f), s = [s sprintf('%s x%.3g  ', f{j}, H(i).S.(f{j}))]; end %#ok<AGROW>
    fprintf(fid, '%-14s %s\n', H(i).name, s);
end
fclose(fid);
type(fullfile(resH, [TAG '.txt']));
diary off;
