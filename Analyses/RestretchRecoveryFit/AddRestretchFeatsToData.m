% AddRestretchFeatsToData.m
% Add the new post-restretch recovery features (rsK, rsA, rsT0, rsR2, rsK63)
% to the stored features_data of every slack protocol .mat file.
%
% SAFETY: this re-runs extractSlackAttributes and then compares EVERY
% pre-existing field against what is stored. Existing fields are NEVER
% overwritten, and any field that would have changed is reported loudly — the
% stored features are live fit targets and must not move silently.
%
% Set DRYRUN=false to actually write.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('DRYRUN', 'var'); DRYRUN = true; end
NEWFIELDS = {'rsK', 'rsA', 'rsT0', 'rsR2', 'rsK63'};
RELTOL    = 1e-6;

DATA_FILES = { ...
    'data/protocol_03_27_2026_8mM_slack.mat', ...
    'data/protocol_03_27_2026_2mM_slack.mat', ...
    'data/protocol_04_10_2026_8mM_slack.mat', ...
    'data/protocol_04_10_2026_2mM_slack.mat', ...
    'data/protocol_04_10_2026_ActivePNBMava_slack.mat', ...
    'data/protocol_04_03_2026_8mM_slack.mat', ...
    'data/protocol_04_03_2026_2mM_slack.mat'};

fprintf('DRYRUN = %d\n', DRYRUN);

for i_file = 1:numel(DATA_FILES)
    F = DATA_FILES{i_file};
    if ~isfile(F)
        fprintf('\n--- %-52s  SKIP (no such file)\n', F);
        continue;
    end
    ds = load(F);
    if ~isfield(ds, 'features_data') || isempty(fieldnames(ds.features_data))
        fprintf('\n--- %-52s  SKIP (no features_data)\n', F);
        continue;
    end

    fprintf('\n=== %s ===\n', F);
    fd_old = ds.features_data;

    fresh = extractSlackAttributes(ds.datatable(:,1), ds.datatable(:,3), ...
                                   ds.datatable(:,2), ds.velocitytable, ...
                                   [], [], false, true);

    % --- regression check on every pre-existing numeric field ---
    moved = {};
    for f = fieldnames(fd_old)'
        nm = f{1};
        if ~isfield(fresh, nm); continue; end
        a = fd_old.(nm); b = fresh.(nm);
        if ~isnumeric(a) || ~isnumeric(b) || numel(a) ~= numel(b); continue; end
        d = abs(a(:) - b(:));
        s = max(abs(a(:)), 1e-12);
        fin = isfinite(d);
        if any(d(fin) > RELTOL * s(fin))
            moved{end+1} = sprintf('%s(%.1f%%)', nm, 100*max(d(fin)./s(fin))); %#ok<SAGROW>
        end
    end
    if isempty(moved)
        fprintf('  regression: OK — all %d shared fields reproduce to %g rel.\n', ...
                numel(fieldnames(fd_old)), RELTOL);
    else
        fprintf(2, ['  regression: %d field(s) DIFFER, max rel. change in ()\n' ...
                    '              -> stored values KEPT; these were extracted by an\n' ...
                    '                 older extractSlackAttributes. %s\n'], ...
                numel(moved), strjoin(moved, ', '));
    end

    % --- add only the new fields ---
    fd_new = fd_old;
    for f = NEWFIELDS
        nm = f{1};
        if ~isfield(fresh, nm)
            fprintf(2, '  MISSING new field %s in fresh extraction!\n', nm); continue;
        end
        fd_new.(nm) = fresh.(nm);
        fprintf('  + %-6s = %s\n', nm, mat2str(round(fresh.(nm), 4, 'significant')));
    end

    if ~DRYRUN
        ds.features_data = fd_new;
        save(F, '-struct', 'ds');
        fprintf('  WROTE %s\n', F);
    end
end

fprintf('\nDone (DRYRUN=%d).\n', DRYRUN);
