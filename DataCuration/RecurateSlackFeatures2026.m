% RecurateSlackFeatures2026.m -- re-extract features_data for the three 8 mM
% slack files with ONE consistent call, after the peak-prominence fix.
%
% WHY THIS AND NOT UpdateSlackFeatures.m
%   That script's DATA_FILES list does not contain 04/03 at all, and its
%   featRunIdx defaults to 3:5 (the 04/10 rows), so the three preparations'
%   stored features were not all produced by the same code path. Since the
%   pooled target averages them, that is a silent inconsistency.
%
% WHAT CHANGED UPSTREAM
%   extractSlackAttributes' peak prominence was a fixed 0.1 kPa = 0.07-0.11
%   sigma of these recordings' own noise, so findpeaks returned 11-15 ripples
%   per re-stretch and peak1_t(1) took the earliest. It is now
%   max(3*sigma, 5% of the window's force range). See Analyses/RestretchPeakAudit.
%
% SAFETY: data/ is not under version control, so every file is copied to
% data/_backup_recurate_<timestamp>/ before it is touched.
%
%   matlab -batch "run('DataCuration/RecurateSlackFeatures2026.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));

FILES = { 'protocol_03_27_2026_8mM_slack.mat'
          'protocol_04_03_2026_8mM_slack.mat'
          'protocol_04_10_2026_8mM_slack.mat' };
SHOW  = {'peak1_dSL','vall_t','peak1_y','vall_y','A','steady','peak2','ktr','rsK80'};

bdir = fullfile(root,'data',sprintf('_backup_recurate_%s', datestr(now,'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
mkdir(bdir);
fprintf('backup -> %s\n', bdir);

OLD = struct(); NEW = struct();
for i = 1:numel(FILES)
    f = fullfile(root,'data',FILES{i});
    copyfile(f, fullfile(bdir, FILES{i}));
    D = load(f);
    OLD(i).fd = D.features_data;

    % EXACTLY the call runSlackExperiment makes for the data path
    % (runSlackExperiment.m:390), minus the plotting.
    fd = extractSlackAttributes(D.datatable(:,1), D.datatable(:,3), D.datatable(:,2), ...
                                D.velocitytable, struct(), [], false, true);
    NEW(i).fd = fd;

    D.features_data = fd;
    save(f, '-struct', 'D');
    fprintf('\n=== %s ===\n', FILES{i});
    for q = SHOW
        if ~isfield(OLD(i).fd,q{1}) && ~isfield(fd,q{1}); continue; end
        o = []; n = [];
        if isfield(OLD(i).fd,q{1}); o = OLD(i).fd.(q{1})(:)'; end
        if isfield(fd,q{1});        n = fd.(q{1})(:)';        end
        fprintf('  %-11s old %s\n', q{1}, mat2str(round(o,4)));
        fprintf('  %-11s NEW %s\n', '',   mat2str(round(n,4)));
    end
end

%% ---- what the between-preparation spread did -------------------------
fprintf('\n=============== BETWEEN-PREPARATION SPREAD ===============\n');
fprintf('%-12s %10s %10s %10s\n','feature','old','new','change');
for q = SHOW
    vo = nan(3,5); vn = nan(3,5);
    ok = true;
    for i = 1:3
        if ~isfield(OLD(i).fd,q{1}) || ~isfield(NEW(i).fd,q{1}); ok = false; break; end
        a = OLD(i).fd.(q{1})(:)'; b = NEW(i).fd.(q{1})(:)';
        m = min([5 numel(a) numel(b)]);
        vo(i,1:m) = a(1:m); vn(i,1:m) = b(1:m);
    end
    if ~ok; continue; end
    so = mean(std(vo,0,1)./abs(mean(vo,1)),'omitnan');
    sn = mean(std(vn,0,1)./abs(mean(vn,1)),'omitnan');
    fprintf('%-12s %10.3f %10.3f %9.2fx\n', q{1}, so, sn, sn/so);
end
fprintf('\nNext: rebuild the pooled target (Analyses/DataPooling/WritePooledTarget.m),\n');
fprintf('then re-derive the objective floors from the new spreads.\n');
fprintf('\nDONE\n');
