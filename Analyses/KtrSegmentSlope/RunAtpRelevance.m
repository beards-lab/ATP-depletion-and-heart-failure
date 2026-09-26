% RunAtpRelevance.m -- which scored metrics actually change between 8 and 2 mM,
% measured on ONE consistent curation of all six slack files and all six ktr
% files, with the robust rate estimators (ktr80, ktrP_k80).
%
% WHY. The deliverable is the low-ATP fit, and a neighbouring campaign found
% the 2 mM fit is wrong exactly where the 8 mM fit is wrong. The 2 mM data
% enter as RATIOS against 8 mM, so what the 8 mM fit must get right first is
% the baseline of every metric the ATP effect acts on. Earlier measurement
% (Analyses/FeatureWeighting, 2026-09-01): the ATP effect is one force scale
% plus three rates; the force SHAPE features are flat (z 0.07-0.21). This
% re-measures it on today's curation, with the new estimators, before the
% objective is rebuilt around it.
%
% WHAT IT DOES
%   1. re-extracts features_data for all six slack files through the exact
%      call runSlackExperiment makes (backs data/ up first -- not under git)
%   2. extracts the ktr-protocol features for all six ktr files
%   3. per metric: 8 mM between-prep spread (level and shape), and the paired
%      2mM/8mM log-ratio across the three preparations -> mean ratio and a
%      paired t = mean/(SD/sqrt(3)). |t| > 4.3 is p < 0.05 two-sided at 2 dof.
%   NOTE 04/10 ran 2 mM FIRST, so its force ratio carries rundown the other
%   way; rates are not rundown-corrected (it makes their CV worse).
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunAtpRelevance.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
DAYS = {'03_27','04_03','04_10'};
ATP  = {'8mM','2mM'};

bdir = fullfile(root,'data',sprintf('_backup_atprel_%s', datestr(now,'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
mkdir(bdir); fprintf('backup -> %s\n', bdir);

%% ---- 1. slack: re-curate all six ----------------------------------------
FS = cell(3,2);
for d = 1:3
    for a = 1:2
        fn = sprintf('protocol_%s_2026_%s_slack.mat', DAYS{d}, ATP{a});
        f = fullfile(root,'data',fn);
        copyfile(f, fullfile(bdir, fn));
        D = load(f);
        fd = extractSlackAttributes(D.datatable(:,1), D.datatable(:,3), D.datatable(:,2), ...
                                    D.velocitytable, struct(), [], false, true);
        D.features_data = fd; save(f, '-struct', 'D');
        FS{d,a} = fd;
        fprintf('%s  ktr %s | ktr80 %s\n', fn, mat2str(round(fd.ktr,1)), mat2str(round(fd.ktr80,1)));
    end
end

%% ---- 2. ktr protocol: all six -------------------------------------------
FK = cell(3,2);
for d = 1:3
    for a = 1:2
        if a == 1
            fn = sprintf('protocol_%s_2026_velocitytable_ktr.mat', DAYS{d});
        else
            fn = sprintf('protocol_%s_2026_velocitytable_ktr_2mM.mat', DAYS{d});
        end
        K = load(fullfile(root,'data',fn));
        FK{d,a} = extractKtrAttributes(K.datatable(:,1), K.datatable(:,3), K.velocitytable);
        fprintf('%-45s k63 %6.2f  k80 %6.2f  F20 %.3f\n', fn, FK{d,a}.ktrP_k63, ...
                FK{d,a}.ktrP_k80, FK{d,a}.ktrP_F20);
    end
end

%% ---- 3. the table -------------------------------------------------------
% {name, source 'S'/'K', isForce (shape = divide by steady)}
M = { 'steady','S',0; 'A','S',0; 'peak1_y','S',0; 'vall_y','S',0; 'peak2','S',0; ...
      'vall2_dy','S',0; 'restretchSlopeStart','S',0; 'A_over_steady','S',0; ...
      'peak1_y/steady','S',1; 'vall_y/steady','S',1; 'peak2/steady','S',1; ...
      'peak1_dSL','S',0; 'vall_t','S',0; 't0','S',0; ...
      'ktr','S',0; 'ktr80','S',0; 'ktr80_shape','S',0; ...
      'rsK80','S',0; 'rsT20','S',0; 'rsKfit','S',0; 'rsTfit','S',0; 'ktrfit','S',0; ...
      'ktrP_k63','K',0; 'ktrP_k80','K',0; 'ktrP_F20','K',0 };

getv = @(fd, nm) localGet(fd, nm);
fprintf('\n%-20s %8s %8s %8s | %8s %8s %8s %7s  %s\n', 'metric', 's8', 's8lev', 's8shp', ...
        'x2/8', 'sd_log', 't', 'n', 'verdict');
rows = {};
for i = 1:size(M,1)
    nm = M{i,1}; src = M{i,2};
    V8 = []; V2 = [];
    for d = 1:3
        if src == 'S'; a8 = getv(FS{d,1}, nm); a2 = getv(FS{d,2}, nm);
        else;          a8 = getv(FK{d,1}, nm); a2 = getv(FK{d,2}, nm); end
        V8(d,1:numel(a8)) = a8; V2(d,1:numel(a2)) = a2; %#ok<AGROW>
    end
    % 8 mM between-prep spread
    mu = mean(V8,1,'omitnan');
    s_raw = mean(std(V8,0,1,'omitnan') ./ abs(mu), 'omitnan');
    lv = mean(V8,2,'omitnan');
    s_lev = std(lv)/abs(mean(lv));
    Xs = V8 ./ lv; s_shp = mean(std(Xs,0,1,'omitnan')./abs(mean(Xs,1,'omitnan')),'omitnan');
    if size(V8,2) == 1; s_shp = NaN; end
    % paired ATP effect: per prep, mean over elements of log(v2/v8)
    ok = V8 > 0 & V2 > 0;
    L = nan(3,1);
    for d = 1:3
        if any(ok(d,:)); L(d) = mean(log(V2(d,ok(d,:)) ./ V8(d,ok(d,:)))); end
    end
    n = nnz(isfinite(L));
    mL = mean(L,'omitnan'); sL = std(L,'omitnan');
    tt = mL / (sL/sqrt(max(n,1)));
    verdict = 'flat';
    if abs(tt) > 4.3 && abs(exp(mL)-1) > 0.05; verdict = 'ATP-CARRYING';
    elseif abs(tt) > 2.9 && abs(exp(mL)-1) > 0.05; verdict = 'likely'; end
    fprintf('%-20s %8.3f %8.3f %8.3f | %8.3f %8.3f %8.2f %7d  %s\n', nm, s_raw, s_lev, s_shp, ...
            exp(mL), sL, tt, n, verdict);
    rows(end+1,:) = {nm, s_raw, s_lev, s_shp, exp(mL), sL, tt, verdict, V8, V2}; %#ok<SAGROW>
end

%% ---- 4. the ktr ladder itself: does its SLOPE change with ATP? ----------
fprintf('\n=== ktr80 ladder per prep (slope per segment, relative to the mean) ===\n');
for d = 1:3
    for a = 1:2
        k = FS{d,a}.ktr80(:)'; b = polyfit(1:5, k, 1);
        fprintf('  %s %s  %s   slope %+6.2f  (%+5.1f %%/seg)\n', DAYS{d}, ATP{a}, ...
                sprintf('%6.1f', k), b(1), 100*b(1)/mean(k));
    end
end

save(fullfile(resdir,'atp_relevance.mat'),'rows','FS','FK','M','DAYS');
fprintf('\nDONE\n');

function v = localGet(fd, nm)
    if contains(nm,'/')
        p = split(nm,'/'); v = fd.(p{1})(:)' ./ fd.(p{2})(:)';
    elseif isfield(fd, nm)
        v = fd.(nm)(:)';
    else
        v = NaN;
    end
end
