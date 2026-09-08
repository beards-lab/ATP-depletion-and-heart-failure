% PredictionMetrics.m -- put numbers on the two UNSCORED protocols (ktr, staircase)
% from a saved RunFitReport run, so "the prediction held" is a measurement and not
% an impression formed from a figure. Uses RunFitReport's own normalisation: the
% data column is ALREADY F/F_iso, and the model is divided by its force at the time
% of the first velocity command.
%
% Requires Analyses/OptimResultReport/results/fit_report_<TAG>.mat (RunFitReport).
clearvars -except TAG;                % keep a caller-set TAG, drop everything else
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); addpath(genpath('.'));
if ~exist('TAG','var'); TAG = 'pooled8_opt'; end
% RunFitReport names its .mat after the snapshot (<tag>_opt), but callers tend to
% hold the campaign tag (<tag>). Accept either rather than failing on a stale one.
MATF = sprintf('Analyses/OptimResultReport/results/fit_report_%s.mat', TAG);
if ~isfile(MATF)
    alt = sprintf('Analyses/OptimResultReport/results/fit_report_%s_opt.mat', TAG);
    if isfile(alt); MATF = alt; else
        error('PredictionMetrics:noRun', ['neither %s nor %s exists — run ' ...
            'Analyses/OptimResultReport/RunFitReport.m with TAG_NEW first.'], MATF, alt);
    end
end
fprintf('reading %s\n', MATF);
S = load(MATF);
R = S.R; NAMES = S.NAMES;

for which = {'stairs','ktr'}
    w = which{1};
    if strcmp(w,'ktr')
        Sd = load(['data/' R(2).p.ktr_velocitytableonfile]); fld = 'out_ktr';
    else
        Sd = load(['data/' R(2).p.stairs_velocitytableonfile]); fld = 'out_stairs';
    end
    dt2 = Sd.datatable; vt2 = Sd.velocitytable;
    td = dt2(:,1); yd = dt2(:,3);                     % already normalised
    if strcmp(w,'ktr'); tref = vt2(2,1); else; tref = vt2(find(vt2(:,2)~=0,1),1); end

    T = cell(1,3); Y = cell(1,3); NM = {'data', NAMES{1}, NAMES{2}};
    T{1} = td; Y{1} = yd;
    for k = 1:2
        o = R(k).(fld);
        F0 = interp1(o.t, o.Force, tref, 'linear','extrap'); if F0 == 0; F0 = 1; end
        T{k+1} = o.t; Y{k+1} = o.Force/F0;
    end

    fprintf('\n===== %s protocol (NEVER scored) =====\n', upper(w));

    if strcmp(w,'stairs')
        % steps = commanded stretch onsets from the velocity table
        iv = find(vt2(:,2) > 1e-6);
        st = vt2(iv,1); st = st([true; diff(st) > 0.02]);
        fprintf('%-14s %8s %8s %8s   %8s %8s\n','step peak/iso','data',NM{2},NM{3},'d base','d new');
        P = nan(numel(st),3);
        for s = 1:numel(st)
            for k = 1:3
                m = T{k} >= st(s) & T{k} <= st(s)+0.05;
                if any(m); P(s,k) = max(Y{k}(m)); end
            end
            fprintf('  step %-8d %8.3f %8.3f %8.3f   %+7.1f%% %+7.1f%%\n', s, P(s,1), P(s,2), P(s,3), ...
                100*(P(s,2)/P(s,1)-1), 100*(P(s,3)/P(s,1)-1));
        end
        g = all(isfinite(P),2);
        fprintf('  %-12s %8.3f %8.3f %8.3f   %+7.1f%% %+7.1f%%\n','MEAN', ...
            mean(P(g,1)), mean(P(g,2)), mean(P(g,3)), ...
            100*(mean(P(g,2)./P(g,1))-1), 100*(mean(P(g,3)./P(g,1))-1));

        % final plateau after the last step = the force-length end point
        t1 = max(td); pl = nan(1,3);
        for k = 1:3
            m = T{k} >= t1-0.12 & T{k} <= t1;
            pl(k) = mean(Y{k}(m),'omitnan');
        end
        fprintf('\n  %-12s %8.3f %8.3f %8.3f   %+7.1f%% %+7.1f%%\n','plateau', pl(1),pl(2),pl(3), ...
            100*(pl(2)/pl(1)-1), 100*(pl(3)/pl(1)-1));
        fprintf('  (F at SL 2.20 over F at SL 2.00 -- the force-length end point)\n');
    else
        for k = 1:3
            t = T{k}; y = Y{k};
            [ovs, iM] = max(y);
            A  = mean(y(t > max(td)-0.006 & t <= max(td)),'omitnan');
            % redevelopment: 20-80 % of A after the post-stretch release to zero
            iz = find(t > t(iM) & y < 0.1*A);
            if isempty(iz); iz = iM; else; iz = iz(end); end
            tw = t(t >= t(iz)); yw = y(t >= t(iz));
            i20 = find(yw >= 0.2*A, 1); i80 = find(yw >= 0.8*A, 1);
            kr = NaN; if ~isempty(i20) && ~isempty(i80) && i80 > i20; kr = log(4)/(tw(i80)-tw(i20)); end
            fprintf('%-16s overshoot %6.3f   ktr(20-80) %6.1f s^-1   plateau %6.3f\n', NM{k}, ovs, kr, A);
        end
    end
end
