% RunFVConsistency.m -- does the in-preparation force-velocity point make the
% force-velocity data consistent with everything else?
%
% THE PROBLEM
% The F-V curve in the objective comes from two Baker-lab isovelocity datasets
% recorded on OTHER preparations, in another era, and it is scored as
% F(v)/F(0) -- normalised by the MODEL's own isometric force. So F-V carries no
% information about how strong the fibre is, and nothing ties it to the same
% preparations the slack protocol was run on. At FV_SHARE = 0.70 that makes most
% of the objective scale-blind (FeatureWeighting, the anchor warning in S4).
%
% WHAT IS ACTUALLY IN THE DATA
% The 04-series slack recordings open with a load-clamp manoeuvre whose second
% phase is a 34 ms ISOVELOCITY RAMP at -2 ML/s, ending just before the first
% slack release. It is a genuine sub-Vmax force-velocity point measured in the
% SAME activation, on the SAME preparation, against the SAME isometric
% reference (Fss22) as every slack feature. Call it FV2.
%
% NOTE: DataCuration/CreateFVRampProtocolVelocityTable.m builds FV2 files for
% 04/10 only. The ramp is present and identical (-2.00 ML/s, 34.6 ms, same
% commanded length) in 04/03 as well, in all three of its conditions. 03/27
% predates the manoeuvre and has none. That is 2 preparations x (8 mM, 2 mM,
% PNB+Mava) = six FV2 measurements, not three.
%
% THREE QUESTIONS
%   1. Does FV2 agree with the Baker curves at 2 ML/s, in normalised units?
%      If yes, the Baker F-V is transferable and may be pooled; if no, it is a
%      separate population and the 70 % pin is resting on borrowed data.
%   2. How much of the ramp force is PASSIVE? The PNB+Mava run runs the same
%      ramp, so this is measured, not assumed.
%   3. Does low ATP move the F-V point, and by how much relative to its spread?
%
% Output -> results/fv_consistency.mat, results/fv_consistency.png

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

%% ---- every slack burst, scanned for the ramp --------------------------
SRC = { '03/27','8mM','03 27 2026 M/8mM_slack.txt'
        '03/27','2mM','03 27 2026 M/2mM_slack.txt'
        '03/27','PNB','03 27 2026 M/8mM_slack_PNB_Mava.txt'
        '04/03','8mM','04 03 2026 F/8mM_slack.txt'
        '04/03','2mM','04 03 2026 F/2mM_slack.txt'
        '04/03','PNB','04 03 2026 F/8mM_slack_PNB_Mava.txt'
        '04/10','8mM','04 10 2026 Male 2-8/8mM_slack.txt'
        '04/10','2mM','04 10 2026 Male 2-8/2mM_slack.txt'
        '04/10','PNB','04 10 2026 Male 2-8/8mM_slack_PNB_Mava.txt' };

FV = struct('day',{},'cond',{},'v',{},'dur',{},'Fss',{},'Fpl',{},'Fnorm',{}, ...
            'L1',{},'L2',{},'t',{},'F',{},'L',{},'t1',{},'t2',{});
for i = 1:size(SRC,1)
    M = readmatrix(fullfile('data', SRC{i,3}));
    i0 = find(M(:,1) >= -0.1, 1);
    t = M(i0:end,1)/1000; L = M(i0:end,2); F = M(i0:end,3);
    dt = median(diff(t));
    v  = gradient(movmean(L, max(1,round(5e-4/dt))), t);
    isRamp = t >= 0 & t <= 0.45 & v < -0.5 & v > -5;
    d = diff([0; isRamp; 0]); rs = find(d==1); re = find(d==-1)-1;
    if isempty(rs); continue; end
    [dur, iL] = max(t(re)-t(rs));
    if dur < 0.010; continue; end
    t1 = t(rs(iL)); t2 = t(re(iL));
    Fss = median(F(t >= 0.002 & t < 0.040));
    w   = t >= t1+0.005 & t <= t2;
    FV(end+1) = struct('day',SRC{i,1},'cond',SRC{i,2}, ...
        'v',(L(re(iL))-L(rs(iL)))/(t2-t1), 'dur',dur, 'Fss',Fss, ...
        'Fpl',mean(F(w)), 'Fnorm',mean(F(w))/Fss, ...
        'L1',L(rs(iL)),'L2',L(re(iL)), ...
        't',t(t<0.12),'F',F(t<0.12),'L',L(t<0.12),'t1',t1,'t2',t2); %#ok<SAGROW>
end

fprintf('=== FV2 ramps found ===\n');
fprintf('%-7s %-5s %8s %8s %9s %9s %8s %14s\n','day','cond','v(ML/s)','dur(ms)', ...
        'Fss(kPa)','Fpl(kPa)','P/P0','L span (ML)');
for i = 1:numel(FV)
    fprintf('%-7s %-5s %8.2f %8.1f %9.2f %9.2f %8.3f  %.3f -> %.3f\n', FV(i).day, ...
        FV(i).cond, FV(i).v, 1000*FV(i).dur, FV(i).Fss, FV(i).Fpl, FV(i).Fnorm, ...
        FV(i).L1, FV(i).L2);
end
fprintf('03/27 has no ramp: the manoeuvre was added after that session.\n');

%% ---- 2. how much of it is passive? ------------------------------------
% The PNB+Mava run drives the identical ramp, so the passive/viscous force
% borne at -2 ML/s is measured on the same preparation, not modelled.
days = unique({FV.day},'stable');
fprintf('\n=== passive share of the ramp, measured on the same preparation ===\n');
fprintf('%-7s %-5s %10s %10s %12s %12s\n','day','cond','Fpl','Fpl_PNB','P/P0 raw','P/P0 active');
ACT = struct('day',{},'cond',{},'raw',{},'act',{});
for p = 1:numel(days)
    ip = find(strcmp({FV.day},days{p}) & strcmp({FV.cond},'PNB'));
    if isempty(ip); continue; end
    for c = {'8mM','2mM'}
        i = find(strcmp({FV.day},days{p}) & strcmp({FV.cond},c{1}));
        if isempty(i); continue; end
        act = (FV(i).Fpl - FV(ip).Fpl) / (FV(i).Fss - FV(ip).Fss);
        fprintf('%-7s %-5s %10.2f %10.2f %12.3f %12.3f\n', days{p}, c{1}, ...
                FV(i).Fpl, FV(ip).Fpl, FV(i).Fnorm, act);
        ACT(end+1) = struct('day',days{p},'cond',c{1},'raw',FV(i).Fnorm,'act',act); %#ok<SAGROW>
    end
end
fprintf(['  The PNB+Mava ramp bears %.0f-%.0f %% of its own isometric force -- passive\n' ...
         '  and viscous force are NOT negligible at 2 ML/s and must be handled\n' ...
         '  explicitly, whichever way the comparison is made.\n'], ...
         100*min([FV(strcmp({FV.cond},'PNB')).Fnorm]), ...
         100*max([FV(strcmp({FV.cond},'PNB')).Fnorm]));

%% ---- 1. against the Baker curves --------------------------------------
AllV = [0 0.5 1 2 3 4 5 6];
BK = { 'Baker2022 8 mM',      [56.4048 51.8120 37.4459 17.8025 11.4430 6.2643 3.2759 2.2120]
       'AB2021 8 mM',         [67.3942 60.9885 42.7059 18.8539 12.5956 6.4426 3.5136 1.8556]
       'Baker2022 2 mM',      [61.3192 47.4794 31.3870 10.2112  6.3895 2.5781 1.5940 1.2117]
       'AB2021 2 mM',         [65.4207 43.9473 27.8404  8.0630  5.1020 2.9306 1.8704 1.0416] };
iv2 = find(AllV == 2);
fprintf('\n=== the 2 ML/s point: in-preparation FV2 vs the borrowed Baker curves ===\n');
fprintf('%-24s %10s %10s\n','source','P/P0 @2','8mM basis');
for j = 1:size(BK,1)
    fprintf('%-24s %10.3f %10.1f kPa\n', BK{j,1}, BK{j,2}(iv2)/BK{j,2}(1), BK{j,2}(1));
end
for i = 1:numel(ACT)
    fprintf('%-24s %10.3f %10s\n', sprintf('FV2 %s %s (raw)',ACT(i).day,ACT(i).cond), ACT(i).raw, '-');
end
b8 = [BK{1,2}(iv2)/BK{1,2}(1), BK{2,2}(iv2)/BK{2,2}(1)];
b2 = [BK{3,2}(iv2)/BK{3,2}(1), BK{4,2}(iv2)/BK{4,2}(1)];
f8 = [ACT(strcmp({ACT.cond},'8mM')).raw];  f8a = [ACT(strcmp({ACT.cond},'8mM')).act];
f2 = [ACT(strcmp({ACT.cond},'2mM')).raw];  f2a = [ACT(strcmp({ACT.cond},'2mM')).act];
fprintf('\n%-14s %14s %14s %12s\n','','Baker (n=2 sets)','FV2 (n=2 preps)','agreement');
fprintf('%-14s %6.3f +- %5.3f %6.3f +- %5.3f %11.2f sd\n','8 mM raw', ...
        mean(b8),std(b8), mean(f8),std(f8), abs(mean(b8)-mean(f8))/sqrt(var(b8)+var(f8)));
fprintf('%-14s %6.3f +- %5.3f %6.3f +- %5.3f %11.2f sd\n','8 mM active', ...
        mean(b8),std(b8), mean(f8a),std(f8a), abs(mean(b8)-mean(f8a))/sqrt(var(b8)+var(f8a)));
fprintf('%-14s %6.3f +- %5.3f %6.3f +- %5.3f %11.2f sd\n','2 mM raw', ...
        mean(b2),std(b2), mean(f2),std(f2), abs(mean(b2)-mean(f2))/sqrt(var(b2)+var(f2)));

%% ---- is the ramp force actually at a plateau? --------------------------
% A 34 ms ramp entered through a step is only a force-velocity point if force
% has settled. Fit a line across the ramp window and express the drift as a
% fraction of the plateau.
fprintf('\n=== is the ramp force settled? (drift across the 34 ms window) ===\n');
fprintf('%-7s %-5s %12s %14s\n','day','cond','slope(kPa/s)','drift / Fpl');
for i = 1:numel(FV)
    w = FV(i).t >= FV(i).t1+0.005 & FV(i).t <= FV(i).t2;
    pf = polyfit(FV(i).t(w), FV(i).F(w), 1);
    fprintf('%-7s %-5s %12.1f %13.1f%%\n', FV(i).day, FV(i).cond, pf(1), ...
            100*pf(1)*(FV(i).t2-FV(i).t1-0.005)/FV(i).Fpl);
end

%% ---- the honest spread of a normalised F-V point ----------------------
% FeatureWeighting takes FV's spread BETWEEN THE TWO BAKER SOURCE SETS (11.3 %
% aggregated over velocities). FV2 adds two independent PREPARATIONS at the
% same velocity, which is the spread that actually applies when the model is
% asked to predict a new fibre.
pool2 = [b8 f8];
fprintf('\n=== spread of P/P0 at 2 ML/s, 8 mM ===\n');
fprintf('  between the two Baker source sets  : %.1f%%  (n = 2 curve fits)\n', ...
        100*std(b8)/mean(b8));
fprintf('  between the two 2026 preparations  : %.1f%%  (n = 2 preparations)\n', ...
        100*std(f8)/mean(f8));
fprintf('  all four together                  : %.1f%%  (mean %.3f)\n', ...
        100*std(pool2)/mean(pool2), mean(pool2));
fprintf('  -> a normalised F-V point is reproducible to ~%.0f%%, not to the ~8%%\n', ...
        100*std(pool2)/mean(pool2));
fprintf('     that two fits of the same lab''s data suggest.\n');

%% ---- 3. the ATP effect on the F-V point -------------------------------
fprintf('\n=== ATP effect on P/P0 at 2 ML/s ===\n');
fprintf('%-7s %10s %10s %10s\n','day','8 mM','2 mM','ratio');
rr = [];
for p = 1:numel(days)
    i8 = find(strcmp({ACT.day},days{p}) & strcmp({ACT.cond},'8mM'));
    i2 = find(strcmp({ACT.day},days{p}) & strcmp({ACT.cond},'2mM'));
    if isempty(i8) || isempty(i2); continue; end
    fprintf('%-7s %10.3f %10.3f %10.3f\n', days{p}, ACT(i8).raw, ACT(i2).raw, ACT(i2).raw/ACT(i8).raw);
    rr(end+1) = ACT(i2).raw/ACT(i8).raw; %#ok<SAGROW>
end
fprintf('%-7s %10s %10s %10.3f  (Baker sets: %.3f, %.3f)\n','mean','','', mean(rr), ...
        b2(1)/b8(1), b2(2)/b8(2));

%% ---- figure ------------------------------------------------------------
f = figure('Position',[80 80 1180 460],'Color','w'); tiledlayout(f,1,3,'Padding','compact');

nexttile; hold on; box on
cols = lines(4); k = 0;
for i = 1:numel(FV)
    if strcmp(FV(i).cond,'PNB'); sty = ':'; else; sty = '-'; end
    k = k+1;
    plot(1000*FV(i).t, movmean(FV(i).F,41)/FV(i).Fss, sty, 'LineWidth',1.3, ...
         'Color', cols(1+mod(k-1,4),:), 'DisplayName', sprintf('%s %s',FV(i).day,FV(i).cond));
end
xline(52.8,'k--','HandleVisibility','off'); xline(87.4,'k--','HandleVisibility','off');
xlabel('time (ms)'); ylabel('F / F_{ss22}'); xlim([0 120]); ylim([-0.1 1.25]);
title('the FV2 manoeuvre (10 ms smoothing; dashed = the -2 ML/s ramp)');
legend('Location','southwest','FontSize',7);

nexttile; hold on; box on
for j = 1:2
    plot(BK{j,2}/BK{j,2}(1), AllV, 'o-', 'LineWidth',1.4, 'DisplayName',BK{j,1});
end
for j = 3:4
    plot(BK{j,2}/BK{j,2}(1), AllV, 's--', 'LineWidth',1.2, 'DisplayName',BK{j,1});
end
plot(f8, 2*ones(size(f8)), 'k^','MarkerFaceColor','k','MarkerSize',9,'DisplayName','FV2 8 mM');
plot(f2, 2*ones(size(f2)), 'kv','MarkerFaceColor',[.6 .6 .6],'MarkerSize',9,'DisplayName','FV2 2 mM');
set(gca,'YDir','reverse'); xlabel('P / P_0'); ylabel('shortening velocity (ML/s)');
title('in-prep point vs the borrowed curves'); legend('Location','southeast','FontSize',7);

nexttile; hold on; box on
bar([mean(b8) mean(f8) mean(f8a); mean(b2) mean(f2) mean(f2a)]');
set(gca,'XTick',1:3,'XTickLabel',{'Baker','FV2 raw','FV2 active'});
ylabel('P / P_0 at 2 ML/s'); legend({'8 mM','2 mM'},'Location','northeast');
title('agreement, and what the passive correction costs');

exportgraphics(f, fullfile(resDir,'fv_consistency.png'), 'Resolution',150);
save(fullfile(resDir,'fv_consistency.mat'), 'FV','ACT','BK','AllV','b8','b2','f8','f2','f8a','f2a');
fprintf('\nsaved -> %s\n', fullfile(resDir,'fv_consistency.mat'));
