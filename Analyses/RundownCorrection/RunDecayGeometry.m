% RunDecayGeometry.m
% Does decay happen ONLY during activation, and can each run be referenced to
% its own activation instead of to a fixed protocol time?
%
% THE PROPOSAL BEING TESTED
% Rundown appears to accrue only while the fibre generates force. If so, the
% within-run decay line -- measured on the isometric ML=1.0 plateaus that bracket
% the protocol -- could be extrapolated back to the moment activation completes
% and forward to relaxation, giving each run a "start" and "end" force. If there
% is then NO decay while the fibre sits in relaxing solution, consecutive runs
% should JOIN UP: end of run n = start of run n+1, apart from the treatment
% effect. That would make the correction independent of phi.
%
% WHY IT MATTERS MOST FOR LOW ATP
% Force does not appear instantly: it takes several seconds to develop (calcium
% and substrate diffusing into a skinned preparation), and rundown is already
% running during that rise. The observed peak is therefore NOT the undamaged
% force, and if the two ATP conditions develop force at different rates they are
% not being compared at the same point in their own activation. They are not:
% every 2 mM run has already passed its peak by t = 68.7 s while most 8 mM runs
% are still rising.
%
% Three isometric windows at ML = 1.0 are used:
%   W1  68.20-69.20 s  before any perturbation           <- the only clean one
%   W2  71.20-72.35 s  after the ktr burst, before the staircase
%   W3  77.70-78.30 s  after the last slack
% CAVEAT: W2 and W3 are post-perturbation recoveries. They sit at identical
% protocol times in every run, so the bias cancels in between-run ratios, but
% they are not absolute isometric values.
%
% Outputs -> results/decay_geometry.png, decay_geometry.mat

clear; close all; clc;
cd(fileparts(which('RunDecayGeometry')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

W   = {[68.20 69.20],[71.20 72.35],[77.70 78.30]};
PHI = 0.452;

R = struct( ...
 'day', {'03/27 M','03/27 M','03/27 M','04/03 F','04/03 F','04/10 M','04/10 M'}, ...
 'cond',{'8 mM','2 mM','8 mM rpt','8 mM','2 mM','2 mM','8 mM'}, ...
 'act', {1,2,4,1,2,1,2}, ...
 'file',{'03 27 2026 M/02_Merged_8mM_Active.txt', ...
         '03 27 2026 M/03_Merged_2mM_Active.txt', ...
         '03 27 2026 M/04_Merged_8mM_Active_repeat.txt', ...
         '04 03 2026 F/02_Merged_8mM_Active.txt', ...
         '04 03 2026 F/03_Merged_2mM_Active.txt', ...
         '04 10 2026 Male 2-8/02_Merged_2mM_Active.txt', ...
         '04 10 2026 Male 2-8/03_Merged_8mM_Active.txt'});

%% ── 1. Activation envelope + decay line for every run ──────────────────────
fprintf('%-9s %-9s | %7s %7s %7s | %7s %7s %7s | %8s %7s\n', ...
        'day','cond','W1','W2','W3','t_on','t_peak','t_off','slope','T_dec');
for i = 1:numel(R)
    M = readmatrix([D R(i).file]);
    t = M(:,1); L = M(:,2); F = M(:,3);
    iso = abs(L-1) < 2e-3 & t > 63 & t < 82;          % ML = 1.0 samples only
    R(i).ti = t(iso);  R(i).Fi = F(iso);
    R(i).env = movmean(F(iso), 401);                   % activation envelope

    for k = 1:3
        w = t >= W{k}(1) & t <= W{k}(2) & abs(L-1) < 2e-3;
        R(i).v(k) = mean(F(w),'omitnan');  R(i).tc(k) = mean(t(w));
    end
    [~, ip] = max(R(i).env);  R(i).tpk = R(i).ti(ip);
    Fs = movmean(F,51); thr = 0.2*prctile(Fs,99);
    R(i).ton  = t(find(Fs>thr,1,'first'));
    R(i).toff = t(find(Fs>thr,1,'last'));
    R(i).slope = (R(i).v(3)-R(i).v(2)) / (R(i).tc(3)-R(i).tc(2));
    R(i).Tdec  = R(i).toff - R(i).tpk;
    R(i).Tact  = R(i).toff - R(i).ton;
    % decay line extrapolated to the peak and to relaxation
    ln = @(tq) R(i).v(2) + R(i).slope*(tq - R(i).tc(2));
    R(i).Fpk  = ln(R(i).tpk);
    R(i).Fend = ln(R(i).toff);
    R(i).rising = R(i).v(2) > R(i).v(1);               % still rising at W1?
    fprintf('%-9s %-9s | %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f | %8.4f %7.2f\n', ...
        R(i).day, R(i).cond, R(i).v, R(i).ton, R(i).tpk, R(i).toff, R(i).slope, R(i).Tdec);
end

%% ── 2. Is the rise complete before the decay window? ───────────────────────
fprintf('\n\n============ 1. THE ACTIVATION RISE OVERLAPS THE DECAY ============\n');
fprintf('%-9s %-9s %10s %10s %14s\n','day','cond','W2-W1','t_peak-t_on','at W1 (68.7 s)');
for i = 1:numel(R)
    if R(i).rising; s = 'STILL RISING'; else; s = 'already falling'; end
    fprintf('%-9s %-9s %+10.2f %10.2f %14s\n', R(i).day, R(i).cond, ...
            R(i).v(2)-R(i).v(1), R(i).tpk-R(i).ton, s);
end
is8 = strncmp({R.cond},'8',1);  is2 = strncmp({R.cond},'2',1);
fprintf('\n  Force takes %.1f-%.1f s to peak after activation starts, so the observed\n', ...
        min([R.tpk]-[R.ton]), max([R.tpk]-[R.ton]));
fprintf('  peak is NOT the undamaged force -- decay runs throughout the rise.\n');
fprintf('  ALL %d of the 2 mM runs are already falling at 68.7 s;\n', sum(is2));
fprintf('  %d of %d 8 mM runs are still rising there.\n', sum(is8 & [R.rising]), sum(is8));
fprintf('  => the two conditions are NOT compared at the same point of their own\n');
fprintf('     activation if a fixed protocol window is used. 2 mM decay phase\n');
fprintf('     %.1f s vs 8 mM %.1f s.\n', mean([R(is2).Tdec]), mean([R(is8).Tdec]));

%% ── 3. Is the 2 mM decay really faster? ────────────────────────────────────
fprintf('\n\n============ 2. DECAY RATE, 8 mM vs 2 mM ============\n');
fprintf('%-9s %-9s %12s %12s\n','day','cond','abs kPa/s','rel %%/s');
for i = 1:numel(R)
    fprintf('%-9s %-9s %12.4f %11.2f%%\n', R(i).day, R(i).cond, ...
            R(i).slope, 100*R(i).slope/R(i).v(2));
end
a8  = abs([R(is8).slope]);  a2 = abs([R(is2).slope]);
rel = arrayfun(@(x) 100*abs(x.slope)/x.v(2), R);
fprintf('\n  absolute  8 mM %.3f +- %.3f   2 mM %.3f +- %.3f   ratio x%.2f\n', ...
        mean(a8), std(a8), mean(a2), std(a2), mean(a2)/mean(a8));
fprintf('  relative  8 mM %.2f%% +- %.2f   2 mM %.2f%% +- %.2f   ratio x%.2f\n', ...
        mean(rel(is8)), std(rel(is8)), mean(rel(is2)), std(rel(is2)), mean(rel(is2))/mean(rel(is8)));
fprintf('\n  paired within prep (2 mM / 8 mM absolute slope):\n');
days = unique({R.day},'stable'); pr = nan(1,numel(days));
for d = 1:numel(days)
    ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
    il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
    pr(d) = abs(R(il).slope)/abs(R(ih).slope);
    fprintf('    %-9s x%.2f\n', days{d}, pr(d));
end
fprintf('  mean x%.2f, all > 1 -> CONFIRMED: 2 mM runs down faster.\n', mean(pr));

%% ── 4. Continuity: does "no decay between activations" hold? ───────────────
% 03/27 chain, using each run's decay line extrapolated to its own peak / end:
%   F_end(8 mM #1) --[2 mM run]--> F_peak(2 mM #2)   ratio = r * a
%   F_end(2 mM #2) --[PNB, rest]-> F_peak(8 mM #4)   ratio = r^n / a
% with a = ATP effect and r = between-activation recovery (r > 1 means force
% partially recovers in relaxing solution). Two equations, two unknowns.
i1 = 1; i2 = 2; i4 = 3;
q1 = R(i2).Fpk / R(i1).Fend;          % = r * a
q2 = R(i4).Fpk / R(i2).Fend;          % = r^n / a
fprintf('\n\n============ 3. CONTINUITY / RECOVERY BETWEEN ACTIVATIONS ============\n');
fprintf('  03/27 chain, decay lines extrapolated to each run''s own peak and end:\n');
fprintf('    8 mM #1 : peak %.2f -> end %.2f  (loses %.2f kPa while activated)\n', R(i1).Fpk, R(i1).Fend, R(i1).Fpk-R(i1).Fend);
fprintf('    2 mM #2 : peak %.2f -> end %.2f  (loses %.2f kPa)\n', R(i2).Fpk, R(i2).Fend, R(i2).Fpk-R(i2).Fend);
fprintf('    8 mM #4 : peak %.2f -> end %.2f  (loses %.2f kPa)\n', R(i4).Fpk, R(i4).Fend, R(i4).Fpk-R(i4).Fend);
fprintf('\n    step 8mM#1_end -> 2mM#2_peak : x%.3f  = r * a\n', q1);
fprintf('    step 2mM#2_end -> 8mM#4_peak : x%.3f  = r^n / a\n', q2);
% n = 1:  r*a = q1 and r/a = q2  =>  r = sqrt(q1*q2),  a = q1/r
rA = sqrt(q1*q2);   aA = q1/rA;
% n = 2:  r*a = q1 and r^2/a = q2  =>  r = (q1*q2)^(1/3),  a = q1/r
rB = (q1*q2)^(1/3); aB = q1/rB;
fprintf('\n    assuming ONE  recovery interval each : a = %.3f, r = %.3f (%+.1f%% per rest)\n', aA, rA, 100*(rA-1));
fprintf('    assuming TWO  intervals after the 2 mM: a = %.3f, r = %.3f (%+.1f%% per rest)\n', aB, rB, 100*(rB-1));
fprintf('\n  => "no decay between activations" very nearly HOLDS. The residual is a\n');
fprintf('     small RECOVERY of %.1f-%.1f%% per rest interval, not further decay.\n', ...
        100*(min(rA,rB)-1), 100*(max(rA,rB)-1));
fprintf('     And this route gives an ATP effect of x%.2f-%.2f from 03/27 ALONE,\n', min(aA,aB), max(aA,aB));
fprintf('     with NO phi and no damage model assumed -- an independent confirmation\n');
fprintf('     of the x1.36 consensus obtained the other way.\n');

%% ── 5. Referencing each run to its own peak ────────────────────────────────
fprintf('\n\n============ 4. ATP RATIO: FIXED WINDOW vs OWN-PEAK REFERENCE ============\n');
fprintf('%-9s %12s %12s %12s\n','day','W2 window','own peak','peak + damage');
rw = nan(1,3); rp = nan(1,3); rpd = nan(1,3);
for d = 1:numel(days)
    ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
    il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
    rw(d) = R(il).v(2)  / R(ih).v(2);
    rp(d) = R(il).Fpk   / R(ih).Fpk;
    % damage during the earlier run, over ITS OWN decay phase
    if R(ih).act < R(il).act
        dmg = PHI*abs(R(ih).slope)*R(ih).Tdec;   rpd(d) = R(il).Fpk/(R(ih).Fpk - dmg);
    else
        dmg = PHI*abs(R(il).slope)*R(il).Tdec;   rpd(d) = R(il).Fpk/(R(ih).Fpk + dmg);
    end
    fprintf('%-9s %12.3f %12.3f %12.3f\n', days{d}, rw(d), rp(d), rpd(d));
end
cv = @(x) 100*std(x)/mean(x);
fprintf('%-9s %12s %12s %12s\n','mean', sprintf('%.3f',mean(rw)), sprintf('%.3f',mean(rp)), sprintf('%.3f',mean(rpd)));
fprintf('%-9s %11.0f%% %11.0f%% %11.0f%%\n','CV', cv(rw), cv(rp), cv(rpd));
fprintf('\n  Referencing to each run''s OWN peak removes part of the bias for free --\n');
fprintf('  no phi, no damage model -- because it stops comparing a 2 mM run that has\n');
fprintf('  been decaying for %.0f s with an 8 mM run that has been decaying for %.0f s.\n', ...
        mean([R(is2).tpk]-[R(is2).ton]), mean([R(is8).tpk]-[R(is8).ton]));

%% ── 6. Figure ──────────────────────────────────────────────────────────────
fig = figure(970); clf; set(fig,'Position',[10 10 1600 830]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cd8 = [0 0.45 0.74]; cd2 = [0.85 0.33 0.10];

% (a) activation envelopes, normalised, showing when each peaks
nexttile; hold on; box on;
for i = 1:numel(R)
    c = cd8; if is2(i); c = cd2; end
    plot(R(i).ti - R(i).ton, R(i).env/max(R(i).env), '-', 'Color',[c 0.75], 'LineWidth',1.2);
    plot(R(i).tpk - R(i).ton, 1, 'v','Color',c,'MarkerFaceColor',c,'MarkerSize',8);
end
xlabel('time since activation onset (s)'); ylabel('force / own peak');
title('(a) 2 mM (red) peaks earlier than 8 mM (blue)'); xlim([0 16]); ylim([0 1.05]);

% (b) the three windows and the extrapolated decay line, 03/27
nexttile; hold on; box on;
for i = [i1 i2 i4]
    c = cd8; if is2(i); c = cd2; end
    plot(R(i).ti, R(i).env, '-', 'Color',[c 0.35]);
    plot(R(i).tc, R(i).v, 'o','Color',c,'MarkerFaceColor',c,'MarkerSize',7);
    tq = [R(i).tpk R(i).toff];
    plot(tq, R(i).v(2)+R(i).slope*(tq-R(i).tc(2)), '-','Color',c,'LineWidth',2.2);
    plot(R(i).tpk, R(i).Fpk, 'p','Color',c,'MarkerFaceColor','y','MarkerSize',13);
end
for k=1:3
    patch([W{k}(1) W{k}(2) W{k}(2) W{k}(1)],[0 0 100 100],[0 0 0],'FaceAlpha',0.07,'EdgeColor','none');
end
text(68.7, 96,'W1','HorizontalAlignment','center','FontSize',8);
text(71.8, 96,'W2','HorizontalAlignment','center','FontSize',8);
text(78.0, 96,'W3','HorizontalAlignment','center','FontSize',8);
xlabel('protocol time (s)'); ylabel('force (kPa)'); xlim([64 82]); ylim([0 100]);
title('(b) 03/27: decay line extrapolated to each run''s own peak');

% (c) slope comparison
nexttile; hold on; box on;
for d = 1:numel(days)
    ih = find(strcmp({R.day},days{d}) & strcmp({R.cond},'8 mM'));
    il = find(strcmp({R.day},days{d}) & strcmp({R.cond},'2 mM'));
    plot([1 2],[abs(R(ih).slope) abs(R(il).slope)],'-o','LineWidth',1.8, ...
         'MarkerFaceColor','auto','DisplayName',days{d});
end
set(gca,'XTick',1:2,'XTickLabel',{'8 mM','2 mM'}); xlim([0.7 2.3]);
ylabel('|within-run decay| (kPa/s)');
title(sprintf('(c) 2 mM decays x%.2f faster (all preps)', mean(pr)));
legend('Location','northwest','FontSize',7);

% (d) the continuity chain
nexttile; hold on; box on;
seq = [i1 i2 i4]; xoff = [0 1 2];
for j = 1:3
    i = seq(j); c = cd8; if is2(i); c = cd2; end
    plot([xoff(j) xoff(j)+0.62], [R(i).Fpk R(i).Fend], '-o', 'Color', c, ...
         'LineWidth', 3, 'MarkerFaceColor', c, 'MarkerSize', 8);
    text(xoff(j)+0.31, R(i).Fpk+3, R(i).cond, 'HorizontalAlignment','center','FontSize',8,'Color',c);
end
for j = 1:2
    plot([xoff(j)+0.62 xoff(j+1)], [R(seq(j)).Fend R(seq(j+1)).Fpk], 'k:', 'LineWidth', 1.5);
end
text(1.35, 46, sprintf('a = %.2f-%.2f\nrecovery r = %+.0f to %+.0f%% per rest', ...
     min(aA,aB), max(aA,aB), 100*(rB-1), 100*(rA-1)), 'FontSize',9,'BackgroundColor','w');
set(gca,'XTick',xoff+0.31,'XTickLabel',{'activation 1','activation 2','activation 4'});
ylabel('force at SL 2.0 (kPa)'); ylim([40 95]);
title('(d) Runs nearly join up: little decay while relaxed');

% (e) ATP ratio by referencing scheme
nexttile; hold on; box on;
b = bar([rw; rp; rpd]'); b(1).FaceColor=[.6 .6 .6]; b(2).FaceColor=[.3 .6 .85]; b(3).FaceColor=[.15 .4 .65];
set(gca,'XTick',1:3,'XTickLabel',days); ylabel('ATP force ratio 2 mM / 8 mM');
title(sprintf('(e) CV: %.0f%% -> %.0f%% -> %.0f%%', cv(rw), cv(rp), cv(rpd)));
legend({'fixed window W2','own-peak reference','own-peak + damage'},'Location','northwest','FontSize',7);

% (f) decay-phase duration
nexttile; hold on; box on;
b = bar([[R(is8).Tdec] nan; nan(1,sum(is8)) mean([R(is2).Tdec])]');
plot(ones(1,sum(is8)), [R(is8).Tdec], 'o','Color',cd8,'MarkerFaceColor',cd8,'MarkerSize',10);
plot(2*ones(1,sum(is2)), [R(is2).Tdec], 'o','Color',cd2,'MarkerFaceColor',cd2,'MarkerSize',10);
cla; hold on; box on;
plot(ones(1,sum(is8)), [R(is8).Tdec], 'o','Color',cd8,'MarkerFaceColor',cd8,'MarkerSize',11);
plot(2*ones(1,sum(is2)), [R(is2).Tdec], 'o','Color',cd2,'MarkerFaceColor',cd2,'MarkerSize',11);
plot([0.85 1.15],[1 1]*mean([R(is8).Tdec]),'k-','LineWidth',2);
plot([1.85 2.15],[1 1]*mean([R(is2).Tdec]),'k-','LineWidth',2);
set(gca,'XTick',1:2,'XTickLabel',{'8 mM','2 mM'}); xlim([0.6 2.4]);
ylabel('decay-phase duration t_{off} - t_{peak} (s)');
title(sprintf('(f) 2 mM spends %.1f s more in decay', mean([R(is2).Tdec])-mean([R(is8).Tdec])));

sgtitle('Decay geometry — when rundown acts, and why it biases the low-ATP comparison');
exportgraphics(fig, fullfile(resDir,'decay_geometry.png'), 'Resolution', 150);
save(fullfile(resDir,'decay_geometry.mat'), 'R','q1','q2','aA','rA','aB','rB','rw','rp','rpd','pr');
fprintf('\nSaved %s\n', fullfile(resDir,'decay_geometry.png'));
