% CompareRecoveryRates.m
% =========================================================================
% Does the ktr protocol's force redevelopment really have the same rate as
% the post-slack redevelopment and the post-restretch recovery inside the
% slack protocol?  Three protocol days, 8 mM ATP only.
%
% THREE WINDOWS per dataset
%   (S) post-SLACK      window = [end of release ramp .. start of restretch ramp]
%                       time origin = START of the release ramp (project
%                       convention, see extractSlackAttributes) so the fitted
%                       t0 is the physical slack take-up dead time.
%                       Force starts at ZERO -> baseline B = 0.
%   (R) post-RESTRETCH  window = [deepest valley after the restretch ramp .. next slack]
%                       Force NEVER reaches zero -> baseline B = that valley.
%   (K) ktr protocol    window = [end of the final 1.05->1.00 release .. end]
%                       time origin = start of that release ramp.
%                       Force starts at ZERO -> baseline B = 0.
%
% WHICH MINIMUM (this matters -- there are two valleys)
%   The restretch ramp produces:  peak1 (during the ramp, ~1.20 F_iso)
%                              -> vall_y  (first dip, still during the ramp)
%                              -> peak2   (at the ramp end, ~1.00 F_iso)
%                              -> vall2   (deepest, ~7 ms AFTER the ramp ends,
%                                          ~0.81 F_iso)  <-- THIS is the baseline
%                              -> monotone recovery to steady.
%   extractSlackAttributes' `vall_y` is the FIRST dip and is NOT the start of
%   the recovery. Using it would understate the amplitude and start the fit in
%   the middle of the ringing. This script uses vall2 = min over the 80 ms
%   after the ramp ends, and reports where it was found.
%
% RATE ESTIMATORS -- four, because THE ESTIMATOR ITSELF BIASES THE ANSWER
%   All use the plateau P, which is MEASURED, never fitted:
%     postSlack / postRestretch : median of the last 15% of the window
%     ktr                       : 1.0, i.e. the pre-release isometric force
%                                 (the trace is normalised to it by
%                                 CreateKtrProtocolVelocityTable, and the
%                                 muscle is back at the same length) -- the
%                                 58.5 ms window ends before the plateau, so
%                                 taking the last samples would understate A.
%   k_63    model free: 1/(t63 - t_on), first crossings of B+0.05*(P-B) and
%           B+(1-1/e)*(P-B) on a 2 ms median-smoothed trace. This is the
%           "2/3 of (max-min)" estimator, and it is the most robust one here.
%   k_fixA  F = B + A*(1-exp(-k*max(t-t0,0))) with A FIXED at P-B, [k t0] free.
%   k_fit   same but with A FREE -- the definition as usually written.
%   k_fit58 k_fit restricted to the first 58.5 ms after onset, i.e. exactly
%           what a ktr recording of this window would have been able to report.
%           Comparing k_fit58 with k_fit on the (S)/(R) windows CALIBRATES the
%           bias that the ktr protocol's short window necessarily carries.
%
% COMPARABILITY CONTROLS -- the reason this script exists
%   * A/F_iso is reported for every window. The three rises do NOT have the
%     same amplitude, so a normalised overlay is not a like-for-like test.
%   * (K) is observed for only 58.5 ms vs ~275 ms for (S)/(R). Every window is
%     therefore ALSO fitted over the first 58.5 ms AFTER ITS OWN ONSET
%     (k_fit58) so all three are estimated over an equal span of the rise.
%   * (S) happens at a SHORTER length than (R)/(K) -- overlap differs, so even
%     the rate is not a pure comparison. The hold length is printed per window.
%   * cycle 5 restretches back to 1.00 ML, not 1.10 ML, and is held ~840 ms.
%     Reported, but excluded from the means.
%
% Run: cd(root); addpath(genpath('.')); Analyses/RestretchVsKtrRecovery/CompareRecoveryRates
% =========================================================================

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','RestretchVsKtrRecovery','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

DS = { '03 27', 'protocol_03_27_2026_8mM_slack.mat', 'protocol_03_27_2026_velocitytable_ktr.mat'
       '04 03', 'protocol_04_03_2026_8mM_slack.mat', 'protocol_04_03_2026_velocitytable_ktr.mat'
       '04 10', 'protocol_04_10_2026_8mM_slack.mat', 'protocol_04_10_2026_velocitytable_ktr.mat' };

KTR_SPAN = 0.0585;     % the ktr observation window (s) -- equal-span control
SMOOTH_S = 0.002;      % median-smoothing span for crossing detection (s)
VALL_WIN = 0.080;      % search span for the post-restretch valley (s)

W = struct('name',{},'win',{},'Fiso',{},'ev',{});

for d = 1:size(DS,1)
    S = load(fullfile(root,'data',DS{d,2}));
    K = load(fullfile(root,'data',DS{d,3}));
    svt = S.velocitytable; sdt = S.datatable;
    kvt = K.velocitytable; kdt = K.datatable;

    iS = find(svt(:,2) < -1);  iR = find(svt(:,2) > 1);  nC = min(numel(iS),numel(iR));

    % slack trace is in kPa -> express in units of the pre-protocol isometric force
    t1   = svt(iS(1),1);
    Fiso = median(sdt(sdt(:,1) > t1-0.30 & sdt(:,1) < t1-0.01, 3));
    ts   = sdt(:,1);  Fs = sdt(:,3)/Fiso;
    tk   = kdt(:,1);  Fk = kdt(:,3);          % ktr trace is already F/F_iso

    win = {}; ev = {};
    for c = 1:nC
        % ---------- (S) post-SLACK ----------
        t_org = svt(iS(c),1);                 % time origin: release ramp START
        ta    = svt(iS(c)+1,1);               % window start: release ramp END
        tb    = svt(iR(c),1);
        m = ts >= ta & ts <= tb;
        win{end+1} = mkwin(sprintf('S%d',c),'postSlack',c, ts(m)-t_org, Fs(m), 0, ...
                           svt(iS(c)+1,4), KTR_SPAN, SMOOTH_S); %#ok

        % ---------- (R) post-RESTRETCH ----------
        t_ramp0 = svt(iR(c),1);  t_ramp1 = svt(iR(c)+1,1);
        if c < nC; td = svt(iS(c+1),1); else; td = max(ts); end
        m  = ts >= t_ramp1 & ts <= td;  tt = ts(m); FF = Fs(m);
        mv = tt <= tt(1)+VALL_WIN;  [Bv,iv] = min(FF(mv));      % vall2
        win{end+1} = mkwin(sprintf('R%d',c),'postRestretch',c, tt(iv:end)-tt(iv), FF(iv:end), Bv, ...
                           svt(iR(c)+1,4), KTR_SPAN, SMOOTH_S); %#ok
        % keep the whole restretch event for the "which minimum" panel
        me = ts >= t_ramp0-0.005 & ts <= t_ramp1+0.150;
        ev{end+1} = struct('cyc',c,'t',ts(me)-t_ramp1,'F',Fs(me), ...
                           'tv',tt(iv)-t_ramp1,'Bv',Bv,'tramp',t_ramp0-t_ramp1); %#ok
    end

    % ---------- (K) ktr protocol ----------
    t_org = kvt(end-2,1);  ta = kvt(end-1,1);
    m = tk >= ta & tk <= kvt(end,1);
    win{end+1} = mkwin('K','ktr',0, tk(m)-t_org, Fk(m), 0, kvt(end-1,4), KTR_SPAN, SMOOTH_S);

    W(d).name = DS{d,1};  W(d).win = [win{:}];  W(d).Fiso = Fiso;  W(d).ev = [ev{:}];
end

%% ======================= console report =======================
for d = 1:numel(W)
    fprintf('\n================ %s   (slack F_iso = %.1f kPa) ================\n', W(d).name, W(d).Fiso);
    fprintf('%-16s %6s | %6s %6s %7s | %7s %7s %7s %8s | %6s %6s\n', ...
        'window','ML','B','P','A/Fiso','k_63','k_fixA','k_fit','k_fit58','t0_ms','r2');
    for j = 1:numel(W(d).win)
        w = W(d).win(j);
        fprintf('%-16s %6.3f | %6.3f %6.3f %7.3f | %7.1f %7.1f %7.1f %8.1f | %6.1f %6.3f\n', ...
            [w.type ' ' w.id], w.ML, w.B, w.P, w.A, w.k63, w.kfixA, w.kfit, w.kfit58, 1000*w.t0, w.r2A);
    end
end

fprintf('\n============ MEANS over cycles 1-4 (cycle 5 excluded: returns to 1.00 ML) ============\n');
fprintf('%-8s | %-31s | %-31s | %-31s\n','day','post-SLACK','post-RESTRETCH','ktr protocol');
fprintf('%-8s | %7s %7s %7s %7s | %7s %7s %7s %7s | %7s %7s %7s %7s\n','', ...
    'k63','kfixA','kfit58','A/Fiso','k63','kfixA','kfit58','A/Fiso','k63','kfixA','kfit58','A/Fiso');
SUM = zeros(numel(W),12);
for d = 1:numel(W)
    sel = @(ty) W(d).win(strcmp({W(d).win.type},ty) & [W(d).win.cyc]>=1 & [W(d).win.cyc]<=4);
    a=sel('postSlack'); b=sel('postRestretch'); k=W(d).win(strcmp({W(d).win.type},'ktr'));
    SUM(d,:) = [mean([a.k63]) mean([a.kfixA]) mean([a.kfit58]) mean([a.A]) ...
                mean([b.k63]) mean([b.kfixA]) mean([b.kfit58]) mean([b.A]) ...
                k.k63        k.kfixA        k.kfit58        k.A];
    fprintf('%-8s | %7.1f %7.1f %7.1f %7.3f | %7.1f %7.1f %7.1f %7.3f | %7.1f %7.1f %7.1f %7.3f\n', W(d).name, SUM(d,:));
end
m = mean(SUM,1);
fprintf('%-8s | %7.1f %7.1f %7.1f %7.3f | %7.1f %7.1f %7.1f %7.3f | %7.1f %7.1f %7.1f %7.3f\n','MEAN', m);

fprintf('\n--- DOES ktr OVERLAY THE OTHER TWO?  (k_63, the robust model-free estimator) ---\n');
fprintf('  post-slack     k = %5.1f /s  (tau %4.1f ms)\n', m(1), 1000/m(1));
fprintf('  post-restretch k = %5.1f /s  (tau %4.1f ms)   ratio vs post-slack %.2f\n', m(5), 1000/m(5), m(5)/m(1));
fprintf('  ktr protocol   k = %5.1f /s  (tau %4.1f ms)   ratio vs post-slack %.2f\n', m(9), 1000/m(9), m(9)/m(1));

fprintf('\n--- SHORT-WINDOW BIAS: what the 58.5 ms ktr window does to a KNOWN rate ---\n');
rs = m(3)/m(1); rr = m(7)/m(5);
fprintf('  post-slack     : k_fit58 %5.1f vs k_63 %5.1f  -> the 58.5 ms window reports %.0f%% of the true rate\n', m(3), m(1), 100*rs);
fprintf('  post-restretch : k_fit58 %5.1f vs k_63 %5.1f  -> %.0f%%\n', m(7), m(5), 100*rr);
fprintf('  So the ktr protocol''s own number (%.1f /s) is biased low by the same mechanism;\n', m(11));
fprintf('  bias-corrected it is ~%.1f /s.\n', m(11)/mean([rs rr]));

fprintf('\n--- AMPLITUDE: why a normalised overlay is NOT a like-for-like test ---\n');
fprintf('  post-slack     A = %.3f F_iso   0 -> a SHORT-length plateau (%.0f%% of F_iso)\n', m(4), 100*m(4));
fprintf('  ktr protocol   A = %.3f F_iso   0 -> the FULL isometric force\n', m(12));
fprintf('  post-restretch A = %.3f F_iso   -> %.1fx smaller than post-slack, %.1fx smaller than ktr\n', ...
        m(8), m(4)/m(8), m(12)/m(8));

save(fullfile(resdir,'recovery_rates.mat'),'W','SUM','DS');

%% ======================= FIGURE 1 =======================
gS=[0.10 0.60 0.25]; gR=[0.85 0.15 0.15]; gK=[0 0 0];
fig = figure('Position',[10 10 1860 980],'Color','w');
tl  = tiledlayout(fig,numel(W),4,'TileSpacing','compact','Padding','compact');
for d = 1:numel(W)
    % ---- col 1 : the restretch event, and WHICH minimum is the baseline ----
    nexttile(tl); hold on; box on; grid on;
    for e = 1:numel(W(d).ev)
        if W(d).ev(e).cyc==5, continue; end
        E = W(d).ev(e);
        plot(1000*E.t, E.F, '-', 'Color', [gR*(0.5+0.12*E.cyc) 0.85], 'LineWidth',1.0, ...
             'HandleVisibility', tern(E.cyc==1,'on','off'), 'DisplayName','restretch event');
        plot(1000*E.tv, E.Bv, 'ko','MarkerFaceColor','y','MarkerSize',6, ...
             'HandleVisibility', tern(E.cyc==1,'on','off'), 'DisplayName','vall2 = baseline B');
    end
    xline(0,'k--','HandleVisibility','off'); yline(1,'k:','HandleVisibility','off');
    xlabel('ms rel. end of restretch ramp'); ylabel('F / F_{iso}'); xlim([-25 100]);
    title(sprintf('%s — which minimum?', W(d).name));
    if d==1, legend('Location','southeast','FontSize',8); end

    % ---- col 2 : ABSOLUTE ----
    nexttile(tl); hold on; box on; grid on;
    plotset(W(d).win, gS,gR,gK, false);
    yline(1,'k:','HandleVisibility','off');
    xlim([0 150]); ylim([-0.05 1.15]); xlabel('ms since event'); ylabel('F / F_{iso}');
    title(sprintf('%s — ABSOLUTE', W(d).name));
    if d==1, legend('Location','southeast','FontSize',8); end

    % ---- col 3 : NORMALISED ----
    nexttile(tl); hold on; box on; grid on;
    plotset(W(d).win, gS,gR,gK, true);
    yline(1-exp(-1),'k:','HandleVisibility','off'); yline(1,'k:','HandleVisibility','off');
    xlim([0 150]); ylim([-0.15 1.35]); xlabel('ms since event'); ylabel('(F-B)/(P-B)');
    title(sprintf('%s — NORMALISED', W(d).name));

    % ---- col 4 : rates per cycle ----
    nexttile(tl); hold on; box on; grid on;
    a=W(d).win(strcmp({W(d).win.type},'postSlack'));
    b=W(d).win(strcmp({W(d).win.type},'postRestretch'));
    k=W(d).win(strcmp({W(d).win.type},'ktr'));
    plot([a.cyc],[a.k63],'o-','Color',gS,'LineWidth',1.8,'MarkerFaceColor',gS,'DisplayName','post-slack k_{63}');
    plot([b.cyc],[b.k63],'s-','Color',gR,'LineWidth',1.8,'MarkerFaceColor',gR,'DisplayName','post-restretch k_{63}');
    plot([a.cyc],[a.kfit58],'o:','Color',gS,'LineWidth',1.2,'DisplayName','post-slack, 58 ms window');
    yline(k.k63,'k-','LineWidth',2,'DisplayName','ktr k_{63}');
    xlabel('slack cycle (deeper \rightarrow)'); ylabel('k (s^{-1})');
    xlim([0.5 5.5]); ylim([0 90]); title(sprintf('%s — rate constants', W(d).name));
    if d==1, legend('Location','northwest','FontSize',8); end
end
sgtitle('ktr vs post-slack vs post-restretch — 8 mM, three protocol days');
exportgraphics(fig, fullfile(resdir,'recovery_rates_traces.png'),'Resolution',120);

%% ======================= FIGURE 2 =======================
fig2 = figure('Position',[40 40 1450 520],'Color','w');
tiledlayout(fig2,1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on; grid on;
bar([SUM(:,1) SUM(:,5) SUM(:,9)]); xticks(1:numel(W)); xticklabels({W.name});
ylabel('k_{63} (s^{-1})'); ylim([0 60]);
legend({'post-slack','post-restretch','ktr'},'Location','northwest','FontSize',9);
title('Rate constant (model-free k_{63})');
nexttile; hold on; box on; grid on;
bar([SUM(:,4) SUM(:,8) SUM(:,12)]); xticks(1:numel(W)); xticklabels({W.name});
ylabel('amplitude A / F_{iso}');
legend({'post-slack','post-restretch','ktr'},'Location','northwest','FontSize',9);
title('Amplitude — the three rises are NOT the same size');
nexttile; hold on; box on; grid on;
for d=1:numel(W); for j=1:numel(W(d).win)
    w=W(d).win(j); if w.cyc==5, continue; end
    switch w.type
        case 'postSlack',     scatter(w.k63,w.kfit58,46,gS,'o','filled','HandleVisibility',tern(d==1&&w.cyc==1,'on','off'),'DisplayName','post-slack');
        case 'postRestretch', scatter(w.k63,w.kfit58,46,gR,'s','filled','HandleVisibility',tern(d==1&&w.cyc==1,'on','off'),'DisplayName','post-restretch');
        otherwise,            scatter(w.k63,w.kfit58,90,gK,'p','filled','HandleVisibility',tern(d==1,'on','off'),'DisplayName','ktr');
    end
end; end
lim=[0 70]; plot(lim,lim,'k--','HandleVisibility','off'); xlim(lim); ylim(lim);
xlabel('k_{63}, full window (s^{-1})'); ylabel('k from a 58.5 ms window (s^{-1})');
title('Short-window bias (points below the line = underestimate)');
legend('Location','northwest','FontSize',9);
sgtitle('Summary across the three protocol days (8 mM)');
exportgraphics(fig2, fullfile(resdir,'recovery_rates_summary.png'),'Resolution',120);
disp('DONE');

%% ======================= helpers =======================
function plotset(win, gS,gR,gK, normalise)
    for j = 1:numel(win)
        w = win(j);
        if w.cyc == 5, continue; end
        switch w.type
            case 'postSlack',     co = gS*(0.5+0.12*w.cyc); nm='post-slack';     lw=1.0;
            case 'postRestretch', co = gR*(0.5+0.12*w.cyc); nm='post-restretch'; lw=1.0;
            otherwise,            co = gK;                  nm='ktr protocol';   lw=2.0;
        end
        y = w.F;
        if normalise; y = (w.F - w.B)/max(w.P - w.B, eps); end
        plot(1000*(w.t - w.t(1)), y, '-', 'Color', co, 'LineWidth', lw, ...
             'HandleVisibility', tern(w.cyc<=1,'on','off'), 'DisplayName', nm);
    end
end

function w = mkwin(id,type,cyc,t,F,B,ML,SPAN,SM)
    t=t(:); F=F(:); ok=isfinite(t)&isfinite(F); t=t(ok); F=F(ok);
    w.id=id; w.type=type; w.cyc=cyc; w.t=t; w.F=F; w.B=B; w.ML=ML;
    % --- plateau: MEASURED, never fitted (see header) ---
    if strcmp(type,'ktr')
        w.P = 1.0;                       % pre-release isometric force
        w.Pnote = 'F_iso (window ends before the plateau)';
    else
        w.P = median(F(t >= t(end)-0.15*(t(end)-t(1))));
        w.Pnote = 'measured end-of-window plateau';
    end
    Afix = w.P - B;
    [w.Afree, w.kfit,  w.t0,  w.r2   ] = fitKtr(t,F,B,[],  []);          % A free, full
    [~,       w.kfixA, w.t0A, w.r2A  ] = fitKtr(t,F,B,Afix,[]);          % A fixed, full
    m58 = t <= (t(1)+w.t0) + SPAN;
    if nnz(m58) < 20; m58 = t <= t(1)+SPAN; end
    [w.Afree58, w.kfit58, ~, w.r258] = fitKtr(t(m58),F(m58),B,[],[w.Afree w.kfit w.t0]);
    w.A   = Afix;                        % reported amplitude = measured
    w.k63 = crossRate(t,F,B,w.P,SM);
end

function [A,k,t0,r2] = fitKtr(t,F,B,Afix,q0)
    % F = B + A*(1-exp(-k*max(t-t0,0))). B always FIXED at the baseline.
    % Afix = [] -> A free (3 params);  Afix = value -> A fixed ([k t0] free).
    tt = t - t(1);  y = F - B;
    A0 = max(median(y(max(1,end-round(0.15*numel(y))):end)), 1e-6);
    opt = optimset('Display','off','MaxFunEvals',8e3,'MaxIter',8e3);
    if isempty(Afix)
        if isempty(q0); q0 = [A0, 45, 0.004]; end
        obj = @(q) sum((y - abs(q(1)).*(1-exp(-abs(q(2)).*max(tt-max(q(3),0),0)))).^2);
        q = fminsearch(obj, q0, opt);
        A = abs(q(1)); k = abs(q(2)); t0 = max(q(3),0);
    else
        A = Afix;
        obj = @(q) sum((y - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
        q = fminsearch(obj, [45, 0.004], opt);
        k = abs(q(1)); t0 = max(q(2),0);
    end
    pred = A.*(1-exp(-k.*max(tt-t0,0)));
    r2 = 1 - sum((y-pred).^2)/max(sum((y-mean(y)).^2),eps);
end

function k = crossRate(t,F,B,P,SM)
    n  = max(3, round(SM/median(diff(t))));
    y  = (movmedian(F,n) - B)/max(P-B,eps);
    % onset: LAST time the trace is still below 5% (skips any pre-event samples)
    below = find(y < 0.05);
    if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(y >= 1-exp(-1) & (1:numel(y))' > i_on, 1);
    if isempty(i63); k = NaN; else; k = 1/(t(i63)-t(i_on)); end
end

function o = tern(c,a,b); if c; o=a; else; o=b; end; end
