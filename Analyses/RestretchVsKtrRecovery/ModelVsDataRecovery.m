% ModelVsDataRecovery.m
% =========================================================================
% The data says (CompareRecoveryRates.m / conclusions.md): at 8 mM the ktr
% protocol, the post-slack redevelopment and the post-restretch recovery all
% run at ~41-47 /s. The cleanest pair is ktr vs post-slack CYCLE 2, because
% cycle 2 holds at exactly SL 2.0 um -- the same length the ktr protocol
% redevelops at. Same length, same starting force (zero), same everything.
%
% This script asks: does the model reproduce that pair, and if not, why.
%
% Variants:
%   V1 TNF               baseline (params/ThursdayNightFever.m)
%   V2 +TensionOnly      UseMaxwellTensionOnly = 1  (removes the titin-dashpot
%                        force floor left by the 1.05 -> 1.00 release)
%   V3 +long dwell       V2, but the ktr protocol's hold at 0.80 ML is
%                        stretched from 5 ms to 275 ms so it matches the slack
%                        protocol's low-force dwell. DIAGNOSTIC: if the model's
%                        ktr slows down, its rate is set by how long it sat at
%                        zero force, not by cross-bridge kinetics. The data
%                        excludes any such dependence (10 ms and 275 ms dwells
%                        give the same rate).
%   V4 +A_reg pinned     V2 with v_ref_reg -> huge, so the shortening boost of
%                        the registration-availability gate is removed and
%                        A_reg stays at A0.
%   V5 +no SRX entry     V2 with ksr0 = 0, so heads cannot be parked in SRX
%                        during the low-force dwell.
%
% Run: cd(root); addpath(genpath('.')); Analyses/RestretchVsKtrRecovery/ModelVsDataRecovery
% =========================================================================

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','RestretchVsKtrRecovery','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

SLACK_F = 'protocol_03_27_2026_8mM_slack.mat';
KTR_F   = 'protocol_03_27_2026_velocitytable_ktr.mat';
S = load(fullfile(root,'data',SLACK_F));  svt = S.velocitytable;
K = load(fullfile(root,'data',KTR_F));    kvt = K.velocitytable;

%% ---------------- DATA reference ----------------
iS = find(svt(:,2) < -1);
FisoD = median(S.datatable(S.datatable(:,1) > svt(iS(1),1)-0.30 & ...
                           S.datatable(:,1) < svt(iS(1),1)-0.01, 3));
Wd_s = recoveryWindows(S.datatable(:,1), S.datatable(:,3), svt, 'slack', FisoD);
Wd_k = recoveryWindows(K.datatable(:,1), K.datatable(:,3), kvt, 'ktr',   1.0);

%% ---------------- MODEL variants ----------------
base = getParams(loadParams('ThursdayNightFever'), [], true, false);
base.PlotEachSeparately = 0; base.ghostSave = ''; base.ghostLoad = '';
base.EvalFeatures = 0; base.RunSlackSegments = 'All';

V = { 'V1 TNF',           @(p) p
      'V2 +TensionOnly',  @(p) setf(p,'UseMaxwellTensionOnly',1)
      'V3 +long dwell',   @(p) setf(p,'UseMaxwellTensionOnly',1)
      'V4 +A_reg pinned', @(p) setf(setf(p,'UseMaxwellTensionOnly',1),'v_ref_reg',1e6)
      'V5 +no SRX entry', @(p) setf(setf(p,'UseMaxwellTensionOnly',1),'ksr0',0)
      'V6 +tau_reg=2ms',  @(p) setf(setf(p,'UseMaxwellTensionOnly',1),'tau_reg',0.002) };
LONG_DWELL = 3;                       % which variant gets the stretched ktr hold
DWELL_EXTRA = 0.270;                  % s, to match the slack protocol's dwell

M = struct('name',{},'ws',{},'wk',{},'Fiso',{},'probe',{},'E',{});
for v = 1:size(V,1)
    p = V{v,2}(base);
    p = getParams(p, p.g, true);

    % ---- slack protocol ----
    [~, os] = runSlackExperiment(p);
    FisoM = interp1(os.t, os.Force, svt(iS(1),1)-0.02, 'linear','extrap');
    ws = recoveryWindows(os.t, os.Force, svt, 'slack', FisoM);

    % ---- ktr protocol (optionally with a stretched low-force dwell) ----
    vt = kvt;
    if v == LONG_DWELL; vt(4:end,1) = vt(4:end,1) + DWELL_EXTRA; end
    pk = p; pk.SL0 = 2*vt(1,4); pk = getParams(pk, pk.g, true);
    pk.Velocity = vt(:,2);
    [~, ok] = evaluateModel(str2func(pk.modelFcn), vt(:,1), pk);
    FisoK = interp1(ok.t, ok.Force, vt(2,1), 'linear','extrap');
    wk = recoveryWindows(ok.t, ok.Force, vt, 'ktr', FisoK);

    % ---- full experiment battery, so a "fix" can be judged against the rest ----
    [~, Ev] = evaluateBakersExp([], p, false, false);   % [FV ktr stairs slack]

    M(v).name = V{v,1};  M(v).ws = ws;  M(v).wk = wk;  M(v).Fiso = [FisoM FisoK];
    M(v).probe = probeStates(os, ok, svt, vt, iS);
    M(v).E = Ev;
    M(v).os = os; M(v).ok = ok; M(v).vt = vt;
end

%% ---------------- report ----------------
gsel = @(w,ty,c) w(strcmp({w.type},ty) & [w.cyc]==c);
d_s2 = gsel(Wd_s,'postSlack',2);  d_k = Wd_k;

fprintf('\n=================== THE SL-MATCHED PAIR (both at SL 2.0 um) ===================\n');
fprintf('%-18s | %-24s | %-24s | %s\n','variant','post-slack cycle 2','ktr protocol','ktr/postSlack2');
fprintf('%-18s | %7s %7s %8s | %7s %7s %8s |\n','','k63','kfixA','A/Fiso','k63','kfixA','A/Fiso');
fprintf('%-18s | %7.1f %7.1f %8.3f | %7.1f %7.1f %8.3f | %.2f\n','DATA 03/27', ...
    d_s2.k63,d_s2.kfixA,d_s2.A, d_k.k63,d_k.kfixA,d_k.A, d_k.k63/d_s2.k63);
for v = 1:numel(M)
    s2 = gsel(M(v).ws,'postSlack',2);  kk = M(v).wk;
    fprintf('%-18s | %7.1f %7.1f %8.3f | %7.1f %7.1f %8.3f | %.2f\n', M(v).name, ...
        s2.k63,s2.kfixA,s2.A, kk.k63,kk.kfixA,kk.A, kk.k63/s2.k63);
end

fprintf('\n=================== post-slack k63 per cycle (the SL dependence) ===================\n');
fprintf('%-18s | %s\n','source','cycle 1    2    3    4    5   (SL 2.04 2.00 1.96 1.92 1.88 um)');
pc = @(w) sprintf('%6.1f', arrayfun(@(c) getfield(gsel(w,'postSlack',c),'k63'), 1:5)); %#ok
fprintf('%-18s | %s\n','DATA 03/27', pc(Wd_s));
for v = 1:numel(M); fprintf('%-18s | %s\n', M(v).name, pc(M(v).ws)); end

fprintf('\n=================== post-slack amplitude A/Fiso per cycle ===================\n');
pa = @(w) sprintf('%6.3f', arrayfun(@(c) getfield(gsel(w,'postSlack',c),'A'), 1:5)); %#ok
fprintf('%-18s | %s\n','DATA 03/27', pa(Wd_s));
for v = 1:numel(M); fprintf('%-18s | %s\n', M(v).name, pa(M(v).ws)); end

fprintf('\n=================== post-restretch k63 per cycle (the SEPARATE, unfixed defect) ===================\n');
pr = @(w) sprintf('%6.1f', arrayfun(@(c) getfield(gsel(w,'postRestretch',c),'k63'), 1:5)); %#ok
fprintf('%-18s | %s\n','DATA 03/27', pr(Wd_s));
for v = 1:numel(M); fprintf('%-18s | %s\n', M(v).name, pr(M(v).ws)); end

fprintf('\n=================== COST BATTERY — does the fix break anything else? ===================\n');
fprintf('%-18s | %8s %8s %8s %10s\n','variant','E_FV','E_ktr','E_stairs','E_slack');
for v = 1:numel(M)
    if v == LONG_DWELL; continue; end          % diagnostic only, ktr VT is altered
    fprintf('%-18s | %8.3f %8.3f %8.3f %10.2f\n', M(v).name, M(v).E(1), M(v).E(2), M(v).E(3), M(v).E(4));
end

fprintf('\n=================== WHY: state of the head pool at the start of each rise ===================\n');
fprintf('(fractions at the moment force starts redeveloping)\n');
fprintf('%-18s %-12s | %6s %6s %6s %6s %6s %6s %7s\n','variant','window','PD','SR+SRD','PT','p1_0','p2_0','A_reg','F/Fiso');
for v = 1:numel(M)
    for q = 1:2
        pr = M(v).probe(q);
        fprintf('%-18s %-12s | %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %7.3f\n', ...
            M(v).name, pr.what, pr.PD, pr.SRX, pr.PT, pr.p1, pr.p2, pr.A, pr.F);
    end
end

save(fullfile(resdir,'model_vs_data_recovery.mat'),'M','Wd_s','Wd_k','FisoD');

%% ---------------- figure ----------------
fig = figure('Position',[20 20 1560 860],'Color','w');
tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
cols = lines(numel(M));

nexttile(tl); hold on; box on; grid on;
plot(1000*(d_k.t-d_k.t(1)), d_k.F, 'k-','LineWidth',1.4,'DisplayName','DATA ktr');
for v=1:numel(M)
    if v==LONG_DWELL; continue; end
    plot(1000*(M(v).wk.t-M(v).wk.t(1)), M(v).wk.F, '-','Color',cols(v,:),'LineWidth',1.6,'DisplayName',M(v).name);
end
yline(1,'k:','HandleVisibility','off'); xlabel('ms'); ylabel('F / F_{iso}'); ylim([-0.1 1.3]);
title('ktr protocol'); legend('Location','southeast','FontSize',8);

nexttile(tl); hold on; box on; grid on;
plot(1000*(d_s2.t-d_s2.t(1)), d_s2.F, 'k-','LineWidth',1.4,'DisplayName','DATA post-slack #2');
for v=1:numel(M)
    s2 = gsel(M(v).ws,'postSlack',2);
    plot(1000*(s2.t-s2.t(1)), s2.F, '-','Color',cols(v,:),'LineWidth',1.6,'DisplayName',M(v).name);
end
xlabel('ms'); ylabel('F / F_{iso}'); xlim([0 150]); ylim([-0.1 1.1]);
title('post-slack cycle 2 (same SL = 2.0 \mum)'); legend('Location','southeast','FontSize',8);

nexttile(tl); hold on; box on; grid on;
kk = [d_k.k63 arrayfun(@(v) M(v).wk.k63, 1:numel(M))];
ss = [d_s2.k63 arrayfun(@(v) getfield(gsel(M(v).ws,'postSlack',2),'k63'), 1:numel(M))]; %#ok
b = bar([ss(:) kk(:)]); b(1).FaceColor=[0.1 0.6 0.25]; b(2).FaceColor=[0 0 0];
xticks(1:numel(kk)); xticklabels([{'DATA'} {M.name}]); xtickangle(30);
ylabel('k_{63} (s^{-1})'); legend({'post-slack #2','ktr'},'Location','northwest','FontSize',8);
title('The pair the model has to match');

nexttile(tl); hold on; box on; grid on;
for v=1:numel(M)
    pcv = arrayfun(@(c) getfield(gsel(M(v).ws,'postSlack',c),'k63'), 1:5); %#ok
    plot(1:5, pcv, 'o-','Color',cols(v,:),'LineWidth',1.5,'DisplayName',M(v).name);
end
plot(1:5, arrayfun(@(c) getfield(gsel(Wd_s,'postSlack',c),'k63'), 1:5), 'ks--','LineWidth',2,'DisplayName','DATA'); %#ok
xlabel('slack cycle (SL 2.04 \rightarrow 1.88 \mum)'); ylabel('k_{63} (s^{-1})');
title('post-slack rate vs SL'); legend('Location','southwest','FontSize',8);

nexttile(tl); hold on; box on; grid on;
os = M(2).os; ok = M(2).ok; g=[true;diff(os.t(:))>0]; gk=[true;diff(ok.t(:))>0];
t_s2 = svt(iS(2),1); t_k = M(2).vt(end-1,1);
T=os.t(g); Tk=ok.t(gk);
plot(1000*(T-t_s2), os.SR(g)+os.SRD(g), '-','Color',[0.1 0.6 0.25],'LineWidth',1.8,'DisplayName','SRX, post-slack #2');
plot(1000*(Tk-t_k), ok.SR(gk)+ok.SRD(gk), 'k-','LineWidth',1.8,'DisplayName','SRX, ktr');
plot(1000*(T-t_s2), os.PuR(g), '--','Color',[0.1 0.6 0.25],'LineWidth',1.5,'DisplayName','PD, post-slack #2');
plot(1000*(Tk-t_k), ok.PuR(gk), 'k--','LineWidth',1.5,'DisplayName','PD, ktr');
xlim([-30 150]); xlabel('ms rel. start of redevelopment'); ylabel('fraction');
title('V2: the head pool the model recruits from'); legend('Location','east','FontSize',8);

nexttile(tl); hold on; box on; grid on;
lb = {'DATA', M.name};
vals = [d_k.k63/d_s2.k63, arrayfun(@(v) M(v).wk.k63/getfield(gsel(M(v).ws,'postSlack',2),'k63'), 1:numel(M))]; %#ok
bar(vals); yline(1,'k--','LineWidth',1.5);
xticks(1:numel(vals)); xticklabels(lb); xtickangle(30);
ylabel('ktr / post-slack #2'); title('Ratio — data says 1.0');
sgtitle('Model vs data: the SL-matched ktr / post-slack pair');
exportgraphics(fig, fullfile(resdir,'model_vs_data_recovery.png'),'Resolution',125);
disp('DONE');

%% ---------------- helpers ----------------
function p = setf(p,f,v); p.(f) = v; end

function pr = probeStates(os, ok, svt, kvt, iS)
    g=[true;diff(os.t(:))>0]; gk=[true;diff(ok.t(:))>0];
    ss = os.params.ss;
    t_s2 = svt(iS(2)+1,1);           % end of slack ramp, cycle 2
    t_k  = kvt(end-1,1);             % start of the ktr recovery hold
    pr(1) = one(os, g, ss, t_s2, 'postSlack#2');
    pr(2) = one(ok, gk, ok.params.ss, t_k, 'ktr');
end

function s = one(o, g, ss, tq, what)
    T=o.t(g); f=@(x) interp1(T, x(g), tq, 'linear','extrap');
    Fiso = max(o.Force);
    s.what=what; s.PD=f(o.PuR); s.SRX=f(o.SR)+f(o.SRD); s.PT=f(o.PuATP);
    s.p1=f(o.p1_0); s.p2=f(o.p2_0);
    A = o.PU(:,2*ss+8); s.A = interp1(T, A(g), tq, 'linear','extrap');
    s.F = f(o.Force)/Fiso;
end
