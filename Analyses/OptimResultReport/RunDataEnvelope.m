% RunDataEnvelope.m — build the THREE-PREPARATION data envelope, features and
% timecourse, so the model can be judged against "the data" rather than against
% one preparation.
%
% WHY THIS IS POSSIBLE AT ALL (checked, not assumed)
%   All three 2026 protocol days share:
%     sampling      4000 Hz, dt = 0.00025 s, identical
%     slack cycles  5, spaced 0.600 s, in all three
%     release vel   -36 / -45 / -53 / -62 / -70 ML/s
%     restretch     4, 5, 6, 7, 3 ML/s
%     length        SL 1.879 -> 2.201 um (datatable col2 is SL, not ML)
%   Only the absolute session-time offset differs (74.435 / 74.883 / 74.925 s).
%   So aligning on each slack release and cutting a common relative window is a
%   pure index shift -- no resampling, no protocol reconciliation.
%
% WHAT IS AVERAGED, AND WHY NOT THE RAW FORCE
%   Absolute force is NOT transferable between preparations: mean `steady` is
%   77.3 / 64.6 / 47.9 kPa, a 1.6x span that is cross-section, damage and
%   rundown, not model error (SlackDataAnalysis Sec 4). Averaging raw traces would
%   produce a mean nobody measured and an SD band dominated by that scale.
%   Every trace is therefore divided by ITS OWN mean `steady` first. The band is
%   then the spread of the SHAPE, which is the reproducible quantity (CV 4-7 %).
%   The model is normalised the same way, by its own steady -- so the comparison
%   is scale-free on both sides, and the absolute level is judged separately.
%
% Two products in results/data_envelope.mat:
%   ENV.feat   per-feature, per-cycle mean/SD across the 3 preps (normalised)
%   ENV.tc     per-cycle timecourse mean/SD on a common relative-time grid
%              + ENV.tcAll, the whole 5-cycle protocol aligned on release 1
%
% Figures: results/env_features.png, results/env_timecourse.png

clearvars -except IREF; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir), mkdir(resDir); end

DS = { '03/27 M','protocol_03_27_2026_8mM_slack.mat'
       '04/03 F','protocol_04_03_2026_8mM_slack.mat'
       '04/10 M','protocol_04_10_2026_8mM_slack.mat' };
nP = size(DS,1);

% relative window around each slack release (cycle spacing is 0.600 s)
TPRE = 0.05; TPOST = 0.545; DT = 0.00025;
tg = (-TPRE:DT:TPOST)';

%% ---- load, extract, align ----------------------------------------------
P = struct();
for i = 1:nP
    S = load(['data/' DS{i,2}]);
    P(i).t   = S.datatable(:,1);
    P(i).ML  = S.datatable(:,2);
    P(i).F   = S.datatable(:,3);
    P(i).vt  = S.velocitytable;
    P(i).rel = S.velocitytable(S.velocitytable(:,2) < -1, 1);   % release times
    P(i).ft  = extractSlackAttributes(P(i).t, P(i).F, P(i).ML, P(i).vt, [], [], false, true);
    P(i).sc  = mean(double(P(i).ft.steady), 'omitnan');          % per-prep force scale
    fprintf('%-9s  scale(steady) = %5.1f kPa | %d releases | dt = %.5f s\n', ...
            DS{i,1}, P(i).sc, numel(P(i).rel), median(diff(P(i).t)));
end
nC = min(arrayfun(@(x) numel(x.rel), P));

% The preparation the model is actually fitted to. Marked separately in every
% figure: the envelope says what "the data" is, but the model was only ever
% shown THIS one, so a deviation from the mean is only a defect if it is also a
% deviation from 03/27.
IREF = find(strcmp(DS(:,1), '03/27 M'), 1);
ENVREF = DS{IREF,1};

%% ---- 1. FEATURE envelope ----------------------------------------------
% class F = force-dimension (divide by own steady); S = absolute scale anchor;
% everything else is already scale-free.
FEAT = { 'steady','S'; 'A','F'; 'Am','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'rsA','F'; 'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         't0','T'; 'peak1_t','T'; 'peak1_dSL','T'; 'vall_t','T'; 'vall2_t','T'; 'rsT0','T'
         'rsR2','Q'; 'SLslack','Q' };
ENV = struct(); ENV.prep = DS(:,1)'; ENV.scale = [P.sc];
for j = 1:size(FEAT,1)
    nm = FEAT{j,1}; V = nan(nP, nC);
    ok = true;
    for i = 1:nP
        if ~isfield(P(i).ft, nm); ok = false; break; end
        v = double(P(i).ft.(nm)(:))';
        if numel(v) < nC; ok = false; break; end
        v = v(1:nC);
        if FEAT{j,2} == 'F'; v = v / P(i).sc; end
        V(i,:) = v;
    end
    if ~ok; continue; end
    ENV.feat.(nm) = struct('cls',FEAT{j,2},'raw',V,'mean',mean(V,1,'omitnan'), ...
                           'sd',std(V,0,1,'omitnan'));
end
ENV.featNames = fieldnames(ENV.feat)';

%% ---- 2. TIMECOURSE envelope -------------------------------------------
% per cycle: every prep on the same relative grid, force scaled by own steady
ENV.tc.t = tg; ENV.tc.nC = nC;
Fc = nan(nP, numel(tg), nC); Lc = nan(nP, numel(tg), nC);
for c = 1:nC
    for i = 1:nP
        tt = P(i).t - P(i).rel(c);
        m  = tt >= tg(1)-2*DT & tt <= tg(end)+2*DT;
        Fc(i,:,c) = interp1(tt(m), P(i).F(m) / P(i).sc, tg, 'linear', NaN);
        Lc(i,:,c) = interp1(tt(m), P(i).ML(m),          tg, 'linear', NaN);
    end
end
ENV.tc.F     = Fc;
ENV.tc.Fmean = squeeze(mean(Fc,1,'omitnan'));
ENV.tc.Fsd   = squeeze(std (Fc,0,1,'omitnan'));
ENV.tc.Lmean = squeeze(mean(Lc,1,'omitnan'));

% whole protocol, aligned on release 1
tA = (-TPRE : DT : (nC-1)*0.600 + TPOST)';
FA = nan(nP, numel(tA)); LA = nan(nP, numel(tA));
for i = 1:nP
    tt = P(i).t - P(i).rel(1);
    m  = tt >= tA(1)-2*DT & tt <= tA(end)+2*DT;
    FA(i,:) = interp1(tt(m), P(i).F(m) / P(i).sc, tA, 'linear', NaN);
    LA(i,:) = interp1(tt(m), P(i).ML(m),          tA, 'linear', NaN);
end
ENV.tcAll = struct('t',tA,'F',FA,'Fmean',mean(FA,1,'omitnan'),'Fsd',std(FA,0,1,'omitnan'), ...
                   'Lmean',mean(LA,1,'omitnan'));

fprintf('\nband width (SD as %% of mean) over the cycle-1 window: median %.1f %%, p90 %.1f %%\n', ...
    100*median(ENV.tc.Fsd(:,1)./abs(ENV.tc.Fmean(:,1)),'omitnan'), ...
    100*prctile(ENV.tc.Fsd(:,1)./abs(ENV.tc.Fmean(:,1)),90));

ENV.iRef = IREF; ENV.refName = ENVREF;
save(fullfile(resDir,'data_envelope.mat'), 'ENV', 'DS', 'FEAT', 'IREF');

%% ---- figures -----------------------------------------------------------
cD = [0.20 0.25 0.30]; cB = [0.20 0.25 0.30 0.18];
pc = lines(nP);

% --- timecourse: whole protocol + per-cycle detail
fig = figure('Units','pixels','Position',[20 20 1750 900],'Color','w');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

nexttile([1 3]); hold on; box on; grid on
fill([tA; flipud(tA)], [ENV.tcAll.Fmean(:)+ENV.tcAll.Fsd(:); flipud(ENV.tcAll.Fmean(:)-ENV.tcAll.Fsd(:))], ...
     cD, 'FaceAlpha',0.20, 'EdgeColor','none', 'DisplayName','\pm SD (3 preps)');
for i = setdiff(1:nP, IREF)
    plot(tA, FA(i,:), '-', 'Color',[pc(i,:) 0.55], 'LineWidth',0.6, 'DisplayName',DS{i,1});
end
plot(tA, FA(IREF,:), '-', 'Color',pc(IREF,:), 'LineWidth',1.4, ...
     'DisplayName',[ENVREF '  (the fitted prep)']);
plot(tA, ENV.tcAll.Fmean, '-', 'Color',cD, 'LineWidth',2.2, 'DisplayName','mean');
xlabel('time from first slack release (s)'); ylabel('F / own steady');
legend('Location','southeast','FontSize',9); xlim([tA(1) tA(end)]);
title({'Three-preparation envelope, whole slack protocol', ...
       'each trace divided by its OWN mean steady force (77.3 / 64.6 / 47.9 kPa) — the band is SHAPE spread'}, ...
      'FontWeight','bold');

for c = 1:min(6,nC)
    nexttile; hold on; box on; grid on
    mu = ENV.tc.Fmean(:,c); sd = ENV.tc.Fsd(:,c);
    fill([tg; flipud(tg)], [mu+sd; flipud(mu-sd)], cD, 'FaceAlpha',0.20, 'EdgeColor','none');
    for i = setdiff(1:nP, IREF)
        plot(tg, squeeze(Fc(i,:,c)), '-', 'Color',[pc(i,:) 0.45], 'LineWidth',0.6);
    end
    plot(tg, squeeze(Fc(IREF,:,c)), '-', 'Color',pc(IREF,:), 'LineWidth',1.4);   % the fitted prep
    plot(tg, mu, '-', 'Color',cD, 'LineWidth',1.8);
    xlim([-0.02 0.36]); xlabel('t - release (s)');
    if c == 1; ylabel('F / own steady'); end
    title(sprintf('cycle %d  (restretch %.0f ML/s)', c, P(1).ft.v_restretch(c)));
end
exportgraphics(fig, fullfile(resDir,'env_timecourse.png'), 'Resolution',150);

% --- features with between-prep error bars
fn2 = ENV.featNames;
fig = figure('Units','pixels','Position',[20 20 1750 950],'Color','w');
tl = tiledlayout(4,5,'Padding','compact','TileSpacing','compact');
for j = 1:numel(fn2)
    E = ENV.feat.(fn2{j});
    nexttile; hold on; box on; grid on
    errorbar(1:nC, E.mean, E.sd, 'o-', 'Color',cD, 'MarkerFaceColor',cD, ...
             'LineWidth',1.6, 'CapSize',7, 'DisplayName','3-prep mean \pm SD');
    for i = setdiff(1:nP, IREF)
        plot(1:nC, E.raw(i,:), '.', 'Color',[pc(i,:) 0.75], 'MarkerSize',11, 'DisplayName',DS{i,1});
    end
    plot(1:nC, E.raw(IREF,:), 'o-', 'Color',pc(IREF,:), 'LineWidth',1.3, 'MarkerSize',5, ...
         'DisplayName',[ENVREF ' (fitted)']);
    xlim([0.5 nC+0.5]); set(gca,'XTick',1:nC);
    cv = 100*mean(E.sd./abs(E.mean),'omitnan');
    title(sprintf('%s   (CV %.0f %%)', fn2{j}, cv), 'Interpreter','none','FontSize',10);
    if j > numel(fn2)-5; xlabel('slack cycle'); end
    if j == 1; legend('Location','best','FontSize',7); end
end
title(tl, {'Every scored feature with its BETWEEN-PREPARATION spread (3 preps, 8 mM)', ...
           sprintf(['force-dimension features divided by each prep''s own steady;  ' ...
                    'coloured line = %s, the preparation the model is fitted to'], ENVREF)}, ...
      'FontWeight','bold','FontSize',13);
exportgraphics(fig, fullfile(resDir,'env_features.png'), 'Resolution',150);

fprintf('wrote %s\n', fullfile(resDir,'data_envelope.mat'));
fprintf('wrote env_timecourse.png, env_features.png\n');
