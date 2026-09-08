% RunTargetComparison.m -- the old target, the new target, and where the
% current model sits between them.
%
% THE OLD TARGET is one preparation: params/spreadw_opt.m fits
% protocol_03_27_2026_8mM_slack.mat, in raw kPa, with no error bar anywhere.
% Every residual is therefore unfalsifiable -- a 15 % miss and a 2 % miss look
% the same to the objective, because nothing says how big 15 % is.
%
% THE NEW TARGET is buildPooledTarget(): three preparations, each divided by its
% own pre-slack isometric at SL 2.2, the run recorded second in each session
% lifted for rundown, features extracted per preparation and then averaged,
% spreads moderated (Smyth 2004) because n = 3 gives 2 dof. It carries
%   * a value and a 95 % CI for every feature,
%   * ONE absolute level term (F_ref = 69.3 kPa, CV 15.5 %), because the
%     normalisation removes the scale and an objective without it is
%     scale-blind,
%   * the 2 mM / 8 mM ratios, with the four that are actually resolved marked,
%   * the audited re-stretch rate rsK80 in place of rsK
%     (Analyses/RestretchMetricAudit).
%
% This script puts the three side by side, in the SAME units, and reports for
% every feature whether the model's residual is inside or outside the error bar
% the data can actually support. It also builds the pooled TIMECOURSE with a
% between-preparation band and the model on top -- which is the thing to look
% at when a feature table says "fine" and the eye says otherwise.
%
% Output -> results/target_comparison.mat
%           results/target_table.csv          the numbers
%           results/target_features.png       old vs new vs model, per feature
%           results/target_timecourse.png     the pooled trace + band + model
%           results/target_atp.png            the ATP ratio block
%
%   TAG   which snapshot to score (default 'spreadw_opt')

clearvars -except TAG REUSE; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('TAG','var'); TAG = 'spreadw_opt'; end

DAYS = {'03/27','04/03','04/10'};
REF  = 1;                          % the preparation the old target is built on
nS   = 5;

%% ======================= the two targets ==============================
T = buildPooledTarget();
R = poolLoadRuns();

% the old target: 03/27's own features. Force features are put on the same
% normalised footing (divided by 03/27's own Fss22) so the comparison is about
% SHAPE; the scale difference is reported separately as the level term.
iRef = find(strcmp({R.day},DAYS{REF}) & strcmp({R.cond},'8mM'),1);
aRef = R(iRef).A.Fss22;
old  = nan(numel(T.name), nS);
for k = 1:numel(T.name)
    v = double(R(iRef).feat.(T.name{k})(:))'; m = min(nS,numel(v));
    old(k,1:m) = v(1:m);
    if strcmp(T.cls{k},'F'); old(k,:) = old(k,:)/aRef; end
end

%% ======================= the model ====================================
MATF = fullfile(resDir, sprintf('target_model_%s.mat', TAG));
if exist('REUSE','var') && REUSE && isfile(MATF)
    S = load(MATF); fm = S.fm; out_slack = S.out_slack;
    fprintf('REUSE: replaying %s\n', MATF);
else
    if isempty(gcp('nocreate')); parpool('local',4); end
    params0 = getParams(loadParams(TAG), [], true, false);
    params0.mods = {}; params0.g = ones(1,0);
    params0.RunSlack = true;  params0.RunSlackPassive = false;
    params0.RunForceVelocity = false; params0.RunKtr = false; params0.RunStairs = false;
    params0.RunForceVelocityTime = false;
    params0.EvalFeatures = true; params0.OptimizeOn = 'Feats';
    params0.RunSlackSegments = 'AllPar'; params0.BreakOnODEUnstable = false;
    params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
    params0.MaxRunTime = 600;
    LoadData; RunBakersExp;                                     %#ok<*NASGU>
    fm = features_model;
    save(MATF, 'fm', 'out_slack', 'params0');
end
% The model's own anchor, defined exactly as Fss22 is in the data: the isometric
% plateau at ML 1.10 (SL 2.2) immediately BEFORE the first release. The slack
% simulation is initialised to that steady state and its first sample sits at
% SL 2.2000, so the anchor is that first sample -- taken directly rather than by
% averaging a window, because the solver takes long steps through an isometric
% hold and a 40 ms window can contain no sample at all.
Sd0  = load(fullfile('data','protocol_03_27_2026_8mM_slack.mat'));
relM = Sd0.velocitytable(find(Sd0.velocitytable(:,2) < -20),1);
iA   = find(out_slack.t >= relM(1), 1);
aMod = out_slack.Force(iA);
if abs(out_slack.SL(iA) - 2.2) > 0.01
    warning('model anchor taken at SL %.4f, not 2.2 -- check the slack start', ...
            out_slack.SL(iA));
end
fprintf('anchor (pre-slack isometric at ML 1.10):  data 03/27 %.1f kPa   pooled %.1f kPa   model %.1f kPa\n', ...
        aRef, T.scale.mean, aMod);

mod = nan(numel(T.name), nS);
for k = 1:numel(T.name)
    if ~isfield(fm, T.name{k}); continue; end
    v = double(fm.(T.name{k})(:))'; m = min(nS,numel(v));
    mod(k,1:m) = v(1:m);
    if strcmp(T.cls{k},'F'); mod(k,:) = mod(k,:)/aMod; end
end

%% ======================= the table ====================================
% CI on the pooled mean. The moderated variance carries d0 + d degrees of
% freedom, not n-1, which is the whole point of moderating it.
dof  = T.d0 + (T.nP-1);
tcri = tinv(0.975, dof);
halfCI = tcri * T.s8 / sqrt(T.nP);        % relative half-width of the 95 % CI

% Segments 1-4 are the comparable ladder: a common final length (ML 1.10) and a
% re-stretch amplitude rising 0.080 -> 0.140 ML. Segment 5 re-stretches only to
% ML 1.00, by the SMALLEST amount (0.061 ML), at the slowest ramp, and its
% window is 3x longer -- see Analyses/RestretchMetricAudit/RunSegment5Check.
% Averages here are over segments 1-4; the per-segment values are all kept.
gg = 1:4;
mOld = mean(old(:,gg),2,'omitnan');
mNew = mean(T.m8(:,gg),2,'omitnan');
mMod = mean(mod(:,gg),2,'omitnan');

fprintf('\n=== the target, old vs new, and the model against it ===\n');
fprintf('force features are per-preparation-normalised (/Fss22); segments 1-4\n\n');
fprintf('%-20s %3s %11s %11s %9s %11s %9s %8s\n', ...
        'feature','cls','OLD (03/27)','NEW pooled','95% CI','model','resid','z');
zres = nan(numel(T.name),1); dOldNew = zres;
for k = 1:numel(T.name)
    if ~isfinite(mNew(k)); continue; end
    dOldNew(k) = mOld(k)/mNew(k) - 1;
    if isfinite(mMod(k))
        zres(k) = (mMod(k)-mNew(k)) / (T.s8(k)*abs(mNew(k)));
    end
    fprintf('%-20s %3s %11.4g %11.4g %8.1f%% %11.4g %8.1f%% %8.2f\n', ...
        T.name{k}, T.cls{k}, mOld(k), mNew(k), 100*halfCI(k), mMod(k), ...
        100*(mMod(k)/mNew(k)-1), zres(k));
end
fprintf('\n  resid = model/target - 1;  z = resid / s (the moderated spread).\n');
fprintf('  |z| < 1 : inside what three preparations can distinguish.\n');
fprintf('  |z| > 2 : a real miss -- the data resolves it.\n');

nIn = nnz(abs(zres) < 1); nOut = nnz(abs(zres) > 2);
fprintf('\n  model: %d features inside 1 sd, %d beyond 2 sd, of %d scored\n', ...
        nIn, nOut, nnz(isfinite(zres)));

fprintf('\n=== how far the old target was from the pool it stands for ===\n');
[~, iw] = sort(abs(dOldNew),'descend');
fprintf('%-20s %11s %11s %9s %9s\n','feature','03/27','pooled','differ','in sd');
for j = iw(1:min(8,numel(iw)))'
    if ~isfinite(dOldNew(j)); continue; end
    fprintf('%-20s %11.4g %11.4g %+8.1f%% %8.2f\n', T.name{j}, mOld(j), mNew(j), ...
            100*dOldNew(j), dOldNew(j)/T.s8(j));
end
fprintf('  chi-2 of 03/27 against the pool = %.2f (expected %.2f for n=3)\n', ...
        T.chi2Prep(REF), (T.nP-1)/T.nP);

fprintf('\n=== the level term the normalisation removes ===\n');
fprintf('  F_ref pooled       %.1f kPa,  CV %.1f%%  (95%% CI %.1f - %.1f kPa)\n', ...
        T.scale.mean, 100*T.scale.cv, ...
        T.scale.mean*(1-tinv(0.975,T.nP-1)*T.scale.cv/sqrt(T.nP)), ...
        T.scale.mean*(1+tinv(0.975,T.nP-1)*T.scale.cv/sqrt(T.nP)));
fprintf('  03/27 alone        %.1f kPa   (%+.0f%% of the pooled level)\n', ...
        aRef, 100*(aRef/T.scale.mean-1));
fprintf('  model              %.1f kPa   (%+.0f%% of the pooled level, z = %.2f)\n', ...
        aMod, 100*(aMod/T.scale.mean-1), (aMod/T.scale.mean-1)/T.scale.cv);

%% ======================= the ATP block ================================
fprintf('\n=== the ATP block: score these as RATIOS, not as absolute 2 mM values ===\n');
fprintf('%-20s %10s %9s %7s %s\n','feature','2mM/8mM','s_ratio','z','');
[~, iz] = sort(T.z,'descend');
for j = iz
    mark = '';
    if T.z(j) >= 2; mark = '  <-- resolved'; end
    if T.z(j) < 1;  mark = '      (flat)';   end
    fprintf('%-20s %10.3f %8.1f%% %7.2f%s\n', T.name{j}, T.r(j), 100*T.sr(j), T.z(j), mark);
end
fprintf('\n  the force LEVEL carries the ATP effect for every force feature\n');
fprintf('  (x%.3f +- %.1f%%); their SHAPE ratios are all 1.00 +- 0.05, so nine\n', ...
        T.level, 100*T.sLevel);
fprintf('  separate force ratios would be nine copies of one number.\n');

%% ======================= the pooled timecourse ========================
% Per segment: each preparation's post-release trace divided by its own Fss22,
% interpolated onto a common clock referenced to the RELEASE, then the mean and
% the between-preparation SD. The band is a genuine +-1 SD of the normalised
% traces -- but read it knowing what RestretchMetricAudit found: it is honest
% where the dynamics are slow, and in the first ~20 ms of the re-stretch
% recovery the preparations' onsets differ, so the band there is wider than any
% single preparation's own transient. Features, not bands, decide that window.
TG = (-0.010 : 2e-4 : 0.560)';   % a whole slack cycle: release, hold, re-stretch, recovery
% the actual manoeuvre of each segment, for the panel titles: segment 5 is not
% a continuation of the ladder (it re-stretches to ML 1.00, not 1.10)
sg0  = find(Sd0.velocitytable(:,2) < -1)';
SLSL = arrayfun(@(g) Sd0.velocitytable(sg0(g)+1,4), 1:nS);
SLBK = arrayfun(@(g) Sd0.velocitytable(sg0(g)+3,4), 1:nS);
BAND = nan(nS, numel(TG), 3);                  % mean, -sd, +sd
for g = 1:nS
    V = nan(numel(DAYS), numel(TG));
    for p = 1:numel(DAYS)
        i = find(strcmp({R.day},DAYS{p}) & strcmp({R.cond},'8mM'),1);
        S = load(fullfile('data', slackFile(DAYS{p},'8mM')));
        t = S.datatable(:,1); F = S.datatable(:,3); v = S.velocitytable;
        rel = find(v(:,2) < -20);
        V(p,:) = interp1(t - v(rel(g),1), F/R(i).A.Fss22, TG, 'linear', NaN);
    end
    BAND(g,:,1) = mean(V,1,'omitnan');
    BAND(g,:,2) = BAND(g,:,1) - std(V,0,1,'omitnan');
    BAND(g,:,3) = BAND(g,:,1) + std(V,0,1,'omitnan');
end

f = figure('Position',[40 40 1400 720],'Color','w');
tl = tiledlayout(f,2,5,'Padding','compact','TileSpacing','compact');
cB = [0.15 0.18 0.22]; cM = [0.72 0.14 0.18];
for g = 1:nS
    nexttile(g); hold on; box on
    fill([1000*TG; flipud(1000*TG)], [BAND(g,:,2)'; flipud(BAND(g,:,3)')], cB, ...
         'FaceAlpha',0.18, 'EdgeColor','none', 'HandleVisibility','off');
    plot(1000*TG, BAND(g,:,1), '-','Color',cB,'LineWidth',1.5,'DisplayName','pooled mean');
    m = out_slack.t-relM(g) >= TG(1) & out_slack.t-relM(g) <= TG(end);
    plot(1000*(out_slack.t(m)-relM(g)), out_slack.Force(m)/aMod, '-', ...
         'Color',cM,'LineWidth',1.5,'DisplayName',TAG);
    xlabel('time from release (ms)');
    if g == 1; ylabel('F / F_{ref}'); legend('Location','south','FontSize',7); end
    title(sprintf('segment %d: slack to ML %.3f, back to %.3f', g, SLSL(g), SLBK(g)), ...
          'FontSize',8);
    ylim([-0.05 1.25]);
end
for g = 1:nS
    nexttile(5+g); hold on; box on
    fill([1000*TG; flipud(1000*TG)], [BAND(g,:,2)'; flipud(BAND(g,:,3)')], cB, ...
         'FaceAlpha',0.18,'EdgeColor','none');
    plot(1000*TG, BAND(g,:,1), '-','Color',cB,'LineWidth',1.5);
    m = out_slack.t-relM(g) >= TG(1) & out_slack.t-relM(g) <= TG(end);
    plot(1000*(out_slack.t(m)-relM(g)), out_slack.Force(m)/aMod, '-','Color',cM,'LineWidth',1.5);
    xlim([310 420]); ylim([0.4 1.15]); xlabel('time from release (ms)');
    if g == 1; ylabel('F / F_{ref}'); end
    title(sprintf('seg %d: the re-stretch recovery', g),'FontSize',9);
end
title(tl, sprintf('pooled 8 mM timecourse, mean \\pm between-preparation SD (n=3), vs %s', TAG));
exportgraphics(f, fullfile(resDir,'target_timecourse.png'), 'Resolution',140);

%% ---- per-feature figure ----------------------------------------------
sel = find(isfinite(mNew) & isfinite(mMod));
f2 = figure('Position',[40 40 1300 560],'Color','w');
tiledlayout(f2,1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on; box on
x = 1:numel(sel);
errorbar(x, ones(size(x)), halfCI(sel), 'o','Color',cB,'LineWidth',1.4, ...
         'MarkerFaceColor',cB,'DisplayName','new target (95 % CI)');
plot(x, mOld(sel)./mNew(sel), 's','Color',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74], ...
     'DisplayName','old target = 03/27 alone');
plot(x, mMod(sel)./mNew(sel), 'd','Color',cM,'MarkerFaceColor',cM,'DisplayName',TAG);
yline(1,'k:','HandleVisibility','off');
set(gca,'XTick',x,'XTickLabel',T.name(sel),'XTickLabelRotation',50,'FontSize',8);
ylabel('value / pooled target'); legend('Location','best','FontSize',8);
title('(a) everything relative to the new target');

nexttile; hold on; box on
zz = zres(sel); [~, o] = sort(abs(zz),'descend');
barh(abs(zz(o)));
set(gca,'YTick',1:numel(o),'YTickLabel',T.name(sel(o)),'FontSize',8,'YDir','reverse');
xline(1,'k--','1 sd'); xline(2,'r--','2 sd');
xlabel('|model - target| / s'); title('(b) which residuals the data resolves');
exportgraphics(f2, fullfile(resDir,'target_features.png'), 'Resolution',140);

%% ---- ATP figure -------------------------------------------------------
f3 = figure('Position',[60 60 900 480],'Color','w'); hold on; box on
[~, io] = sort(T.z,'descend'); io = io(:)';
xr = 1:numel(io);
errorbar(xr, T.r(io), T.r(io).*T.sr(io)*tinv(0.975,T.nP-1)/sqrt(T.nP), 'o', ...
         'Color',cB,'MarkerFaceColor',cB,'LineWidth',1.3);
yline(1,'k:');
for j = 1:numel(io)
    if T.z(io(j)) >= 2
        plot(xr(j), T.r(io(j)), 'o','MarkerSize',12,'Color',cM,'LineWidth',1.6);
    end
end
set(gca,'XTick',xr,'XTickLabel',T.name(io),'XTickLabelRotation',55,'FontSize',8);
ylabel('2 mM / 8 mM'); title('the ATP block: circled = z \geq 2, the four the data resolves');
exportgraphics(f3, fullfile(resDir,'target_atp.png'), 'Resolution',140);

%% ---- csv + mat --------------------------------------------------------
Tab = table(T.name(:), T.cls(:), mOld, mNew, halfCI(:), mMod, zres, T.s8(:), ...
            T.s8raw(:), T.r(:), T.sr(:), T.z(:), ...
    'VariableNames',{'feature','cls','old_0327','new_pooled','ci95_rel','model', ...
                     'z_resid','s_used','s_raw','atp_ratio','atp_sd','atp_z'});
writetable(Tab, fullfile(resDir,'target_table.csv'));
save(fullfile(resDir,'target_comparison.mat'),'T','old','mod','mOld','mNew','mMod', ...
     'zres','halfCI','BAND','TG','aRef','aMod','TAG','Tab');
fprintf('\nwrote results/target_table.csv, target_features.png, target_timecourse.png, target_atp.png\n');

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
