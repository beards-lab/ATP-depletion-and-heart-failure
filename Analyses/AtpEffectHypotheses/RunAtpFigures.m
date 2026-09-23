% RunAtpFigures.m -- the 8 mM and the 2 mM fit, data and model, on shared axes.
%
%   matlab -batch "run('Analyses/AtpEffectHypotheses/RunAtpFigures.m')"
%
% ONE parameter file gives both conditions. Every ATP gate is normalised to 1 at
% MgATP_ref = 8 mM, so the 8 mM leg is bit-identical with the gates on or off;
% the two legs differ by [MgATP], the protocol file and nothing else. The gate
% exponents are READ from the search's own output file rather than retyped, so
% these figures cannot drift away from the number that was actually optimised.
%
% NORMALISATION -- this is what makes the two conditions comparable:
%   data   each preparation's trace is rundown-lifted and then divided by ONE
%          divisor, its own 8 mM pre-slack isometric at SL 2.2 (buildPooledTarget
%          steps 1-2). Dividing the 2 mM trace by its own 2 mM isometric instead
%          would erase the ATP force effect, which is the effect being plotted.
%   model  divided by the 8 mM model's force at the same instant -- the first
%          sample of the first release, where SL is still 2.2. Same rule.
% So the vertical gap between the blue and the red family IS the ATP effect, in
% the data and in the model, on the same scale.
%
% Output -> results/atp_slack.png        the five slack cycles, both conditions
%           results/atp_restretch.png    the re-stretch, zoomed
%           results/atp_features.png     every scored feature + the 2/8 ratios
%           results/atp_fv.png           force-velocity and power
%           results/atp_figures.mat
%
% CAVEAT on the F-V panel: the objective scores F-V at 8 mM only. The 2 mM curve
% there is a PREDICTION -- runFVExperiment hardcodes a = 1, i.e. it always runs
% the F-V protocol at ATP_c(1), so the 2 mM F-V is obtained by calling it
% directly with the concentration list rotated. Baker's own 2 mM column is drawn
% for reference and is not in any cost.

clearvars -except CORE GATEFILE; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('CORE','var')     || isempty(CORE);     CORE = 'lean5_a_opt'; end
if ~exist('GATEFILE','var') || isempty(GATEFILE); GATEFILE = fullfile(resDir,'lowatp5_best.txt'); end
PHI  = 0.452;                      % rundown factor, buildPooledTarget's own
nS   = 5;
DAYS = {'03/27','04/03','04/10'};
cH = [0.10 0.32 0.65];   cL = [0.78 0.20 0.12];    % 8 mM, 2 mM
if isempty(gcp('nocreate')); parpool('local',5); end

%% ---- 1. the gate the search actually found -----------------------------
gate = struct();
txt = strsplit(fileread(GATEFILE), newline);
for i = 1:numel(txt)
    tk = regexp(strtrim(txt{i}), '^params0\.(pAtp_\w+)\s*=\s*([-\d\.eE+]+);', 'tokens','once');
    if ~isempty(tk); gate.(tk{1}) = str2double(tk{2}); end
end
assert(~isempty(fieldnames(gate)), 'no pAtp_* assignments parsed from %s', GATEFILE);
gn = fieldnames(gate);
fprintf('gate from %s:\n', GATEFILE);
for i = 1:numel(gn); fprintf('   %-16s %+.4f\n', gn{i}, gate.(gn{i})); end

%% ---- 2. the two runs ---------------------------------------------------
% Built exactly as RunOptimLowAtp built them, so the costs printed here are the
% same function that was minimised -- not a lookalike assembled for a figure.
p0 = getParams(loadParams(CORE), [], true, false);
p0.EvalFeatures=true; p0.OptimizeOn='Feats'; p0.BreakOnODEUnstable=false;
p0.RunKtr=false; p0.RunStairs=false; p0.RunForceVelocityTime=false;
p0.RunSlackPassive=false; p0.RunSlack=true; p0.RunForceVelocity=true;
p0.MaxRunTime=400; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.UsePooledTarget=true; p0.mods={}; p0.g=ones(1,0);
p0.UseAtpPowerGates=true; p0.MgATP_ref=8;
p0.UseAtpDetachR2=0; p0.UseAtpKmsrd=0; p0.UseStrainDepAtpGate=0;
p0.recalculateDataFeats=false;            % see the header of RunOptimLowAtp
p0.RunSlackSegments='AllPar';

prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn8 = p0.fn; [~,iu]=unique(fn8,'stable'); fn8=fn8(iu);
fn8 = fn8(~startsWith(prim(fn8),'PS_'));
fn2 = rescaleFnToCondition(fn8, 2);       % same terms, 2 mM chi-2 scale

p8 = p0; p8.MgATP=8; p8.Pi=0; p8.fn=fn8;
p8.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
for i=1:numel(gn); p8.(gn{i}) = gate.(gn{i}); end     % inert at 8 mM, by construction

p2 = p8; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8; p2.fn=fn2;
p2.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';

rep = checkFeatureCoverage(p2);
if ~rep.ok
    bad = regexprep(unique([rep.missingModel(:);rep.missingData(:)])', '\s*\(in\s.*$', '');
    fprintf(2,'dropping unscorable at 2 mM: %s\n', strjoin(bad,', '));
    fn2 = fn2(~ismember(prim(fn2), bad)); p2.fn = fn2;
end

[fm8, fd8, os8] = iRun(p8);
[fm2, fd2, os2] = iRun(p2);
c8 = sum(evalFeatureCost(fd8, fm8, fn8));
c2 = sum(evalFeatureCost(fd2, fm2, fn2));
fprintf('\nobjective   8 mM %.4f  (%d terms)\n            2 mM %.4f  (%d terms)\n', ...
        c8, numel(fn8), c2, numel(fn2));

% Gate inertness at 8 mM, asserted on this very run rather than assumed.
pOff = p8; for i=1:numel(gn); pOff.(gn{i}) = 0; end
[fmOff, fdOff] = iRun(pOff);
dInert = abs(sum(evalFeatureCost(fdOff, fmOff, fn8)) - c8);
fprintf('8 mM gate inertness |d| = %.3e  (MUST be 0)\n', dInert);
assert(dInert < 1e-9, 'the low-ATP gate leaks into the 8 mM fit');

%% ---- 3. the model anchor, and the pooled data bands --------------------
Sd8  = load(fullfile('data','protocol_03_27_2026_8mM_slack.mat'));
Sd2  = load(fullfile('data','protocol_03_27_2026_2mM_slack.mat'));
rel8 = Sd8.velocitytable(find(Sd8.velocitytable(:,2) < -20),1);
rel2 = Sd2.velocitytable(find(Sd2.velocitytable(:,2) < -20),1);
iA   = find(os8.t >= rel8(1), 1);
aMod = os8.Force(iA);                       % ONE divisor for both model legs
if abs(os8.SL(iA) - 2.2) > 0.01
    warning('model anchor taken at SL %.4f, not 2.2', os8.SL(iA));
end

R = poolLoadRuns();
CONDS = {'8mM','2mM'};
anc = nan(3,2); ord = anc; slp = anc; tac = anc;
for p = 1:3
    for c = 1:2
        i = find(strcmp({R.day},DAYS{p}) & strcmp({R.cond},CONDS{c}),1);
        anc(p,c)=R(i).A.Fss22; ord(p,c)=R(i).order; slp(p,c)=R(i).slope; tac(p,c)=R(i).Tact;
    end
end
lift = ones(3,2);
for p = 1:3
    e = ord(p,:)==1; l = ord(p,:)==2;
    lift(p,l) = 1 + PHI*abs(slp(p,e))*tac(p,e)/anc(p,l);
end
alpha = anc(:,1).*lift(:,1);                % the one divisor, per preparation
fprintf('anchors (kPa, rundown-lifted): %s  | model %.1f\n', ...
        strtrim(sprintf('%.1f ', alpha)), aMod);

TG   = (-0.010 : 2e-4 : 0.560)';
BAND = nan(nS, numel(TG), 3, 2);            % seg x time x [mean -sd +sd] x cond
for c = 1:2
    for g = 1:nS
        V = nan(3, numel(TG));
        for p = 1:3
            S = load(fullfile('data', iSlackFile(DAYS{p}, CONDS{c})));
            t = S.datatable(:,1); F = S.datatable(:,3);
            rl = S.velocitytable(find(S.velocitytable(:,2) < -20),1);
            V(p,:) = interp1(t - rl(g), F*lift(p,c)/alpha(p), TG, 'linear', NaN);
        end
        BAND(g,:,1,c) = mean(V,1,'omitnan');
        BAND(g,:,2,c) = BAND(g,:,1,c) - std(V,0,1,'omitnan');
        BAND(g,:,3,c) = BAND(g,:,1,c) + std(V,0,1,'omitnan');
    end
end

%% ---- 4. the slack cycles ----------------------------------------------
MOD = {os8, os2}; REL = {rel8, rel2}; CC = {cH, cL};
f = figure('Position',[30 30 1500 620],'Color','w');
tl = tiledlayout(f,1,nS,'Padding','compact','TileSpacing','compact');
for g = 1:nS
    nexttile(g); hold on; box on
    for c = 1:2
        fill([1000*TG; flipud(1000*TG)], [BAND(g,:,2,c)'; flipud(BAND(g,:,3,c)')], ...
             CC{c}, 'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
        plot(1000*TG, BAND(g,:,1,c), '-','Color',CC{c},'LineWidth',1.6, ...
             'DisplayName',sprintf('%s data (n=3)', CONDS{c}));
        o = MOD{c}; r = REL{c}(g);
        m = o.t-r >= TG(1) & o.t-r <= TG(end);
        plot(1000*(o.t(m)-r), o.Force(m)/aMod, '--','Color',CC{c},'LineWidth',1.5, ...
             'DisplayName',sprintf('%s model', CONDS{c}));
    end
    xlabel('time from release (ms)');
    if g == 1
        ylabel('F / F_{ref} (8 mM isometric, SL 2.2)');
        legend('Location','south','FontSize',7);
    end
    title(sprintf('slack cycle %d', g),'FontSize',9);
    ylim([-0.05 1.35]);
end
title(tl, sprintf(['8 mM and 2 mM, pooled data (mean \\pm between-preparation SD, n=3) vs %s' ...
                   '   |   objective  8 mM %.3f   2 mM %.3f'], CORE, c8, c2), 'Interpreter','tex');
exportgraphics(f, fullfile(resDir,'atp_slack.png'), 'Resolution',140);

%% ---- 5. the re-stretch -------------------------------------------------
% 270-420 ms: starts 40 ms BEFORE the re-stretch so the first peak and the
% approach to it are both on screen, and the y-range clears the MODEL's
% overshoot rather than only the data's.
f = figure('Position',[30 30 1500 620],'Color','w');
tl = tiledlayout(f,1,nS,'Padding','compact','TileSpacing','compact');
for g = 1:nS
    nexttile(g); hold on; box on
    for c = 1:2
        fill([1000*TG; flipud(1000*TG)], [BAND(g,:,2,c)'; flipud(BAND(g,:,3,c)')], ...
             CC{c}, 'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
        plot(1000*TG, BAND(g,:,1,c), '-','Color',CC{c},'LineWidth',1.6, ...
             'DisplayName',sprintf('%s data', CONDS{c}));
        o = MOD{c}; r = REL{c}(g);
        m = o.t-r >= TG(1) & o.t-r <= TG(end);
        plot(1000*(o.t(m)-r), o.Force(m)/aMod, '--','Color',CC{c},'LineWidth',1.5, ...
             'DisplayName',sprintf('%s model', CONDS{c}));
    end
    xlim([270 420]); ylim([0.5 1.35]); xlabel('time from release (ms)');
    if g == 1; ylabel('F / F_{ref}'); legend('Location','southeast','FontSize',7); end
    title(sprintf('cycle %d: the re-stretch', g),'FontSize',9);
end
title(tl, 'the re-stretch recovery -- the transient that carries the ATP effect (rsK80 x0.38 in the data)');
exportgraphics(f, fullfile(resDir,'atp_restretch.png'), 'Resolution',140);

%% ---- 6. every scored feature, both conditions --------------------------
% chi-2 per element in units of the measured between-preparation spread:
% 1.0 = one preparation-width = exactly the deadband floor.
P8 = pooledTarget8mM(); P2 = pooledTarget2mM();
FEAT = P8.name(:)'; nF = numel(FEAT);
X  = nan(nF,2); MOfD = nan(nF,2); RAT = nan(nF,2);
FMS = {fm8, fm2}; PTS = {P8, P2};
for c = 1:2
    for k = 1:nF
        nm = FEAT{k}; PT = PTS{c}; fmc = FMS{c};
        if ~isfield(fmc, nm); continue; end
        m = double([fmc.(nm)]); m = m(:)'; d = PT.val(k,:);
        n = min(numel(m), numel(d)); if n == 0; continue; end
        mu = mean(d(1:n),'omitnan'); if ~isfinite(mu) || mu == 0; continue; end
        X(k,c)    = mean(((m(1:n)-d(1:n))/mu).^2,'omitnan') / PT.s(k)^2;
        MOfD(k,c) = mean(m(1:n),'omitnan') / mu;
    end
end
for k = 1:nF
    nm = FEAT{k};
    if ~isfield(fm8,nm) || ~isfield(fm2,nm); continue; end
    a = mean(double([fm2.(nm)]),'omitnan'); b = mean(double([fm8.(nm)]),'omitnan');
    if isfinite(b) && b ~= 0; RAT(k,1) = a/b; end
    RAT(k,2) = P2.r(k);                       % the measured 2/8 ratio
end

f = figure('Position',[30 30 1500 820],'Color','w');
tiledlayout(f,3,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on
b = bar(1:nF, X, 'grouped'); b(1).FaceColor=cH; b(2).FaceColor=cL;
% Explicit limits. On auto-scale the log axis clipped the three LARGEST terms
% flat against the top edge and dropped the smallest below the bottom -- i.e.
% it hid exactly the terms the panel exists to show. Range covers 5e-3 to 11.
set(gca,'YScale','log'); ylim([3e-3 30]);
yline(1,'k--','LineWidth',1.3);
text(nF+0.4, 1.35, 'one preparation-width = the deadband floor', ...
     'HorizontalAlignment','right','FontSize',8,'FontWeight','bold');
set(gca,'XTick',1:nF,'XTickLabel',FEAT,'XTickLabelRotation',35,'FontSize',8);
ylabel('\chi^2 per element'); xlim([0.4 nF+0.6]);
legend({'8 mM','2 mM'},'Location','northwest','FontSize',8);
title(sprintf(['every pooled feature, in units of its OWN between-preparation spread   ' ...
               '|   same 19 terms both sides:  8 mM %.3f,  2 mM %.3f'], c8, c2));

nexttile; hold on; box on
b = bar(1:nF, MOfD, 'grouped'); b(1).FaceColor=cH; b(2).FaceColor=cL;
yline(1,'k-','LineWidth',1);
for c = 1:2
    s = PTS{c}.s(:);
    errorbar((1:nF)+0.15*(2*c-3), ones(nF,1), s, 'k', 'LineStyle','none', ...
             'CapSize',3, 'HandleVisibility','off');
end
set(gca,'XTick',1:nF,'XTickLabel',FEAT,'XTickLabelRotation',35,'FontSize',8);
ylabel('model / pooled data'); xlim([0.4 nF+0.6]); ylim([0 3]);
title(['model over data, with the between-preparation SD drawn on the unity line ' ...
       '(a bar inside the whisker is inside the noise)']);

nexttile; hold on; box on
b = bar(1:nF, RAT, 'grouped'); b(1).FaceColor=[0.45 0.45 0.50]; b(2).FaceColor=[0.15 0.55 0.35];
errorbar((1:nF)+0.15, P2.r(:), P2.sr(:).*P2.r(:), 'k', 'LineStyle','none', ...
         'CapSize',3, 'HandleVisibility','off');   % on the DATA bar, not between them
yline(1,'k-','LineWidth',1);
set(gca,'XTick',1:nF,'XTickLabel',FEAT,'XTickLabelRotation',35,'FontSize',8);
ylabel('2 mM / 8 mM'); xlim([0.4 nF+0.6]);
legend({'model','data (\pm s_r)'},'Location','northwest','FontSize',8);
title('the ATP effect itself: the 2/8 ratio the model produces against the one measured');
exportgraphics(f, fullfile(resDir,'atp_features.png'), 'Resolution',140);

%% ---- 7. force-velocity -------------------------------------------------
% The 2 mM leg is a prediction: runFVExperiment pins a = 1, so the only way to
% get the protocol at 2 mM is to hand it a concentration list starting at 2.
LoadData;                                      %#ok<*NASGU>  gives ATP_c, Data_ATP
% The data features come out of this same call, not out of the merged
% features_data of the battery -- that one is a struct ARRAY once the slack
% block has been merged into it, and dotting a field off it is a different
% operation from reading a plain struct.
[~,~,fv8,fvD] = runFVExperiment(p8, [8 4 2], Data_ATP);
[~,~,fv2]     = runFVExperiment(p2, [2 4 8], Data_ATP);

f = figure('Position',[30 30 1400 440],'Color','w');
tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on
plot(Data_ATP(:,1), Data_ATP(:,2), 'o','Color',cH,'MarkerFaceColor',cH, ...
     'MarkerSize',6,'DisplayName','Baker 8 mM');
plot(Data_ATP(:,1), Data_ATP(:,4), 'o','Color',cL,'MarkerFaceColor','w', ...
     'MarkerSize',6,'DisplayName','Baker 2 mM (not scored)');
plot(fv8.FV_v, fv8.FV_f, '-s','Color',cH,'LineWidth',1.7,'MarkerSize',5,'DisplayName','model 8 mM');
plot(fv2.FV_v, fv2.FV_f, '--s','Color',cL,'LineWidth',1.7,'MarkerSize',5, ...
     'DisplayName','model 2 mM (prediction)');
xlabel('shortening velocity (ML/s)'); ylabel('force (kPa)');
legend('Location','northeast','FontSize',7); title('force-velocity, absolute');

nexttile; hold on; box on
plot(fvD.FV_v, fvD.FV_fnorm, 'ko','MarkerFaceColor','k','MarkerSize',6, ...
     'DisplayName','data 8 mM (scored)');
plot(Data_ATP(:,1), Data_ATP(:,4)/Data_ATP(1,4), 'o','Color',cL,'MarkerFaceColor','w', ...
     'MarkerSize',6,'DisplayName','Baker 2 mM');
plot(fv8.FV_v, fv8.FV_fnorm, '-s','Color',cH,'LineWidth',1.7,'MarkerSize',5,'DisplayName','model 8 mM');
plot(fv2.FV_v, fv2.FV_fnorm, '--s','Color',cL,'LineWidth',1.7,'MarkerSize',5,'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)'); ylabel('F / F_0');
legend('Location','northeast','FontSize',7); title('force-velocity, each normalised by its OWN isometric');

nexttile; hold on; box on
plot(fvD.FV_powerv, fvD.FV_normpowerAvg, 'ko','MarkerFaceColor','k','MarkerSize',6, ...
     'DisplayName','data 8 mM (scored)');
plot(fv8.FV_powerv, fv8.FV_normpowerAvg, '-s','Color',cH,'LineWidth',1.7,'MarkerSize',5, ...
     'DisplayName','model 8 mM');
plot(fv2.FV_powerv, fv2.FV_normpowerAvg, '--s','Color',cL,'LineWidth',1.7,'MarkerSize',5, ...
     'DisplayName','model 2 mM');
xlabel('shortening velocity (ML/s)'); ylabel('normalised power (F/F_0 \cdot v)');
legend('Location','northeast','FontSize',7); title('power');
exportgraphics(f, fullfile(resDir,'atp_fv.png'), 'Resolution',140);

%% ---- 8. the headline numbers ------------------------------------------
fprintf('\n=== the ATP effect ===\n%-22s %9s %9s %9s\n','2 mM / 8 mM','model','data','s_r');
for nm = {'steady','A','ktr','rsK80','peak1_y','vall_y','vall2_dy','t0'}
    k = find(strcmp(FEAT,nm{1}),1); if isempty(k); continue; end
    fprintf('%-22s %9.3f %9.3f %8.1f%%\n', nm{1}, RAT(k,1), RAT(k,2), 100*P2.sr(k));
end
fprintf('\nsaved 4 figures to %s\n', resDir);
save(fullfile(resDir,'atp_figures.mat'), 'X','MOfD','RAT','FEAT','BAND','TG', ...
     'c8','c2','gate','CORE','aMod','alpha');

%% ---- local functions ---------------------------------------------------
function [fm, fd, out_slack] = iRun(p)
    % RunBakersExp is a SCRIPT -- it runs in the caller's workspace and
    % clobbers whatever it uses, including loop counters. A function workspace
    % of its own is the only reason the surrounding loops survive.
    params0 = p;                                               %#ok<NASGU>
    LoadData; RunBakersExp;                                    %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end

function f = iSlackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
