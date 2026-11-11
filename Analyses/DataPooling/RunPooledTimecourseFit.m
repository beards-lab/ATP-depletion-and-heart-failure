% RunPooledTimecourseFit.m -- the model against the pooled TIMECOURSE, in kPa.
%
% Everything else in this project compares the model to the pool in one of two
% scale-free ways: per-feature against the pooled error bars
% (RunTargetComparison) or trace-shape with each side over its own steady
% (RunDataEnvelope / RunFitReport's envelope overlay). Both deliberately divide
% the absolute level out, and both say so. The consequence is that no figure
% has ever shown the model, the preparation it was fitted to, and the pooled
% average on the SAME axis in kPa -- which is the one view where a level error
% and a shape error cannot cancel each other.
%
% They do cancel: pooled8_opt sits +14 % above the pooled level (z = 0.94) and
% its first re-stretch peak is +33 % in kPa but only +15 % once normalised.
% Both numbers are true; this figure is where you can see which is which.
%
% NORMALISATION -- matched to the TARGET, not to the envelope. The pooled
% target divides each preparation's force by its own RAW `Fss22` (the lift
% cancels within the 8 mM condition, so the rundown correction does not enter
% here) and rescales by the ONE pooled level Fref. This script does exactly
% that to the traces, so the curve drawn IS the timecourse whose features the
% model was fitted to. Using `steady` instead -- what RunDataEnvelope uses for
% its scale-free band -- would draw a slightly different pooled trace and the
% features would no longer sit on it.
%
% REUSE, so there is one alignment in the project and not two:
%   ENV   Analyses/OptimResultReport/results/data_envelope.mat -- the 3 preps
%         aligned on each release (pure index shift: same 4 kHz, same 0.600 s
%         spacing, same ladder). Stored over each prep's own steady, so
%         F/Fss22 = ENV.tc.F * ENV.scale(p) / Fss22(p).
%   R     poolLoadRuns() -- the Fss22 anchors (cached in pool_raw.mat).
%   PT    pooledTarget8mM() -- Fref. Committed, no data needed.
%   MODEL Analyses/OptimResultReport/results/fit_report_<TAG>.mat if present,
%         otherwise simulated here.
%
% DATA REQUIRED: all THREE 8 mM preparations, not the two files the campaign
% needs -- see checkCampaignData. A pooled average cannot be drawn from one
% preparation.
%
% Output -> results/pooled_timecourse_fit.png, results/pooled_timecourse_fit.txt

clearvars -except TAG TAG_BASE; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('TAG','var');      TAG      = 'pooled8_opt';    end
if ~exist('TAG_BASE','var'); TAG_BASE = 'estiff125_opt';  end

ENVF = fullfile(root,'Analyses','OptimResultReport','results','data_envelope.mat');
if ~isfile(ENVF)
    error('RunPooledTimecourseFit:noEnv', ...
        ['%s not found -- run Analyses/OptimResultReport/RunDataEnvelope.m first.\n' ...
         'It needs all three 8 mM protocol files in data/.'], ENVF);
end
S = load(ENVF); ENV = S.ENV;
PT = pooledTarget8mM();
R  = poolLoadRuns();

%% ---- per-preparation anchor, in the target's own terms ------------------
% ENV.prep is '03/27 M' style; PT.days and R(i).day are '03/27'.
nP  = numel(ENV.prep);
a   = nan(1,nP);
for i = 1:nP
    d  = strtrim(erase(erase(ENV.prep{i},' M'),' F'));
    k  = find(strcmp({R.day}, d) & strcmp({R.cond},'8mM'), 1);
    if isempty(k)
        error('RunPooledTimecourseFit:noAnchor','no 8 mM run for %s in poolLoadRuns', d);
    end
    a(i) = R(k).A.Fss22;                       % RAW Fss22 -- see header
end
Fref = PT.Fref;
fprintf('anchors (raw Fss22, kPa): %s   Fref %.2f\n', strtrim(sprintf('%.1f ', a)), Fref);

% F/steady -> F/Fss22 -> xFref.  ENV.scale(p) is that prep's mean steady.
cvt = reshape(ENV.scale(:) ./ a(:), [], 1) * Fref;      % one factor per prep

%% ---- pooled timecourse in kPa ------------------------------------------
% TWO different absolute traces, and the distinction is the point of the figure:
%   FA   each prep rescaled to the POOLED level (F * Fref/Fss22_p) -- this is
%        what the pooled mean and its band are made of, and it is the target.
%   FR   each prep in the kPa it was actually RECORDED at (F). 03/27's raw
%        trace is what the model was fitted to before pooling, and it sits at
%        a 79.8 kPa anchor against the pool's 69.3 -- drawing it rescaled
%        would hide exactly the level question this figure exists to show.
tA = ENV.tcAll.t;
FA = ENV.tcAll.F .* cvt;                                % nP x nt, at pooled level
FR = ENV.tcAll.F .* reshape(ENV.scale(:),[],1);         % nP x nt, as recorded
pooledMean = mean(FA,1,'omitnan');
pooledSD   = std (FA,0,1,'omitnan');
IREF = find(contains(ENV.prep,'03/27'),1);              % the fitted preparation

nC = ENV.tc.nC; tc = ENV.tc.t;
Fc  = ENV.tc.F .* reshape(cvt,[],1,1);                  % pooled level
FcR = ENV.tc.F .* reshape(ENV.scale(:),[],1,1);         % as recorded
cMean = squeeze(mean(Fc,1,'omitnan')); cSD = squeeze(std(Fc,0,1,'omitnan'));

%% ---- the model ---------------------------------------------------------
MODF = fullfile(root,'Analyses','OptimResultReport','results', ...
                sprintf('fit_report_%s.mat', TAG));
if isfile(MODF)
    M = load(MODF); Rm = M.R; NAMES = M.NAMES;
    fprintf('model traces from %s (no re-simulation)\n', MODF);
else
    fprintf('no saved run for %s -- simulating\n', TAG);
    NAMES = {TAG_BASE, TAG};
    for k = 1:2; Rm(k) = iRun(NAMES{k}); end                       %#ok<SAGROW>
end

% align the model on ITS first release, the same way the data was aligned
mdl = struct();
for k = 1:numel(Rm)
    o  = Rm(k).out_slack;
    vt = o.datatable;                                    %#ok<NASGU>
    rel = Rm(k).p.velocitytableonfile;                   %#ok<NASGU>
    Sd  = load(fullfile('data', Rm(k).p.velocitytableonfile));
    rl  = Sd.velocitytable(Sd.velocitytable(:,2) < -1, 1);
    % 'AllPar' stitches the five segments, so the time vector repeats at the
    % joins. makeMonotonous is what extractSlackAttributes uses before its own
    % interp1 -- same fix, same helper.
    nr = makeMonotonous(o.t);
    mdl(k).t = o.t(nr) - rl(1);  mdl(k).F = o.Force(nr);  mdl(k).rl = rl - rl(1);
end

%% ---- numbers -----------------------------------------------------------
L = {};
L{end+1} = sprintf('MODEL vs POOLED TIMECOURSE, absolute kPa   (TAG %s)', TAG);
L{end+1} = sprintf('pooled level Fref = %.2f kPa from %d preps %s', Fref, PT.nP, strjoin(PT.days,'/'));
L{end+1} = '';
L{end+1} = sprintf('%-18s %10s %10s %10s %10s', 'reference','RMS kPa','RMS %','bias kPa','bias %');
L{end+1} = repmat('-',1,64);
ref = {sprintf("pooled mean @%.1f", Fref), pooledMean; sprintf("%s raw", ENV.prep{IREF}), FR(IREF,:)};
for k = numel(Rm):-1:1
    for r = 1:size(ref,1)
        y = reshape(interp1(mdl(k).t, mdl(k).F, tA, 'linear', NaN), 1, []);
        rr = reshape(ref{r,2}, 1, []);
        d = y - rr; g = isfinite(d);
        mu = mean(rr(g),'omitnan');
        L{end+1} = sprintf('%-18s %10.2f %10.1f %10.2f %10.1f   <- %s', ...
            ref{r,1}, sqrt(mean(d(g).^2)), 100*sqrt(mean(d(g).^2))/mu, ...
            mean(d(g)), 100*mean(d(g))/mu, NAMES{k});                  %#ok<SAGROW>
    end
end
L{end+1} = '';
L{end+1} = sprintf('band width: mean +-SD is %.1f %% of the pooled mean', ...
    100*mean(pooledSD(isfinite(pooledSD)))/mean(pooledMean(isfinite(pooledMean))));
txt = strjoin(L,newline); disp(txt);
fid = fopen(fullfile(resDir,'pooled_timecourse_fit.txt'),'w'); fprintf(fid,'%s\n',txt); fclose(fid);

%% ---- figure ------------------------------------------------------------
cBand=[0.55 0.60 0.66]; cDat=[0.10 0.12 0.16]; cRef=[0.10 0.45 0.75];
CC={[0.55 0.60 0.66],[0.72 0.14 0.18]}; LW={1.4,1.8};

fig = figure('Units','pixels','Position',[20 20 1850 950],'Color','w');
tiledlayout(3,nC,'Padding','compact','TileSpacing','compact');

nexttile([1 nC]); hold on; box on; grid on
g = isfinite(pooledMean);
fill([tA(g); flipud(tA(g))], [pooledMean(g)+pooledSD(g), fliplr(pooledMean(g)-pooledSD(g))]', ...
     cBand,'FaceAlpha',0.25,'EdgeColor','none','DisplayName','3 preps, \pm SD');
plot(tA, pooledMean,'-','Color',cDat,'LineWidth',2.6,'DisplayName','pooled mean (kPa)');
plot(tA, FR(IREF,:),'--','Color',cRef,'LineWidth',1.2, ...
     'DisplayName',[ENV.prep{IREF} ' as recorded (the prep fitted)']);
for k=1:numel(Rm); plot(mdl(k).t, mdl(k).F,'-','Color',CC{k},'LineWidth',LW{k},'DisplayName',NAMES{k}); end
xlim([tA(1) tA(end)]); ylabel('force (kPa)'); xlabel('time from first release (s)');
legend('Location','southeast','FontSize',9,'Interpreter','none');
title({'Model against the POOLED timecourse in absolute kPa — nothing normalised away', ...
       sprintf('pooled level %.1f kPa; %s alone is %.1f kPa', Fref, ENV.prep{IREF}, a(IREF))}, ...
      'FontWeight','bold');

for c = 1:nC
    nexttile(nC + c); hold on; box on; grid on
    gg = isfinite(cMean(:,c))';
    fill([tc(gg); flipud(tc(gg))], [cMean(gg,c)+cSD(gg,c); flipud(cMean(gg,c)-cSD(gg,c))], ...
         cBand,'FaceAlpha',0.25,'EdgeColor','none');
    plot(tc, cMean(:,c),'-','Color',cDat,'LineWidth',2.2);
    plot(tc, squeeze(FcR(IREF,:,c)),'--','Color',cRef,'LineWidth',1.1);
    for k=1:numel(Rm)
        tt = mdl(k).t - mdl(k).rl(c);
        plot(tt, mdl(k).F,'-','Color',CC{k},'LineWidth',LW{k});
    end
    xlim([-0.02 0.30]); title(sprintf('cycle %d — release + recovery',c),'FontSize',10);
    if c==1; ylabel('force (kPa)'); end

    nexttile(2*nC + c); hold on; box on; grid on
    fill([tc(gg); flipud(tc(gg))], [cMean(gg,c)+cSD(gg,c); flipud(cMean(gg,c)-cSD(gg,c))], ...
         cBand,'FaceAlpha',0.25,'EdgeColor','none');
    plot(tc, cMean(:,c),'-','Color',cDat,'LineWidth',2.2);
    plot(tc, squeeze(FcR(IREF,:,c)),'--','Color',cRef,'LineWidth',1.1);
    for k=1:numel(Rm)
        tt = mdl(k).t - mdl(k).rl(c);
        plot(tt, mdl(k).F,'-','Color',CC{k},'LineWidth',LW{k});
    end
    xlim([0.295 0.40]); title(sprintf('cycle %d — the re-stretch',c),'FontSize',10);
    xlabel('t - release (s)'); if c==1; ylabel('force (kPa)'); end
end
exportgraphics(fig, fullfile(resDir,'pooled_timecourse_fit.png'),'Resolution',150);
fprintf('\nwrote results/pooled_timecourse_fit.{png,txt}\n');

%% ======================= helpers ========================================
function Rk = iRun(tag)
    p = getParams(loadParams(tag), [], true, false);
    p.mods = {}; p.g = ones(1,0);
    p.RunForceVelocity = false; p.RunSlack = true; p.RunSlackPassive = false;
    p.RunKtr = false; p.RunStairs = false; p.RunForceVelocityTime = false;
    p.EvalFeatures = true; p.OptimizeOn = 'Feats';
    p.RunSlackSegments = 'AllPar'; p.BreakOnODEUnstable = false;
    p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0; p.MaxRunTime = 300;
    [~, out_slack] = runSlackExperiment(p);
    Rk = struct('tag',tag,'p',p,'out_slack',out_slack);
end
