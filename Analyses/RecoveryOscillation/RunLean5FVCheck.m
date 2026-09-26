% RunLean5FVCheck.m -- did force-velocity HOLD while it was paying nothing?
%
% THE QUESTION. F-V is the term the floor switched off: at the lean5 seed it
% sits at chi-2 ~0.07 and contributes 0.0000 to the objective, and it
% contributed 0.0 % of the campaign's 48 % drop. That is the intended
% behaviour -- but "contributes nothing" and "is free to drift" are different
% claims, and only one of them is true. `fl:X` is a DEADBAND: the term pays
% max(0, c - X), so it is silent inside its band and starts paying again the
% moment it leaves. This figure is the check that the distinction is real in
% practice and not just in the source comment.
%
% If the model's F-V at the campaign's best point still lies on the data, the
% floor did what it was built to do: it stopped the search BUYING a solved
% term without letting it be SOLD. If it has visibly drifted, the deadband is
% too wide and FLOOR_CHI2 needs lowering for F-V specifically.
%
% WHAT IS PLOTTED, and why two panels rather than one:
%   (a) P/P0 vs v. The conventional picture, but it is dominated by the
%       isometric point (fnorm = 1 at v = 0 by construction), which carries no
%       kinetic information at all.
%   (b) NORMALISED POWER vs v -- fnorm*v, the quantity actually scored. The
%       v = 0 point is excluded by construction (power is identically zero
%       there in every dataset and every model). This is where the ATP effect
%       lives and where the two source curves genuinely disagree, so it is the
%       panel to judge the fit on.
%   (c) per-point residual in units of the DATA's own spread, with the +-1 sd
%       band drawn. The floor is a whole-term deadband, not per-point, so this
%       shows how much of the term's budget is actually near its own noise.
%
% The data's spread here is FV_normpowerVar, which is the disagreement between
% the two source curves -- between LABS, not between preparations. Treat it as
% the honest error bar on a borrowed curve, and note DataPooling sec. 6: the
% in-preparation FV2 ramp lands 0.21 sd from this at 8 mM, so it is a fair
% target at this concentration.
%
%   matlab -batch "run('Analyses/RecoveryOscillation/RunLean5FVCheck.m')"
%
% Output -> results/lean5_fv.png

clearvars -except CAMPAIGN TAGS; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('CAMPAIGN','var') || isempty(CAMPAIGN); CAMPAIGN = 'lean5'; end
SPEC = campaignSpec(CAMPAIGN);
if ~exist('TAGS','var') || isempty(TAGS); TAGS = SPEC.tags; end
TAGS = TAGS(cellfun(@(t) isfile(fullfile('params',[t '_opt.m'])), TAGS));
snaps = [{SPEC.seedTag}, cellfun(@(t) [t '_opt'], TAGS, 'UniformOutput', false)];

if isempty(gcp('nocreate')); parpool('local',5); end

R = struct('name',{},'fv',{},'cost',{});
for k = 1:numel(snaps)
    pk = iLoadFV(snaps{k});
    [fm, fd] = iRunModelFV(pk);
    fc = evalFeatureCost(fd, fm, pk.fn);
    iFV = find(contains(pk.fn, 'FV_normpowerAvg'), 1);
    R(k).name = snaps{k};
    R(k).fv   = fm(1);
    R(k).cost = fc(iFV);

    % How much of the deadband has been eaten. The scored cost is 0 whether
    % the term is 1 % or 99 % of the way to its floor, so "pays nothing" says
    % nothing about whether it has drifted -- this is the number that does.
    % NOTE the floor is subtracted from the UNWEIGHTED cost (third output),
    % before the scalar weight; comparing it to the weighted cost is a
    % category error and reads ~1000x too high.
    tk   = strsplit(pk.fn{iFV}, '|');
    flTk = tk(startsWith(tk,'fl:'));
    R(k).floor = 0;
    if ~isempty(flTk); R(k).floor = str2double(flTk{1}(4:end)); end
    [~,~,R(k).raw] = evalFeatureCost(fd, fm, {strjoin(tk(~startsWith(tk,'fl:')), '|')});
    R(k).band = 100 * R(k).raw / max(R(k).floor, eps);

    if k == 1; D = fd(1); end
    fprintf('  %-14s  pays %.6f | raw %.6f | %.1f %% of its %.4f deadband\n', ...
            snaps{k}, R(k).cost, R(k).raw, R(k).band, R(k).floor);
end

%% ---- draw --------------------------------------------------------------
fig = figure('Position',[40 40 1500 460],'Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
co = lines(4); mk = {'-o','-s','-^'};

% --- (a) P/P0
nexttile; hold on; box on
plot(D.FV_v, D.FV_fnorm, 'k o','MarkerFaceColor','k','MarkerSize',7, ...
     'LineStyle','none','DisplayName','data (2 source curves, avg)');
if isfield(D,'FV_expFit_v') && isfield(D,'FV_expFit_f')
    plot(D.FV_expFit_v, D.FV_expFit_f, 'k--','LineWidth',1.1,'DisplayName','Hill fit to data');
end
for k = 1:numel(R)
    plot(R(k).fv.FV_v, R(k).fv.FV_fnorm, mk{min(k,3)}, 'LineWidth',1.7, ...
         'MarkerSize',5,'Color',co(k,:),'DisplayName',iLabel(R(k).name,k));
end
xlabel('shortening velocity  (ML s^{-1})'); ylabel('P / P_0');
title('(a)  force-velocity'); legend('Box','off','Location','northeast','Interpreter','none');

% --- (b) normalised power -- the quantity actually scored
nexttile; hold on; box on
% Forced to ROWS. features_* fields are sometimes columns, and a row minus a
% column is implicit expansion into a matrix -- which silently turns each
% residual line in panel (c) into n lines and quadruples the legend.
pv = D.FV_powerv(:)'; pa = D.FV_normpowerAvg(:)';
ps = sqrt(max(D.FV_normpowerVar(:)',0));
errorbar(pv, pa, ps, 'k o','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.4, ...
         'LineStyle','none','DisplayName','data \pm 1 sd (between source curves)');
for k = 1:numel(R)
    plot(R(k).fv.FV_powerv, R(k).fv.FV_normpowerAvg, mk{min(k,3)}, 'LineWidth',1.7, ...
         'MarkerSize',5,'Color',co(k,:),'DisplayName',iLabel(R(k).name,k));
end
xlabel('shortening velocity  (ML s^{-1})'); ylabel('normalised power  (P/P_0 \cdot v)');
title('(b)  normalised power  --  THIS is the scored term');
legend('Box','off','Location','southeast','Interpreter','none');

% --- (c) residual in units of the spread THE COST ACTUALLY USES.
% Not sqrt(FV_normpowerVar) raw: at v = 0.5 the two source curves happen to
% agree almost exactly, so the raw sd is ~0 there and the residual would read
% +-40 sd on a point panel (b) shows the model sitting on. evalFeatureCost
% never sees that -- it floors the variance at VAR_FLOOR_FRAC = 0.1 of its own
% mean for exactly this reason. Use the same floor here or the panel disagrees
% with the objective it is supposed to explain.
VAR_FLOOR_FRAC = 0.1;
vRaw = ps.^2;
vEff = max(vRaw, VAR_FLOOR_FRAC * mean(vRaw(vRaw > 0), 'omitnan'));
psEff = sqrt(vEff);
nexttile; hold on; box on
for k = 1:numel(R)
    m = R(k).fv.FV_normpowerAvg(:)'; n = min(numel(m), numel(pa));
    z = (m(1:n) - pa(1:n)) ./ psEff(1:n);
    plot(pv(1:n), z, mk{min(k,3)}, 'LineWidth',1.7,'MarkerSize',5,'Color',co(k,:), ...
         'DisplayName',iLabel(R(k).name,k));
end
yline(0,'k-','HandleVisibility','off');
yline( 1,'r--','HandleVisibility','off'); yline(-1,'r--','DisplayName','\pm1 sd (variance-floored)');
xlabel('shortening velocity  (ML s^{-1})'); ylabel('(model - data) / sd');
title(sprintf('(c)  residual  --  pays 0 throughout, but deadband use %.0f%% -> %.0f%%', ...
      R(1).band, R(2).band));
legend('Box','off','Location','best','Interpreter','none');

outPng = fullfile(resDir,'lean5_fv.png');
exportgraphics(fig, outPng, 'Resolution',150);
fprintf('\nwrote %s\n', outPng);

%% ---- helpers -----------------------------------------------------------
function s = iLabel(name, k)
if k == 1; s = ['seed: ' name]; else; s = ['best: ' name]; end
end

function [fm, fd] = iRunModelFV(p)
% LoadData/RunBakersExp are scripts; give them their own workspace.
params0 = p;                                                    %#ok<NASGU>
LoadData; RunBakersExp;
fm = features_model; fd = features_data;
end

function p = iLoadFV(tag)
p = getParams(); s = loadParams(tag);
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunKtr=false; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
p.UsePooledTarget=true; p.mods={}; p.g=ones(1,0);
end
