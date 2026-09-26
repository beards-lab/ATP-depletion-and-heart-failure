% RunLean5Overview.m -- the four pictures that answer "did this campaign work,
% and did it work for the right reason?"
%
% The total is the least interesting number here. Two campaigns in a row
% lowered it while the features they existed to fix went backwards
% (conclusions.md sec. 7), so an overview that shows only a descent curve is
% the overview that hid the last failure. Hence panels 2-4.
%
%   1  THE DESCENT, per profile, against the seed.
%   2  rsK80 PER SLACK CYCLE -- data, seed, best. This is the campaign's own
%      question: cycle 5 was put back because it is the only re-stretch to
%      ~2.0 um SL, so the panel that matters is whether the model now tracks
%      the rate across all five cycles or only the four it used to be scored on.
%   3  PER-FEATURE chi-2, seed vs best, with the one-preparation-width floor
%      drawn. Above the line a term pays; below it the deadband holds it at
%      zero. This is where "what did it buy" is actually legible.
%   4  OBJECTIVE COMPOSITION, seed vs best. The F-V bar is the control: two
%      campaigns ago it was 40 % of the objective and absorbed 110 % of the
%      gain. If it is still ~0 here, the floor is doing its job.
%
%   matlab -batch "run('Analyses/RecoveryOscillation/RunLean5Overview.m')"
%   matlab -batch "CAMPAIGN='lean5'; TAGS={'lean5_a','lean5_c'}; run('...')"
%
% Output -> results/lean5_overview.png
%           results/lean5_overview.mat

clearvars -except CAMPAIGN TAGS; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('CAMPAIGN','var') || isempty(CAMPAIGN); CAMPAIGN = 'lean5'; end
SPEC = campaignSpec(CAMPAIGN);
if ~exist('TAGS','var') || isempty(TAGS); TAGS = SPEC.tags; end
TAGS = TAGS(cellfun(@(t) isfile(fullfile('params',[t '_opt.m'])), TAGS));
if isempty(TAGS); error('RunLean5Overview:noSnaps','no <tag>_opt.m for campaign %s', CAMPAIGN); end

PT  = pooledTarget8mM();
kR  = find(strcmp(PT.name,'rsK80'),1);
fprintf('campaign %s | seed %.4f | snapshots: %s\n', SPEC.name, SPEC.seedExpect, strjoin(TAGS,', '));

if isempty(gcp('nocreate')); parpool('local',5); end

%% ---- 1. evaluate the seed and every snapshot ---------------------------
% Re-evaluated from disk, never read off a snapshot header: a header is what
% the search BELIEVED it had, and this project has ranked a campaign on a fake
% one before (Analyses/StairsPrediction).
% ONE model run per snapshot. evaluateBakersExp would give the cost but not
% features_model, and calling both doubles a ~40 s evaluation for nothing:
% RunBakersExp leaves features_model AND features_data in the workspace, and in
% OptimizeOn='Feats' the objective IS sum(evalFeatureCost(...)) -- so the cost
% is recoverable from the same run. Verified against evaluateBakersExp below.
% features_model/features_data are STRUCT ARRAYS -- one element per slack
% segment -- so they are assigned field-by-field, never through struct(),
% which would expand the whole record into a 1x5 array. Per-cycle values are
% read back as [EV(i).fm.<name>] for the same reason.
% RunBakersExp is a SCRIPT: it executes in the caller's workspace and clobbers
% whatever it happens to use -- including a loop counter named `i`. Calling it
% through iRunModel gives it a function workspace of its own, so the loop
% survives. (This is why the surrounding code cannot simply inline it.)
EV = struct('name',{},'cost',{},'featCost',{},'featNames',{},'fm',{});
% TAGS are campaign tags (lean5_a); the snapshots they wrote are <tag>_opt.
snaps = [{SPEC.seedTag}, cellfun(@(t) [t '_opt'], TAGS, 'UniformOutput', false)];
for kSnap = 1:numel(snaps)
    pk = iLoad(snaps{kSnap});
    [fm, fd] = iRunModel(pk);
    fc = evalFeatureCost(fd, fm, pk.fn);
    EV(kSnap).name      = snaps{kSnap};
    EV(kSnap).cost      = sum(fc);
    EV(kSnap).featCost  = fc(:)';
    EV(kSnap).featNames = iTrimNames(pk.fn);
    EV(kSnap).fm        = fm;
    fprintf('  %-14s %8.4f\n', snaps{kSnap}, EV(kSnap).cost);
end

% Cross-check the shortcut on the seed only (one extra evaluation, not three).
% If this ever disagrees, the objective is no longer sum(evalFeatureCost) and
% every number on the figure is suspect.
cChk = evaluateBakersExp(1, iLoad(SPEC.seedTag), true, false);
if abs(cChk - EV(1).cost) > 1e-6
    warning('RunLean5Overview:costPath', ...
        ['seed cost via evalFeatureCost (%.6f) != evaluateBakersExp (%.6f). ' ...
         'The objective is not sum(evalFeatureCost) any more -- fix this ' ...
         'script before reading the figure.'], EV(1).cost, cChk);
else
    fprintf('  cost path verified against evaluateBakersExp (%.6f)\n', cChk);
end
iSeed = 1;
[~, iBest] = min([EV(2:end).cost]); iBest = iBest + 1;
fprintf('  best: %s at %.4f (%.1f %% below seed)\n', EV(iBest).name, EV(iBest).cost, ...
        100*(1 - EV(iBest).cost/EV(iSeed).cost));

%% ---- 2. the descents ---------------------------------------------------
H = struct('tag',{},'round',{},'cost',{});
for t = 1:numel(TAGS)
    f = fullfile('params',[TAGS{t} '_state.mat']);
    if ~isfile(f); continue; end
    L = load(f,'state');
    if ~isfield(L,'state') || ~isfield(L.state,'optIter') || isempty(L.state.optIter); continue; end
    it = L.state.optIter;
    imp = it([it.isImproved]);
    H(end+1) = struct('tag',TAGS{t}, 'round',[0 [imp.round]], ...
                      'cost',[SPEC.seedExpect [imp.cost]]);      %#ok<SAGROW>
end

%% ---- 3. per-element chi-2, computed from the features directly ---------
% chi2 = mean over scored cycles of ((model-data)/mean(data))^2 / s^2, i.e.
% the residual in units of the measured between-preparation spread. 1.0 = one
% preparation-width = exactly the floor.
FEAT = PT.name(:)';
X = nan(numel(EV), numel(FEAT));
for i = 1:numel(EV)
    for k = 1:numel(FEAT)
        nm = FEAT{k};
        if ~isfield(EV(i).fm, nm); continue; end
        m = double([EV(i).fm.(nm)]); m = m(:)'; d = PT.val(k,:);
        n = min(numel(m), numel(d)); if n == 0; continue; end
        mu = mean(d(1:n),'omitnan'); if ~isfinite(mu) || mu == 0; continue; end
        X(i,k) = mean(((m(1:n)-d(1:n))/mu).^2,'omitnan') / PT.s(k)^2;
    end
end

%% ---- 4. draw -----------------------------------------------------------
fig = figure('Position',[40 40 1500 900],'Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% --- panel 1: descents
nexttile; hold on; box on
co = lines(max(1,numel(H)));
for i = 1:numel(H)
    stairs(H(i).round, H(i).cost, '-o','LineWidth',1.8,'MarkerSize',4,'Color',co(i,:));
end
yline(SPEC.seedExpect,'k--','seed','LabelHorizontalAlignment','left');
xlabel('round'); ylabel('objective (floored)');
title(sprintf('1  the descent  --  %.0f %% off the seed', ...
      100*(1-EV(iBest).cost/EV(iSeed).cost)));
if ~isempty(H); legend({H.tag},'Box','off','Interpreter','none','Location','northeast'); end

% --- panel 2: rsK80 per cycle -- the campaign's own question
nexttile; hold on; box on
d = PT.val(kR,:); sd = PT.s(kR)*abs(d);
errorbar(1:5, d, sd, 'k-o','LineWidth',2,'MarkerFaceColor','k','MarkerSize',7, ...
         'DisplayName','pooled data \pm 1 prep-width');
mk = {'-s','-^','-v'};
for i = [iSeed iBest]
    if ~isfield(EV(i).fm,'rsK80'); continue; end
    m = double([EV(i).fm.rsK80]); m = m(:)'; m = m(1:min(5,numel(m)));
    lbl = 'seed'; if i ~= iSeed; lbl = ['best: ' EV(i).name]; end
    plot(1:numel(m), m, mk{1+(i~=iSeed)}, 'LineWidth',1.8,'MarkerSize',7, ...
         'DisplayName',lbl);
end
xlim([0.5 5.5]); set(gca,'XTick',1:5);
xline(4.5,':','Color',[.5 .5 .5],'HandleVisibility','off');
text(4.52, double(min(ylim)+0.06*range(ylim)), ...
     sprintf('cycle 5\nre-stretch to ~2.0 \\mum SL\n(was UNSCORED)'), ...
     'FontSize',8,'Color',[.3 .3 .3],'VerticalAlignment','bottom');
xlabel('slack cycle'); ylabel('rsK80  (s^{-1})');
title('2  post-re-stretch rate, per cycle');
legend('Box','off','Location','northwest','Interpreter','none');

% --- panel 3: per-element chi2, seed vs best
nexttile; hold on; box on
scored = find(any(isfinite(X([iSeed iBest],:)),1));
[~, o] = sort(X(iSeed,scored),'descend','MissingPlacement','last');
scored = scored(o);
xs = 1:numel(scored);
bar(xs-0.2, X(iSeed,scored), 0.38, 'FaceColor',[.75 .75 .78],'DisplayName','seed');
bar(xs+0.2, X(iBest,scored), 0.38, 'FaceColor',co(1,:),'DisplayName','best');
yline(1,'r--','LineWidth',1.5,'DisplayName','floor = 1 prep-width');
set(gca,'XTick',xs,'XTickLabel',FEAT(scored),'XTickLabelRotation',55, ...
        'TickLabelInterpreter','none','YScale','log');
ylabel('\chi^2 per scored cycle'); title('3  residual in preparation-widths');
legend('Box','off','Location','northeast');

% --- panel 4: objective composition
nexttile; hold on; box on
nm  = EV(iSeed).featNames; cs = EV(iSeed).featCost(:)'; cb = EV(iBest).featCost(:)';
n   = min([numel(nm) numel(cs) numel(cb)]);
nm  = nm(1:n); cs = cs(1:n); cb = cb(1:n);
keep = (cs > 1e-4) | (cb > 1e-4);
nm = nm(keep); cs = cs(keep); cb = cb(keep);
[~, o] = sort(cs,'descend'); nm = nm(o); cs = cs(o); cb = cb(o);
xs = 1:numel(nm);
bar(xs-0.2, cs, 0.38, 'FaceColor',[.75 .75 .78],'DisplayName','seed');
bar(xs+0.2, cb, 0.38, 'FaceColor',co(1,:),'DisplayName','best');
iFV = find(contains(nm,'FV_'),1);
if ~isempty(iFV)
    text(iFV, max([cs cb])*0.5, sprintf('F-V\n%.3f \\rightarrow %.3f', cs(iFV), cb(iFV)), ...
         'FontSize',8,'Color',[.6 0 0],'HorizontalAlignment','center');
end
set(gca,'XTick',xs,'XTickLabel',nm,'XTickLabelRotation',55,'TickLabelInterpreter','none');
ylabel('weighted cost');
title('4  where the objective sits   (no bar = term is at its floor)');
legend('Box','off','Location','northeast');

outPng = fullfile(resDir,'lean5_overview.png');
exportgraphics(fig, outPng, 'Resolution',150);
save(fullfile(resDir,'lean5_overview.mat'),'EV','H','X','FEAT','SPEC','PT');
fprintf('\nwrote %s\n', outPng);

%% ---- helpers -----------------------------------------------------------
function [fm, fd] = iRunModel(p)
% Run one snapshot and hand back the feature structs. Exists purely to give
% LoadData/RunBakersExp -- both scripts -- a workspace that is not the
% caller's, so they cannot overwrite loop counters or accumulators.
params0 = p;                                                    %#ok<NASGU>
LoadData; RunBakersExp;
fm = features_model; fd = features_data;
end

function nm = iTrimNames(fn)
% Mirrors evaluateBakersExp's private iTrimFn: keep the primary field, drop
% selectors, weights, ranges and any grouped tail.
nm = cellfun(@(s) regexprep(s, '[\[|,].*$', ''), fn, 'UniformOutput', false);
end

function p = iLoad(tag)
p = getParams(); s = loadParams(tag);
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunKtr=false; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=300; p.OptimizeOn='Feats';
p.UsePooledTarget=true; p.mods={}; p.g=ones(1,0);
end
