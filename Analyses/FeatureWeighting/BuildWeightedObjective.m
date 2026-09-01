% BuildWeightedObjective.m — turn the measured data spread into objective weights.
%
% Reads  results/feature_spread.mat  (from RunFeatureSpread.m) and the current
% best parametrisation, and produces a REWEIGHTED params0.fn plus a seed
% snapshot that scores EXACTLY the same total as the current objective, so a
% new optimisation starts from a comparable number.
%
% THE RULE
%
%   w_f  =  lambda * R_f / s_f^2
%
%   s_f   between-preparation relative SD of the feature at 8 mM (3 preps,
%         Baker excluded). evalFeatureCost forms  sum_i |fd_i-fs_i|^2/mean(fd)^2,
%         i.e. a SUM of squared RELATIVE residuals, so dividing by s_f^2 makes
%         every term a chi-square: "1.0 = the model is one preparation-width
%         away from the data". Terms become directly comparable, which they are
%         not today (weights 0.2 ... 50 were set by hand).
%
%   R_f   ATP relevance = 1 + KAPPA*min(z_f, Z_CAP), z = |ratio-1|/sqrt(s8^2+sr^2).
%         The fit runs at 8 mM only, but its purpose is to be the platform the
%         ATP perturbation is applied to. A baseline error in a feature that
%         MOVES with ATP is later indistinguishable from an ATP mechanism; an
%         error in a feature that is ATP-flat is not. So buy accuracy where the
%         ATP signal will be read. With the defaults R spans 1.0 (ATP-blind)
%         to 3.0 (ktr, rsK).
%
%   lambda ONE global constant, chosen so the sum over the reweighted
%         (data-fit) terms equals what those terms contribute today.
%         Regularisers/guardrails are not data and are left untouched, so the
%         grand total is preserved exactly.
%
% Output -> results/weighted_objective.mat, results/weight_table.txt
%           params/<SEEDTAG>.m  (the reweighted seed snapshot)

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results');

%% ---- knobs -------------------------------------------------------------
BASE     = 'realign4_opt';     % the parametrisation to be cost-neutral at
SEEDTAG  = 'spreadw_seed';     % params/<SEEDTAG>.m
% KAPPA = 0 => R == 1, i.e. the weights are PURE PRECISION and nothing else.
% This is deliberate. A weight that also encodes "relevance" is no longer a
% variance, so `cost/s^2` stops being a chi-square and the whole table stops
% being readable as a goodness-of-fit -- which is the most useful thing this
% analysis produces. ATP relevance (`z`) is still measured, still reported, and
% still drawn (results/feature_weighting.png); it just does not touch w.
% The honest place for the ATP effect is a term that SCORES the 2 mM condition
% (conclusions.md Sec 8), not a tilt smuggled into the 8 mM weights.
% Set KAPPA = 0.5 to restore the tilt; R spans 1.0 .. 1+KAPPA*Z_CAP.
KAPPA    = 0;
Z_CAP    = 4;
S_FLOOR  = 0.05;               % no feature is better determined than 5 % across preps
S_CEIL   = 1.00;               % ... nor worse than 100 %
DEDUPE   = true;               % drop the duplicated rsR2 boundary entry

% --- objective CONDITIONING guards (engineering, not data) ---------------
% FV_SHARE -- *** ARBITRARY. THIS NUMBER IS A JUDGEMENT CALL, NOT A
%   MEASUREMENT. *** Force-velocity is a separate, independently acquired
%   experiment, and the whole muscle work profile is in it, so it is held at a
%   fixed share of the objective by fiat rather than allowed to find its own
%   level. Its measured spread (11.3 %, between two LABS -- isometric 56.4 vs
%   67.4 kPa) would put it at 1.4 %: the model already fits FV to well inside
%   the disagreement between its own two source curves, so 1/s^2 sees nothing
%   left to buy there. That is a statement about the FV data's precision, not
%   about FV's importance, and this pin is where the difference is asserted.
%   Consequence to keep in view: FV_SHARE = 0.70 leaves ~19 % of the total for
%   the fourteen other data terms, so the slack SHAPE is along for the ride.
%   0.40-0.50 keeps FV dominant while the slack still has a vote.
% SHARE_CAP: no single OTHER term may exceed this share. Pure 1/s^2 hands rsK
%   70 % (5.8 preparation-widths off -- a model rejection, not a fitting
%   signal), and an objective one feature can buy is how realign2 spent FV
%   monotonically. Applied AFTER the FV pin, so at FV_SHARE >= ~0.55 it is
%   inert (there is no longer 30 % of the total available to any one term).
FV_SHARE  = 0.70;
SHARE_CAP = 0.30;
% Candidate features NOT currently scored; reported, added only if listed here.
ADD_FEATS    = {};             % e.g. {'rsT0','vall2_t','Am'}
CAND         = {'rsT0','vall2_t','Am','peak1_t','rsA'};

% ANCHOR_GUARD -- optional insurance for the ABSOLUTE FORCE LEVEL, which the
%   FV pin quietly weakens. extractForceVelocityAttributes normalises the
%   model's power by the MODEL's OWN isometric force (FV_fnorm = FV_f/FV_f(1)),
%   so FV carries ZERO information about how strong the fibre is: 70 % of the
%   objective is scale-blind by construction. The weight backing the absolute
%   level (A + peak2 + peak1_y + vall_y + steady) falls 125 -> 9.4, i.e. 13x.
%   It is not a free ride -- a uniform +20 % force drift still costs ~1.9,
%   about half the objective -- but it is 13x cheaper than it was.
%   Turning this on adds ONE boundary entry on steady[1]. A boundary term costs
%   exactly 0 inside its range, so cost-neutrality is untouched; it only bites
%   if the force level leaves the range the three preps actually span
%   (50.6 - 80.6 kPa at 8 mM; 03/27, the fitted prep, is 80.6).
ANCHOR_GUARD = false;
ANCHOR_ENTRY = 'steady[1]|50.6-88.7|5';   % ub = 80.6 x 1.10, headroom for 03/27

SP = load(fullfile(resDir,'feature_spread.mat'));

%% ---- fn entry -> spread source ----------------------------------------
% name in fn            source feature in feature_spread   ('' = regulariser,
%                                                            left untouched)
MAP = { 'ktr'                , 'ktr'
        'A'                  , 'A'
        't0_crossing'        , 't0'                  % data value IS t0
        'restretchSlopeStart', 'restretchSlopeStart'
        'peak1_y'            , 'peak1_y'
        'peak1_dSL'          , 'peak1_dSL'
        'vall_y'             , 'vall_y'
        'vall_t'             , 'vall_t'
        'peak2'              , 'peak2'
        'steady'             , 'steady'
        'vall2_dy'           , 'vall2_dy'
        'rsK'                , 'rsK'
        'PS_restretchPeak'   , '#PS_restretchPeak'
        'PS_steady22'        , '#PS_steady22'
        'PS_steady20'        , '#PS_steady20'
        'FV_normpowerAvg'    , '#FV' };

%% ---- evaluate the base point, get UNWEIGHTED per-feature costs ---------
params0 = getParams(loadParams(BASE), [], true, false);
params0.mods = {}; params0.g = ones(1,0);
params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false; params0.RunSlackPassive = true;
params0.EvalFeatures = true; params0.OptimizeOn = 'Feats';
params0.PlotEachSeparately = 0; params0.PlotFeatureFitting = 0;
params0.RunSlackSegments = 'AllPar'; params0.BreakOnODEUnstable = false;
params0.MaxRunTime = 300;

fn0 = params0.fn;
if DEDUPE
    [~, iu] = unique(fn0, 'stable');
    if numel(iu) < numel(fn0)
        dropped = fn0(setdiff(1:numel(fn0), sort(iu(:))'));
        fprintf('DEDUPE: dropping %d duplicated fn entr(y/ies): %s\n', ...
            numel(dropped), strjoin(dropped, ' , '));
        fn0 = fn0(sort(iu(:))');
    end
end
params0.fn = fn0;

if isempty(gcp('nocreate')); parpool('local', 4); end
fprintf('evaluating base point %s ...\n', BASE);
LoadData; RunBakersExp;                                     %#ok<*NASGU>
[ctot, wOld, cRaw] = evalFeatureCost(features_data, features_model, fn0);
TOT0 = sum(ctot);
fprintf('base total = %.4f over %d fn entries\n', TOT0, numel(fn0));

prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), fn0, 'uni', 0);

%% ---- assemble the new weights -----------------------------------------
n   = numel(fn0);
sF  = nan(1,n); zF = nan(1,n); isData = false(1,n);
for i = 1:n
    k = find(strcmp(MAP(:,1), prim{i}), 1);
    if isempty(k); continue; end
    src = MAP{k,2};
    switch src
        case '#FV'
            sF(i) = SP.s8fv; zF(i) = 0;
        case {'#PS_restretchPeak','#PS_steady22','#PS_steady20'}
            j = find(strcmp(SP.PSN(:,1), src(2:end)), 1);
            sF(i) = SP.s8p(j); zF(i) = 0;
        otherwise
            j = find(strcmp(SP.FEAT(:,1), src), 1);
            sF(i) = SP.s8(j);  zF(i) = SP.zsc(j);
    end
    if ~isfinite(sF(i)); continue; end
    isData(i) = true;
end
s   = min(max(sF, S_FLOOR), S_CEIL);
R   = 1 + KAPPA * min(zF, Z_CAP);
wNewRaw       = nan(1,n);
wNewRaw(isData) = R(isData) ./ s(isData).^2;

% lambda: preserve the DATA-term subtotal exactly
num = sum(wOld(isData)    .* cRaw(isData));
den = sum(wNewRaw(isData) .* cRaw(isData));
lambda = num / den;
wNew = wOld; wNew(isData) = lambda * wNewRaw(isData);
fprintf('lambda = %.6g   (data subtotal %.4f preserved; regularisers untouched)\n', lambda, num);

% --- conditioning guards ------------------------------------------------
% Order matters. The FV pin is an editorial decision and goes FIRST; the share
% cap then polices whatever budget is left. Doing both in one pass can pin more
% mass than `num` and drive the free terms negative.
TOT_SEED = sum(wOld .* cRaw);
pinned   = false(1,n);
iFV      = find(strcmp(prim,'FV_normpowerAvg') & isData, 1);

if FV_SHARE > 0 && ~isempty(iFV) && cRaw(iFV) > 0
    wNew(iFV) = FV_SHARE * TOT_SEED / cRaw(iFV);
    pinned(iFV) = true;
    wNew = iRenorm(wNew, cRaw, isData, pinned, num, FV_SHARE, SHARE_CAP);
    fprintf('  FV pinned at %.0f %% of the seed total (ARBITRARY — see header)\n', 100*FV_SHARE);
end

nCap = 0;
for pass = 1:20
    v    = wNew .* cRaw;
    hitC = isData & ~pinned & v > SHARE_CAP*TOT_SEED & cRaw > 0;
    if ~any(hitC); break; end
    wNew(hitC) = SHARE_CAP * TOT_SEED ./ cRaw(hitC);
    pinned = pinned | hitC;
    wNew = iRenorm(wNew, cRaw, isData, pinned, num, FV_SHARE, SHARE_CAP);
    nCap = nCap + nnz(hitC);
    fprintf('  share cap pass %d: pinned %s\n', pass, strjoin(prim(hitC), ', '));
end
if nCap == 0
    fprintf('  share cap inert at FV_SHARE = %.2f (no other term reaches %.0f %%)\n', ...
            FV_SHARE, 100*SHARE_CAP);
end

%% ---- rewrite fn with the new scalar weights ---------------------------
fnNew = fn0;
for i = find(isData)
    tk = strsplit(fn0{i}, '|');
    % strip a trailing numeric weight if present, then append the new one
    if numel(tk) > 1 && ~isnan(str2double(tk{end})) && ~isRangeTok(tk{end})
        tk(end) = [];
    end
    tk{end+1} = num2str(wNew(i), '%.6g');                               %#ok<SAGROW>
    fnNew{i} = strjoin(tk, '|');
end

%% ---- report -----------------------------------------------------------
lines = {};
lines{end+1} = sprintf('%-20s %8s %8s %9s %9s %9s %9s %9s %9s', ...
    'feature','s8','z','R','w_old','w_new','w x','cost_old','cost_new');
lines{end+1} = repmat('-', 1, 100);
[~, ord] = sort(wNew.*cRaw, 'descend');
for i = ord
    if isData(i)
        tagp = ''; if pinned(i); tagp = ' *'; end
        lines{end+1} = sprintf('%-20s %7.1f%% %8.2f %9.2f %9.4g %9.4g %9.2f %9.4f %9.4f%s', ...
            prim{i}, 100*sF(i), zF(i), R(i), wOld(i), wNew(i), wNew(i)/wOld(i), ...
            wOld(i)*cRaw(i), wNew(i)*cRaw(i), tagp);                    %#ok<SAGROW>
    else
        lines{end+1} = sprintf('%-20s %8s %8s %9s %9.4g %9.4g %9.2f %9.4f %9.4f', ...
            prim{i}, 'reg','-','-', wOld(i), wNew(i), 1, ...
            wOld(i)*cRaw(i), wNew(i)*cRaw(i));                          %#ok<SAGROW>
    end
end
lines{end+1} = repmat('-', 1, 100);
lines{end+1} = sprintf('%-20s %8s %8s %9s %9s %9s %9s %9.4f %9.4f', ...
    'TOTAL','','','','','','', sum(wOld.*cRaw), sum(wNew.*cRaw));
lines{end+1} = sprintf('%-20s %8s %8s %9s %9s %9s %9s %9.4f %9.4f', ...
    ' data terms','','','','','','', num, sum(wNew(isData).*cRaw(isData)));
lines{end+1} = sprintf('%-20s %8s %8s %9s %9s %9s %9s %9.4f %9.4f', ...
    ' regularisers','','','','','','', sum(wOld(~isData).*cRaw(~isData)), sum(wNew(~isData).*cRaw(~isData)));
lines{end+1} = '';
lines{end+1} = 'chi-square view: model residual in units of the measured between-prep spread';
lines{end+1} = sprintf('%-20s %12s %12s', 'feature', 'chi2 = c/s^2', 'per element');
for i = ord
    if ~isData(i); continue; end
    ne = 5; if strcmp(prim{i},'FV_normpowerAvg'); ne = SP.nFV; end
    if startsWith(prim{i},'PS_steady'); ne = 1; end
    lines{end+1} = sprintf('%-20s %12.2f %12.2f', prim{i}, cRaw(i)/s(i)^2, cRaw(i)/s(i)^2/ne); %#ok<SAGROW>
end
%% ---- candidates: data-backed features the objective does NOT score -----
lines{end+1} = '';
lines{end+1} = 'CANDIDATES — measured features not currently in fn (weight/cost if added):';
lines{end+1} = sprintf('%-20s %8s %8s %9s %9s %9s', 'feature','s8','z','w_would','cost_would','share');
for k = 1:numel(CAND)
    j = find(strcmp(SP.FEAT(:,1), CAND{k}), 1);
    if isempty(j) || ~isfinite(SP.s8(j)); continue; end
    sk = min(max(SP.s8(j), S_FLOOR), S_CEIL);
    Rk = 1 + KAPPA*min(SP.zsc(j), Z_CAP);
    wk = lambda * Rk / sk^2;
    ck = NaN;
    if isfield(features_data, CAND{k}) && isfield(features_model, CAND{k})
        [ctk,~,~] = evalFeatureCost(features_data, features_model, {sprintf('%s|%.6g',CAND{k},wk)});
        ck = ctk;
    end
    lines{end+1} = sprintf('%-20s %7.1f%% %8.2f %9.4g %9.4f %8.1f%%', ...
        CAND{k}, 100*SP.s8(j), SP.zsc(j), wk, ck, 100*ck/TOT_SEED);       %#ok<SAGROW>
end

txt = strjoin(lines, newline);
disp(txt);
fid = fopen(fullfile(resDir,'weight_table.txt'),'w'); fprintf(fid,'%s\n',txt); fclose(fid);

%% ---- optional absolute-force anchor guard ------------------------------
if ANCHOR_GUARD
    [ca,~,~] = evalFeatureCost(features_data, features_model, {ANCHOR_ENTRY});
    fnNew{end+1} = ANCHOR_ENTRY;
    fprintf('\nANCHOR_GUARD on: appended %s  (costs %.4g at the seed)\n', ANCHOR_ENTRY, ca);
    if ca > 1e-9
        warning('BuildWeightedObjective:anchorBites', ...
                'anchor guard is NOT zero at the seed (%.4g) — cost-neutrality is broken by that much', ca);
    end
end

%% ---- figure: where the ATP signal is, and where the model is wrong -----
% ATP relevance no longer enters the weights (KAPPA = 0), so it has to be
% VISIBLE instead. Every panel marks the ATP-responsive features in red.
Z_HI = 1.0;                                   % "carries the ATP effect"
fi   = find(isData);
[~,o] = sort(sF(fi));  fi = fi(o);            % best-determined first
lbl  = prim(fi);
isATP = zF(fi) >= Z_HI;
cA = [0.72 0.14 0.18]; cN = [0.40 0.45 0.50];
col = repmat(cN, numel(fi), 1); col(isATP,:) = repmat(cA, nnz(isATP), 1);

fig = figure('Units','pixels','Position',[20 20 1650 900],'Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on
b = bar(100*sF(fi),'FaceColor','flat'); b.CData = col;
set(gca,'XTick',1:numel(fi),'XTickLabel',lbl,'XTickLabelRotation',45,'FontSize',11,'TickLabelInterpreter','none');
ylabel('between-prep spread  s  (%)','FontSize',12);
title({'How well is each feature determined?','3 preps, 8 mM, Baker excluded — red = carries the ATP effect (z \geq 1)'}, ...
      'FontSize',12,'FontWeight','bold');

nexttile; hold on; box on; grid on
b = bar(zF(fi),'FaceColor','flat'); b.CData = col;
yline(Z_HI,'k--','z = 1','LineWidth',1.3,'FontSize',11);
set(gca,'XTick',1:numel(fi),'XTickLabel',lbl,'XTickLabelRotation',45,'FontSize',11,'TickLabelInterpreter','none');
ylabel('ATP discriminability  z','FontSize',12);
title({'Where the ATP effect actually is','|ratio-1| in noise widths — every force SHAPE feature is flat'}, ...
      'FontSize',12,'FontWeight','bold');

nexttile; hold on; box on; grid on
chi = cRaw(fi)./s(fi).^2;
b = bar(chi,'FaceColor','flat'); b.CData = col;
yline(1,'k--','one preparation-width','LineWidth',1.3,'FontSize',11);
set(gca,'XTick',1:numel(fi),'XTickLabel',lbl,'XTickLabelRotation',45,'FontSize',11, ...
        'YScale','log','TickLabelInterpreter','none');
ylabel('\chi^2 = cost / s^2','FontSize',12);
title({'How far is the model from the data, in noise units?', ...
       sprintf('%d of %d terms at or below the noise floor — rsK is the defect', nnz(chi<4), numel(chi))}, ...
      'FontSize',12,'FontWeight','bold');

nexttile; hold on; box on; grid on
% log scale: with FV pinned at 70 % a linear axis hides everything else
Sh = max([wOld(fi).*cRaw(fi); wNew(fi).*cRaw(fi)]' / TOT_SEED * 100, 1e-3);
b = bar(Sh,'grouped'); b(1).FaceColor = [0.75 0.75 0.75]; b(2).FaceColor = [0.15 0.35 0.55];
set(gca,'YScale','log'); ylim([1e-3 200]);
for k = find(isATP(:)')      % red tick under the ATP-relevant features
    plot(k, 1.5e-3, '^', 'MarkerFaceColor', cA, 'MarkerEdgeColor', cA, 'MarkerSize', 7);
end
set(gca,'XTick',1:numel(fi),'XTickLabel',lbl,'XTickLabelRotation',45,'FontSize',11,'TickLabelInterpreter','none');
ylabel('share of the objective (%)','FontSize',12);
legend({'realign4 (hand-set)','spread-weighted'},'Location','northeast','FontSize',11);
title({sprintf('Objective composition — total held at %.4f', TOT_SEED), ...
       sprintf('FV pinned at %.0f %% (arbitrary); everything else = \\lambda/s^2', 100*FV_SHARE)}, ...
      'FontSize',12,'FontWeight','bold');

exportgraphics(fig, fullfile(resDir,'feature_weighting.png'), 'Resolution', 160);
fprintf('wrote %s\n', fullfile(resDir,'feature_weighting.png'));

%% ---- write the reweighted seed ----------------------------------------
params0.fn = fnNew; params0.mods = {}; params0.g = [];
snap = fullfile('params', [SEEDTAG '.m']);
writeParamsToMFile(snap, params0, [], ...
    sprintf('%s reweighted by measured 8 mM between-prep spread (kappa=%g, lambda=%.6g); cost-neutral at %.4f', ...
            BASE, KAPPA, lambda, TOT0));
fprintf('\nwrote %s\n', snap);

save(fullfile(resDir,'weighted_objective.mat'), ...
     'fn0','fnNew','wOld','wNew','cRaw','sF','zF','R','isData','lambda','TOT0','prim','BASE','SEEDTAG');

%% ---- verify neutrality -------------------------------------------------
pv = params0; pv.fn = fnNew;
[cv,~,~] = evalFeatureCost(features_data, features_model, fnNew);
fprintf('\nNEUTRALITY CHECK: old total %.4f  ->  new total %.4f  (delta %.2e)\n', ...
        TOT0, sum(cv), sum(cv)-TOT0);

function tf = isRangeTok(t)
    tf = ~isempty(regexp(t, '^-?[\d.]+\s*-\s*-?[\d.]+$', 'once'));
end

function w = iRenorm(w, cRaw, isData, pinned, num, FV_SHARE, SHARE_CAP)
% Free (unpinned) data terms carry whatever budget the pins left over, so the
% data subtotal stays exactly `num` and the grand total stays cost-neutral.
    free  = isData & ~pinned;
    resid = num - sum(w(pinned & isData) .* cRaw(pinned & isData));
    if resid <= 0 || ~any(free)
        error('BuildWeightedObjective:infeasible', ...
            ['pinned terms claim %.4f of the %.4f data budget — lower ' ...
             'FV_SHARE (now %.2f) or SHARE_CAP (now %.2f).'], ...
            num - resid, num, FV_SHARE, SHARE_CAP);
    end
    w(free) = w(free) * (resid / sum(w(free) .* cRaw(free)));
end
