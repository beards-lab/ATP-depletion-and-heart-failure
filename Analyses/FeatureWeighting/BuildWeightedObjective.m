% BuildWeightedObjective.m — turn the measured data spread into objective weights.
%
% Reads  results/feature_spread.mat  (from RunFeatureSpread.m) and the current
% best parametrisation, and produces a REWEIGHTED params0.fn plus a seed
% snapshot that scores EXACTLY the same total as the current objective, so a
% new optimisation starts from a comparable number.
%
% Caller-settable knobs (see the block below): BASE, SEEDTAG, SPREAD_SOURCE,
% POOLED_TARGET, FN_RENAME, ADD_FEATS, FV_SHARE, SHARE_CAP, KAPPA,
% ANCHOR_GUARD. SPREAD_SOURCE / POOLED_TARGET are what switch this script from
% the single-preparation objective to the 3-preparation pooled one; see
% Analyses/DataPooling/RunOptimPooled.m for the campaign that uses them.
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

clearvars -except BASE SEEDTAG SPREAD_SOURCE POOLED_TARGET FN_RENAME ADD_FEATS FLOOR_CHI2 ...
                  FV_SHARE SHARE_CAP KAPPA ANCHOR_GUARD ANCHOR_ENTRY FN_DROP W_BOOST; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results');

%% ---- knobs -------------------------------------------------------------
% All of these can be set by the caller before running this script; the
% `clearvars -except` above keeps them. A bare `clear` here would silently
% ignore them, which is how two earlier runs were mis-attributed.
if ~exist('BASE','var');    BASE    = 'realign4_opt'; end  % be cost-neutral at this parametrisation
if ~exist('SEEDTAG','var'); SEEDTAG = 'spreadw_seed'; end  % params/<SEEDTAG>.m
% KAPPA = 0 => R == 1, i.e. the weights are PURE PRECISION and nothing else.
% This is deliberate. A weight that also encodes "relevance" is no longer a
% variance, so `cost/s^2` stops being a chi-square and the whole table stops
% being readable as a goodness-of-fit -- which is the most useful thing this
% analysis produces. ATP relevance (`z`) is still measured, still reported, and
% still drawn (results/feature_weighting.png); it just does not touch w.
% The honest place for the ATP effect is a term that SCORES the 2 mM condition
% (conclusions.md Sec 8), not a tilt smuggled into the 8 mM weights.
% Set KAPPA = 0.5 to restore the tilt; R spans 1.0 .. 1+KAPPA*Z_CAP.
if ~exist('KAPPA','var'); KAPPA = 0; end
Z_CAP    = 4;
S_FLOOR  = 0.05;               % no feature is better determined than 5 % across preps
S_CEIL   = 1.00;               % ... nor worse than 100 %
DEDUPE   = true;               % drop the duplicated rsR2 boundary entry

% --- WHICH TARGET, AND WHOSE SPREAD ------------------------------------
% POOLED_TARGET  goes straight into params0.UsePooledTarget, so the base point
%   below is evaluated against the 3-preparation POOLED slack target rather
%   than the single curated preparation (Model/pooledTargetFeatures). It is
%   also written into the seed snapshot, so the optimisation inherits it.
% SPREAD_SOURCE  where s_f comes from:
%   'spread'  results/feature_spread.mat — RAW between-prep SD, force features
%             normalised by their own run's mean `steady`, no rundown lift.
%   'pooled'  Model/pooledTarget8mM — the SAME three preparations, but anchored
%             on Fss22, rundown-corrected and variance-moderated (Smyth 2004,
%             d0 = 2.74 against d = 2). Moderation is the substantive change:
%             it compresses the 1/s^2 dynamic range 446x -> 13x, which stops
%             `steady` (4 % raw, a within-run repeatability, not a between-prep
%             spread) from taking the largest weight in the objective, and
%             stops `vall_t` (84 % raw) from being effectively unscored.
%             FV and the passive PS_* terms are NOT pooled yet, so they always
%             fall back to feature_spread.mat.
% The two knobs are independent on purpose: 'pooled' spreads with the old
% target isolates the reweighting, and vice versa.
if ~exist('POOLED_TARGET','var'); POOLED_TARGET = false;    end
if ~exist('SPREAD_SOURCE','var'); SPREAD_SOURCE = 'spread'; end
% FLOOR_CHI2 -- deadband, in preparation-widths, below which a data term stops
% paying. 0 = off (previous behaviour exactly). 1.0 = "being inside the
% measurement's own resolution is free". See the block where it is applied.
if ~exist('FLOOR_CHI2','var');    FLOOR_CHI2    = 0;        end
SPREAD_SOURCE = validatestring(SPREAD_SOURCE, {'spread','pooled'});

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
if ~exist('FV_SHARE','var');  FV_SHARE  = 0.70; end
if ~exist('SHARE_CAP','var'); SHARE_CAP = 0.30; end
% Candidate features NOT currently scored; reported, added only if listed here.
if ~exist('ADD_FEATS','var'); ADD_FEATS = {}; end   % e.g. {'rsT0','vall2_t','Am'}
CAND         = {'rsT0','rsT20','vall2_t','Am','peak1_t','rsA'};

% FN_RENAME — {old, new} pairs applied to the fn entries before anything else.
% Its use is swapping a feature for a BETTER ESTIMATOR OF THE SAME QUANTITY,
% e.g. {'rsK','rsK80'}: Analyses/RestretchMetricAudit showed rsK reads its
% amplitude off the tail of the fit window, so a 60 ms window shift moved it by
% +77 %; rsK80 = ln4/(t80-t20) takes the amplitude from `steady` instead. The
% pooled target is built on rsK80, so scoring rsK against a pooled number would
% mix the audited target with the unaudited estimator.
% This is a rename, not an addition: the weight is re-derived from the new
% feature's own spread and the term count does not change.
%
% A source containing '[' matches the WHOLE first token, selector included, so
% a segment restriction can be lifted: {'rsK80[1:4]','rsK80'} puts slack cycle
% 5 back. Without the bracket the match is on the bare name and any selector is
% preserved, which is the original behaviour.
if ~exist('FN_RENAME','var'); FN_RENAME = cell(0,2); end

% FN_DROP -- primary feature names to REMOVE from the objective entirely.
% Unlike a rename this changes the term count, so the dropped term's old
% contribution leaves the cost-neutrality budget with it and lambda re-spreads
% that mass over what remains. Used to retire `coolDownLS`, an unbounded raw
% MSE that is not a preparation-width chi-2 like every other data term and
% that dominated the movers of nearly every MechanismAblation arm
% (Analyses/MechanismAblation/conclusions.md sec. 8).
if ~exist('FN_DROP','var'); FN_DROP = {}; end

% W_BOOST -- {name, factor} pairs multiplying a data term's derived weight,
% applied AFTER lambda and after the FV pin, then pinned so the renormalisation
% does not undo them. A boost is a deliberate departure from w = lambda*R/s^2,
% so it breaks that term's chi-2 reading and must be justified in the caller.
% The one in use: rsK80 x2, inheriting the post-re-stretch role that dropping
% coolDownLS vacated (MechanismAblation sec. 8).
if ~exist('W_BOOST','var'); W_BOOST = cell(0,2); end

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
if ~exist('ANCHOR_GUARD','var'); ANCHOR_GUARD = false; end
% Range = what the three preparations actually span. Under the POOLED target
% the force block is already at the pooled level, so the guard's job changes
% from "do not drift away from 03/27" to "stay inside the measured spread".
if ~exist('ANCHOR_ENTRY','var')
    ANCHOR_ENTRY = 'steady[1]|50.6-88.7|5';   % ub = 80.6 x 1.10, headroom for 03/27
end

SP = load(fullfile(resDir,'feature_spread.mat'));

% Pooled spreads override SP's per-feature entries when SPREAD_SOURCE is
% 'pooled'. FV / PS_* have no pooled counterpart and keep SP's numbers.
if strcmp(SPREAD_SOURCE,'pooled')
    PT = pooledTarget8mM();
    for k = 1:numel(PT.name)
        j = find(strcmp(SP.FEAT(:,1), PT.name{k}), 1);
        if isempty(j)
            SP.FEAT(end+1,1) = PT.name(k);  j = size(SP.FEAT,1);
        end
        SP.s8(j)  = PT.s(k);
        SP.zsc(j) = PT.z(k);
    end
    fprintf(['SPREAD_SOURCE = pooled : s from Model/pooledTarget8mM (%d preps %s, ' ...
             'moderated dof %.2f)\n'], PT.nP, strjoin(PT.days,'/'), PT.dof);
    fprintf('  FV and PS_* keep feature_spread.mat (not pooled yet)\n');
end

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
        'rsK80'              , 'rsK80'               % audited rate  (RestretchMetricAudit)
        'rsT20'              , 'rsT20'               % audited delay
        'rsT0'               , 'rsT0'
        'rsA'                , 'rsA'
        'Am'                 , 'Am'
        'peak1_t'            , 'peak1_t'
        'vall2_t'            , 'vall2_t'
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
params0.UsePooledTarget = POOLED_TARGET;

fn0 = params0.fn;

%% ---- fn surgery: rename estimators, add unscored features --------------
% Both happen BEFORE the base evaluation, so cRaw and the neutrality budget
% are computed on the entry set that will actually be optimised.
renamed = {};
for r = 1:size(FN_RENAME,1)
    src = FN_RENAME{r,1};
    % A bracketed source addresses the whole token (name AND selector), so the
    % selector can be replaced or lifted. A bare source matches the name only
    % and leaves any selector in place.
    withSel = contains(src, '[');
    for i = 1:numel(fn0)
        tk = strsplit(fn0{i}, '|');
        if withSel
            if ~strcmp(tk{1}, src); continue; end
            tk{1} = FN_RENAME{r,2};
        else
            if ~strcmp(regexprep(tk{1}, '\[.*$', ''), src); continue; end
            tk{1} = regexprep(tk{1}, ['^' src], FN_RENAME{r,2});
        end
        fn0{i} = strjoin(tk, '|');
        renamed{end+1} = sprintf('%s -> %s', src, FN_RENAME{r,2}); %#ok<SAGROW>
    end
end
if ~isempty(renamed)
    fprintf('FN_RENAME: %s\n', strjoin(unique(renamed,'stable'), ' , '));
end

% Drop retired terms. Done here, before the base evaluation, so the dropped
% term is absent from cRaw and therefore from the neutrality budget -- its mass
% is released to the survivors rather than silently preserved.
if ~isempty(FN_DROP)
    primD = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), fn0, 'uni', 0);
    kill  = ismember(primD, FN_DROP);
    if any(kill)
        fprintf('FN_DROP: removing %s\n', strjoin(fn0(kill), ' , '));
        fn0 = fn0(~kill);
    end
    miss = setdiff(FN_DROP, primD);
    if ~isempty(miss)
        fprintf('FN_DROP: not present, nothing to drop: %s\n', strjoin(miss, ', '));
    end
end

% Added entries carry weight 1 for now; the loop below gives them lambda/s^2.
% They are held OUT of the neutrality budget (isAdded) — they were not in the
% old objective, so pretending the total is unchanged would hide them.
nOld  = numel(fn0);
for a = 1:numel(ADD_FEATS)
    if any(strcmp(regexprep(fn0, '[|\[].*$', ''), ADD_FEATS{a})); continue; end
    fn0{end+1} = ADD_FEATS{a};                                          %#ok<SAGROW>
end
if numel(fn0) > nOld
    fprintf('ADD_FEATS: appended %s (excluded from the neutrality budget)\n', ...
        strjoin(fn0(nOld+1:end), ', '));
end

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
isAdded = ismember(prim, ADD_FEATS);

%% ---- assemble the new weights -----------------------------------------
n   = numel(fn0);
sF  = nan(1,n); zF = nan(1,n); isData = false(1,n);
for i = 1:n
    % A BOUNDARY entry is not a data term, whatever it is named. evalFeatureCost
    % never reads feats_data for one -- it scores the simulation against a fixed
    % range -- so 1/s^2 is meaningless for it and rewriting its weight silently
    % destroys a hand-set guardrail. This bit the anchor guard
    % 'steady[1]|52.5-87.7|5': `steady` is in MAP, so it was being treated as a
    % fit term and its deliberate weight of 5 rewritten to 0.46.
    if any(cellfun(@isRangeTok, strsplit(fn0{i}, '|'))); continue; end
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

% lambda: preserve the DATA-term subtotal exactly. Newly ADDED entries are
% held out of the budget — they had no old weight to preserve, so folding them
% in would silently shrink every other term to pay for them.
isBudget = isData & ~isAdded;
num = sum(wOld(isBudget)    .* cRaw(isBudget));
den = sum(wNewRaw(isBudget) .* cRaw(isBudget));
lambda = num / den;
wNew = wOld; wNew(isData) = lambda * wNewRaw(isData);
fprintf('lambda = %.6g   (data subtotal %.4f preserved; regularisers untouched)\n', lambda, num);
if any(isAdded)
    fprintf('  added terms sit ON TOP: %s = %.4f\n', ...
        strjoin(prim(isAdded),'+'), sum(wNew(isAdded).*cRaw(isAdded)));
end

% --- conditioning guards ------------------------------------------------
% Order matters. The FV pin is an editorial decision and goes FIRST; the share
% cap then polices whatever budget is left. Doing both in one pass can pin more
% mass than `num` and drive the free terms negative.
TOT_SEED = sum(wOld(~isAdded) .* cRaw(~isAdded));
pinned   = isAdded;                 % added terms are not renormalised
iFV      = find(strcmp(prim,'FV_normpowerAvg') & isData, 1);

if FV_SHARE > 0 && ~isempty(iFV) && cRaw(iFV) > 0
    wNew(iFV) = FV_SHARE * TOT_SEED / cRaw(iFV);
    pinned(iFV) = true;
    wNew = iRenorm(wNew, cRaw, isBudget, pinned, num, FV_SHARE, SHARE_CAP);
    fprintf('  FV pinned at %.0f %% of the seed total (ARBITRARY — see header)\n', 100*FV_SHARE);
end

% Deliberate weight boosts. Pinned after boosting so iRenorm redistributes the
% cost off the OTHER free terms rather than undoing the boost it was just told
% to make.
for b = 1:size(W_BOOST,1)
    iB = find(strcmp(prim, W_BOOST{b,1}) & isData, 1);
    if isempty(iB) || cRaw(iB) <= 0
        fprintf('  W_BOOST: %s is not a scored data term — skipped\n', W_BOOST{b,1});
        continue;
    end
    wNew(iB)   = wNew(iB) * W_BOOST{b,2};
    pinned(iB) = true;
    wNew = iRenorm(wNew, cRaw, isBudget, pinned, num, FV_SHARE, SHARE_CAP);
    fprintf('  W_BOOST: %s weight x%.3g -> %.1f %% of the seed total\n', ...
            W_BOOST{b,1}, W_BOOST{b,2}, 100*wNew(iB)*cRaw(iB)/TOT_SEED);
end

nCap = 0;
for pass = 1:20
    v    = wNew .* cRaw;
    hitC = isBudget & ~pinned & v > SHARE_CAP*TOT_SEED & cRaw > 0;
    if ~any(hitC); break; end
    wNew(hitC) = SHARE_CAP * TOT_SEED ./ cRaw(hitC);
    pinned = pinned | hitC;
    wNew = iRenorm(wNew, cRaw, isBudget, pinned, num, FV_SHARE, SHARE_CAP);
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

%% ---- noise floors ------------------------------------------------------
% FLOOR_CHI2 > 0 gives every data term a deadband: it stops paying once it is
% within FLOOR_CHI2 preparation-widths of the target, and starts paying again
% if it drifts back out. The cost is a sum of squared relative residuals, so
% one preparation-width over nEl elements is nEl*s^2.
%
% This is the fix for the failure the rigaware campaign demonstrated. A fixed
% objective budget is spent wherever descent is cheapest, and that is whichever
% term is already best determined: F-V held 40 % of the objective at chi-2 0.10
% and absorbed 110 % of a 14-hour gain while every feature the campaign existed
% to fix went backwards. Re-pinning the WEIGHT cannot fix that -- any nonzero
% weight on a smooth, nearly-solved term still buys the cheapest descent. The
% floor removes the reward below the noise instead, without removing the term's
% defence of itself.
floorRaw = zeros(1,n);
if FLOOR_CHI2 > 0
    nFloored = 0;
    for i = find(isData)
        if ~isfinite(s(i)) || s(i) <= 0; continue; end
        nEl = iNumEl(fnNew{i}, features_data);
        if nEl < 1; continue; end
        floorRaw(i) = FLOOR_CHI2 * nEl * s(i)^2;
        tk = strsplit(fnNew{i}, '|');
        tk = [tk(1:end-1), {sprintf('fl:%.6g', floorRaw(i))}, tk(end)];  % before the weight
        fnNew{i} = strjoin(tk, '|');
        nFloored = nFloored + 1;
    end
    fprintf(['FLOOR_CHI2 = %.2f : %d data terms given a deadband at %.2f ' ...
             'preparation-width(s)\n'], FLOOR_CHI2, nFloored, FLOOR_CHI2);
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
    % Drop every existing guard on the same feature+range before adding the
    % canonical one, REGARDLESS of the weight it carries. The BASE snapshot may
    % already have one (it was built with the guard on), and until the isData
    % fix above a boundary entry could have its weight silently rewritten --
    % rigaware_a_opt carries both 'steady[1]|52.5-87.7|5' and a corrupted
    % '...|0.547416'. Two guards on one range double it and neither is the
    % number anyone chose.
    [aF, aR] = iGuardKey(ANCHOR_ENTRY);
    dup = false(1, numel(fnNew));
    for i = 1:numel(fnNew)
        [f2, r2] = iGuardKey(fnNew{i});
        dup(i) = ~isempty(f2) && strcmp(f2, aF) && strcmp(r2, aR);
    end
    if any(dup)
        fprintf('\nANCHOR_GUARD: replacing %d existing guard(s) on %s|%s: %s\n', ...
                nnz(dup), aF, aR, strjoin(fnNew(dup), ' , '));
        fnNew(dup) = [];
    end
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
[~, iWorst] = max(chi);
title({'How far is the model from the data, in noise units?', ...
       sprintf('%d of %d terms at or below the noise floor — worst is %s (%.0f)', ...
               nnz(chi<4), numel(chi), lbl{iWorst}, chi(iWorst))}, ...
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
legend({sprintf('%s (incoming)', BASE), sprintf('%s-weighted', SPREAD_SOURCE)}, ...
       'Location','northeast','FontSize',11,'Interpreter','none');
title({sprintf('Objective composition — total held at %.4f  (target: %s)', TOT_SEED, ...
               ternary(isequal(POOLED_TARGET,false), 'single prep', 'POOLED')), ...
       sprintf('FV pinned at %.0f %% (arbitrary); everything else = \\lambda/s^2', 100*FV_SHARE)}, ...
      'FontSize',12,'FontWeight','bold');

exportgraphics(fig, fullfile(resDir,'feature_weighting.png'), 'Resolution', 160);
fprintf('wrote %s\n', fullfile(resDir,'feature_weighting.png'));

%% ---- write the reweighted seed ----------------------------------------
params0.fn = fnNew; params0.mods = {}; params0.g = [];
snap = fullfile('params', [SEEDTAG '.m']);
writeParamsToMFile(snap, params0, [], ...
    sprintf(['%s reweighted by %s 8 mM between-prep spread (kappa=%g, lambda=%.6g); ' ...
             'target=%s; cost-neutral at %.4f'], ...
            BASE, SPREAD_SOURCE, KAPPA, lambda, ...
            ternary(isequal(POOLED_TARGET,false),'single-prep','pooled'), TOT0));
fprintf('\nwrote %s\n', snap);

save(fullfile(resDir,'weighted_objective.mat'), ...
     'fn0','fnNew','wOld','wNew','cRaw','sF','zF','R','isData','lambda','TOT0','prim','BASE','SEEDTAG');

%% ---- verify neutrality -------------------------------------------------
pv = params0; pv.fn = fnNew;
[cv,~,~] = evalFeatureCost(features_data, features_model, fnNew);
if FLOOR_CHI2 > 0
    % Cost-neutrality is deliberately broken by a floor: the floor's whole job
    % is to stop paying for residual that is inside the noise, so the total MUST
    % fall. Report what it removed rather than flagging it as drift.
    fprintf(['\nFLOORED TOTAL: %.4f  ->  %.4f   (%.4f = %.0f %% of the objective was\n' ...
             '  residual inside one preparation-width, and now earns nothing)\n'], ...
            TOT0, sum(cv), TOT0-sum(cv), 100*(TOT0-sum(cv))/TOT0);
    fprintf('  the seed for the next campaign is %.4f, NOT %.4f\n', sum(cv), TOT0);
else
    fprintf('\nNEUTRALITY CHECK: old total %.4f  ->  new total %.4f  (delta %.2e)\n', ...
            TOT0, sum(cv), sum(cv)-TOT0);
end

function tf = isRangeTok(t)
    tf = ~isempty(regexp(t, '^-?[\d.]+\s*-\s*-?[\d.]+$', 'once'));
end

function o = ternary(c, a, b); if c; o = a; else; o = b; end; end

function [feat, rng] = iGuardKey(entry)
%IGUARDKEY  Identify a boundary entry by feature+selector and range, ignoring
% its weight, so two guards on the same range are recognised as the same guard.
    feat = ''; rng = '';
    if contains(entry, ',') || contains(entry, '@'); return; end   % grouped/custom
    tk = strsplit(entry, '|');
    ir = find(cellfun(@isRangeTok, tk), 1);
    if isempty(ir); return; end
    feat = tk{1}; rng = tk{ir};
end

function n = iNumEl(entry, fd)
%INUMEL  How many elements this fn entry actually scores, selector included.
% The floor is nEl*s^2, so it has to count the elements that survive a
% selector like rsK80[1:4] -- not the raw field length.
    n = 0;
    [nm, sel] = parseFeatSelector(strtok(entry, '|'));
    if ~isfield(fd, nm); return; end
    n = numel(applyFeatSelector(double([fd.(nm)]), sel));
end

function w = iRenorm(w, cRaw, isBudget, pinned, num, FV_SHARE, SHARE_CAP)
% Free (unpinned) data terms carry whatever budget the pins left over, so the
% data subtotal stays exactly `num` and the grand total stays cost-neutral.
% isBudget excludes newly ADDED entries: those are not part of the preserved
% total, so they must not absorb or release any of it.
    free  = isBudget & ~pinned;
    resid = num - sum(w(pinned & isBudget) .* cRaw(pinned & isBudget));
    if resid <= 0 || ~any(free)
        error('BuildWeightedObjective:infeasible', ...
            ['pinned terms claim %.4f of the %.4f data budget — lower ' ...
             'FV_SHARE (now %.2f) or SHARE_CAP (now %.2f).'], ...
            num - resid, num, FV_SHARE, SHARE_CAP);
    end
    w(free) = w(free) * (resid / sum(w(free) .* cRaw(free)));
end
