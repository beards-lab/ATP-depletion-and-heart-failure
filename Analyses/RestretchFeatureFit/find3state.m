% find3state.m
% =========================================================================
% Hand-search for a NON-COLLAPSED 3-state parameter seed (TWO-STAGE, fast).
%
% Diagnosis (see labdiary.md "3-state investigation"): topology is
% p1->p2->p3->detach(PT), so p2's ONLY exit is R2 = k2*PieceWiseStrainDep2(strain).
% That strain shape was tuned as 2-state DETACHMENT (ramps UP at negative/
% shortening strain); as the p2->p3 FEED it drains force-bearing p2 during
% shortening -> FV tail collapses to 0 (no force at velocity).
%
% FIX under test: flatten the negative-strain R2 knots
% (PieceWiseStrainDep2Params indices 5,6,7; currently ~[9.4 16 50]) toward
% small values; if still collapsed, also drop k2 hard so p2 lingers.
%
% STAGE 1 (fast): FV-only screen, 4 velocities, MaxRunTime=20, no slack/ktr.
%                 Pass = FV_fnorm(end) > 0.05 and no NaN.  (~4 evals/config)
% STAGE 2 (slow): top 1-2 stage-1 configs -> full costOfSnap (FV+slack) for
%                 total feature cost; save the best (or least-collapsed).
%
% Output of record: Analyses/RestretchFeatureFit/find3state_results.txt
% (fclose+reopen-append each iteration so it survives a process kill --
%  MATLAB text streams have no fflush, so we physically reopen).
% Base snapshot: params/params_3state_DRAFT.m (NumberOfStates=3; ODE has the
% p3(p3<0)=0 clamp -> should not hang).
% =========================================================================

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

baseSnap = 'params/params_3state_DRAFT.m';
resFile  = 'Analyses/RestretchFeatureFit/find3state_results.txt';

% --- feature list (verbatim from RunOptim3State.m / task spec) ---
USERFN = {'FV_fnorm|FV_v|10', 'ktr|2', 'A|50', 'ktr_rmse|0-.5|.1', ...
    'XTOR[1]|3-15|15,XTOR_vmax[1]|5-25,SRX_ss[1]|0.01-0.40|5,attached_ss[1]|0.10-0.45|10,PT_ss[1]|0.01-0.70',...
    't0_crossing|SLdiff|2', 'restretchSlopeStart|1', 'peak1_y|10', 'peak1_dSL|1', 'vall_y|10', ...
    'vall_t|0.2', 'peak2|5', 'steady|50', 'vall2_dy|1.0', 'ovrsht_dy|1', ...
    'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001'};

% --- crash-safe append helper (defined at end of file) ---
logline = @(s) appendLine(resFile, s);

% start the results file fresh
fid = fopen(resFile, 'w');
fprintf(fid, '# find3state two-stage search  %s\n', datestr(now));
fprintf(fid, '# STAGE 1: FV-only screen (4 vel, MaxRunTime=20). pass = FV_fnorm(end)>0.05, no NaN\n');
fclose(fid);

% --- build an override struct (scalars + __N array knots) ---
mkcfg = @(k2, k3, k3m, kstiff3, ka, r5, r6, r7) struct( ...
    'k2', k2, 'k3', k3, 'k3m', k3m, 'kstiff3', kstiff3, 'ka', ka, ...
    'PieceWiseStrainDep2Params__5', r5, 'PieceWiseStrainDep2Params__6', r6, ...
    'PieceWiseStrainDep2Params__7', r7);

% =====================================================================
% STAGE 1 grid.  k2 in {80,100,120}, k3 in {150,250}, k3m=10,
% kstiff3 in {22000,28000}, ka=160, R2 knots in
% {[1 1.5 3],[0.5 0.8 2],[0.3 0.5 1]}. => 3(R2)*3(k2)*2(k3)*2(ks3)=36 configs.
% =====================================================================
R2sets  = {[1 1.5 3], [0.5 0.8 2], [0.3 0.5 1]};
R2names = {'r113', 'r052', 'r031'};

cfg = {}; names = {};
for ir = 1:numel(R2sets)
    R = R2sets{ir};
    for k2 = [80 100 120]
        for k3 = [150 250]
            for ks3 = [22000 28000]
                cfg{end+1} = mkcfg(k2, k3, 10, ks3, 160, R(1), R(2), R(3)); %#ok<SAGROW>
                names{end+1} = sprintf('%s_k2%d_k3%d_ks%d', R2names{ir}, k2, k3, ks3); %#ok<SAGROW>
            end
        end
    end
end
nStage1 = numel(cfg);

stage1 = struct('name', {}, 'fvEnd', {}, 'FV', {}, 'ok', {}, 'cfg', {});
for i = 1:nStage1
    [FVfnorm, fvEnd, okS, msg] = fvScreen(baseSnap, cfg{i});
    if isempty(msg)
        logline(sprintf('[S1 %2d/%d] %-22s FV=%-20s fvEnd=%.4f ok=%d', ...
            i, nStage1, names{i}, mat2str(round(FVfnorm,3)), fvEnd, okS));
    else
        logline(sprintf('[S1 %2d/%d] %-22s FAIL: %s', i, nStage1, names{i}, msg));
    end
    stage1(end+1) = struct('name', names{i}, 'fvEnd', fvEnd, 'FV', FVfnorm, 'ok', okS, 'cfg', cfg{i}); %#ok<SAGROW>
end

% ---- STAGE 1b rescue: if nothing passed, drop k2 hard so p2 lingers ----
if ~any([stage1.ok])
    logline('# STAGE 1b RESCUE: no config passed -> lowering k2 to {40,60}');
    cfgB = {}; namesB = {};
    for ir = 1:numel(R2sets)
        R = R2sets{ir};
        for k2 = [40 60]
            for k3 = [150 250]
                cfgB{end+1} = mkcfg(k2, k3, 10, 28000, 160, R(1), R(2), R(3)); %#ok<SAGROW>
                namesB{end+1} = sprintf('%s_k2%d_k3%d_ks28000', R2names{ir}, k2, k3); %#ok<SAGROW>
            end
        end
    end
    for i = 1:numel(cfgB)
        [FVfnorm, fvEnd, okS, msg] = fvScreen(baseSnap, cfgB{i});
        if isempty(msg)
            logline(sprintf('[S1b %2d/%d] %-22s FV=%-20s fvEnd=%.4f ok=%d', ...
                i, numel(cfgB), namesB{i}, mat2str(round(FVfnorm,3)), fvEnd, okS));
        else
            logline(sprintf('[S1b %2d/%d] %-22s FAIL: %s', i, numel(cfgB), namesB{i}, msg));
        end
        stage1(end+1) = struct('name', namesB{i}, 'fvEnd', fvEnd, 'FV', FVfnorm, 'ok', okS, 'cfg', cfgB{i}); %#ok<SAGROW>
    end
end

% ---- rank stage-1 by fvEnd (desc), NaN last ----
fvEnds = [stage1.fvEnd];
fvEnds(isnan(fvEnds)) = -Inf;
[~, order] = sort(fvEnds, 'descend');

logline('# STAGE 1 RANKING (by FV tail, best first):');
for j = 1:min(6, numel(order))
    s = stage1(order(j));
    logline(sprintf('#   %-22s fvEnd=%.4f FV=%s', s.name, s.fvEnd, mat2str(round(s.FV,3))));
end

% Candidates for stage 2: passing configs preferred; else least-collapsed.
passOrder = order(arrayfun(@(o) stage1(o).ok, order));
if ~isempty(passOrder)
    stage2idx = passOrder(1:min(2, numel(passOrder)));
    logline(sprintf('# STAGE 2: full costOfSnap on top %d PASSING config(s)', numel(stage2idx)));
    anyPass = true;
else
    validOrder = order(~isinf(fvEnds(order)));
    if isempty(validOrder)
        logline('# ALL STAGE-1 CONFIGS FAILED (NaN/error). Nothing saved.');
        return;
    end
    stage2idx = validOrder(1);
    logline('# STAGE 2: NO config passed FV>0.05; running least-collapsed non-NaN config only.');
    anyPass = false;
end

% ---- STAGE 2: full cost, pick best ----
bestCost = Inf; bestCfg = []; bestName = ''; bestFV = NaN;
bestKtr = NaN; bestSteady = NaN; bestPeak1y = NaN;
for jj = 1:numel(stage2idx)
    s = stage1(stage2idx(jj));
    tc = NaN; ktrMean = NaN; steadyMean = NaN; peak1yMean = NaN; FVfull = s.FV; msg = '';
    try
        [tc, ~, fm, ~] = costOfSnap(baseSnap, USERFN, 30, s.cfg);
        FVfull = fm.FV_fnorm;
        ktrMean = mean(fm.ktr, 'omitnan');
        steadyMean = mean(fm.steady, 'omitnan');
        peak1yMean = mean(fm.peak1_y, 'omitnan');
    catch ME
        msg = ME.message;
    end
    if isempty(msg)
        logline(sprintf('[S2] %-22s cost=%.3f FV=%s ktr=%.1f steady=%.1f peak1_y=%.1f', ...
            s.name, tc, mat2str(round(FVfull,3)), ktrMean, steadyMean, peak1yMean));
    else
        logline(sprintf('[S2] %-22s FAIL: %s', s.name, msg));
    end
    if (~isnan(tc) && tc < bestCost) || isempty(bestCfg)
        bestCost = tc; bestCfg = s.cfg; bestName = s.name;
        bestFV = FVfull; bestKtr = ktrMean; bestSteady = steadyMean; bestPeak1y = peak1yMean;
    end
end

% ---- SAVE winner ----
if isempty(bestCfg)
    logline('# NO WINNER (stage 2 produced nothing). Nothing saved.');
    return;
end

params0 = getParams();
run(baseSnap);
ef = fieldnames(bestCfg);
for k = 1:numel(ef); params0.(ef{k}) = bestCfg.(ef{k}); end
params0 = getParams(params0, [], true, false);
params0.mods = {};
params0.g = [];
writeParamsToMFile('params/params_3state_seed.m', params0);

logline(sprintf('# WINNER: %s  cost=%.3f  FV=%s  ktr=%.1f  steady=%.1f  peak1_y=%.1f', ...
    bestName, bestCost, mat2str(round(bestFV,3)), bestKtr, bestSteady, bestPeak1y));
if ~anyPass
    logline('# NOTE: winner did NOT pass FV>0.05 (least-collapsed fallback) -- FV tail still weak.');
end
logline(sprintf('SAVED %.3f', bestCost));

% =====================================================================
% local helpers (must be at end of a script file)
% =====================================================================
function appendLine(fname, s)
    f = fopen(fname, 'a');
    fprintf(f, '%s\n', s);
    fclose(f);
end

function [FVfnorm, fvEnd, okS, msg] = fvScreen(baseSnap, c)
% Fast FV-only screen for one override config. Returns FV_fnorm, its tail,
% pass flag, and an error message ('' on success).
    FVfnorm = NaN; fvEnd = NaN; okS = false; msg = '';
    try
        params0 = getParams(); %#ok<NASGU>
        run(baseSnap);
        ef = fieldnames(c);
        for k = 1:numel(ef); params0.(ef{k}) = c.(ef{k}); end
        params0.RunForceVelocity = true;
        params0.RunSlack = false;
        params0.RunKtr = false;
        params0.RunStairs = false;
        params0.RunForceVelocityTime = false;
        params0.EvalFeatures = true;
        params0.FV_velocities = -[0 1 2 4];
        params0.MaxRunTime = 20;
        params0.BreakOnODEUnstable = false;
        params0.PlotEachSeparately = 0;
        params0.PlotFeatureFitting = 0;
        params0 = getParams(params0, [], true, false);
        RunBakersExp;                       % produces features_model in this scope
        FVfnorm = features_model.FV_fnorm;
        fvEnd   = FVfnorm(end);
        okS     = (fvEnd > 0.05) && ~any(isnan(FVfnorm));
    catch ME
        msg = ME.message;
    end
end
