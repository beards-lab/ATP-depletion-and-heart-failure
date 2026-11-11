% RunObjectiveFloor.m -- how much of the objective is unreachable?
%
% Three terms in params0.fn score a model trace against a data trace by least
% squares, with the data value pinned at 0 (a perfect model would score 0):
%
%     coolDownLS         0.15 * mean(resid^2) over the 150 ms after the
%                        re-stretch of each active slack cycle
%     PS_rampupRMSE      RMSE over the passive re-stretch ramp
%     PS_holdDecayRMSE   RMSE over the passive post-re-stretch hold
%
% That 0 is not attainable. The data carries white transducer noise of
% sigma ~ 1.2 kPa (RunNoiseBudget: measured in situ, identical with the
% cross-bridges switched off), so even a model sitting exactly on the truth
% scores RMSE = sigma. The part of each term below its floor is not a fit, it
% is the rig, and any weight spent on it is spent on nothing.
%
% This script measures sigma in each of those windows -- from the data alone,
% by the two independent estimators that agreed in RunNoiseBudget -- turns it
% into a floor for each term, and reports what fraction of the objective at
% params/pooled8_opt.m is irreducible.
%
% Output -> results/objective_floor.mat / .txt

clearvars -except TAG; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('TAG','var'); TAG = 'pooled8_opt'; end

% ---- the two windows, straight out of the extractors -----------------------
ACT = 'protocol_03_27_2026_8mM_slack.mat';           % coolDownLS
PAS = 'protocol_03_27_2026_ActivePNBMava_slack.mat'; % PS_*
COOLDOWN_LS_SCALE = 0.15;                            % Model/extractSlackAttributes.m

S = load(fullfile('data',ACT));  dtA = S.datatable;  vtA = S.velocitytable;
P = load(fullfile('data',PAS));  dtP = P.datatable;  vtP = P.velocitytable;

F = struct();

% ---- coolDownLS: the cool-down window of each active slack cycle -----------
% extractSlackAttributes decimates the data 1:10 before the LSQ, which does
% not average -- it subsamples -- so the per-sample noise is unchanged.
segA = find(vtA(:,2) < -1)';
sw = [];
for i = 1:numel(segA)
    s = segA(i);
    vs4 = vtA(min(s+3,size(vtA,1)),1);  vs5 = vtA(min(s+4,size(vtA,1)),1);
    w = dtA(:,1) >= vs4 & dtA(:,1) <= min(vs5, vs4+0.15);
    if nnz(w) < 40; continue; end
    sw(end+1) = noiseSigma(dtA(w,1), dtA(w,3)); %#ok<SAGROW>
end
F.sig_cool = mean(sw);
F.floor_coolDownLS = COOLDOWN_LS_SCALE * F.sig_cool^2;

% ---- PS_rampupRMSE / PS_holdDecayRMSE: the passive windows ----------------
segP = find(vtP(:,2) < -1)';
sr = []; sh = [];
for i = 1:numel(segP)
    s = segP(i);
    t_rs = vtP(s+2,1); t_hold = vtP(s+3,1);
    if s+4 > size(vtP,1); t_end = vtP(end,1); else; t_end = vtP(s+4,1); end
    wr = dtP(:,1) >= t_rs   & dtP(:,1) <= t_hold;
    wh = dtP(:,1) >= t_hold & dtP(:,1) <= t_end;
    if nnz(wr) > 12; sr(end+1) = noiseSigma(dtP(wr,1), dtP(wr,3)); end %#ok<SAGROW>
    if nnz(wh) > 40; sh(end+1) = noiseSigma(dtP(wh,1), dtP(wh,3)); end %#ok<SAGROW>
end
F.sig_ramp = mean(sr);  F.floor_PS_rampupRMSE    = F.sig_ramp;
F.sig_hold = mean(sh);  F.floor_PS_holdDecayRMSE = F.sig_hold;

% the passive signal itself, for scale
F.passiveSpan = range(dtP(dtP(:,1)>=vtP(segP(1),1),3));

% ---- what the model actually scores ---------------------------------------
params0 = getParams(loadParams(TAG), [], true, false);
params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
params0.RunStairs = false; params0.RunForceVelocityTime = false;
params0.RunSlackPassive = true; params0.EvalFeatures = true;
params0.BreakOnODEUnstable = false; params0.PlotEachSeparately = 0;
params0.PlotFeatureFitting = 0; params0.RunSlackSegments = 'AllPar';
params0.OptimizeOn = 'Feats'; params0.mods = {}; params0.g = ones(1,0);
[Et, ~, featCost, featNames] = evaluateBakersExp(1, params0, true, false);

% weights, straight out of fn
wOf = @(nm) fnWeight(params0.fn, nm);

TERMS = {'coolDownLS', F.floor_coolDownLS, COOLDOWN_LS_SCALE
         'PS_rampupRMSE',    F.floor_PS_rampupRMSE,    1
         'PS_holdDecayRMSE', F.floor_PS_holdDecayRMSE, 1};

fid = fopen(fullfile(resDir,'objective_floor.txt'),'w');
p = @(varargin) tee(fid, varargin{:});
p('THE UNREACHABLE PART OF THE OBJECTIVE      (%s, total %.4f)\n', TAG, Et);
p('%s\n', repmat('=',1,92));
p('sigma is the white transducer noise IN THAT WINDOW, from the data alone.\n');
p('A model sitting exactly on the truth still scores the floor.\n\n');
p('%-20s %9s %9s %9s %9s %9s %8s\n','term','sigma','floor','model','cost','floorCost','irred%');
tot = 0; totFloor = 0;
for i = 1:size(TERMS,1)
    nm = TERMS{i,1}; fl = TERMS{i,2};
    j  = find(strcmp(featNames, nm), 1);
    if isempty(j); continue; end
    w  = wOf(nm);
    % cost = w * mean(model_value^2) because the data value is 0
    modelVal = sqrt(featCost(j)/max(w,eps));
    floorCost = w * fl^2;
    tot = tot + featCost(j); totFloor = totFloor + min(floorCost, featCost(j));
    if strcmp(nm,'coolDownLS'); sg = F.sig_cool; else; sg = F.(sprintf('sig_%s', ternary(strcmp(nm,'PS_rampupRMSE'),'ramp','hold'))); end
    p('%-20s %9.3f %9.4f %9.4f %9.4f %9.4f %7.0f %%\n', nm, sg, fl, modelVal, ...
      featCost(j), floorCost, 100*min(floorCost,featCost(j))/featCost(j));
end
p('%s\n', repmat('-',1,92));
p('trace-LSQ terms total       %.4f   = %.1f %% of the objective\n', tot, 100*tot/Et);
p('of which irreducible        %.4f   = %.1f %% of the objective\n', totFloor, 100*totFloor/Et);
p('\npassive force span in the PNB run: %.1f kPa, noise sigma %.2f kPa -> SNR %.1f\n', ...
  F.passiveSpan, F.sig_ramp, F.passiveSpan/F.sig_ramp);
p('\nrsR2 scores %.4f (boundary 0.4-1.0, model 0.999 sits inside). A noise-free\n', ...
  featCost(strcmp(featNames,'rsR2')));
p('model CANNOT reach the data value 0.795 -- its ceiling is 0.834 -- so this\n');
p('term must stay a boundary and must never be turned into a fit-to-data term.\n');
fclose(fid);
F.Et = Et; F.featCost = featCost; F.featNames = {featNames{:}};
save(fullfile(resDir,'objective_floor.mat'),'F');

% ------------------------------------------------------------------------
function s = noiseSigma(t, y)
%NOISESIGMA  White noise in a window, by the two estimators that agreed in
%   RunNoiseBudget: content above 60 Hz, and the first difference. The smaller
%   is taken, so the floor this produces is a LOWER bound on the floor.
    fs = 1/median(diff(t));
    hp = y - movmean(y, max(3, round(fs/60)));
    s1 = std(hp);
    if numel(y) < 3*max(3, round(fs/60)); s1 = Inf; end   % window too short for the high-pass
    s2 = std(diff(y,2))/sqrt(6);   % 2nd difference: immune to any linear ramp
    s  = min(s1, s2);
end

function w = fnWeight(fn, nm)
    w = NaN;
    for i = 1:numel(fn)
        tk = split(fn{i}, '|');
        if strcmp(strtrim(tk{1}), nm)
            w = 1;
            for j = 2:numel(tk)
                v = str2double(tk{j});
                if ~isnan(v); w = v; end
            end
            return
        end
    end
end

function v = ternary(c, a, b); if c; v = a; else; v = b; end; end
function tee(fid, varargin); fprintf(1, varargin{:}); fprintf(fid, varargin{:}); end
