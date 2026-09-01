% RunFeatureSpread.m — measure, from the DATA alone, how well each scored
% feature is determined, and how strongly it responds to ATP.
%
% This is the empirical input to the objective's weights. No model is run.
%
% Two numbers per feature:
%
%   s8  BETWEEN-PREPARATION SPREAD AT 8 mM  (relative SD, 3 preps)
%       = how reproducible the number is. This is the honest noise scale for a
%         residual: asking the model to beat it is asking it to fit one
%         preparation's idiosyncrasy.
%       Force-dimension features are divided by their OWN run's mean `steady`
%       first (SlackDataAnalysis §4: absolute force is not transferable, CV
%       19-28 %; shape is, CV 4-7 %). `steady` itself is the one ABSOLUTE
%       anchor and keeps its raw CV -- it is what sets the scale.
%
%   z   ATP DISCRIMINABILITY = |mean ratio - 1| / sqrt(s8^2 + sr^2)
%       where ratio = 2 mM / 8 mM per prep (rundown-corrected for force-
%       dimension features, raw for rates/times -- ATPEffectReconciliation),
%       and sr is the cross-prep spread OF THE RATIO.
%       = how many noise widths the ATP effect moves this feature. A feature
%         with z ~ 0 cannot carry the ATP story no matter how well it is fit.
%
% The Baker slack recordings are DELIBERATELY EXCLUDED: different protocol,
% recovery windows 2-6x too short, amplitude features truncated ATP-dependently
% (SlackDataAnalysis §1-2). Only the three 2026 protocol days are used.
%
% Features are re-extracted fresh through extractSlackAttributes on every
% dataset rather than read from the stored features_data, so every number comes
% from one code path (SlackDataAnalysis's rule).
%
% Output -> results/feature_spread.mat  + printed table.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results'); if ~isfolder(resDir), mkdir(resDir); end
D = 'data/';

PHI = 0.452;                        % permanent rundown fraction (RundownCorrection)
zA  = [71.50 72.35]; zB = [77.70 78.30];   % quiescent windows for the rundown slope

%% ---- the three 2026 protocol days (Baker excluded on purpose) -----------
% day, earlier cond, earlier merged .txt, earlier slack .mat, later cond, later merged, later .mat
P = { '03/27 M','8 mM','03 27 2026 M/02_Merged_8mM_Active.txt','protocol_03_27_2026_8mM_slack.mat', ...
                '2 mM','03 27 2026 M/03_Merged_2mM_Active.txt','protocol_03_27_2026_2mM_slack.mat'
      '04/03 F','8 mM','04 03 2026 F/02_Merged_8mM_Active.txt','protocol_04_03_2026_8mM_slack.mat', ...
                '2 mM','04 03 2026 F/03_Merged_2mM_Active.txt','protocol_04_03_2026_2mM_slack.mat'
      '04/10 M','2 mM','04 10 2026 Male 2-8/02_Merged_2mM_Active.txt','protocol_04_10_2026_2mM_slack.mat', ...
                '8 mM','04 10 2026 Male 2-8/03_Merged_8mM_Active.txt','protocol_04_10_2026_8mM_slack.mat' };
nP = size(P,1);
PNB = { 'protocol_03_27_2026_ActivePNBMava_slack.mat'
        'protocol_04_03_2026_ActivePNBMava_slack.mat'
        'protocol_04_10_2026_ActivePNBMava_slack.mat' };

%% ---- feature registry: name, class -------------------------------------
%  F = force-dimension (normalise by own steady; rundown-correctable)
%  S = the absolute scale anchor (raw CV, no normalisation)
%  R = rate      T = time / length      Q = dimensionless quality
FEAT = { 'steady','S'; 'A','F'; 'Am','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'rsA','F'
         'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         't0','T'; 'peak1_t','T'; 'peak1_dSL','T'; 'vall_t','T'; 'vall2_t','T'; 'rsT0','T'
         'rsR2','Q'; 'ktr_rmse','Q' };
nF = size(FEAT,1);

%% ---- rundown dose per prep ---------------------------------------------
frac = nan(1,nP);
fprintf('=== rundown dose ===\n%-9s %-12s %10s %9s %8s\n','day','order','loss(kPa)','F_later','frac');
for i = 1:nP
    Me = readmatrix([D P{i,3}]); te = Me(:,1); Fe = Me(:,3);
    m = (te>=zA(1)&te<=zA(2)) | (te>=zB(1)&te<=zB(2)); pf = polyfit(te(m),Fe(m),1);
    Fs = movmean(Fe,51); thr = 0.2*prctile(Fs,99);
    on = find(Fs>thr,1); off = find(Fs>thr,1,'last'); Ta = te(off)-te(on);
    loss = PHI*abs(pf(1))*Ta;
    Ml = readmatrix([D P{i,6}]); tl = Ml(:,1);
    Fl = mean(Ml(tl>=zA(1)&tl<=zA(2),3),'omitnan');
    frac(i) = loss/Fl;
    fprintf('%-9s %-12s %10.2f %9.2f %7.1f%%\n', P{i,1}, sprintf('%s->%s',P{i,2},P{i,5}), loss, Fl, 100*frac(i));
end

%% ---- fresh feature extraction, one code path ---------------------------
F8 = cell(1,nP); F2 = cell(1,nP);
for i = 1:nP
    is8first = strcmp(P{i,2},'8 mM');
    fe = extractFresh([D P{i,4}]);  fl = extractFresh([D P{i,7}]);
    if is8first, F8{i} = fe; F2{i} = fl; else, F8{i} = fl; F2{i} = fe; end
end
FP = cell(1,numel(PNB));
for i = 1:numel(PNB); FP{i} = extractFresh([D PNB{i}]); end

%% ---- s8: between-prep spread at 8 mM -----------------------------------
s8 = nan(1,nF); nEl = nan(1,nF);
for j = 1:nF
    nm = FEAT{j,1}; cls = FEAT{j,2};
    V = collect(F8, nm, cls, F8);        % nP x nEl, normalised per class
    if isempty(V); continue; end
    s8(j) = relSD(V); nEl(j) = size(V,2);
end

%% ---- ATP ratio, its spread, and z --------------------------------------
rat = nan(nP,nF);
for i = 1:nP
    s8i = mean(double(F8{i}.steady),'omitnan');
    s2i = mean(double(F2{i}.steady),'omitnan');
    for j = 1:nF
        nm = FEAT{j,1}; cls = FEAT{j,2};
        if ~isfield(F8{i},nm) || ~isfield(F2{i},nm); continue; end
        v8 = mean(double(F8{i}.(nm)),'omitnan');
        v2 = mean(double(F2{i}.(nm)),'omitnan');
        if ~isfinite(v8) || ~isfinite(v2) || v8 == 0; continue; end
        r = v2/v8;
        switch cls
            case {'F','S'}
                % undo the depression of whichever run was recorded SECOND
                if strcmp(P{i,2},'8 mM'); r = r*(1+frac(i)); else; r = r/(1+frac(i)); end
                if cls == 'F'
                    % shape ratio: each run scaled by its own steady first, so
                    % what is left is ATP's effect on SHAPE, not on level.
                    r = (v2/s2i)/(v8/s8i);
                end
        end
        rat(i,j) = r;
    end
end
rbar = mean(rat,1,'omitnan');
sr   = std(rat,0,1,'omitnan');
zsc  = abs(rbar-1) ./ sqrt(s8.^2 + sr.^2);

%% ---- passive (PNB+Mava) features ---------------------------------------
% PS_restretchPeak / PS_steady22 / PS_steady20 are the passive data features
% that runPassiveExperiment scores. Approximated here by their slack-extraction
% equivalents on the three PNB files: peak2 (post-restretch peak) and steady.
PSN = {'PS_restretchPeak','peak2'; 'PS_steady22','steady'; 'PS_steady20','steady'};
s8p = nan(1,size(PSN,1));
for j = 1:size(PSN,1)
    V = collect(FP, PSN{j,2}, 'S', FP);     % passive force IS absolute here
    if ~isempty(V); s8p(j) = relSD(V); end
end

%% ---- FV: between-SOURCE spread of the normalised power ------------------
% Already consumed per-point by 'w:FV_normpowerVar'. This is the aggregate
% relative spread, so FV can be placed on the same 1/s^2 footing as the rest.
AllV = [0, 0.5, 1, 2, 3, 4, 5, 6];
S1 = [56.40, 51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];
S2 = [67.3942, 60.9885, 42.7059, 18.8539, 12.5956, 6.4426, 3.5136, 1.8556];
ip = AllV > 0;
Pw = [ (S1(ip)/S1(1)).*AllV(ip); (S2(ip)/S2(1)).*AllV(ip) ]';
s8fv = relSD(Pw');  nFV = size(Pw,1);

%% ---- report ------------------------------------------------------------
fprintf('\n=== 8 mM reproducibility and ATP response (3 preps, Baker excluded) ===\n');
fprintf('%-20s %3s %4s %8s %9s %8s %8s %7s\n','feature','cls','nEl','s8(rel)','ATP ratio','sr','z','1/s8^2');
for j = 1:nF
    fprintf('%-20s %3s %4d %7.1f%% %9.3f %7.1f%% %8.2f %7.1f\n', FEAT{j,1}, FEAT{j,2}, nEl(j), ...
        100*s8(j), rbar(j), 100*sr(j), zsc(j), 1/s8(j)^2);
end
fprintf('\n--- passive (PNB+Mava, 3 preps; absolute) ---\n');
for j = 1:size(PSN,1)
    fprintf('%-20s %3s %4d %7.1f%% %9s %7s %8s %7.1f\n', PSN{j,1}, 'S', 1, 100*s8p(j), '-','-','-', 1/s8p(j)^2);
end
fprintf('\n--- FV normalised power (2 sources, %d points) ---\n', nFV);
fprintf('%-20s %3s %4d %7.1f%% %9s %7s %8s %7.1f\n','FV_normpowerAvg','F',nFV,100*s8fv,'-','-','-',1/s8fv^2);

save(fullfile(resDir,'feature_spread.mat'), ...
     'FEAT','s8','nEl','rat','rbar','sr','zsc','frac','PSN','s8p','s8fv','nFV','P','PNB');
fprintf('\nsaved -> %s\n', fullfile(resDir,'feature_spread.mat'));

%% ======================= helpers ========================================
function f = extractFresh(file)
    S = load(file);
    f = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), S.datatable(:,2), ...
                               S.velocitytable, [], [], false, true);
end

function V = collect(FS, nm, cls, FSnorm)
% nP x nEl matrix of feature nm, normalised per prep for class 'F'.
    nP = numel(FS); V = [];
    for i = 1:nP
        if ~isfield(FS{i},nm); return; end
        v = double(FS{i}.(nm)(:))';
        if cls == 'F'
            v = v / mean(double(FSnorm{i}.steady),'omitnan');
        end
        if isempty(V); V = nan(nP, numel(v)); end
        if numel(v) ~= size(V,2); return; end     % length mismatch -> skip feature
        V(i,:) = v;                                                     %#ok<AGROW>
    end
end

function s = relSD(V)
% RMS over elements of the between-row relative SD. Matches the SUM-over-
% elements form of evalFeatureCost, so 1/s^2 turns each term into a chi-square.
    mu = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
    good = isfinite(mu) & abs(mu) > 0 & isfinite(sd);
    if ~any(good); s = NaN; return; end
    s = sqrt(mean((sd(good)./abs(mu(good))).^2));
end
