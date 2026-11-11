% VerifyPooledRun.m -- re-evaluate the pooled8 gallery before trusting it.
%
% WHY THIS EXISTS. `optimizeFeatures` has no repeat-evaluation guard on accept,
% and the objective's own numerical roughness is ~0.076 sd / 0.25 peak-to-peak
% (Docs objective_noise_floor). Twice now a rough evaluation has bought a step
% that was not there: in the spreadw campaign `spreadw_iter_013` carried a
% header of 1.4769 and re-evaluated at 1.6615 (+12.5 %) on a deterministic
% model. So every snapshot header is a claim, not a measurement, until it has
% been re-run here.
%
% Also scores the winner on the OLD single-preparation objective, which is the
% only way to see whether the pooled refit moved the model or just moved the
% yardstick.
%
% Output -> results/pooled8_verify.txt, results/pooled8_verify.mat

clearvars -except TAG NWORKERS; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('TAG','var');      TAG      = 'pooled8'; end
if ~exist('NWORKERS','var'); NWORKERS = 5;         end
SEED_TAG = 'pooled8_seed';
BASE_TAG = 'estiff125_opt';        % what the campaign started from
NOISE    = 0.25;                   % objective peak-to-peak roughness

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

%% ---- one common objective ---------------------------------------------
% Everything below is scored with the WINNER's fn and the WINNER's target, so
% the numbers are comparable across snapshots. iP() applies both.
pW  = getParams(loadParams([TAG '_opt']), [], true, false);
FN  = pW.fn;  TGT = pW.UsePooledTarget;

iter = dir(fullfile('params', [TAG '_iter'], [TAG '_iter_*.m']));
[~, ord] = sort({iter.name}); iter = iter(ord);
tags = [{SEED_TAG}, cellfun(@(f) erase(f,'.m'), {iter.name}, 'uni', 0), {[TAG '_opt']}];

hdr = nan(1, numel(tags));      % the cost each snapshot's provenance claims
for i = 1:numel(iter); hdr(i+1) = iHeaderCost(fullfile(iter(i).folder, iter(i).name)); end

L = {};
L{end+1} = sprintf('%-24s %10s %10s %9s   %s', 'snapshot','header','re-eval','delta','verdict');
L{end+1} = repmat('-', 1, 78);
got = nan(1, numel(tags)); ct = cell(1, numel(tags)); nmAll = [];
for i = 1:numel(tags)
    p = iP(tags{i}, FN, TGT);
    [c, ~, cti, nmi] = evaluateBakersExp(1, p, true, false);
    got(i) = c; ct{i} = cti; if isempty(nmAll); nmAll = nmi; end
    if isfinite(hdr(i))
        d = got(i) - hdr(i);
        if abs(d) <= NOISE; v = 'ok'; else; v = sprintf('*** %+.1f %% ***', 100*d/hdr(i)); end
        L{end+1} = sprintf('%-24s %10.4f %10.4f %+9.4f   %s', tags{i}, hdr(i), got(i), d, v); %#ok<SAGROW>
    else
        L{end+1} = sprintf('%-24s %10s %10.4f %9s   %s', tags{i}, '-', got(i), '-', 'no header'); %#ok<SAGROW>
    end
end

%% ---- did the model move, or only the yardstick? -----------------------
% Four numbers. The diagonal is each parametrisation on its own objective; the
% off-diagonal is the cross-check. If the winner beats the base on the OLD
% objective too, the refit found genuinely better parameters. If it does not,
% the pooled target bought pooled-target performance at 03/27's expense --
% which is legitimate (the pool is the better target) but must be SAID.
%
% EACH COLUMN CARRIES ITS OWN fn, and it has to. The single-preparation column
% cannot use the winner's fn: that list scores `rsK80`, and the curated
% 03/27 file predates the metric audit, so the field is simply absent from its
% features_data and evalFeatureCost charges a flat MISSING_FEATURE_COST = 100.
% Scoring the old target with the new list reports ~64 for BOTH rows -- a
% missing field, not a fit. (The real fix is re-curating the protocol files to
% carry rsK80/rsT20; until then the old column uses the old list.)
L{end+1} = '';
L{end+1} = 'the model, or the yardstick?  (rows = parametrisation, cols = objective)';
L{end+1} = sprintf('%-24s %14s %14s', '', 'POOLED (3 prep)', 'single 03/27');
FN_OLD = getParams(loadParams(BASE_TAG), [], true, false).fn;
cross  = nan(2,2); rows = {BASE_TAG, [TAG '_opt']};
for r = 1:2
    cross(r,1) = evaluateBakersExp(1, iP(rows{r}, FN,     true),  true, false);
    cross(r,2) = evaluateBakersExp(1, iP(rows{r}, FN_OLD, false), true, false);
    L{end+1} = sprintf('%-24s %14.4f %14.4f', rows{r}, cross(r,1), cross(r,2)); %#ok<SAGROW>
end
L{end+1} = sprintf('%-24s %+13.1f%% %+13.1f%%', 'change', ...
    100*(cross(2,1)-cross(1,1))/cross(1,1), 100*(cross(2,2)-cross(1,2))/cross(1,2));

%% ---- per feature, in chi-square units ---------------------------------
% chi2 per slack segment: 1.0 = one preparation-width from the pooled mean.
% This is the only scale on which "is the fit good" has an answer.
PT = pooledTarget8mM();
w  = zeros(1, numel(nmAll)); s = nan(1, numel(nmAll));
for i = 1:numel(nmAll)
    j = find(strcmp(PT.name, nmAll{i}), 1);
    if isempty(j) && strcmp(nmAll{i}, 't0_crossing'); j = find(strcmp(PT.name,'t0'),1); end
    if ~isempty(j); s(i) = PT.s(j); end
    w(i) = iWeight(FN{i});
end
iSeed = 1; iWin = numel(tags);
L{end+1} = '';
L{end+1} = 'per feature (chi2/segment vs the pooled spread; 1.0 = one preparation-width)';
L{end+1} = sprintf('%-22s %9s %9s %9s   %8s %8s', ...
    'feature','seed','winner','delta','chi2 sd','chi2 win');
[~, o] = sort(ct{iWin} - ct{iSeed});
for i = o(:)'
    if ct{iSeed}(i) < 1e-9 && ct{iWin}(i) < 1e-9; continue; end
    if isfinite(s(i)) && w(i) > 0
        x0 = ct{iSeed}(i)/w(i)/s(i)^2/5;  x1 = ct{iWin}(i)/w(i)/s(i)^2/5;
        xs = sprintf('%8.2f %8.2f', x0, x1);
    else
        xs = sprintf('%8s %8s','-','-');
    end
    L{end+1} = sprintf('%-22s %9.4f %9.4f %+9.4f   %s', ...
        nmAll{i}, ct{iSeed}(i), ct{iWin}(i), ct{iWin}(i)-ct{iSeed}(i), xs); %#ok<SAGROW>
end

txt = strjoin(L, newline); disp(txt);
fid = fopen(fullfile(resDir,'pooled8_verify.txt'),'w'); fprintf(fid,'%s\n',txt); fclose(fid);
save(fullfile(resDir,'pooled8_verify.mat'), 'tags','hdr','got','ct','nmAll','cross','w','s','FN');
fprintf('\nwrote results/pooled8_verify.{txt,mat}\n');

%% ======================= helpers ========================================
function p = iP(tag, FN, TGT)
    p = getParams(loadParams(tag), [], true, false);
    p.mods = {}; p.g = ones(1,0);
    p.fn = FN; p.UsePooledTarget = TGT;
    p.RunForceVelocity = true; p.RunKtr = false; p.RunSlack = true;
    p.RunStairs = false; p.RunForceVelocityTime = false; p.RunSlackPassive = true;
    p.EvalFeatures = true; p.OptimizeOn = 'Feats';
    p.RunSlackSegments = 'AllPar'; p.BreakOnODEUnstable = false;
    p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0; p.MaxRunTime = 300;
end

function c = iHeaderCost(f)
% writeParamsToMFile puts the provenance comment on line 2. recordOptIter
% writes it as "% param_opt_iter #3 | cost 3.1939 | round 3 | direct", so the
% number is NOT at the end of the line — match the label, not the position.
    c = NaN;
    fid = fopen(f,'r'); if fid < 0; return; end
    fgetl(fid); l = fgetl(fid); fclose(fid);
    if ~ischar(l); return; end
    t = regexp(l, 'cost\s+([\d.eE+-]+)', 'tokens', 'once');
    if ~isempty(t); c = str2double(t{1}); end
end

function w = iWeight(entry)
% Trailing scalar weight of an fn entry, or 0 for boundary/grouped entries.
    w = 0;
    if contains(entry, ',') || contains(entry, '@'); return; end
    tk = strsplit(entry, '|');
    v  = str2double(tk{end});
    if ~isnan(v) && isempty(regexp(tk{end}, '^-?[\d.]+\s*-\s*-?[\d.]+$', 'once')); w = v; end
end
