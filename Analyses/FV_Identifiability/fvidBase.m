function [pBat, pBeard, pool, recs] = fvidBase()
% fvidBase  The incumbent in the two configurations this analysis evaluates, the
% perturbed parameter pool, and the three Beard2022 records that drive the model.
%
%   pBat    robust1_seed set up as the slack-battery objective (slack + ktr +
%           passive + F-V on the pooled 8 mM target), as RunRefit9Optim does
%   pBeard  robust1_seed + fit A (k2 power gate, R2 x0.54 at 2 mM), as
%           Analyses/FV_Beard2022/RunFVBeard2022 does
%   pool    the optimiser's lean pool (zero-valued parameters dropped)
%   recs    one record per [ATP] (8, 4, 2 mM): t (ms), L (Lo), vt (velocity table).
%           Within one [ATP] all 12 / 6 / 6 records share their timing to < 0.1 ms,
%           so one record per [ATP] gives the model trace for all of them.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'FV_Beard2022'));

% ---- slack battery --------------------------------------------------------
pBat = getParams(); s = loadParams('robust1_seed');
f = fieldnames(s); for i = 1:numel(f), pBat.(f{i}) = s.(f{i}); end
pBat.RunForceVelocity = true; pBat.RunSlack = true; pBat.RunStairs = false; pBat.RunKtr = true;
pBat.RunForceVelocityTime = false; pBat.RunSlackPassive = true;
pBat.EvalFeatures = true; pBat.BreakOnODEUnstable = false;
pBat.PlotEachSeparately = 0; pBat.PlotFeatureFitting = 0;
pBat.RunSlackSegments = 'AllPar'; pBat.MaxRunTime = 300;
pBat.OptimizeOn = 'Feats'; pBat.ghostSave = ''; pBat.ghostLoad = '';
pBat.SaveBest = 0; pBat.mods = {}; pBat.g = ones(1, 0);

% ---- Beard2022 records ------------------------------------------------------
pBeard = getParams(loadParams('robust1_seed'), [], true, false);
pBeard.mods = {}; pBeard.g = ones(1, 0);
pBeard.EvalFeatures = false; pBeard.BreakOnODEUnstable = false; pBeard.MaxRunTime = 1800;
pBeard.PlotEachSeparately = 0; pBeard.PlotFeatureFitting = 0; pBeard.PlotFullscreen = 0;
pBeard.ghostSave = ''; pBeard.ghostLoad = ''; pBeard.SaveBest = 0;
pBeard.UseAtpPowerGates = true; pBeard.MgATP_ref = 8; pBeard.UseAtpDetachR2 = 0; pBeard.UseAtpKmsrd = 0;
pBeard.pAtp_k2 = log(0.54) / log(1/4);                % fit A: R2 x0.54 at 2 mM

% ---- parameter pool ---------------------------------------------------------
pool = pBat.leanPool; keep = true(1, numel(pool));
for i = 1:numel(pool)
    nm = pool{i}; us = strfind(nm, '__');
    if isempty(us)
        keep(i) = isfield(pBat, nm) && isscalar(pBat.(nm)) && pBat.(nm) ~= 0;
    else
        b = nm(1:us(1)-1); ix = str2double(nm(us(1)+2:end));
        keep(i) = isfield(pBat, b) && numel(pBat.(b)) >= ix && pBat.(b)(ix) ~= 0;
    end
end
pool = pool(keep);

% ---- records ----------------------------------------------------------------
if nargout > 3
    cf = fullfile(here, 'results', 'beard_records.mat');
    if isfile(cf)
        S = load(cf); recs = S.recs;
    else
        D = loadBeard2022FV(); recs = struct('atp', {}, 'prep', {}, 't', {}, 'L', {}, 'vt', {});
        for a = [8 4 2]
            k = find([D.atp] == a, 1);
            recs(end+1) = struct('atp', a, 'prep', D(k).prep, 't', D(k).t, 'L', D(k).L, ...
                'vt', buildVelocityTable(D(k).t, D(k).L, 5e-4, 2)); %#ok<AGROW>
        end
        if ~isfolder(fileparts(cf)), mkdir(fileparts(cf)); end
        save(cf, 'recs');
    end
end
end
