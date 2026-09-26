%% RunModelNoise.m -- numerical noise of every slack feature definition, per solver setting.
%
% Same 1 + NREP parameter points in every configuration (the incumbent and
% NREP copies with each lean-pool scalar multiplied by 1 + 1e-6*randn). One
% slack run per point with the objective's own semantics (RunSlackSegments =
% 'AllPar': init hold, then every cycle from the same init state), serial
% inside a parfor worker. The raw merged trace is saved, so every feature
% definition is read from the IDENTICAL trace:
%   pipe : features_model exactly as the objective computes it now
%          (extractSlackAttributes' internal 1 kHz resample, anchored on the
%          trace's first stored time)
%   a1k  : robustSlackFeatures on a 1 kHz grid anchored on absolute time
%          (whole milliseconds) -- current definitions AND replacements
%   a4k  : the same on a 4 kHz grid, the rate the data are recorded at
% plus, on the unperturbed point of every configuration, a PHASE test: the
% trace read at 1 kHz with the grid shifted by 0, 0.1, ... 0.9 ms.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('NREP', 'var'), NREP = 8; end
if ~exist('NW', 'var'), NW = 5; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root)); res = fullfile(here, 'results');
fx = @(varargin) struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, varargin{:});
if ~exist('CFG', 'var')
CFG = { 'A old: 1e-2, step output, hard rectifier', struct('SmoothRectifierEps', 0, 'ClampStrainDepExtrapolation', 0, 'odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 0)
        'B SiLU 1e-3, step output',                 fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 0)
        'C SiLU 1e-3, dense in cycles only',        fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'D SiLU 1e-2, dense in cycles only',        fx('odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'E reference: SiLU 1e-4, dense in cycles',  fx('odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5) };
end
if ~exist('CFGSEL', 'var'), CFGSEL = 1:size(CFG, 1); end

P = getParams(loadParams(PARAMSET), [], true, false);
P.RunSlackSegments = 'AllPar'; P.MaxRunTime = 1200; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.BreakOnODEUnstable = false; P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
D = load(fullfile(root, 'data', P.velocitytableonfile)); vt = D.velocitytable; vt(1,1) = -20;
rng(11); pool = P.leanPool; DR = randn(NREP, numel(pool));
if isempty(gcp('nocreate')), parpool(NW); end

outFile = fullfile(res, ['model_noise_' PARAMSET '.mat']);
if isfile(outFile), load(outFile, 'M'); else, M = struct(); end
for c = CFGSEL
    ov = CFG{c, 2}; Pc = P; f = fieldnames(ov); for k = 1:numel(f), Pc.(f{k}) = ov.(f{k}); end
    Pc = getParams(Pc, Pc.g, true, false);
    TR = cell(NREP+1, 1); FM = cell(NREP+1, 1); FDc = cell(NREP+1, 1); CT = nan(NREP+1, 1); W = nan(NREP+1, 1);
    fnS = Pc.fn;
    parfor r = 1:NREP+1
        pv = Pc;
        if r > 1
            for i = 1:numel(pool)
                nm = pool{i};
                if ~contains(nm, '__') && isfield(pv, nm) && isscalar(pv.(nm)) && pv.(nm) ~= 0
                    pv.(nm) = pv.(nm) * (1 + 1e-6*DR(r-1, i));
                end
            end
            pv = getParams(pv, pv.g, true, false);
        end
        t0 = tic;
        [~, os, fm, fd] = runSlackExperiment(pv);
        W(r) = toc(t0);
        TR{r} = struct('t', os.t(:), 'F', os.Force(:), 'SL', os.SL(:), 'LSE', os.LSE(:));
        FM{r} = fm;
        FDc{r} = fd;
        keep = cellfun(@(e) ~contains(e, ',') && isfield(fm, regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '')), fnS);
        CT(r) = sum(evalFeatureCost(fd, fm, fnS(keep)));
    end
    O = struct('name', CFG{c,1}, 'cfg', ov, 'TR', {TR}, 'FM', {FM}, 'cost', CT, 'wall', W);
    O.FD = FDc{1};
    M.(sprintf('c%d', c)) = O;
    save(outFile, 'M', 'CFG', 'vt', 'PARAMSET', '-v7.3');
    fprintf('%-44s slack cost base %.4f | SD %.4f (%.1f %%) | %.0f s/run\n', CFG{c,1}, CT(1), std(CT), 100*std(CT)/mean(CT), mean(W));
end
fprintf('DONE\n');
