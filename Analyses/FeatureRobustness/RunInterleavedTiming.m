%% RunInterleavedTiming.m -- wall-time RATIOS that survive a shared, fluctuating machine.
% The configurations are run in rotation (A B C F A B C F ...), NROUND rounds,
% one serial slack run each (AllPar semantics, no pool), so a load swing from
% another session hits every configuration alike. Reported: the median over
% rounds of each configuration's time divided by A's time in the same round.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
if ~exist('NROUND', 'var'), NROUND = 4; end
here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
addpath(genpath(root)); res = fullfile(here, 'results');
fx = @(varargin) struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, varargin{:});
CFG = { 'A old: 1e-2, step output, hard rectifier', struct('SmoothRectifierEps', 0, 'ClampStrainDepExtrapolation', 0, 'odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 0)
        'B SiLU 1e-3, step output',                 fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 0)
        'C SiLU 1e-3, dense in cycles only',        fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'F SiLU 1e-3, dense everywhere (campaign)', fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4) };
P = getParams(loadParams(PARAMSET), [], true, false);
P.RunSlackSegments = 'AllPar'; P.MaxRunTime = 3000; P.PlotEachSeparately = 0; P.PlotFeatureFitting = 0;
P.BreakOnODEUnstable = false; P.ghostSave = ''; P.ghostLoad = ''; P.SaveBest = 0; P.mods = {}; P.g = ones(1,0);
delete(gcp('nocreate'));
PC = cell(size(CFG,1), 1);
for c = 1:size(CFG,1)
    Pc = P; ov = CFG{c,2}; f = fieldnames(ov); for k = 1:numel(f), Pc.(f{k}) = ov.(f{k}); end
    PC{c} = getParams(Pc, Pc.g, true, false);
end
runSlackExperiment(PC{1});                                % warm-up, not timed
W = nan(NROUND, size(CFG,1));
for r = 1:NROUND
    for c = 1:size(CFG,1)
        t0 = tic; runSlackExperiment(PC{c}); W(r,c) = toc(t0);
    end
    fprintf('round %d: %s s\n', r, mat2str(round(W(r,:))));
    save(fullfile(res, ['timing_interleaved_' PARAMSET '.mat']), 'W', 'CFG');
end
Rt = W ./ W(:,1);
for c = 1:size(CFG,1)
    fprintf('%-44s median %6.1f s   ratio to A: median %.2f  (range %.2f-%.2f)\n', CFG{c,1}, median(W(:,c)), median(Rt(:,c)), min(Rt(:,c)), max(Rt(:,c)));
end
