%% RunStairsCheck.m -- held-out test: does the yield gate help the staircase?
%
% The staircase is not in any objective here (RunStairs = 0), and every fit so
% far overshoots each 1 % step (Analyses/KtrSegmentSlope sec. 8) -- a re-stretch
% protocol the gate acts on but was never tuned against. For each snapshot in
% SNAPS: E_stairs (MSE*10) and the step features, model vs data, one instrument
% (Analyses/StairsPrediction/extractStairsAttributes).
if ~exist('SNAPS', 'var'), SNAPS = {'refit9_opt'}; end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
d = load('data/protocol_03_27_2026_velocitytable_stairs.mat'); dt = d.datatable; vt = d.velocitytable;
fD = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
fprintf('%-14s %8s | %-34s | %-34s | %7s\n', 'snapshot', 'E_stairs', 'st_peak (step amplitude, F0)', 'st_nadir', 'platEnd');
fprintf('%-14s %8s | %-34s | %-34s | %7.3f\n', 'DATA', '', mat2str(round(fD.st_peak, 3)), mat2str(round(fD.st_nadir, 3)), fD.st_platEnd);
ST = struct();
for k = 1:numel(SNAPS)
    p = getParams(loadParams(SNAPS{k}), [], true, false);
    p.RunStairs = 1; p.PlotEachSeparately = false; p.ghostSave = ''; p.ghostLoad = ''; p.MaxRunTime = 900;
    [E, out] = runStairsExperiment(p);
    i_move = find(vt(:,2) ~= 0, 1);
    Fn = out.Force(:) / interp1(out.t, out.Force, vt(i_move,1), 'linear', 'extrap');
    fM = extractStairsAttributes(out.t(:), Fn, vt, 'st_', struct(), false);
    fprintf('%-14s %8.4f | %-34s | %-34s | %7.3f\n', SNAPS{k}, E, mat2str(round(fM.st_peak, 3)), mat2str(round(fM.st_nadir, 3)), fM.st_platEnd);
    ST(k).snap = SNAPS{k}; ST(k).E = E; ST(k).fM = fM; ST(k).t = out.t(:); ST(k).Fn = Fn;
end
save(fullfile(here, 'results', 'stairs_check.mat'), 'ST', 'fD', 'SNAPS');
