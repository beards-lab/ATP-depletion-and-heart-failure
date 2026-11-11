% RunCrossProtocol.m — two more out-of-sample tests, and the strain-offset lever.
%
% (A) CROSS-PROTOCOL. 04/03 and 04/10 ran a COARSE staircase (1x1% then 3x3%)
%     instead of 03/27's fine one (7x1% then 1x3%). Neither was ever in any
%     objective, so running the incumbent fit on them is a second and third
%     independent prediction test -- and, because the step sizes differ, it
%     probes the amplitude nonlinearity directly rather than the accumulated
%     offset.
%
% (B) dr LEVER. The elastic peak/base ratio of a power-law bridge is
%         peak/base ~ e * dx_eff / x0
%     so it can be brought down either by lowering the exponent e (estiff) or by
%     raising the resting strain x0. dr sets x0 (the isometric population sits at
%     x = s + dr ~ +9.5 nm). If dr does the same job as estiff, then the staircase
%     peak does not identify the exponent -- it identifies the STRAIN SENSITIVITY,
%     and estiff is only one of at least two ways to write that down.
%
% Outputs -> results/crossprotocol.mat

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end
if ~exist('SNAP','var'); SNAP = 'spreadw_opt'; end

base = getParams(loadParams(SNAP), [], true, false);
files = {'protocol_03_27_2026_velocitytable_stairs.mat', ...
         'protocol_04_03_2026_velocitytable_stairs.mat', ...
         'protocol_04_10_2026_velocitytable_stairs.mat'};
days = {'03/27 M (fine)','04/03 F (coarse)','04/10 M (coarse)'};

fprintf('=== (A) cross-protocol prediction, %s ===\n', SNAP);
CP = struct();
for i = 1:numel(files)
    q  = load(['data/' files{i}]); dt = q.datatable; vt = q.velocitytable;
    fDi = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
    p = base; p.RunStairs = 1; p.PlotEachSeparately = false;
    p.ghostSave=''; p.ghostLoad=''; p.stairs_velocitytableonfile = files{i};
    [Elsq, os] = runStairsExperiment(p);
    tMov = vt(find(vt(:,2)~=0,1),1);
    Fm = interp1(os.t, os.Force, tMov, 'linear','extrap');
    [Cf, parts, fMi] = stairsCost(os.t, os.Force/Fm, vt, fDi);
    CP(i).day=days{i}; CP(i).fD=fDi; CP(i).fM=fMi; CP(i).Elsq=Elsq; CP(i).Cfeat=Cf;
    CP(i).parts=parts; CP(i).t=os.t; CP(i).Fn=os.Force/Fm; CP(i).vt=vt; CP(i).dt=dt;

    s = fDi.st_dSL; i1 = find(abs(s-0.010)<0.003); i3 = find(abs(s-0.030)<0.006);
    rD = mean(fDi.st_peak(i3))/mean(fDi.st_peak(i1));
    rM = mean(fMi.st_peak(i3))/mean(fMi.st_peak(i1));
    fprintf('\n%s\n  LSQ=%.5f  featCost=%.3f  [peak %.3f nadir %.3f plat %.3f rec %.3f]\n', ...
        days{i}, Elsq, Cf, parts.st_peak, parts.st_nadir, parts.st_plat, parts.st_rec);
    fprintf('  peak 1%%: data %.4f  model %.4f  (model/data %.2f)\n', ...
        mean(fDi.st_peak(i1)), mean(fMi.st_peak(i1)), mean(fMi.st_peak(i1))/mean(fDi.st_peak(i1)));
    fprintf('  peak 3%%: data %.4f  model %.4f  (model/data %.2f)\n', ...
        mean(fDi.st_peak(i3)), mean(fMi.st_peak(i3)), mean(fMi.st_peak(i3))/mean(fDi.st_peak(i3)));
    fprintf('  AMPLITUDE RATIO peak(3%%)/peak(1%%):  data %.3f   model %.3f\n', rD, rM);
    fprintf('  envelope 2.0->2.2:                 data %.4f   model %.4f\n', ...
        fDi.st_plat(end), fMi.st_plat(end));
    CP(i).rD = rD; CP(i).rM = rM;
end

%% ---- (B) dr lever, on the 03/27 staircase -------------------------------
q  = load(['data/' files{1}]); dt = q.datatable; vt = q.velocitytable;
fD = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
fprintf('\n=== (B) dr (resting strain offset) lever, base dr = %.5g um ===\n', base.dr);
b2 = base; b2.stairs_velocitytableonfile = files{1};
S4 = scanLever(b2, fD, vt, struct('name','dr', ...
     'vals', base.dr*[0.6 0.8 1.0 1.4 2.0 2.8], 'comp','none'), true);

save(fullfile(resdir,'crossprotocol.mat'), 'CP','S4','base','SNAP','days');
fprintf('\nsaved -> %s\n', fullfile(resdir,'crossprotocol.mat'));
