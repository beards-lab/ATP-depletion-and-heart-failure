% RunEstiffAllDays.m — is the estiff verdict a property of the 03/27 staircase,
% or of the staircase experiment?
%
% 03/27 is a fine staircase (7x1% + 1x3%) and 04/03 / 04/10 are coarse ones
% (1x1% + 3x3%). The incumbent model reproduces 03/27's amplitude ratio
% peak(3%)/peak(1%) ~ 1.4 but INVERTS it on the coarse days (0.88 vs a measured
% 1.34-1.48), so the coarse days are the ones where a case for a stiffer
% exponent could plausibly be made. Run the same compensated estiff sweep on all
% three and see whether any of them prefers e > 1.
%
% Outputs -> results/estiff_alldays.mat

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end
if ~exist('SNAP','var'); SNAP = 'spreadw_opt'; end

base  = getParams(loadParams(SNAP), [], true, false);
files = {'protocol_03_27_2026_velocitytable_stairs.mat', ...
         'protocol_04_03_2026_velocitytable_stairs.mat', ...
         'protocol_04_10_2026_velocitytable_stairs.mat'};
days  = {'03/27 fine','04/03 coarse','04/10 coarse'};
EV    = [0.70 0.85 0.965322 1.15 1.40 1.75];

AD = struct();
for i = 1:numel(files)
    q = load(['data/' files{i}]); dt = q.datatable; vt = q.velocitytable;
    fD = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
    b  = base; b.stairs_velocitytableonfile = files{i};
    ref = stairsIsoRef(b, vt);
    fprintf('\n=== %s ===\n', days{i});
    S = scanLever(b, fD, vt, struct('name','estiff','vals',EV,'comp','kstiff2','ref',ref), false);
    % amplitude ratio at each exponent
    s = fD.st_dSL; i1 = find(abs(s-0.010)<0.003); i3 = find(abs(s-0.030)<0.006);
    rD = mean(fD.st_peak(i3))/mean(fD.st_peak(i1));
    fprintf('  amplitude ratio: data %.3f | model:', rD);
    for k = 1:numel(S)
        if isempty(S(k).fM); fprintf('   n/a'); continue; end
        fprintf('  %.2f', mean(S(k).fM.st_peak(i3))/mean(S(k).fM.st_peak(i1)));
    end
    fprintf('\n');
    AD(i).day=days{i}; AD(i).S=S; AD(i).fD=fD; AD(i).rD=rD; AD(i).EV=EV;
end

fprintf('\n%-14s', 'estiff'); fprintf('%9.3f', EV); fprintf('\n');
for i = 1:numel(AD)
    fprintf('%-14s', [AD(i).day ' feat']); fprintf('%9.3f',[AD(i).S.Cfeat]); fprintf('\n');
    fprintf('%-14s', [AD(i).day ' LSQ']);  fprintf('%9.4f',[AD(i).S.Elsq]);  fprintf('\n');
end
save(fullfile(resdir,'estiff_alldays.mat'), 'AD','EV','days','SNAP');
fprintf('\nsaved -> %s\n', fullfile(resdir,'estiff_alldays.mat'));
