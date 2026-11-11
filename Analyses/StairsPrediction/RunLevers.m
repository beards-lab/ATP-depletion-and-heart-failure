% RunLevers.m — three sweeps that together answer the estiff question.
%
%   (1) UNCOMPENSATED estiff, staircase only.  Robustness check: the staircase
%       cost is normalised by the model's own F0, so if (1) and the compensated
%       scan agree in DIRECTION the result is not an artefact of the kstiff2
%       rescale.
%   (2) FINE compensated estiff over [0.55 .. 1.00], with the full paid-for
%       objective, to locate the staircase optimum and price it.
%   (3) SERIES COMPLIANCE (kSE) at the incumbent estiff.  The staircase's worst
%       defect is a 2x too large elastic peak; end compliance is the textbook
%       alternative cure. If kSE fixes the peak as well as estiff does, then the
%       staircase does not identify the exponent at all.
%
% Outputs -> results/levers.mat

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end
if ~exist('SNAP','var'); SNAP = 'spreadw_opt'; end

base = getParams(loadParams(SNAP), [], true, false);
d  = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
dt = d.datatable; vt = d.velocitytable;
fD = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
ref = stairsIsoRef(base, vt);
fprintf('iso ref: F=%.4g kPa, p2 share %.1f%%\n', ref.Fiso, 100*(ref.Fiso-ref.Fp1)/ref.Fiso);

t0 = tic;
fprintf('\n(1) UNCOMPENSATED estiff, staircase only\n');
S1 = scanLever(base, fD, vt, struct('name','estiff', ...
     'vals',[0.5 0.7 0.85 0.965322 1.10 1.25 1.50 1.75 2.00], 'comp','none'), false);

fprintf('\n(2) COMPENSATED estiff, fine grid + full objective\n');
S2 = scanLever(base, fD, vt, struct('name','estiff', ...
     'vals',[0.55 0.60 0.65 0.70 0.75 0.80 0.90 1.00], 'comp','kstiff2','ref',ref), true);

fprintf('\n(3) SERIES COMPLIANCE kSE at estiff=%.4g (base kSE=%.5g)\n', base.estiff, base.kSE);
S3 = scanLever(base, fD, vt, struct('name','kSE', ...
     'vals',base.kSE*[0.2 0.35 0.5 0.75 1 1.5 2.5], 'comp','none'), true);

fprintf('\ntotal %.1f min\n', toc(t0)/60);
save(fullfile(resdir,'levers.mat'), 'S1','S2','S3','fD','base','SNAP','ref');
fprintf('saved -> %s\n', fullfile(resdir,'levers.mat'));
