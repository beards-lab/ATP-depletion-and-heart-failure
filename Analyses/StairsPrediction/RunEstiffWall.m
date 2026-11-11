% RunEstiffWall.m — why a 60-parameter optimiser that HAD estiff in its pool
% never moved it.
%
% estiff was in the spreadw pool (RunOptimSpreadW.m:76) and kstiff2 was
% COMPULSORY, i.e. the two were jointly drawable in every single round. Yet
% across 18 improved iterations estiff went 0.965542 -> 0.964708 -> 0.965322 and
% then never moved again: a total excursion of 0.09 %.
%
% This script measures the reason. At the isometric point the attached
% population is a spike at <x> ~ 0.0102 um, so
%       d(log F2)/d(estiff) = log(<x>) = -4.59 per unit e
% i.e. a 1 % change in estiff is a 4.4 % change in active force unless kstiff2
% follows it exactly. Two curves:
%
%   RAW  : estiff moved alone           -> how steep the wall is
%   RIDGE: estiff + iso-force kstiff2   -> how steep the FLOOR of the valley is
%
% The gap between them is the answer: any optimiser step that does not track
% the ridge to within a fraction of a percent is rejected long before it can
% discover that a different basin exists 0.29 units away.
%
% Usage: matlab -batch "run('Analyses/StairsPrediction/RunEstiffWall.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
if ~exist('NWORKERS','var'); NWORKERS = 3; end
if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

inc = getParams(loadParams('spreadw_opt'), [], true, false);
d   = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
ref = stairsIsoRef(inc, d.velocitytable);
pw0 = ref.pw(inc.estiff);
e0  = inc.estiff;

DE = [-0.05 -0.02 -0.01 -0.005 0 0.005 0.01 0.02 0.05];

base = inc; base.mods = {}; base.g = ones(1,0);
base.RunStairs = false; base.EvalFeatures = true; base.OptimizeOn = 'Feats';
base.PlotEachSeparately = 0; base.PlotFeatureFitting = 0;

raw = nan(size(DE)); rid = nan(size(DE)); Fer = nan(size(DE));
for i = 1:numel(DE)
    e = e0 + DE(i);

    p = base; p.estiff = e;                       % RAW: kstiff2 untouched
    raw(i) = evaluateBakersExp(1, getParams(p,[],true,false), true, false);

    sc = pw0 / ref.pw(e);                         % RIDGE: iso-force kstiff2
    q = base; q.estiff = e; q.kstiff2 = inc.kstiff2*sc; q.kstiff2_n = inc.kstiff2_n*sc;
    rid(i) = evaluateBakersExp(1, getParams(q,[],true,false), true, false);

    Fer(i) = 100*(1/sc - 1);   % % isometric active-force error the RAW point carries
    fprintf('de %+7.3f (e=%.4f)  Fiso err %+7.2f %%   RAW Et %10.4f   RIDGE Et %8.4f\n', ...
            DE(i), e, Fer(i), raw(i), rid(i));
end

% ---- the same wall expressed in the OTHER coordinate ---------------------
% How much does the objective cost for a pure kstiff2 error, no estiff change?
KS = [0.90 0.95 0.98 0.99 1 1.01 1.02 1.05 1.10];
ks = nan(size(KS));
for i = 1:numel(KS)
    p = base; p.kstiff2 = inc.kstiff2*KS(i); p.kstiff2_n = inc.kstiff2_n*KS(i);
    ks(i) = evaluateBakersExp(1, getParams(p,[],true,false), true, false);
    fprintf('kstiff2 x%.3f   Et %10.4f\n', KS(i), ks(i));
end

save(fullfile(here,'results','estiff_wall.mat'), 'DE','raw','rid','Fer','KS','ks','e0','pw0');

fprintf('\n=== WALL (estiff alone) vs RIDGE (estiff + iso kstiff2) ===\n');
fprintf('%8s %8s %10s %10s %10s\n','de','e','Fiso err%','RAW Et','RIDGE Et');
for i = 1:numel(DE)
    fprintf('%+8.3f %8.4f %10.2f %10.4f %10.4f\n', DE(i), e0+DE(i), Fer(i), raw(i), rid(i));
end
fprintf('\n=== kstiff2 alone ===\n');
for i = 1:numel(KS); fprintf('  x%.3f  Et %10.4f\n', KS(i), ks(i)); end
