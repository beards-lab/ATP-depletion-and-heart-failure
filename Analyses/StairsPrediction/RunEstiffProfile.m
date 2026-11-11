% RunEstiffProfile.m — the honest estiff test: a PROFILE, not a slice.
%
% WHY THIS EXISTS. The first sweeps (RunEstiffScan / RunLevers / RunGrid2D) held
% kstiff2 at a value computed to preserve the ISOMETRIC FORCE on the frozen
% isometric p2 distribution. That does hold F_iso to the digit — but F_iso is one
% number, and kstiff2 sets much more than that. estiff and kstiff2 are strongly
% coupled, so a slice that fixes kstiff2 by an iso-force rule measures "kstiff2 is
% no longer optimal here" on top of whatever estiff does. The evidence: one round
% of RunRepairFit at frozen estiff = 1.25 moved essentially kstiff2 ALONE
% (x3.98 where the iso-force rule had said x3.66) and took the objective from
% 5.478 to 1.173 and the staircase from 18.30 to 1.674. A 9 % error in the
% compensation was worth a factor of 11 in the staircase cost.
%
% So: for each estiff, RE-OPTIMISE kstiff2 against the paid-for objective, and
% only then score the staircase (still out of sample). That answers the question
% that matters -- among parameter sets that fit the PAID-FOR data equally well,
% which exponent predicts the staircase best? -- instead of the question a slice
% answers, which is "what happens if I mis-set kstiff2".
%
% kstiff2_n is carried at a fixed ratio to kstiff2 (the optimizer moved them
% together). Everything else is held.
%
% Usage: matlab -batch "run('Analyses/StairsPrediction/RunEstiffProfile.m')"
% Outputs -> results/estiff_profile.mat

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end

if ~exist('SNAP','var');   SNAP  = 'spreadw_opt'; end
if ~exist('EGRID','var');  EGRID = [0.85 0.965322 1.10 1.25 1.40 1.60]; end
if ~exist('NW','var');     NW    = 2; end

if isempty(gcp('nocreate')); parpool('local', NW); end

base = getParams(loadParams(SNAP), [], true, false);
d    = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
vt   = d.velocitytable; dt = d.datatable;
fD   = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);
ref  = stairsIsoRef(base, vt);
tMov = vt(find(vt(:,2)~=0,1),1);
rat_n = base.kstiff2_n / base.kstiff2;

fprintf('=== estiff PROFILE on %s (kstiff2 re-optimised at every node) ===\n', SNAP);
fprintf('incumbent: estiff %.4g  kstiff2 %.5g  Et 0.9786  staircase 5.415\n', ...
        base.estiff, base.kstiff2);

P = struct('estiff',num2cell(EGRID),'k2opt',[],'Et',[],'featCost',[], ...
           'Cfeat',[],'Elsq',[],'parts',[],'nEval',[]);

for i = 1:numel(EGRID)
    e   = EGRID(i);
    k0  = (ref.Fiso - ref.Fp1) / ref.pw(e);        % iso-force rule = the bracket centre
    nEv = 0;
    obj = @(lk) iObj(exp(lk), e, rat_n, base);
    % 1-D golden-section on log(kstiff2) over a factor-of-2.2 bracket either side.
    opts = optimset('TolX', 0.02, 'MaxFunEvals', 14, 'Display','off');
    [lkOpt, EtOpt] = fminbnd(obj, log(k0/2.2), log(k0*2.2), opts);
    kOpt = exp(lkOpt);

    p = base; p.estiff = e; p.kstiff2 = kOpt; p.kstiff2_n = kOpt*rat_n;
    [Et,~,fc] = evaluateBakersExp(p.g, p, true, false);

    ps = p; ps.RunStairs = 1; ps.PlotEachSeparately = false;
    ps.ghostSave=''; ps.ghostLoad='';
    [Elsq, os] = runStairsExperiment(ps);
    Fm = interp1(os.t, os.Force, tMov, 'linear','extrap');
    [Cf, parts] = stairsCost(os.t, os.Force/Fm, vt, fD);

    P(i).k2opt=kOpt; P(i).Et=Et; P(i).featCost=fc; P(i).Cfeat=Cf;
    P(i).Elsq=Elsq; P(i).parts=parts; P(i).Fiso=Fm; P(i).k0=k0;
    fprintf(['e=%5.3f | kstiff2 %8.5g (iso-rule said %8.5g, x%.3f) | Et=%7.4f | ' ...
             'stairs %7.3f LSQ %.5f [pk %.2f nd %.2f pl %.3f rc %.2f]\n'], ...
        e, kOpt, k0, kOpt/k0, Et, Cf, Elsq, parts.st_peak, parts.st_nadir, ...
        parts.st_plat, parts.st_rec);
end

save(fullfile(resdir,'estiff_profile.mat'), 'P','EGRID','base','fD','SNAP');
fprintf('\nsaved -> %s\n', fullfile(resdir,'estiff_profile.mat'));

function Et = iObj(k, e, rat_n, base)
p = base; p.estiff = e; p.kstiff2 = k; p.kstiff2_n = k*rat_n;
try
    Et = evaluateBakersExp(p.g, p, true, false);
catch
    Et = 1e4;
end
if ~isfinite(Et); Et = 1e4; end
end
