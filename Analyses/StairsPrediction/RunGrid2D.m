% RunGrid2D.m — the fair test of "estiff > 1, repaired by other params".
%
% The staircase's dominant defect at the incumbent fit is an elastic peak 2x
% too large. estiff > 1 makes that peak LARGER still, so the hypothesis can
% only survive if some other parameter takes the peak back down. Series
% compliance is the parameter that does that (lever scan (3) in RunLevers:
% kSE x0.2 drops st_peak from 1.25 to 0.12). So sweep the two together and see
% whether any (estiff > 1, kSE) pair reaches a staircase cost the sublinear
% branch cannot, and at what price in the paid-for objective.
%
% estiff is compensated at every node (kstiff2/kstiff2_n rescaled to hold the
% isometric active force) exactly as in RunLevers.
%
% Outputs -> results/grid2d.mat

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
tMov = vt(find(vt(:,2)~=0,1),1);

EV = [0.85 0.965322 1.25 1.50];
KV = [0.35 0.60 1.00 1.80];

nE = numel(EV); nK = numel(KV);
G = struct('estiff',[],'kmult',[],'Elsq',[],'Cfeat',[],'parts',[],'Et',[],'featCost',[],'Fiso',[]);
t0 = tic;
for a = 1:nE
    for b = 1:nK
        p = base;
        p.estiff  = EV(a);
        sc = (ref.Fiso - ref.Fp1) / (base.kstiff2 * ref.pw(EV(a)));
        p.kstiff2 = base.kstiff2*sc; p.kstiff2_n = base.kstiff2_n*sc;
        p.kSE     = base.kSE * KV(b);

        g = struct('estiff',EV(a),'kmult',KV(b),'Elsq',NaN,'Cfeat',NaN, ...
                   'parts',[],'Et',NaN,'featCost',[],'Fiso',NaN);
        ps = p; ps.RunStairs=1; ps.PlotEachSeparately=false; ps.ghostSave=''; ps.ghostLoad='';
        try
            [Elsq, os] = runStairsExperiment(ps);
            Fm = interp1(os.t, os.Force, tMov, 'linear','extrap');
            [Cf, parts] = stairsCost(os.t, os.Force/Fm, vt, fD);
            g.Elsq=Elsq; g.Cfeat=Cf; g.parts=parts; g.Fiso=Fm;
        catch ME
            fprintf('  stairs FAIL e=%.3f k=%.2f: %s\n', EV(a), KV(b), ME.message);
        end
        pf = p; pf.RunStairs = 0;
        try
            [Et,~,fc] = evaluateBakersExp(pf.g, pf, true, false);
            g.Et=Et; g.featCost=fc;
        catch ME
            fprintf('  obj FAIL e=%.3f k=%.2f: %s\n', EV(a), KV(b), ME.message);
        end
        G(a,b) = g;
        fprintf('e=%5.3f kSEx%4.2f | peak=%7.3f nadir=%8.3f rec=%7.3f | stairs=%8.3f LSQ=%.4f | Et=%9.3f\n', ...
            EV(a), KV(b), iget(g.parts,'st_peak'), iget(g.parts,'st_nadir'), iget(g.parts,'st_rec'), ...
            g.Cfeat, g.Elsq, g.Et);
    end
end
fprintf('total %.1f min\n', toc(t0)/60);
save(fullfile(resdir,'grid2d.mat'), 'G','EV','KV','fD','base','SNAP','ref');
fprintf('saved -> %s\n', fullfile(resdir,'grid2d.mat'));

function v = iget(s, f)
if isstruct(s) && isfield(s,f); v = s.(f); else; v = NaN; end
end
