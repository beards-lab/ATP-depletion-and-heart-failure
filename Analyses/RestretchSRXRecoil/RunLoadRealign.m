% RunLoadRealign.m — the two mechanisms §6 asked for.
%
% M1  COMPLIANT REALIGNMENT / load-engaged availability deficit (UseLoadRealign)
%     Bound bridges set the register of neighbouring actin sites through filament
%     compliance. Strip bound mass while the filament stays LOADED and the
%     register is left disturbed; re-equilibration takes tau_cr. One scalar state:
%         dA_def/dt = k_cr * max(0,-d(bound)/dt) * F/(F+F_cr)  -  A_def/tau_cr
%         ka_eff   *= (1 - A_def)
%     Driven by the NET change of bound mass, so it is identically zero at any
%     steady state — it cannot rescale ka, only act on transients. The load
%     weighting is what separates the windows:
%         post-restretch  p2 0.19->0.05 while F holds at 0.84 F_iso -> big deficit
%         post-slack/ktr  p2 collapses but F -> 0                   -> no deficit
%     This is the class §4 argued for: it acts UPSTREAM of attachment (so it
%     escapes the reservoir argument that neuters every return-path mechanism),
%     it is direction-agnostic, and it never touches R2(s).
%
% M2  FIXED-DWELL outflow (SRDOutflowK / D0OutflowK)
%     A first-order rate gives an exponential residence time — a head that just
%     arrived is as likely to leave as one that has waited 50 ms. Saturating the
%     outflow in occupancy makes it a CONSTANT FLUX once the pool is stocked, so
%     the pool drains LINEARLY and a slug of size m takes m/J to clear: a dwell
%     time, not a rate. It also makes a BIGGER slug take PROPORTIONALLY LONGER,
%     which is the sign the amplitude dependence wants.
%       SRDOutflowK  on the shared P_SRD pool  (the literal "SRX dwell")
%       D0OutflowK   on D0, which is fed ONLY by the rupture recoil, so the dwell
%                    is imposed on exactly the stretch-torn heads and nothing else
%
% Run Baseline.m first.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

L = load(fullfile(here,'results','baseline.mat'));
amp = L.amp; rsK_d = L.rsK_d; B = L.B;
pBase = getParams(loadParams(L.SNAP), [], true, false);
% 150 s, not 600: a pathological setting here (tiny dwell K, huge gain) can make
% ode15s crawl, and 15 variants x 600 s is a 2.5 h wall. Anything that needs more
% than 150 s is not a usable operating point anyway.
pBase.MaxRunTime = 150;
pBase.RunForceVelocity = 1; pBase.RunSlack = 1; pBase.RunSlackPassive = 1;
pBase.RunKtr = 0; pBase.RunStairs = 0; pBase.RunForceVelocityTime = 0;

cr = @(k, tau, F) struct('UseLoadRealign', true, 'k_cr', k, 'tau_cr', tau, 'F_cr', F);
% recoil into D0 (the selective private pool established in RunDestControl.m)
rec = @(varargin) struct('SRXFromR2HighStrain',1,'sSRXrip',0.015,'dSRXrip',0.005, ...
                         'SRXRecoilToD0',true,'UseD0State',true,'k_D0T',40,'K_D0T',0.15, ...
                         varargin{:});

V = {};
% --- M1: gain scan at the default timescale
V{end+1} = {'M1 k=1  tau=20 F=25', cr(1,   0.020, 25)};
V{end+1} = {'M1 k=3  tau=20 F=25', cr(3,   0.020, 25)};
V{end+1} = {'M1 k=6  tau=20 F=25', cr(6,   0.020, 25)};
V{end+1} = {'M1 k=12 tau=20 F=25', cr(12,  0.020, 25)};
% --- M1: timescale scan at a working gain
V{end+1} = {'M1 k=6  tau=10 F=25', cr(6,   0.010, 25)};
V{end+1} = {'M1 k=6  tau=40 F=25', cr(6,   0.040, 25)};
% --- M1: is the LOAD weighting doing the work? F_cr->0 removes it (deficit
%         becomes load-blind, so post-slack/ktr should now be hit too)
V{end+1} = {'M1 k=6  tau=20 LOADBLIND', cr(6, 0.020, 1e-6)};
V{end+1} = {'M1 k=6  tau=20 F=60', cr(6,   0.020, 60)};
% --- M2a: fixed dwell on the SHARED SRXD pool
V{end+1} = {'M2a SRD dwell K=.30', struct('SRDOutflowK', 0.30)};
V{end+1} = {'M2a SRD dwell K=.10', struct('SRDOutflowK', 0.10)};
V{end+1} = {'M2a SRD dwell K=.03', struct('SRDOutflowK', 0.03)};
% --- M2b: fixed dwell on D0, which only the rupture recoil fills
V{end+1} = {'M2b recoil+D0 dwell .30', rec('D0OutflowK', 0.30)};
V{end+1} = {'M2b recoil+D0 dwell .10', rec('D0OutflowK', 0.10)};
V{end+1} = {'M2b recoil+D0 dwell .03', rec('D0OutflowK', 0.03)};
V{end+1} = {'M2b recoil+D0 dwell .01', rec('D0OutflowK', 0.01)};

% Chunked: set SEL before running to do a subset (results accumulate in the .mat).
% The whole set is ~15 evaluations at 20-150 s each.
matf = fullfile(here,'results','loadrealign.mat');
if exist(matf, 'file'); S = load(matf); C = S.C; else; C = cell(1, numel(V)); end
if numel(C) < numel(V); C(numel(C)+1:numel(V)) = {[]}; end
if ~exist('SEL', 'var') || isempty(SEL); SEL = 1:numel(V); end

for i = SEL(:)'
    fprintf('[%2d/%2d] %-24s ', i, numel(V), V{i}{1});
    C{i} = srxRecoilCase(pBase, V{i}{1}, V{i}{2}, amp, rsK_d);
    if C{i}.ok
        fprintf('L2 %8.3f  rsK x%.2f  ktr x%.2f  F %.1f  (%.0f s)\n', ...
            C{i}.L2, C{i}.rsK_x, C{i}.ktr_x, C{i}.steady, C{i}.time);
    else
        fprintf('FAILED: %s\n', C{i}.err);
    end
    save(matf, 'C', 'V', 'amp', 'rsK_d');
end
keepC = ~cellfun(@isempty, C);

%% ---- summary -------------------------------------------------------------
fprintf('\n%-24s %9s %8s %8s %8s %9s\n', 'variant', 'L2', 'rsK x', 'ktr x', 'steady', 'slope');
fprintf('%-24s %9.3f %8.2f %8.2f %8.1f %+9.0f\n', 'baseline', B.L2, B.rsK_x, B.ktr_x, ...
        B.steady, B.slope);
for i = find(keepC)
    if isempty(C{i}) || ~C{i}.ok; continue; end
    fprintf('%-24s %9.3f %8.2f %8.2f %8.1f %+9.0f\n', C{i}.name, C{i}.L2, C{i}.rsK_x, ...
            C{i}.ktr_x, C{i}.steady, C{i}.slope);
end
fprintf('%-24s %9s %8.2f %8.2f %8s %+9.0f\n', 'DATA', '-', 1.00, 1.00, '~80', ...
        polyfit(amp(:), rsK_d(:), 1)*[1;0]);

%% ---- per-feature movers for anything that improves rsK -------------------
nm = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), B.fn, 'uni', 0);
for i = find(keepC)
    if isempty(C{i}) || ~C{i}.ok || C{i}.rsK_x >= B.rsK_x; continue; end
    d = C{i}.ct2 - B.ct2;
    [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s   (L2 %.3f -> %.3f)   rsK/cycle %s\n', C{i}.name, B.L2, C{i}.L2, ...
            mat2str(round(C{i}.rsK_m,1)));
    for j = o(1:min(7,numel(o)))
        if abs(d(j)) < 1e-3; break; end
        fprintf('    %-20s %7.3f -> %8.3f  %+8.3f\n', nm{j}, B.ct2(j), C{i}.ct2(j), d(j));
    end
end
fprintf('\ndata rsK/cycle %s\n', mat2str(round(rsK_d,1)));
