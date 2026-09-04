% RunFVFineAllRungs.m -- fine force-velocity sweep at BOTH [ATP], once per rung
% of the ATP hypothesis ladder, exported for the slide figures in the OTHER
% worktree.
%
% Same split of labour as RunBothAtpSlideRuns.m: the ATP power-gate machinery
% (pAtp_*, MgATP_ref) exists only in THIS worktree, so the sweep runs here and a
% plain .mat crosses to Further-optim, where the drawing lives.
%
% WHY A SEPARATE SWEEP AT ALL. The fit evaluates FV at five velocities
% (FV_velocities = [0 -0.5 -1 -2 -4]) because every extra velocity is another
% ODE solve inside every objective evaluation. That is the right economy for
% fitting and the wrong resolution for a figure: four power points cannot show
% the shape of a curve whose peak sits between two of them. Nothing in
% runFVExperiment is tied to the fit grid -- it loops over whatever
% FV_velocities holds and picks up data only where `ismember` finds a match --
% so a finer grid needs no model change. What a fine grid DOES break is the
% feature cost, where model and data vectors must be element-comparable, which
% is why this runs the protocol for PLOTTING and never calls evalFeatureCost.
%
% THE 8 mM SWEEP IS RUN ONCE. Every ATP gate is exactly 1 at MgATP_ref = 8, so
% the 8 mM curve is the ladder's frozen core and is identical for all four
% rungs. Running it per rung would burn four times the ODE budget to produce
% four copies of one curve.
%
% Output -> <figroot>/<rung>/fv_fine.mat, one per rung, each carrying the shared
% 8 mM curve plus that rung's 2 mM curve and the two data sets.
%
% Run:  matlab -batch "run('Analyses/AtpEffectHypotheses/RunFVFineAllRungs.m')"

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

if ~exist('RUNGS','var');  RUNGS  = {'H1','H2','H3','H4'}; end
if ~exist('VMAX','var');   VMAX   = 7;    end   % simulated range
if ~exist('DV','var');     DV     = 0.25; end   % simulated resolution
if ~exist('NWORK','var');  NWORK  = 4;    end
CORE = 'spreadw_iter_018';

FIGROOT = fullfile(root,'..','ATP-depletion-and-heart-failure','Analyses', ...
                   'OptimResultReport','results','bothatp');

vFine = 0:DV:VMAX;
LoadData;
if isempty(gcp('nocreate')); parpool('local', NWORK); end

%% ---- data: the two source curves, at the eight velocities they have ------
% The two FV sets come from different preparations, so their absolute forces
% are not comparable (isometric 56.4 vs 67.4 kPa at 8 mM) but their shapes are.
% Each is normalised by ITS OWN isometric before power is formed. Both sets
% carry a 2 mM column -- Data_ATP columns follow ATP_c = [8 4 2], hence +1.
vD = [0, 0.5, 1, 2, 3, 4, 5, 6];
COL8 = 2; COL2 = 4;
D8 = [Data_ATP(1:8,COL8)'; Data_ATP_AB2021(1:8,COL8)'];
D2 = [Data_ATP(1:8,COL2)'; Data_ATP_AB2021(1:8,COL2)'];
Dn8 = D8 ./ D8(:,1);  Dn2 = D2 ./ D2(:,1);
Dp8 = Dn8 .* vD;      Dp2 = Dn2 .* vD;
mu8 = mean(Dn8,1); sd8 = std(Dn8,0,1); pm8 = mean(Dp8,1); ps8 = std(Dp8,0,1);
mu2 = mean(Dn2,1); sd2 = std(Dn2,0,1); pm2 = mean(Dp2,1); ps2 = std(Dp2,0,1);

%% ---- 8 mM: once, from the core (all gates are 1 there anyway) ------------
fprintf('8 mM sweep, %d velocities (0 : %.2f : %g ML/s)\n', numel(vFine), DV, VMAX);
q8 = iBase(CORE, vFine);
q8.MgATP = 8; q8.MgADP = 0; q8.Pi = 0;
t0 = tic; [~, o8] = runFVExperiment(q8, ATP_c, Data_ATP);
F8 = arrayfun(@(o) o.Force(end), o8);
fprintf('  done in %.0f s\n', toc(t0));
if ~isfinite(F8(1)) || F8(1) <= 0
    error('RunFVFineAllRungs:badIsometric','8 mM isometric is %.4g', F8(1));
end
F8iso = F8(1); Fn8 = F8/F8iso; Pw8 = Fn8 .* vFine;
fprintf('  F_iso %.2f kPa | peak power %.4f at %.2f ML/s\n', ...
        F8iso, max(Pw8), vFine(find(Pw8==max(Pw8),1)));

%% ---- 2 mM: once per rung -------------------------------------------------
for r = 1:numel(RUNGS)
    rung = RUNGS{r}; SNAP = sprintf('atphyp_%s_opt', rung);
    fprintf('\n================ %s (%s) ================\n', rung, SNAP);
    if ~isfile(fullfile('params',[SNAP '.m']))
        fprintf(2,'  MISSING params/%s.m -- skipped\n', SNAP); continue;
    end
    q2 = iBase(SNAP, vFine);
    q2.MgATP = 2; q2.MgADP = 0; q2.Pi = 0.8;
    nz = fieldnames(q2); nz = nz(startsWith(nz,'pAtp_'));
    nz = nz(cellfun(@(f) abs(q2.(f)) > 1e-9, nz));
    fprintf('  active gates: %s\n', strjoin(strrep(nz,'pAtp_',''), ', '));

    t0 = tic; [~, o2] = runFVExperiment(q2, ATP_c, Data_ATP);
    F2 = arrayfun(@(o) o.Force(end), o2);
    fprintf('  done in %.0f s\n', toc(t0));
    if ~isfinite(F2(1)) || F2(1) <= 0
        fprintf(2,'  2 mM isometric is %.4g -- skipped\n', F2(1)); continue;
    end
    F2iso = F2(1); Fn2 = F2/F2iso; Pw2 = Fn2 .* vFine;
    fprintf('  F_iso %.2f kPa (2/8 = %.3f) | peak power %.4f at %.2f ML/s\n', ...
            F2iso, F2iso/F8iso, max(Pw2), vFine(find(Pw2==max(Pw2),1)));

    S = struct('rung',rung,'snap',SNAP,'core',CORE,'DV',DV,'VMAX',VMAX, ...
        'vFine',vFine,'F8',F8,'Fn8',Fn8,'Pw8',Pw8,'F8iso',F8iso, ...
        'F2',F2,'Fn2',Fn2,'Pw2',Pw2,'F2iso',F2iso, ...
        'vD',vD,'mu8',mu8,'sd8',sd8,'pm8',pm8,'ps8',ps8, ...
        'mu2',mu2,'sd2',sd2,'pm2',pm2,'ps2',ps2, ...
        'built',datestr(now,'yyyy-mm-dd HH:MM:SS'));            %#ok<TNOW1,DATST>
    S.gates = struct(); for i = 1:numel(nz); S.gates.(nz{i}) = q2.(nz{i}); end

    OUTS = { fullfile(resDir, sprintf('fv_fine_%s.mat', rung))
             fullfile(FIGROOT, rung, 'fv_fine.mat') };
    for i = 1:numel(OUTS)
        d = fileparts(OUTS{i});
        if ~isfolder(d)
            if contains(d,'ATP-depletion-and-heart-failure') && isfolder(fileparts(d))
                mkdir(d);
            else
                fprintf(2,'  skipped (no folder): %s\n', OUTS{i}); continue;
            end
        end
        save(OUTS{i},'-struct','S'); fprintf('  wrote %s\n', OUTS{i});
    end
end
fprintf('\nDONE\n');

%% ======================= helpers ========================================
% Local functions must be at the END of a script file, not inline where they
% are first used -- MATLAB rejects the latter with a parse error.
function p = iBase(snap, vFine)
% Everything the FV sweep needs and nothing else: gates live, other protocols
% off, and the isometric init reused across velocities (OptimizeFVInit).
    p = getParams(loadParams(snap), [], true, false);
    p.mods = {}; p.g = ones(1,0);
    p.FV_velocities = -vFine;          % model convention: negative = shortening
    p.EvalFeatures = true; p.OptimizeFVInit = true;
    p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0; p.PlotFullscreen = 0;
    p.MaxRunTime = 900; p.UseAtpJointFit = false;
    p.UseAtpPowerGates = true; p.MgATP_ref = 8;
    p.RunSlack = false; p.RunStairs = false; p.RunKtr = false;
    p.RunForceVelocity = true;
end
