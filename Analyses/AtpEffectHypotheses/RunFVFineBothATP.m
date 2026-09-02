% RunFVFineBothATP.m -- force-velocity and power-velocity on a FINE grid, for
% BOTH [ATP], in the style of OptimResultReport/RunFVFine.m in the Further-optim
% worktree, but with the 8 mM and 2 mM curves on the same axes.
%
% WHY A SEPARATE SWEEP. The fit evaluates FV at five velocities only
% (FV_velocities = [0 -0.5 -1 -2 -4]) because every extra velocity is another ODE
% solve inside every objective evaluation. That is right for fitting and far too
% coarse for a figure: four power points cannot show the shape of a curve whose
% peak sits between two of them. Nothing in runFVExperiment is tied to the fit
% grid, so a finer grid needs no model change -- but it DOES break the feature
% cost (model and data vectors must be element-comparable), which is why this
% runs the protocol for PLOTTING only and never calls evalFeatureCost.
%
% THE SD BAND is the spread between the two independent source curves --
% Baker2022 and IsovelocityForFilip2021 -- each normalised by ITS OWN isometric
% force before averaging, because absolute force is not comparable between
% preparations (56.4 vs 67.4 kPa at 8 mM). Both sets carry a 2 mM column
% (ATP_c = [8 4 2]), so the same construction works at both concentrations.
%
% MODELS. 8 mM = the core itself (gates are exactly 1 at MgATP_ref, so the core
% IS the 8 mM model). 2 mM = the best-scoring rung in the hypothesis ladder's
% saved state, chosen automatically and named in the figure.
%
% Cached in results/fv_fine_bothatp.mat; FORCE_RUN=true re-simulates.
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunFVFineBothATP.m')"

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
figDir = fullfile(resDir,'figs_fvfine'); if ~isfolder(figDir); mkdir(figDir); end

if ~exist('CORE','var');      CORE      = 'spreadw_iter_018'; end
if ~exist('VMAX','var');      VMAX      = 7;    end
if ~exist('DV','var');        DV        = 0.25; end
if ~exist('VPLOT_MAX','var'); VPLOT_MAX = 6;    end
% 0.25 ML/s is a protocol artifact, not a model prediction: runFVExperiment takes
% the isometric reference at SL 2.0 but starts every shortening run at SL 2.2, so
% at the slowest speed the 0.4 s ramp spends most of its time at longer, stronger
% lengths and the point lands ABOVE isometric. Kept in the .mat, out of the figure.
if ~exist('VDROP','var');     VDROP     = 0.25; end
if ~exist('NWORKERS','var');  NWORKERS  = 4;    end

MATF = fullfile(resDir,'fv_fine_bothatp.mat');

%% ---- pick the best low-ATP rung from the ladder state ---------------------
ST = fullfile(resDir,'atphyp_state.mat');
if ~isfile(ST)
    error('RunFVFineBothATP:noState', ...
          'no ladder state at %s -- run RunAtpHypothesisOptim.m first', ST);
end
S = load(ST);
[~, ib] = min([S.OUT.c]);
BEST = S.OUT(ib);
fprintf('low-ATP model: %s (%d gates, 2 mM cost %.4f) on core %s\n', ...
        BEST.rung, BEST.dof, BEST.c, S.CORE);
gateNames = fieldnames(BEST.exps);
nzg = gateNames(cellfun(@(f) abs(BEST.exps.(f)) > 1e-9, gateNames));
fprintf('  non-zero gates: %s\n', strjoin(strrep(nzg,'pAtp_',''), ', '));

%% ---- fine sweep, both conditions -----------------------------------------
vFine = 0:DV:VMAX;
haveCache = false;
if isfile(MATF) && ~(exist('FORCE_RUN','var') && FORCE_RUN)
    C = load(MATF);
    haveCache = isfield(C,'DV') && C.DV==DV && isfield(C,'VMAX') && C.VMAX==VMAX ...
             && isfield(C,'rung') && strcmp(C.rung, BEST.rung) && strcmp(C.CORE, S.CORE);
    if haveCache
        fprintf('cached sweep matches -- reusing, no model evaluated\n');
        vFine=C.vFine; Fn8=C.Fn8; Pw8=C.Pw8; Fn2=C.Fn2; Pw2=C.Pw2;
        F8iso=C.F8iso; F2iso=C.F2iso;
    end
end

if ~haveCache
    LoadData;
    if isempty(gcp('nocreate')); parpool('local', NWORKERS); end
    base = getParams(loadParams(S.CORE), [], true, false);
    base.mods = {}; base.g = ones(1,0);
    base.FV_velocities = -vFine;
    base.EvalFeatures = true; base.OptimizeFVInit = true;
    base.PlotEachSeparately = 0; base.PlotFeatureFitting = 0;
    base.MaxRunTime = 900; base.UseAtpJointFit = false;
    base.UseAtpPowerGates = true; base.MgATP_ref = 8;
    for i = 1:numel(gateNames); base.(gateNames{i}) = 0; end   % 8 mM: gates inert anyway

    fprintf('\nFV sweep, %d velocities (0 : %.2f : %g ML/s)\n', numel(vFine), DV, VMAX);

    q8 = base; q8.MgATP = 8; q8.MgADP = 0; q8.Pi = 0;
    t0 = tic; [~, o8] = runFVExperiment(q8, ATP_c, Data_ATP);
    F8 = arrayfun(@(o) o.Force(end), o8);
    fprintf('  8 mM done in %.0f s\n', toc(t0));

    q2 = base; q2.MgATP = 2; q2.MgADP = 0; q2.Pi = 0.8;
    for i = 1:numel(gateNames); q2.(gateNames{i}) = BEST.exps.(gateNames{i}); end
    t0 = tic; [~, o2] = runFVExperiment(q2, ATP_c, Data_ATP);
    F2 = arrayfun(@(o) o.Force(end), o2);
    fprintf('  2 mM done in %.0f s\n', toc(t0));

    for nm = {'F8','F2'}
        v = eval(nm{1});
        if ~isfinite(v(1)) || v(1) <= 0
            error('RunFVFineBothATP:badIsometric','%s isometric is %.4g', nm{1}, v(1));
        end
    end
    F8iso = F8(1); F2iso = F2(1);
    Fn8 = F8/F8iso; Pw8 = Fn8 .* vFine;
    Fn2 = F2/F2iso; Pw2 = Fn2 .* vFine;
    rung = BEST.rung; CORE = S.CORE;
    save(MATF,'vFine','Fn8','Pw8','Fn2','Pw2','F8iso','F2iso','DV','VMAX','rung','CORE');
end

fprintf('\nF_iso: 8 mM %.2f kPa | 2 mM %.2f kPa  (ratio %.3f)\n', F8iso, F2iso, F2iso/F8iso);
fprintf('peak power: 8 mM %.4f at %.2f | 2 mM %.4f at %.2f ML/s\n', ...
        max(Pw8), vFine(find(Pw8==max(Pw8),1)), max(Pw2), vFine(find(Pw2==max(Pw2),1)));

%% ---- data, both conditions, from the two source sets ---------------------
LoadData;
vD = [0, 0.5, 1, 2, 3, 4, 5, 6];
COL8 = 2;    % Data_ATP columns are ATP_c = [8 4 2] -> +1
COL2 = 4;
D8 = [Data_ATP(1:8,COL8)'; Data_ATP_AB2021(1:8,COL8)'];
D2 = [Data_ATP(1:8,COL2)'; Data_ATP_AB2021(1:8,COL2)'];
Dn8 = D8 ./ D8(:,1);  Dn2 = D2 ./ D2(:,1);       % each by its OWN isometric
Dp8 = Dn8 .* vD;      Dp2 = Dn2 .* vD;
mu8 = mean(Dn8,1); sd8 = std(Dn8,0,1); pm8 = mean(Dp8,1); ps8 = std(Dp8,0,1);
mu2 = mean(Dn2,1); sd2 = std(Dn2,0,1); pm2 = mean(Dp2,1); ps2 = std(Dp2,0,1);

%% ---- figure ---------------------------------------------------------------
C8 = [0.00 0.35 0.70]; C2 = [0.80 0.20 0.10];
B8 = [0.62 0.74 0.88]; B2 = [0.92 0.72 0.66];
K8 = [0.06 0.16 0.34]; K2 = [0.42 0.09 0.05];

mM = vFine <= VPLOT_MAX + 1e-9;
for vv = VDROP(:)'; mM = mM & abs(vFine - vv) > 1e-9; end

% Vmax: where each model curve reaches zero force. This is a real prediction --
% low ATP lowers Vmax -- but past it F/F0 goes very slightly negative (~-1e-4),
% which drew a cliff to below zero and a spurious kink in the power panel. Clamp
% at zero for drawing and mark Vmax explicitly instead of hiding the crossing.
Fn8p = max(Fn8, 0); Pw8p = Fn8p .* vFine;
Fn2p = max(Fn2, 0); Pw2p = Fn2p .* vFine;
i8 = find(Fn8 <= 0, 1); v8max = NaN; if ~isempty(i8); v8max = vFine(i8); end
i2 = find(Fn2 <= 0, 1); v2max = NaN; if ~isempty(i2); v2max = vFine(i2); end
% End each curve AT its Vmax. Drawing it flat along zero afterwards implies the
% model still predicts something there, and the near-zero forces just before the
% crossing are numerically noisy (they produced a spurious uptick in power).
m8 = mM; if isfinite(v8max); m8 = m8 & vFine <= v8max + 1e-9; end
m2 = mM; if isfinite(v2max); m2 = m2 & vFine <= v2max + 1e-9; end
fprintf('Vmax (F/F0 -> 0): 8 mM %s | 2 mM %s ML/s\n', ...
        iNum(v8max, VMAX), iNum(v2max, VMAX));

close all force;
fig = figure('Units','pixels','Position',[20 20 1500 640],'Color','w');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

ax = nexttile; hold on; box on; grid on
iBand(vD, mu8, sd8, B8); iBand(vD, mu2, sd2, B2);
p1 = plot(vFine(m8), Fn8p(m8), '-', 'Color',C8,'LineWidth',3.0);
p2 = plot(vFine(m2), Fn2p(m2), '-', 'Color',C2,'LineWidth',3.0);
d1 = plot(vD, mu8, 'o','Color',K8,'MarkerFaceColor',K8,'MarkerSize',9);
d2 = plot(vD, mu2, 'o','Color',K2,'MarkerFaceColor',K2,'MarkerSize',9);
xlabel('shortening velocity (ML/s)','FontSize',15,'FontWeight','bold');
ylabel('F / F_{iso}','FontSize',15,'FontWeight','bold');
title('Force-velocity','FontSize',17,'FontWeight','bold');
if isfinite(v2max) && v2max <= VPLOT_MAX
    xline(v2max,':','Color',C2,'LineWidth',2.0,'HandleVisibility','off');
    text(v2max-0.06, 0.60, sprintf('V_{max} %.2f',v2max),'Color',C2, ...
         'FontSize',12,'FontWeight','bold','Rotation',90,'HorizontalAlignment','left');
end
ylim([0 1.02]);
xlim([0 VPLOT_MAX]); set(ax,'FontSize',14,'LineWidth',1.1);
legend([d1 p1 d2 p2], {'data 8 mM \pm SD','model 8 mM', ...
        'data 2 mM \pm SD',sprintf('model 2 mM (%s)',BEST.rung)}, ...
        'Location','northeast','FontSize',12);

ax = nexttile; hold on; box on; grid on
iBand(vD, pm8, ps8, B8); iBand(vD, pm2, ps2, B2);
plot(vFine(m8), Pw8p(m8), '-', 'Color',C8,'LineWidth',3.0);
plot(vFine(m2), Pw2p(m2), '-', 'Color',C2,'LineWidth',3.0);
plot(vD, pm8, 'o','Color',K8,'MarkerFaceColor',K8,'MarkerSize',9);
plot(vD, pm2, 'o','Color',K2,'MarkerFaceColor',K2,'MarkerSize',9);
xlabel('shortening velocity (ML/s)','FontSize',15,'FontWeight','bold');
ylabel('normalised power','FontSize',15,'FontWeight','bold');
title('Force-velocity (power)','FontSize',17,'FontWeight','bold');
xlim([0 VPLOT_MAX]); set(ax,'FontSize',14,'LineWidth',1.1);

title(tl, sprintf(['8 mM vs 2 mM, fine velocity grid.  core %s;  low-ATP %s, ' ...
        '%d gates, 2 mM cost %.2f.  F_{iso} %.0f \\rightarrow %.0f kPa (\\times%.2f)'], ...
        strrep(S.CORE,'_','\_'), BEST.rung, BEST.dof, BEST.c, F8iso, F2iso, F2iso/F8iso), ...
        'FontSize',14,'FontWeight','bold');
exportgraphics(fig, fullfile(figDir,'fig_fv_fine_bothatp.png'),'Resolution',170);
fprintf('\nsaved fig_fv_fine_bothatp.png\n');

%% ---- the numbers, for the record -----------------------------------------
fprintf('\n%6s | %-24s | %-24s\n','v','F/Fiso  (8 mM)','F/Fiso  (2 mM)');
fprintf('%6s | %8s %8s %6s | %8s %8s %6s\n','ML/s','data','model','SD','data','model','SD');
for i = 1:numel(vD)
    j = find(abs(vFine - vD(i)) < 1e-9, 1);
    if isempty(j); continue; end
    fprintf('%6.1f | %8.3f %8.3f %6.3f | %8.3f %8.3f %6.3f\n', ...
        vD(i), mu8(i), Fn8(j), sd8(i), mu2(i), Fn2(j), sd2(i));
end
fprintf('\nDONE -> %s\n', figDir);

%% ==========================================================================
function iBand(x, mu, sd, col)
g = ~isnan(mu) & ~isnan(sd);
fill([x(g) fliplr(x(g))], [mu(g)-sd(g) fliplr(mu(g)+sd(g))], col, ...
     'FaceAlpha',0.45,'EdgeColor','none','HandleVisibility','off');
end

function s = iNum(v, vmax)
% NB: written via the editor, not shell printf -- a bash printf here ate the
% MATLAB format specifiers (%g -> 0, %.2f -> 0.00) and produced a silently
% meaningless diagnostic while the underlying Vmax detection was correct.
if isnan(v)
    s = sprintf('>%g', vmax);      % never crossed zero within the swept range
else
    s = sprintf('%.2f', v);
end
end
