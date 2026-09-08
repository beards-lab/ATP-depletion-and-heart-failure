% MakeProfileFigure.m — the crux figure: a clean profile that is the wrong experiment.
%
% Left/middle: profiling estiff with kstiff2 re-optimised at every node puts BOTH the
% objective and the staircase minimum at the incumbent exponent. Right: the coupled
% refit — which also moved the FV strain-dependence shape terms — sits at e = 1.25 with
% a nearly-incumbent objective and a 3.2x better staircase, far below the profile curve.
% The gap between the star and the curve is the three-way coupling.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results');
Pf = load(fullfile(resdir,'estiff_profile.mat'));
P  = Pf.P; E = [P.estiff];

ET_INC = 0.9751; ST_INC = 5.415;      % incumbent, re-measured this session
ET_REP = 1.1631; ST_REP = 1.674;      % coupled refit at frozen e = 1.25

f = figure('Position',[60 60 1500 460],'Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on;
plot(E, [P.Et], 'o-','LineWidth',1.8,'Color',[0 0.35 0.75],'MarkerFaceColor',[0 0.35 0.75]);
plot(0.965322, ET_INC,'kp','MarkerSize',16,'MarkerFaceColor','y');
plot(1.25, ET_REP,'rp','MarkerSize',18,'MarkerFaceColor','r');
set(gca,'YScale','log'); xlabel('estiff'); ylabel('paid-for objective E_t');
legend({'profile (kstiff2 re-optimised)','incumbent','coupled refit'},'Location','north');
title('(a) objective');

nexttile; hold on; box on;
plot(E, [P.Cfeat], 'o-','LineWidth',1.8,'Color',[0.85 0.33 0.1],'MarkerFaceColor',[0.85 0.33 0.1]);
plot(0.965322, ST_INC,'kp','MarkerSize',16,'MarkerFaceColor','y');
plot(1.25, ST_REP,'rp','MarkerSize',18,'MarkerFaceColor','r');
set(gca,'YScale','log'); xlabel('estiff'); ylabel('staircase feature cost');
legend({'profile','incumbent','coupled refit'},'Location','north');
title('(b) staircase — the profile says one thing, the refit another');

% (c) the decomposition: which halves of the refit matter
nexttile; hold on; box on;
lbl = {'repair result','minus FV shape','minus estiff+kstiff2','minus k2,v_{ref}','incumbent'};
et  = [1.163 6.420 39.10 1.280 ET_INC];
st  = [1.674 16.24  3.55 5.314 ST_INC];
scatter(st, et, 90, [0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19; 0.49 0.18 0.56; 0 0 0], 'filled');
for i = 1:numel(lbl)
    text(st(i)*1.08, et(i), lbl{i}, 'FontSize',9);
end
set(gca,'XScale','log','YScale','log'); xlim([1 60]); ylim([0.5 80]);
xlabel('staircase feature cost'); ylabel('paid-for objective E_t');
title('(c) remove either half and it collapses');

exportgraphics(f, fullfile(resdir,'profile_vs_refit.png'),'Resolution',140);
fprintf('saved -> %s\n', fullfile(resdir,'profile_vs_refit.png'));
