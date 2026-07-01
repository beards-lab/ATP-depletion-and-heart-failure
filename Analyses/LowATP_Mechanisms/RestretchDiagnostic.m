%% RestretchDiagnostic.m
% Visualises the restretch-peak trade-off diagnosed in Part A/B (kSE vs peak1).
% Values are measured sweeps on the rigor 8 mM baseline (see labdiary); this
% script plots them against the 8 mM data targets. Re-run the sweeps in the
% labdiary to regenerate the numbers.
% Out: results/restretch_tradeoff.png
clear; clc;
here=fileparts(mfilename('fullpath')); resdir=fullfile(here,'results');
if ~exist(resdir,'dir'); mkdir(resdir); end

% measured kSE sweep on the rigor 8 mM baseline (kSE x{0.5,0.75,1,1.5,2})
kSE   = [1601 2402 3203 4804 6406];
slope = [ 793 1063 1285 1611 1860];   % model restretchSlopeStart
peak1 = [108.6 111.8 114.0 116.8 118.5];
% 8 mM data targets
slope_data = 1588; peak1_data = 89;

fig=figure('Position',[100 100 760 460]); hold on; box on;
yyaxis left;
plot(kSE, slope, '-o','LineWidth',1.8,'Color',[0.1 0.35 0.85]);
yline(slope_data,'--','Color',[0.1 0.35 0.85],'LineWidth',1.3);
ylabel('restretchSlopeStart (kSE proxy)','Color',[0.1 0.35 0.85]);
text(4900, slope_data+60,'data slope 1588','Color',[0.1 0.35 0.85],'FontSize',9);
yyaxis right;
plot(kSE, peak1, '-s','LineWidth',1.8,'Color',[0.85 0.3 0.1]);
yline(peak1_data,'--','Color',[0.85 0.3 0.1],'LineWidth',1.3);
ylabel('peak1 (restretch overshoot)','Color',[0.85 0.3 0.1]);
text(1700, peak1_data-4,'data peak1 89','Color',[0.85 0.3 0.1],'FontSize',9);
xline(4804,':k'); text(4804,low(),'');
xlabel('kSE (serial-element stiffness)');
title({'Restretch trade-off: kSE that matches the data slope (~4800) overshoots peak1.',...
       'Levers — peak1\downarrow needs kstiff2\downarrow (drops steady); c\_SE\_visc damps peak1 but crashes slope; A2-shift AMPLIFIES.'});
exportgraphics(fig, fullfile(resdir,'restretch_tradeoff.png'),'Resolution',130);
fprintf('Saved results/restretch_tradeoff.png\n');
function y=low(); yl=ylim; y=yl(1); end
