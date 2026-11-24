% RunLowATPBest.m
% =====================================================================================
%  >> RunLowATPBest
%
%  SINGLE ENTRY POINT to run the current best low-ATP parametrization and show the fit.
%
%  Loads  params/lowATP_TNF_best.m  = ThursdayNightFever high-ATP baseline (merged from
%  Further-optim) + the ATP-driven low-ATP mechanism (two levers, both inert at 8 mM):
%     UseAtpDetachR2 / K_T2  -> R2 detachment slowed by low [ATP]  (force amplitude)
%     UseAtpKmsrd   / K_srx  -> SRX-ADP return slowed by low [ATP] (ktr magnitude+slope)
%
%  It is a ONE-file, TWO-condition model:
%     * as saved (MgATP=2) -> LOW-ATP (2 mM); compared here to the 2 mM slack data
%     * set p.MgATP = 8    -> HIGH-ATP baseline, byte-identical to ThursdayNightFever
%
%  Renders: (1) 2 mM slack-restretch trace, model vs data; (2) per-feature comparison.
% =====================================================================================
clear
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;

% ---- load the best parametrization ----
p = getParams(loadParams('lowATP_TNF_best'), [], true, false);

% ---- slack-only, feature evaluation, driven by the 2 mM protocol ----
p.RunForceVelocity=false; p.RunKtr=false; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.EvalFeatures=true; p.BreakOnODEUnstable=false;
p.RunSlackSegments='All'; p.MaxRunTime=120;
p.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';

[~,out,fm,fd] = runSlackExperiment(p);

% ---- report (absolute high-ATP-style feature cost) ----
fnabs = {'steady|50','A|50','peak1_y|10','vall_y|10','peak2|5','t0|2','ktr|30'};
absCost = sum(evalFeatureCost(fd, fm, fnabs, 2));
fprintf('\n===== BEST low-ATP fit (2 mM) =====\n');
fprintf('K_T2=%.3f  K_srx=%.3f   ->  absolute feature cost %.2f\n', p.K_T2, p.K_srx, absCost);
show={'steady','A','peak1_y','vall_y','peak2','ktr','t0'};
fprintf('feature      model                       data\n');
for i=1:numel(show); f=show{i}; if ~isfield(fm,f)||~isfield(fd,f); continue; end
  fprintf('  %-10s %-26s | %s\n', f, num2str(double(fm.(f)(:))','%6.3g '), num2str(double(fd.(f)(:))','%6.3g ')); end

% ---- (1) trace ----
figure('Position',[60 60 1400 500],'Color','w');
plot(out.datatable(:,1), out.datatable(:,3),'k-','LineWidth',1.2,'DisplayName','data (2 mM)'); hold on; box on;
plot(out.t, out.Force,'r-','LineWidth',1.6,'DisplayName','model (2 mM)');
tmax=max(out.datatable(:,1)); xlim([tmax-4.6, tmax+0.3]);   % crop to the 5 slack-restretch cycles
xlabel('time (s)'); ylabel('Force (rel.)'); legend('Location','best'); grid on;
title(sprintf('Best low-ATP fit (2 mM slack-restretch): K_{T2}=%.2f  K_{srx}=%.1f   [abs cost %.2f]', p.K_T2, p.K_srx, absCost));

% ---- (2) features ----
figure('Position',[60 60 1400 800],'Color','w');
plotFeatures(fd, fm, [], {'steady','A','Am','ktr','t0','peak1_y','vall_y','peak2'});
sgtitle(sprintf('Best low-ATP fit (2 mM): data (o) vs model (x)   [abs cost %.2f]', absCost));
%% (optional) all fitted features at a glance
figure();
plotFeatures(fd, fm, [], p.fn)
