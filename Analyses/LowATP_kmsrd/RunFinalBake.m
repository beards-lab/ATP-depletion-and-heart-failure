% RunFinalBake.m
% =========================================================================
% Bakes and verifies the FINAL one-file low-ATP mechanism: params/params_2state_final.m.
% Fully [ATP]-driven (MgADP=0, physiological CK/CrP-buffered fiber): UseAtpDetachR2 (R2 detachment,
% force engine) + UseAtpKmsrd (SRX-ADP return, ktr weld-break + per-segment slope). Both gates
% normalized to MgATP_ref=8 mM => exactly inert at the 8 mM reference (verified max|d|=0).
% Best (K_T2, K_srx) from a 20-point grid: K_T2=2.0, K_srx=7.0 -> ratio cost 0.1188.
% Produces the final trace + features figures in results/.
% Run: cd(root); addpath(genpath('.')); Analyses/LowATP_kmsrd/RunFinalBake
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

dst=fullfile(root,'params','params_2state_final.m');
if ~exist(dst,'file')
  copyfile(fullfile(root,'params','optfull_opt.m'), dst);
  fid=fopen(dst,'a');
  fprintf(fid,'\r\n%%%% ---- Low-ATP: ATP-ONLY mechanism (no ADP-trap; physiological CK/CrP-buffered fiber) ----\r\n');
  fprintf(fid,'params0.UseAdpTrap = 0;\r\nparams0.UseAtpDetachR2 = 1;\r\nparams0.K_T2 = 2.0;\r\n');
  fprintf(fid,'params0.UseAtpKmsrd = 1;\r\nparams0.K_srx = 7.0;\r\nparams0.MgATP_ref = 8;\r\n');
  fprintf(fid,'params0.UsePiReversal = 0;\r\nparams0.UsePiForce = 0;\r\n');
  fprintf(fid,'params0.MgATP = 2;\r\nparams0.MgADP = 0;\r\nparams0.Pi = 0.8;\r\n'); fclose(fid);
end

setrun=@(p) setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield( ...
  p,'RunForceVelocity',false),'RunKtr',false),'RunSlack',true),'RunStairs',false),'RunForceVelocityTime',false), ...
  'EvalFeatures',true),'BreakOnODEUnstable',false),'RunSlackSegments','All'),'MaxRunTime',120),'velocitytableonfile','protocol_03_27_2026_8mM_slack.mat');

params0 = setrun(getParams(loadParams('params_2state_final'), [], true, false));
opt     = setrun(getParams(loadParams('optfull_opt'),         [], true, false));

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
ff={'steady','A','Am','peak2','peak1_y','vall_y'}; kf={'ktr','t0'}; target=struct('scalar',struct(),'vector',struct());
for i=1:numel(ff); target.scalar.(ff{i})=mean(double(f2d.(ff{i}))./double(f8d.(ff{i})),'omitnan'); end
for i=1:numel(kf); target.vector.(kf{i})=double(f2d.(kf{i}))./double(f8d.(kf{i})); end
target.weights=struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);

opt8=opt; opt8.MgATP=8;opt8.MgADP=0;opt8.Pi=0; [~,~,fmO,~]=runSlackExperiment(opt8);
p8=params0; p8.MgATP=8;p8.MgADP=0;p8.Pi=0; [~,~,fm8,~]=runSlackExperiment(p8);
fprintf('8mM invariance: max|d steady| = %.4g\n', max(abs(double(fm8.steady)-double(fmO.steady))));

[~,out2,fm2,fd2] = runSlackExperiment(params0);
[c,d] = atpRatioCost(fm2,fm8,target);
fprintf('FINAL 2mM ratio cost = %.4f | force %.3f peak1 %.3f | ktr %s (data %s)\n', ...
  c, d.steady.ratio, d.peak1_y.ratio, num2str(d.ktr.ratio(:)','%.2f '), num2str(target.vector.ktr(:)','%.2f '));

% --- trace figure (cropped to the actual restretch window) ---
tmin=min(out2.datatable(:,1)); tmax=max(out2.datatable(:,1));
fig=figure('Position',[50 50 1400 500],'Color','w');
plot(out2.datatable(:,1), out2.datatable(:,3), 'k-', 'LineWidth', 1.3, 'DisplayName','data (2 mM)'); hold on; box on;
plot(out2.t, out2.Force, 'r-', 'LineWidth', 1.6, 'DisplayName','model (2 mM)');
xlim([tmin-1, tmax+1]); ylabel('Force'); xlabel('time (s)');
title(sprintf('Final low-ATP model: 2 mM slack-restretch trace (model vs data, ratio cost %.3f)', c));
legend('Location','best'); grid on;
exportgraphics(fig, fullfile(resdir,'FINAL_2mM_trace.png'), 'Resolution', 140);

% --- features figure ---
figure('Color','w');
plotFeatures(fd2, fm2, [], {'steady','A','Am','ktr','peak1_y','vall_y','peak2','t0'});
sgtitle(sprintf('2 mM: data vs model, all features (ratio cost %.3f)', c));
exportgraphics(gcf, fullfile(resdir,'FINAL_2mM_features.png'), 'Resolution', 140);

save(fullfile(resdir,'FINAL_result.mat'),'fm2','fd2','fm8','target','d','c','out2');
