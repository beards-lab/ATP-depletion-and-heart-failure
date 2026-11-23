% RunFeatureAudit.m
% Honest per-segment comparison: model vs data at BOTH 8 mM and 2 mM, absolute values.
% Clears up whether "A and steady are terrible" is real or a plot artifact
% (RunFinalBake plotted model@2mM against 8mM data because velocitytableonfile='...8mM...').
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); LoadData;

setrun=@(p) setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield(setfield( ...
  p,'RunForceVelocity',false),'RunKtr',false),'RunSlack',true),'RunStairs',false),'RunForceVelocityTime',false), ...
  'EvalFeatures',true),'BreakOnODEUnstable',false),'RunSlackSegments','All'),'MaxRunTime',120), ...
  'velocitytableonfile','protocol_03_27_2026_8mM_slack.mat');

base = setrun(getParams(loadParams('params_2state_final'),[],true,false));

f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;

p8=base; p8.MgATP=8; p8.MgADP=0; p8.Pi=0; [~,~,fm8,~]=runSlackExperiment(p8);
p2=base; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8; [~,~,fm2,~]=runSlackExperiment(p2);

show={'steady','A','Am','peak1_y','vall_y','peak2','ktr','t0'};
fprintf('\n==== 8 mM: model vs data (absolute) ====\n');
for i=1:numel(show); f=show{i};
  fprintf('%-10s model %s | data %s\n', f, num2str(double(fm8.(f)(:))','%7.3g'), num2str(double(f8d.(f)(:))','%7.3g'));
end
fprintf('\n==== 2 mM: model vs data (absolute) ====\n');
for i=1:numel(show); f=show{i};
  fprintf('%-10s model %s | data %s\n', f, num2str(double(fm2.(f)(:))','%7.3g'), num2str(double(f2d.(f)(:))','%7.3g'));
end

% absolute per-condition feature cost (same selectors, per-segment)
fnabs={'A|50','Am|10','steady|50','peak1_y|10','vall_y|10','peak2|5','ktr|2','t0|2','peak1_dSL|1','restretchSlopeStart|1','vall2_dy|1','ovrsht_dy|1'};
c8=sum(evalFeatureCost(f8d,fm8,fnabs,1));
c2=sum(evalFeatureCost(f2d,fm2,fnabs,1));
fprintf('\nAbsolute feature cost:  8mM = %.3f   2mM = %.3f\n', c8, c2);

% ---- figure: model vs data at BOTH conditions, absolute ----
fig=figure('Position',[40 40 1500 820],'Color','w');
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
for i=1:numel(show); f=show{i}; nexttile(tl); hold on; box on;
  plot(1:5,double(f8d.(f)),'bo-','LineWidth',1.2,'DisplayName','data 8mM');
  plot(1:5,double(fm8.(f)),'bx--','LineWidth',1.2,'DisplayName','model 8mM');
  plot(1:5,double(f2d.(f)),'ro-','LineWidth',1.2,'DisplayName','data 2mM');
  plot(1:5,double(fm2.(f)),'rx--','LineWidth',1.2,'DisplayName','model 2mM');
  title(f,'Interpreter','none'); xlabel('segment');
  if i==1; legend('Location','best','FontSize',7); end
end
sgtitle(sprintf('Absolute model vs data, both ATP (8mM abs cost %.2f, 2mM abs cost %.2f)',c8,c2));
resdir=fullfile(root,'Analyses','LowATP_kmsrd','results');
exportgraphics(fig, fullfile(resdir,'AUDIT_abs_both.png'),'Resolution',140);
disp('AUDIT DONE');
