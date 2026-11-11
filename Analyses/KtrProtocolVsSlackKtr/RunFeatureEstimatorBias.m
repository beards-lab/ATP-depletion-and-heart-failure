% RunFeatureEstimatorBias.m
% The objective's own slack-ktr estimator (extractSlackAttributes) drops every
% sample below FORCE_THRESHOLD_KPA = 10 kPa and then fits A(1-exp(-k(t-t0)))
% with F0 pinned at 0 and t0 free. On a clean single exponential that cut is
% harmless. Is it harmless on the MODEL's trace?
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
M=load(fullfile(resdir,'estimator_table.mat')).M;
base=getParams(loadParams('estiff125_opt'),[],true,false);
svt=load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
iS=find(svt(:,2)<-1); iR=find(svt(:,2)>1);
Sd=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));

fprintf('\n=== post-slack #2: what the 10 kPa cut does to the fitted rate ===\n');
fprintf('%-26s %9s %9s %9s\n','trace','thr = 0','thr = 10','change');
rows={}; ta=svt(iS(2)+1,1); tb=svt(iR(2),1); t_org=svt(iS(2),1);
m=Sd.datatable(:,1)>=ta & Sd.datatable(:,1)<=tb;
rows{1}={'DATA', Sd.datatable(m,1)-t_org, Sd.datatable(m,3)};
for v=1:2
    t=M(v).os.t(:); F=M(v).os.Force(:); m=t>=ta & t<=tb;
    rows{end+1}={M(v).name, t(m)-t_org, F(m)}; %#ok
end
B=nan(numel(rows),2);
for i=1:numel(rows)
    B(i,1)=iFit(rows{i}{2},rows{i}{3},0);
    B(i,2)=iFit(rows{i}{2},rows{i}{3},10);
    fprintf('%-26s %9.1f %9.1f %+8.1f%%\n',rows{i}{1},B(i,1),B(i,2),100*(B(i,2)/B(i,1)-1));
end
fprintf(['\nThe cut is worth +1.6%% on the data and -9 to -11%% on the model: the model''s\n' ...
         'rise has a fast early foot that the data does not have, and removing the first\n' ...
         '13%% of the amplitude while holding F0 = 0 makes the remainder fit a slower\n' ...
         'exponential. The feature is therefore biased LOW on model traces only -- a\n' ...
         'one-sided bias, so it pushes the fit toward models that are genuinely too fast.\n']);
save(fullfile(resdir,'feature_bias.mat'),'B','rows');

function k = iFit(t,y,thr)
    w = y > thr; t=t(w); y=y(w);
    f = fit(t(:),y(:),@(A,k,t0,x) A.*(1-exp(-k.*max(x-t0,0))), ...
            'StartPoint',[100,50,0.01],'Lower',[10,0.01,0],'Upper',[200,500,0.1]);
    k = f.k;
end
