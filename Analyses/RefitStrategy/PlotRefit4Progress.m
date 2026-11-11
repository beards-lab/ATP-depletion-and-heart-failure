% PlotRefit4Progress.m -- where the running refit4 campaign stands, on every
% protocol, against the target each one is actually scored on.
%
% Two figures:
%   results/refit4_protocols.png  slack (vs the POOLED trace), F-V, ktr on all
%                                 three days, the STAIRCASE (held out), passive
%   results/refit4_features.png   every scored feature + a paired cost bar
%
% Reads params/refit4_opt.m live, so it can be re-run while the campaign is
% campaign is still going. refit4_seed = refit3_opt PARAMETERS re-scored on the
% re-curated pooled target (peak-prominence fix), so the seed line here is the
% same model as refit3_opt, only the objective underneath it changed.
%
%   matlab -batch "run('Analyses/RefitStrategy/PlotRefit4Progress.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_pg_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 3); end

TAGS = {'refit4_seed','refit4_opt'};
DAYS = {'03_27','04_03','04_10'};
PT   = pooledTarget8mM();
LoadData;

R = struct();
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true);
    params0 = p; RunBakersExp;                     %#ok<*NASGU>
    R(k).p=p; R(k).fm=features_model; R(k).fd=features_data;
    [R(k).cost,~,R(k).fc,R(k).fn] = evaluateBakersExp([], p, false, false);
    R(k).os=out_slack; R(k).ok=out_ktr; R(k).opas=out_pas; R(k).outs=outs;
    % staircase -- NOT in the objective (RunStairs = 0), a pure prediction
    ps = p; ps.RunStairs = 1; ps = getParams(ps, ps.g, true);
    [R(k).Estairs, R(k).ostairs] = runStairsExperiment(ps);
    fprintf('%-14s cost %.4f   stairs(held out) %.4f\n', TAGS{k}, R(k).cost, R(k).Estairs);
end

%% ================= FIGURE 1: protocols =================
f1=figure('Position',[20 20 1700 940],'Color','w');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
CS=[0.85 0.33 0.10]; CO=[0 0.45 0.74]; CD=[0.1 0.1 0.1]; CB=[0.78 0.78 0.78];

% --- pooled slack trace, built as in PlotPooledTimecourse ---
DTs=2.5e-4; TREL=(-0.40:DTs:3.55)'; Y=nan(numel(TREL),3);
for i=1:3
    S=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{i})));
    v=S.velocitytable; iS=find(v(:,2)<-1); t0=v(iS(1),1); d=S.datatable;
    a=median(d(d(:,1)>t0-0.30 & d(:,1)<t0-0.01,3));
    Y(:,i)=interp1(d(:,1)-t0,d(:,3),TREL,'linear',NaN)*(PT.Fref/a);
    if i==1; VT1=v; VT1(:,1)=VT1(:,1)-t0; end
end
POOL=mean(Y,2,'omitnan');
svt=load(fullfile(root,'data',R(1).p.velocitytableonfile)).velocitytable; iSv=find(svt(:,2)<-1);

nexttile([1 2]); hold on; box on; grid on
for i=1:3; plot(TREL,Y(:,i),'-','Color',CB,'LineWidth',0.5); end
plot(TREL,POOL,'-','Color',CD,'LineWidth',1.8);
for k=1:2
    plot(R(k).os.t-svt(iSv(1),1), R(k).os.Force,'-','Color',(k==1)*CS+(k==2)*CO, ...
         'LineWidth',1.0+0.5*(k==2));
end
xlim([-0.3 3.4]); xlabel('t since slack #1 (s)'); ylabel('F (kPa)');
title('A   slack, against the POOLED trace (mean of 3 preps, own-Fss22 normalised)');
legend({'','','','pooled mean','seed','opt'},'Location','southeast','Box','off','NumColumns',2);

% --- FV ---
AllV=[0 0.5 1 2 3 4 5 6];
D1=[56.40,51.8120,37.4459,17.8025,11.4430,6.2643,3.2759,2.2120];
D2=[67.3942,60.9885,42.7059,18.8539,12.5956,6.4426,3.5136,1.8556];
nexttile; hold on; box on; grid on
plot(AllV,D1/D1(1),'ks-','MarkerFaceColor','k','LineWidth',1.2,'MarkerSize',5);
plot(AllV,D2/D2(1),'k^--','LineWidth',1.0,'MarkerSize',5);
for k=1:2
    Fa=arrayfun(@(o) o.Force(end), R(k).outs); vv=-R(k).p.FV_velocities;
    [vv,ix]=sort(vv); Fa=Fa(ix);
    plot(vv,Fa/Fa(1),'o-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.4);
end
xlabel('v (ML/s)'); ylabel('F / F_{iso}');
title(sprintf('B   F-V (cost %.4f -> %.4f)', ...
      R(1).fc(strcmp(R(1).fn,'FV_normpowerAvg')), R(2).fc(strcmp(R(2).fn,'FV_normpowerAvg'))));
legend({'Baker22','Isovel21','seed','opt'},'Box','off');

% --- ktr, three days ---
for d=1:3
    K=load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
    V=K.velocitytable; DTk=K.datatable; tK=V(end-1,1);
    nexttile; hold on; box on; grid on
    m=DTk(:,1)>=tK-0.020;
    plot((DTk(m,1)-tK)*1e3, DTk(m,3),'-','Color',CD,'LineWidth',1.4);
    kk=zeros(1,2);
    for k=1:2
        pk=R(k).p; pk.SL0=2*V(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=V(:,2);
        [~,o]=evaluateModel(str2func(pk.modelFcn),V(:,1),pk);
        Fi=interp1(o.t,o.Force,V(2,1),'linear','extrap');
        mm=o.t>=tK-0.020;
        plot((o.t(mm)-tK)*1e3, o.Force(mm)/Fi,'-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.0+0.4*(k==2));
        f=extractKtrAttributes(o.t,o.Force/Fi,V); kk(k)=f.ktrP_k63;
    end
    fd=extractKtrAttributes(DTk(:,1),DTk(:,3),V);
    xlim([-20 59]); ylim([-0.15 1.15]); xlabel('t since release (ms)'); ylabel('F/F_{iso}');
    title(sprintf('%s   ktr %s  k_{63} %.0f | %.0f | %.0f', char('C'+d), strrep(DAYS{d},'_','/'), ...
          fd.ktrP_k63, kk(1), kk(2)));
    if d==1; legend({'data','seed','opt'},'Location','southeast','Box','off'); end
end

% --- staircase, HELD OUT ---
nexttile; hold on; box on; grid on
St=load(fullfile(root,'data',R(1).p.stairs_velocitytableonfile));
plot(St.datatable(:,1), St.datatable(:,3),'-','Color',CD,'LineWidth',1.3);
for k=1:2
    o=R(k).ostairs;
    Fn=o.Force/interp1(o.t,o.Force,St.velocitytable(2,1),'linear','extrap');
    plot(o.t, Fn,'-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.0+0.4*(k==2));
end
% the table opens with 5.5 s of isometric warm-up (row 1 at t = 67, first
% command at row 2) and the recording itself only starts at 72.3 -- plotting
% from 67 spends 80 % of the axis on a flat line. Zoom to the manoeuvre.
Svt = St.velocitytable; iMove = find(abs(Svt(:,2)) > 1e-6, 1);
tA = Svt(iMove,1) - 0.08; tB = min(Svt(end,1), St.datatable(end,1));
xlim([tA tB]);
yyaxis right
stairs(Svt(:,1), Svt(:,4)*2,'-','Color',[0.6 0.6 0.6],'LineWidth',0.8);
ylabel('SL commanded (\mum)'); ylim([2.0 2.30]);
set(gca,'YColor',[0.6 0.6 0.6]); yyaxis left
xlabel('t (s)'); ylabel('F (norm.)');
title(sprintf('F   STAIRCASE — held out (E %.4f -> %.4f)', R(1).Estairs, R(2).Estairs));
legend({'data','seed','opt'},'Location','northwest','Box','off');

% --- passive ---
nexttile; hold on; box on; grid on
for k=1:2
    if isempty(R(k).opas); continue; end
    plot(R(k).opas.t, R(k).opas.Force,'-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.0+0.4*(k==2));
end
xlabel('t (s)'); ylabel('F (kPa)');
title('G   passive — now floored at seed x1.2 (a guard, not a target)');
legend({'seed','opt'},'Box','off');

% --- A/steady ---
nexttile; hold on; box on; grid on
jA=find(strcmp(PT.name,'A'),1); jS=find(strcmp(PT.name,'steady'),1);
plot(1:5, PT.val(jA,:)./PT.val(jS,:),'o-','Color',CD,'LineWidth',2,'MarkerFaceColor',CD);
for k=1:2
    plot(1:5, R(k).fm.A(:)'./R(k).fm.steady(:)','s-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.4);
end
xlabel('slack segment'); ylabel('A / steady'); xticks(1:5); xlim([0.8 5.2]);
title(sprintf('H   A/steady  (cost %.4f -> %.4f)', ...
      R(1).fc(strcmp(R(1).fn,'A_over_steady')), R(2).fc(strcmp(R(2).fn,'A_over_steady'))));
legend({'pooled target','seed','opt'},'Location','southwest','Box','off');

sgtitle(sprintf('refit4 progress: seed %.4f -> current best %.4f', R(1).cost, R(2).cost), ...
        'FontSize',12,'Interpreter','none');
exportgraphics(f1, fullfile(resdir,'refit4_protocols.png'),'Resolution',110);

%% ================= FIGURE 2: features =================
FN=R(1).fn; lab=cellfun(@(s) regexprep(s,'\|.*$',''), FN,'UniformOutput',false);
f2=figure('Position',[20 20 1700 950],'Color','w');
cmap=containers.Map(FN,num2cell(R(2).fc));
plotMultipleFeatures({R(1).fd, R(1).fm, R(2).fm}, {'data','seed','opt'}, ...
    [0 0 0; CS; CO], {'o','x','s'}, FN, {':','-','-'}, [false false true], cmap);
for ax=findobj(f2,'Type','axes')'
    ln=findobj(ax,'Type','line'); if isempty(ln); continue; end
    ln=ln(strcmp(get(ln,'DisplayName'),'data')); if isempty(ln); continue; end
    set(ln,'MarkerSize',8,'LineWidth',1.2);
    ax.Children=[ln(:); setdiff(ax.Children,ln,'stable')];
end
sgtitle('every scored feature: data (black o) vs seed (orange x) vs current opt (blue s)', ...
        'FontSize',11,'Interpreter','none');
exportgraphics(f2, fullfile(resdir,'refit4_features.png'),'Resolution',110);

f3=figure('Position',[40 40 1200 720],'Color','w');
c1=R(1).fc(:); c2=R(2).fc(:); d=c2-c1;
[~,ord]=sort(max(c1,c2),'descend');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; b=barh(1:numel(ord),[c1(ord) c2(ord)],0.82);
b(1).FaceColor=CS; b(2).FaceColor=CO; b(1).EdgeColor='none'; b(2).EdgeColor='none';
set(gca,'YTick',1:numel(ord),'YTickLabel',lab(ord),'YDir','reverse','FontSize',9, ...
        'TickLabelInterpreter','none'); grid on; xlabel('weighted cost');
legend({sprintf('seed %.4f',sum(c1)),sprintf('opt %.4f',sum(c2))},'Location','southeast','Box','off');
title('per-feature cost');
nexttile; [~,od]=sort(d,'descend'); hold on
for i=1:numel(od)
    col=[0.8 0.2 0.2]; if d(od(i))<0; col=[0.2 0.6 0.3]; end
    barh(i,d(od(i)),0.82,'FaceColor',col,'EdgeColor','none');
end
xline(0,'k-'); set(gca,'YTick',1:numel(od),'YTickLabel',lab(od),'YDir','reverse','FontSize',9, ...
        'TickLabelInterpreter','none'); grid on; xlabel('\Delta cost (opt - seed)');
title(sprintf('what has moved so far (net %+.4f)', sum(d)));
exportgraphics(f3, fullfile(resdir,'refit4_bars.png'),'Resolution',120);

fprintf('\n%-24s %10s %10s %10s\n','feature','seed','opt','delta');
for i=od(:)'
    if abs(d(i))<1e-6 && c1(i)<1e-6; continue; end
    fprintf('%-24s %10.4f %10.4f %+10.4f\n', lab{i}, c1(i), c2(i), d(i));
end
fprintf('%-24s %10.4f %10.4f %+10.4f\n','TOTAL',sum(c1),sum(c2),sum(d));
fprintf('%-24s %10.4f %10.4f %+10.4f   <- HELD OUT\n','stairs',R(1).Estairs,R(2).Estairs,R(2).Estairs-R(1).Estairs);
fprintf('\nwrote refit4_protocols.png / refit4_features.png / refit4_bars.png\n');
fprintf('\nDONE\n');
