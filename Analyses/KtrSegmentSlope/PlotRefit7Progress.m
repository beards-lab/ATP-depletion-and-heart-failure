% PlotRefit7Progress.m -- the 20-minute render: every protocol, with the
% ATP-carrying metrics up front, for the seed and the campaign's current best.
%
% FAST BY CONSTRUCTION. The previous renders took ~30 min because each model
% was simulated three times (RunBakersExp, then evaluateBakersExp again, then
% the staircase). Here each model is simulated ONCE -- the cost comes from
% evalFeatureCost on the features RunBakersExp already produced -- and every
% result is cached against the snapshot file's timestamp, so the seed is
% simulated once per campaign and the best only when it has changed.
%
% STAIRCASE is drawn as a visual check only: its data are a single
% preparation, not pooled, so its normalised force at 1.0 and 1.1 ML does not
% match the pooled level the slack block is scored on. No number is reported.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/PlotRefit7Progress.m')"
%   matlab -batch "TAGS={'refit7_seed','refit7_opt'}; run(...)"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('TAGS','var'); TAGS = {'refit8_seed','refit8_opt'}; end
if ~exist('FIGTAG','var'); FIGTAG = 'refit8'; end
cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_p7_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 2); end

DAYS = {'03_27','04_03','04_10'};
PT = pooledTarget8mM(); PK = pooledKtrTarget();
ATPC = {'ktr80','ktr80_shape','ktrfit','ktrfit_shape','rsK80','rsT20','ktrP_k80','steady','t0_crossing'};
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady','vall_t','peak1_dSL'};

R = struct();
for k = 1:numel(TAGS)
    pf = fullfile(root,'params',[TAGS{k} '.m']);
    di = dir(pf); key = sprintf('%s_%s', TAGS{k}, datestr(di.datenum,'yyyymmddHHMMSS'));
    cf = fullfile(resdir, ['render_cache_' key '.mat']);
    if isfile(cf)
        L = load(cf); R(k).S = L.S; fprintf('%-12s cached (%s)\n', TAGS{k}, key); continue;
    end
    t0 = tic;
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.PlotFeatureFitting=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600; p.RunStairs=0;
    p.RunSlackSegments='AllPar'; p.g = 1;
    params0 = p; LoadData; RunBakersExp;                       %#ok<*NASGU>
    S.fn = p.fn;
    [S.fc, ~, ~] = evalFeatureCost(features_data, features_model, p.fn);
    S.cost = sum(S.fc); S.fm = features_model; S.fd = features_data;
    S.slackT = out_slack.t; S.slackF = out_slack.Force;
    S.fvF = arrayfun(@(o) o.Force(end), outs); S.fvV = -p.FV_velocities;
    S.ktr = struct('t',{},'F',{},'k80',{});
    for d = 1:3
        K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
        V = K.velocitytable;
        pk = p; pk.SL0 = 2*V(1,4); pk = getParams(pk,pk.g,true); pk.Velocity = V(:,2);
        [~,o] = evaluateModel(str2func(pk.modelFcn), V(:,1), pk);
        Fi = interp1(o.t, o.Force, V(2,1), 'linear', 'extrap');
        f = extractKtrAttributes(o.t, o.Force/Fi, V);
        S.ktr(d) = struct('t', o.t, 'F', o.Force/Fi, 'k80', f.ktrP_k80);
    end
    ps = p; ps.RunStairs = 1; ps = getParams(ps, ps.g, true);
    [~, os] = runStairsExperiment(ps);
    S.stT = os.t; S.stF = os.Force;
    S.wall = toc(t0);
    save(cf, 'S');
    R(k).S = S;
    fprintf('%-12s cost %.4f   (simulated in %.0f s)\n', TAGS{k}, S.cost, S.wall);
end

%% ---- the ATP-carrying numbers, model vs data -----------------------------
jk = @(nm) find(strcmp(PT.name, nm), 1);
kd  = PT.val(jk('ktr80'),:);   bd = polyfit(1:5,kd,1);
rd  = PT.val(jk('rsK80'),1:4); td = PT.val(jk('rsT20'),1:4);
kpd = PK.val(strcmp(PK.name,'ktrP_k80'));
M = struct();
for k = 1:numel(TAGS)
    fm = R(k).S.fm;
    b = polyfit(1:5, fm.ktr80(:)', 1);
    M(k).tag = TAGS{k};
    M(k).cost = R(k).S.cost;
    M(k).ktr80_mean = mean(fm.ktr80);
    M(k).ktr80_slope_pct = 100*b(1)/mean(fm.ktr80);
    M(k).rsK80_mean = mean(fm.rsK80(1:4));
    M(k).rsT20_mean_ms = 1e3*mean(fm.rsT20(1:4));
    M(k).ktrP_k80 = mean([R(k).S.ktr.k80]);
    lab = cellfun(@(s) regexprep(regexprep(s,'\|.*$',''),'\[.*$',''), R(k).S.fn, 'uni', 0);
    M(k).atp_cost = sum(R(k).S.fc(ismember(lab, ATPC)));
    M(k).guard_cost = sum(R(k).S.fc(ismember(lab, GUARD)));
end
DATA = struct('ktr80_mean', mean(kd), 'ktr80_slope_pct', 100*bd(1)/mean(kd), ...
              'rsK80_mean', mean(rd), 'rsT20_mean_ms', 1e3*mean(td), 'ktrP_k80', kpd);
J = struct('updated', datestr(now,'yyyy-mm-dd HH:MM'), 'data', DATA, 'models', M);
fid = fopen(fullfile(resdir,['render_' FIGTAG '.json']),'w'); fprintf(fid,'%s',jsonencode(J)); fclose(fid);

%% ================= FIGURE: protocols =================
f1 = figure('Position',[20 20 1700 1250],'Color','w');
tiledlayout(4,3,'Padding','compact','TileSpacing','compact');
CS=[0.85 0.33 0.10]; CO=[0 0.45 0.74]; CD=[0.1 0.1 0.1]; CB=[0.78 0.78 0.78];
col = @(k) (k==1)*CS + (k>1)*CO;
LBL = strrep(TAGS,'_','\_');

% --- A: pooled slack trace ---
DTs=2.5e-4; TREL=(-0.40:DTs:3.55)'; Y=nan(numel(TREL),3);
for i=1:3
    Sd=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{i})));
    v=Sd.velocitytable; iS=find(v(:,2)<-1); tS=v(iS(1),1); dd=Sd.datatable;
    a=median(dd(dd(:,1)>tS-0.30 & dd(:,1)<tS-0.01,3));
    Y(:,i)=interp1(dd(:,1)-tS,dd(:,3),TREL,'linear',NaN)*(PT.Fref/a);
end
p1 = getParams(loadParams(TAGS{1}), [], true, false);
svt=load(fullfile(root,'data',p1.velocitytableonfile)).velocitytable; iSv=find(svt(:,2)<-1);
nexttile([1 2]); hold on; box on; grid on
for i=1:3; plot(TREL,Y(:,i),'-','Color',CB,'LineWidth',0.5); end
plot(TREL,mean(Y,2,'omitnan'),'-','Color',CD,'LineWidth',1.8);
for k=1:numel(TAGS); plot(R(k).S.slackT-svt(iSv(1),1), R(k).S.slackF,'-','Color',col(k),'LineWidth',1+0.5*(k>1)); end
xlim([-0.3 3.4]); xlabel('t since slack #1 (s)'); ylabel('F (kPa)');
title('A   slack against the POOLED trace (3 preps, own-Fss22 normalised)');
legend([{'','','','pooled mean'}, LBL],'Location','southeast','Box','off','NumColumns',2);

% --- B: F-V ---
AllV=[0 0.5 1 2 3 4 5 6];
D1=[56.40,51.8120,37.4459,17.8025,11.4430,6.2643,3.2759,2.2120];
D2=[67.3942,60.9885,42.7059,18.8539,12.5956,6.4426,3.5136,1.8556];
nexttile; hold on; box on; grid on
plot(AllV,D1/D1(1),'ks-','MarkerFaceColor','k','LineWidth',1.2,'MarkerSize',5);
plot(AllV,D2/D2(1),'k^--','LineWidth',1.0,'MarkerSize',5);
for k=1:numel(TAGS)
    [vv,ix]=sort(R(k).S.fvV); Fa=R(k).S.fvF(ix);
    plot(vv,Fa/Fa(1),'o-','Color',col(k),'LineWidth',1.4);
end
xlabel('v (ML/s)'); ylabel('F / F_{iso}'); title('B   F-V');
legend([{'Baker22','Isovel21'}, LBL],'Box','off');

% --- C-E: ktr protocol, three days ---
for d=1:3
    K=load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
    V=K.velocitytable; DTk=K.datatable; tK=V(end-1,1);
    nexttile; hold on; box on; grid on
    m=DTk(:,1)>=tK-0.020; plot((DTk(m,1)-tK)*1e3, DTk(m,3),'-','Color',CD,'LineWidth',1.4);
    for k=1:numel(TAGS)
        s=R(k).S.ktr(d); mm=s.t>=tK-0.020;
        plot((s.t(mm)-tK)*1e3, s.F(mm),'-','Color',col(k),'LineWidth',1+0.4*(k>1));
    end
    fd=extractKtrAttributes(DTk(:,1),DTk(:,3),V);
    xlim([-20 59]); ylim([-0.15 1.15]); xlabel('t since release (ms)'); ylabel('F/F_{iso}');
    ks = arrayfun(@(k) sprintf(' | %.0f', R(k).S.ktr(d).k80), 1:numel(TAGS), 'uni', 0);
    title(sprintf('%s   ktr protocol %s   k_{80} %.0f%s', char('B'+d), strrep(DAYS{d},'_','/'), ...
          fd.ktrP_k80, [ks{:}]));
    if d==1; legend([{'data'}, LBL],'Location','southeast','Box','off'); end
end

% --- F: ktr80 ladder ---
nexttile; hold on; box on; grid on
for d3=1:3
    Sd=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d3})));
    plot(1:5, Sd.features_data.ktr80(:)','-','Color',CB,'LineWidth',0.8);
end
plot(1:5, kd,'ko-','LineWidth',2,'MarkerFaceColor','k');
ss = '';
for k=1:numel(TAGS)
    km=R(k).S.fm.ktr80(:)'; plot(1:5,km,'s-','Color',col(k),'LineWidth',1.4+0.4*(k>1));
    ss = [ss sprintf(' | %+.1f', M(k).ktr80_slope_pct)]; %#ok<AGROW>
end
xlabel('slack segment  (SL 2.04 \rightarrow 1.88 \mum)'); ylabel('ktr_{80} (s^{-1})'); xticks(1:5);
title(sprintf('F   ktr80 LADDER, %%/seg: data %+.1f%s', DATA.ktr80_slope_pct, ss));
legend([{'','','','pooled data'}, LBL],'Location','best','Box','off','NumColumns',2);

% --- G: rsK80 (strongest ATP carrier) ---
nexttile; hold on; box on; grid on
plot(1:5, PT.val(jk('rsK80'),:),'ko-','LineWidth',2,'MarkerFaceColor','k');
for k=1:numel(TAGS); plot(1:5, R(k).S.fm.rsK80(:)','s-','Color',col(k),'LineWidth',1.4+0.4*(k>1)); end
xline(4.5,':','Color',[0.5 0.5 0.5]); text(4.6, min(ylim)+0.05*range(ylim), 'seg 5 not scored','FontSize',8,'Color',[0.4 0.4 0.4]);
xlabel('slack segment'); ylabel('rsK_{80} (s^{-1})'); xticks(1:5);
title(sprintf('G   rsK80 re-stretch recovery rate (x0.38 at 2 mM): data %.1f | %s', ...
      DATA.rsK80_mean, strjoin(arrayfun(@(k) sprintf('%.1f', M(k).rsK80_mean), 1:numel(TAGS), 'uni', 0), ' | ')));
legend([{'pooled data'}, LBL],'Location','best','Box','off');

% --- I: where the objective sits ---
nexttile; hold on; box on; grid on
cats = categorical({'ATP-carrying','force-shape guards','everything else'});
cats = reordercats(cats, {'ATP-carrying','force-shape guards','everything else'});
Bm = zeros(numel(TAGS),3);
for k=1:numel(TAGS); Bm(k,:) = [M(k).atp_cost, M(k).guard_cost, M(k).cost-M(k).atp_cost-M(k).guard_cost]; end
bb = bar(cats, Bm', 0.8); for k=1:numel(TAGS); bb(k).FaceColor = col(k); bb(k).EdgeColor='none'; end
ylabel('weighted cost'); legend(LBL,'Box','off');
title(sprintf('I   objective: %s', strjoin(arrayfun(@(k) sprintf('%.4f', M(k).cost), 1:numel(TAGS), 'uni', 0), ' -> ')));

% --- H: staircase, VISUAL ONLY ---
nexttile([1 3]); hold on; box on; grid on
St=load(fullfile(root,'data',p1.stairs_velocitytableonfile)); Svt=St.velocitytable;
plot(St.datatable(:,1), St.datatable(:,3),'-','Color',CD,'LineWidth',1.3);
for k=1:numel(TAGS)
    Fn=R(k).S.stF/interp1(R(k).S.stT,R(k).S.stF,Svt(2,1),'linear','extrap');
    plot(R(k).S.stT, Fn,'-','Color',col(k),'LineWidth',1+0.4*(k>1));
end
iMove=find(abs(Svt(:,2))>1e-6,1); xlim([Svt(iMove,1)-0.08, min(Svt(end,1),St.datatable(end,1))]);
xlabel('t (s)'); ylabel('F (norm.)');
title('H   staircase -- visual check only (single prep, not pooled)');
legend([{'data'}, LBL],'Location','northwest','Box','off');

sgtitle(sprintf('%s   %s', FIGTAG, strjoin(arrayfun(@(k) sprintf('%s %.4f', TAGS{k}, M(k).cost), 1:numel(TAGS), 'uni', 0), '  ->  ')), ...
        'FontSize',12,'Interpreter','none');
exportgraphics(f1, fullfile(resdir,[FIGTAG '_protocols.png']),'Resolution',110);

%% ================= FIGURE: per-feature bars =================
lab = cellfun(@(s) regexprep(regexprep(s,'\|.*$',''),'\[.*$',''), R(1).S.fn, 'uni', 0);
c1 = R(1).S.fc(:); c2 = R(end).S.fc(:); d = c2 - c1;
f3=figure('Position',[40 40 1200 720],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
[~,ord]=sort(max(c1,c2),'descend'); ord = ord(max(c1(ord),c2(ord)) > 1e-6);
ylab = lab(ord);
for i=1:numel(ord)
    if ismember(ylab{i}, ATPC); ylab{i} = ['\bf' strrep(ylab{i},'_','\_') ' (ATP)'];
    elseif ismember(ylab{i}, GUARD); ylab{i} = [strrep(ylab{i},'_','\_') ' (guard)'];
    else; ylab{i} = strrep(ylab{i},'_','\_'); end
end
nexttile; b=barh(1:numel(ord),[c1(ord) c2(ord)],0.82);
b(1).FaceColor=CS; b(2).FaceColor=CO; b(1).EdgeColor='none'; b(2).EdgeColor='none';
set(gca,'YTick',1:numel(ord),'YTickLabel',ylab,'YDir','reverse','FontSize',9); grid on; xlabel('weighted cost');
legend({sprintf('%s %.4f',LBL{1},sum(c1)),sprintf('%s %.4f',LBL{end},sum(c2))},'Location','southeast','Box','off');
title('per-feature cost (bold = ATP-carrying)');
nexttile; [~,od]=sort(d(ord),'ascend'); hold on
for i=1:numel(od)
    j = ord(od(i)); cc=[0.8 0.2 0.2]; if d(j)<0; cc=[0.2 0.6 0.3]; end
    barh(i,d(j),0.82,'FaceColor',cc,'EdgeColor','none');
end
xline(0,'k-'); set(gca,'YTick',1:numel(od),'YTickLabel',ylab(od),'YDir','reverse','FontSize',9); grid on;
xlabel('\Delta cost (best - seed)'); title(sprintf('what moved (net %+.4f)', sum(d)));
exportgraphics(f3, fullfile(resdir,[FIGTAG '_bars.png']),'Resolution',120);

fprintf('\n%-14s %9s %9s %9s %9s %9s %9s\n','','cost','ktr80','ladder%','rsK80','rsT20ms','ktrP_k80');
fprintf('%-14s %9s %9.2f %9.2f %9.2f %9.2f %9.2f\n','DATA','',DATA.ktr80_mean,DATA.ktr80_slope_pct,DATA.rsK80_mean,DATA.rsT20_mean_ms,DATA.ktrP_k80);
for k=1:numel(TAGS)
    fprintf('%-14s %9.4f %9.2f %9.2f %9.2f %9.2f %9.2f\n', TAGS{k}, M(k).cost, M(k).ktr80_mean, ...
            M(k).ktr80_slope_pct, M(k).rsK80_mean, M(k).rsT20_mean_ms, M(k).ktrP_k80);
end
fprintf('\nwrote %s_protocols.png / %s_bars.png / render_%s.json\nDONE\n', FIGTAG, FIGTAG, FIGTAG);
