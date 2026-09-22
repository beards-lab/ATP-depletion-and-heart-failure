% PlotFVandKtr.m -- the two protocols that are NOT pooled.
%
% pooledTargetFeatures.m pools the SLACK block only. It says so:
%   "Everything outside the slack block (FV, ktr protocol, passive PS_*) --
%    those are not pooled yet, so they keep coming from their own sources."
%
% So:
%   FV   is scored against two HARDCODED LAB CURVES (Baker2022,
%        IsovelocityForFilip2021) from different preparations than the slack
%        data -- while in-preparation FV2 ramps exist on 04/03 and 04/10 and
%        are not used.
%   ktr  is fenced against 03/27 ALONE, which is both the strongest prep
%        (Fss22 80.1) and the fastest (k63 49.4 vs a 3-day mean of 43.2).
%
% This draws both, for the seed and the campaign optimum, with every source
% that exists so the pooling gap is visible.
%
%   matlab -batch "run('Analyses/RefitStrategy/PlotFVandKtr.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_fk_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, 5); end

TAGS = {'ktrfence_seed','ktrfence_opt'};
DAYS = {'03_27','04_03','04_10'};
AllV = [0 0.5 1 2 3 4 5 6];
D1 = [56.40,   51.8120, 37.4459, 17.8025, 11.4430, 6.2643, 3.2759, 2.2120];  % Baker2022
D2 = [67.3942, 60.9885, 42.7059, 18.8539, 12.5956, 6.4426, 3.5136, 1.8556];  % Isovel2021
SCORED = [0 0.5 1 2 4];
VFINE  = -[0, 0.25:0.25:4, 5, 6];
LoadData;

P = cell(1,2); FVf = cell(1,2); KT = cell(1,2);
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=0; p.BreakOnODEUnstable=0; p.MaxRunTime=600;
    p = getParams(p, p.g, true); P{k} = p;

    q = p; q.FV_velocities = VFINE; q = getParams(q,q.g,true);
    [~, outs] = runFVExperiment(q, ATP_c, Data_ATP);
    FVf{k} = arrayfun(@(o) o.Force(end), outs);

    for d = 1:3
        K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
        pk = p; pk.SL0 = 2*K.velocitytable(1,4); pk = getParams(pk,pk.g,true);
        pk.Velocity = K.velocitytable(:,2);
        [~, o] = evaluateModel(str2func(pk.modelFcn), K.velocitytable(:,1), pk);
        KT{k}(d).o = o; KT{k}(d).K = K;
        KT{k}(d).Fiso = interp1(o.t, o.Force, K.velocitytable(2,1),'linear','extrap');
    end
end

%% ---- numbers ----
fprintf('\n================ FORCE-VELOCITY (NOT pooled: two lab curves) ================\n');
fprintf('%6s %2s | %9s %9s | %9s %9s | %9s\n','v','sc','Baker22','Isovel21','seed','opt','avg tgt');
for i = 1:numel(VFINE)
    v = -VFINE(i); j = find(AllV==v,1);
    s1 = FVf{1}(i)/FVf{1}(1); s2 = FVf{2}(i)/FVf{2}(1);
    if isempty(j)
        fprintf('%6.2f %2s | %9s %9s | %9.3f %9.3f | %9s\n', v, '', '-', '-', s1, s2, '-');
    else
        a = mean([D1(j)/D1(1), D2(j)/D2(1)]);
        fprintf('%6.2f %2s | %9.3f %9.3f | %9.3f %9.3f | %9.3f\n', v, ...
            repmat('*',1,any(abs(SCORED-v)<1e-9)), D1(j)/D1(1), D2(j)/D2(1), s1, s2, a);
    end
end
fprintf('\nF_iso: seed %.2f  opt %.2f kPa   (Baker22 56.40, Isovel21 67.39, pooled slack Fref 69.35)\n', ...
        FVf{1}(1), FVf{2}(1));

fprintf('\n================ ktr PROTOCOL (NOT pooled: fenced on 03/27 alone) ================\n');
fprintf('%-8s %10s %10s %10s %12s %12s\n','day','k63 data','k63 seed','k63 opt','peak data','peak opt');
kd = zeros(1,3); ks = zeros(1,3); ko = zeros(1,3);
for d = 1:3
    K = KT{1}(d).K; V = K.velocitytable; DT = K.datatable;
    tA = V(end-4,1); tB = V(end-1,1);
    fd = extractKtrAttributes(DT(:,1), DT(:,3), V);
    pk_d = max(DT(DT(:,1)>=tA & DT(:,1)<=tB, 3));
    fs = extractKtrAttributes(KT{1}(d).o.t, KT{1}(d).o.Force/KT{1}(d).Fiso, V);
    fo = extractKtrAttributes(KT{2}(d).o.t, KT{2}(d).o.Force/KT{2}(d).Fiso, V);
    oo = KT{2}(d).o; m = oo.t>=tA & oo.t<=tB;
    pk_o = max(oo.Force(m))/KT{2}(d).Fiso;
    kd(d)=fd.ktrP_k63; ks(d)=fs.ktrP_k63; ko(d)=fo.ktrP_k63;
    fprintf('%-8s %10.1f %10.1f %10.1f %12.3f %12.3f\n', DAYS{d}, kd(d), ks(d), ko(d), pk_d, pk_o);
end
fprintf('%-8s %10.1f %10.1f %10.1f\n','MEAN', mean(kd), mean(ks), mean(ko));
fprintf('\nthe fence targets 03/27 alone (k63 %.1f). Against the 3-day MEAN (%.1f) the\n', kd(1), mean(kd));
fprintf('opt is +%.0f %%, not the +%.0f %% the fence sees.\n', ...
        100*(mean(ko)/mean(kd)-1), 100*(ko(1)/kd(1)-1));

%% ---- figure ----
f1=figure('Position',[20 20 1600 860],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
CS=[0.85 0.33 0.10]; CO=[0 0.45 0.74];

nexttile; hold on; box on; grid on
plot(AllV, D1/D1(1),'ks-','MarkerFaceColor','k','LineWidth',1.3,'MarkerSize',6);
plot(AllV, D2/D2(1),'k^--','LineWidth',1.1,'MarkerSize',6);
plot(-VFINE, FVf{1}/FVf{1}(1),'-','Color',CS,'LineWidth',1.6);
plot(-VFINE, FVf{2}/FVf{2}(1),'-','Color',CO,'LineWidth',1.8);
plot(SCORED, interp1(-VFINE,FVf{2}/FVf{2}(1),SCORED),'o','Color',[.3 .3 .3],'MarkerSize',11,'LineWidth',1.3);
xlabel('shortening velocity (ML/s)'); ylabel('F / F_{iso}'); xlim([0 6.2]); ylim([0 1.3]);
title('A   F-V, 0.25 ML/s grid — NOT pooled');
legend({'Baker2022','Isovel2021','seed','opt','the 5 scored v'},'Box','off');

nexttile; hold on; box on; grid on
ip = AllV>0;
plot(AllV(ip), D1(ip)/D1(1).*AllV(ip),'ks-','MarkerFaceColor','k','LineWidth',1.3);
plot(AllV(ip), D2(ip)/D2(1).*AllV(ip),'k^--','LineWidth',1.1);
plot(-VFINE, (FVf{1}/FVf{1}(1)).*(-VFINE),'-','Color',CS,'LineWidth',1.6);
plot(-VFINE, (FVf{2}/FVf{2}(1)).*(-VFINE),'-','Color',CO,'LineWidth',1.8);
xlabel('velocity (ML/s)'); ylabel('normalised power'); xlim([0 6.2]);
title('B   normalised power');

nexttile; hold on; box on; grid on
b = bar([kd(:) ks(:) ko(:)]);
b(1).FaceColor=[0.2 0.2 0.2]; b(2).FaceColor=CS; b(3).FaceColor=CO;
set(gca,'XTickLabel',{'03/27','04/03','04/10'});
ylabel('ktr protocol k_{63} (s^{-1})');
title('C   ktr rate by day — fence sees 03/27 only');
legend({'data','seed','opt'},'Box','off','Location','northeast');
yline(mean(kd),':k'); text(2.6, mean(kd)+1.5, '3-day mean', 'FontSize', 8);

for d = 1:3
    nexttile; hold on; box on; grid on
    K = KT{1}(d).K; V = K.velocitytable; DT = K.datatable; tK = V(end-1,1);
    m = DT(:,1) >= tK-0.020;
    plot((DT(m,1)-tK)*1e3, DT(m,3),'k-','LineWidth',1.5);
    for k = 1:2
        o = KT{k}(d).o; mm = o.t >= tK-0.020;
        plot((o.t(mm)-tK)*1e3, o.Force(mm)/KT{k}(d).Fiso,'-', ...
             'Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.1+0.4*(k==2));
    end
    xlim([-20 59]); ylim([-0.15 1.15]);
    xlabel('t since final release (ms)'); ylabel('F / F_{iso}');
    title(sprintf('%s   ktr %s   k_{63} data %.1f, opt %.1f', char('D'+d-1), ...
          strrep(DAYS{d},'_','/'), kd(d), ko(d)));
    if d==1; legend({'data','seed','opt'},'Location','southeast','Box','off'); end
end
exportgraphics(f1, fullfile(resdir,'fv_and_ktr.png'),'Resolution',120);
fprintf('\nwrote %s\n', fullfile(resdir,'fv_and_ktr.png'));
fprintf('\nDONE\n');
