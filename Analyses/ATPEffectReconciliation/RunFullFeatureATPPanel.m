% RunFullFeatureATPPanel.m
% ONE view of the ATP effect across EVERY data-backed scored feature.
%
% RunATPReconciliation.m corrects six force features. This extends that to the
% whole feature list the objective actually scores, and shows three readings of
% the same data side by side:
%
%   RAW              ratio 2 mM / 8 mM, straight from the curated .mat files
%   +RUNDOWN         the same ratio after the ../RundownCorrection correction
%                    (force-dimension features only -- correcting rates makes
%                    the cross-prep spread WORSE, see that analysis)
%   SHAPE (rescaled) each force feature divided by ITS OWN run's steady force
%                    before the ratio is taken. This is the "normalise per prep"
%                    recommendation from ../SlackDataAnalysis: it removes the
%                    overall force change so what is left is ATP's effect on the
%                    SHAPE of the transient.
%
% Output -> results/atp_full_feature_panel.png  +  .mat

clear; close all; clc;
cd(fileparts(which('RunFullFeatureATPPanel')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end
D = '../../data/';

PHI = 0.452;                      % permanent fraction (../RundownCorrection)
zA  = [71.50 72.35];  zB = [77.70 78.30];

P = {  % day, earlier cond, earlier merged, earlier .mat, later cond, later .mat
 '03/27 M','8 mM','03 27 2026 M/02_Merged_8mM_Active.txt','protocol_03_27_2026_8mM_slack.mat', ...
           '2 mM','03 27 2026 M/03_Merged_2mM_Active.txt','protocol_03_27_2026_2mM_slack.mat'
 '04/03 F','8 mM','04 03 2026 F/02_Merged_8mM_Active.txt','protocol_04_03_2026_8mM_slack.mat', ...
           '2 mM','04 03 2026 F/03_Merged_2mM_Active.txt','protocol_04_03_2026_2mM_slack.mat'
 '04/10 M','2 mM','04 10 2026 Male 2-8/02_Merged_2mM_Active.txt','protocol_04_10_2026_2mM_slack.mat', ...
           '8 mM','04 10 2026 Male 2-8/03_Merged_8mM_Active.txt','protocol_04_10_2026_8mM_slack.mat'};
nP = size(P,1);

%% ---- feature registry --------------------------------------------------
% name, pretty label, class:  F = force-dimension (correctable + rescalable)
%                             R = rate            T = time / length
%                             Q = dimensionless quality
FEAT = {
 'A'                  ,'A'                 ,'F'
 'Am'                 ,'Am'                ,'F'
 'steady'             ,'steady'            ,'F'
 'peak1_y'            ,'peak1_y'           ,'F'
 'vall_y'             ,'vall_y'            ,'F'
 'peak2'              ,'peak2'             ,'F'
 'vall2_dy'           ,'vall2_dy'          ,'F'
 'rsA'                ,'rsA'               ,'F'
 'ktr'                ,'k_{tr}'            ,'R'
 'ktr2_k'             ,'ktr2\_k'           ,'R'
 'rsK'                ,'rsK'               ,'R'
 'restretchSlopeStart','restretchSlope'    ,'R'
 't0'                 ,'t0'                ,'T'
 'peak1_t'            ,'peak1\_t'          ,'T'
 'peak1_dSL'          ,'peak1\_dSL'        ,'T'
 'vall_t'             ,'vall\_t'           ,'T'
 'vall2_t'            ,'vall2\_t'          ,'T'
 'rsT0'               ,'rsT0'              ,'T'
 'ovrsht_dy'          ,'ovrsht\_dy'        ,'Q'
 'ktr_rmse'           ,'ktr\_rmse'         ,'Q'
 'rsR2'               ,'rsR2'              ,'Q' };
nF = size(FEAT,1);

%% ---- rundown dose per prep --------------------------------------------
frac = nan(1,nP); ord = cell(1,nP);
fprintf('%-9s %-14s | %10s %9s %11s\n','day','order','loss (kPa)','F_later','frac');
for i = 1:nP
    Me = readmatrix([D P{i,3}]);  te = Me(:,1); Fe = Me(:,3);
    m  = (te>=zA(1)&te<=zA(2)) | (te>=zB(1)&te<=zB(2));  p = polyfit(te(m),Fe(m),1);
    Fs = movmean(Fe,51); thr = 0.2*prctile(Fs,99);
    on = find(Fs>thr,1); off = find(Fs>thr,1,'last');  Ta = te(off)-te(on);
    loss = PHI*abs(p(1))*Ta;
    Ml = readmatrix([D P{i,6}]);  tl = Ml(:,1);
    Fl = mean(Ml(tl>=zA(1)&tl<=zA(2), 3),'omitnan');
    frac(i) = loss/Fl;  ord{i} = sprintf('%s->%s', P{i,2}, P{i,5});
    fprintf('%-9s %-14s | %10.2f %9.2f %10.1f%%\n', P{i,1}, ord{i}, loss, Fl, 100*frac(i));
end

%% ---- ratios ------------------------------------------------------------
raw = nan(nP,nF); cor = nan(nP,nF); shp = nan(nP,nF);
for i = 1:nP
    Se = load([D P{i,4}]);  Sl = load([D P{i,7}]);
    fe = Se.features_data;  fl = Sl.features_data;
    is8first = strcmp(P{i,2},'8 mM');
    se = mean(fe.steady,'omitnan');  sl_ = mean(fl.steady,'omitnan');
    for j = 1:nF
        nm = FEAT{j,1}; cls = FEAT{j,3};
        if ~isfield(fe,nm) || ~isfield(fl,nm), continue; end
        ve = mean(double(fe.(nm)),'omitnan');
        vl = mean(double(fl.(nm)),'omitnan');
        if ~isfinite(ve) || ~isfinite(vl) || ve == 0, continue; end
        if is8first, r = vl/ve; else, r = ve/vl; end
        raw(i,j) = r;

        if cls == 'F'
            % undo the depression of whichever run was recorded SECOND
            if is8first, cor(i,j) = r*(1+frac(i)); else, cor(i,j) = r/(1+frac(i)); end
            % shape: normalise each run by its own steady force first
            if is8first, shp(i,j) = (vl/sl_)/(ve/se); else, shp(i,j) = (ve/se)/(vl/sl_); end
        else
            cor(i,j) = r;         % rates/times are NOT rundown-corrected
            shp(i,j) = r;         % and are already scale-free
        end
    end
end

%% ---- plot --------------------------------------------------------------
cls = string(FEAT(:,3))';
[~,ordF] = sort(categorical(cls,{'F','R','T','Q'}));   % group by class
lab = string(FEAT(ordF,2));
RAW = raw(:,ordF); COR = cor(:,ordF); SHP = shp(:,ordF);

fig = figure('Units','pixels','Position',[20 20 1900 950],'Color','w');
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% --- top: per-prep raw vs corrected vs shape
nexttile; hold on; box on; grid on
x = 1:nF;
cRaw = [0.55 0.55 0.55]; cCor = [0.72 0.14 0.18]; cShp = [0.12 0.42 0.49];
hR = plot(x, mean(RAW,1,'omitnan'), 'o-','Color',cRaw,'MarkerFaceColor',cRaw,'LineWidth',2,'MarkerSize',8);
hC = plot(x, mean(COR,1,'omitnan'), 's-','Color',cCor,'MarkerFaceColor',cCor,'LineWidth',2.4,'MarkerSize',9);
hS = plot(x, mean(SHP,1,'omitnan'), '^-','Color',cShp,'MarkerFaceColor',cShp,'LineWidth',2,'MarkerSize',8);
plot(x, RAW','.','Color',[cRaw 0.9],'MarkerSize',12,'HandleVisibility','off');
plot(x, COR','.','Color',[cCor 0.6],'MarkerSize',12,'HandleVisibility','off');
yline(1,'k--','no ATP effect','LineWidth',1.4,'FontSize',12,'HandleVisibility','off', ...
      'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
set(gca,'YScale','log','XTick',x,'XTickLabel',lab,'XTickLabelRotation',45,'FontSize',12);
xlim([0.5 nF+0.5]); ylabel('2 mM / 8 mM','FontSize',13);
legend([hR hC hS],{'raw','+ rundown correction','shape (rescaled by own steady)'}, ...
       'Location','northwest','FontSize',12);
title({'The ATP effect on every data-backed scored feature (mean of 3 preps; dots = individual preps)', ...
       'Once each run is rescaled by its own steady force, EVERY force feature returns to 1.0 — low ATP scales the trace, it does not reshape it'}, ...
      'FontSize',14,'FontWeight','bold');
% class dividers
nb = cumsum([sum(cls=='F') sum(cls=='R') sum(cls=='T')]);
for b = nb, xline(b+0.5,'-','Color',[0.75 0.75 0.75],'LineWidth',1.2,'HandleVisibility','off'); end
text(sum(cls=='F')/2+0.5, 3.4,'force-dimension','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
text(nb(1)+sum(cls=='R')/2+0.5, 3.4,'rates','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
text(nb(2)+sum(cls=='T')/2+0.5, 3.4,'times / length','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
text(nb(3)+sum(cls=='Q')/2+0.5, 3.4,'quality','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');

% --- bottom: cross-prep spread, raw vs corrected
nexttile; hold on; box on; grid on
cvR = 100*std(RAW,0,1,'omitnan')./abs(mean(RAW,1,'omitnan'));
cvC = 100*std(COR,0,1,'omitnan')./abs(mean(COR,1,'omitnan'));
b = bar(x, [cvR; cvC]', 'grouped');
b(1).FaceColor = cRaw; b(2).FaceColor = cCor;
set(gca,'XTick',x,'XTickLabel',lab,'XTickLabelRotation',45,'FontSize',12);
xlim([0.5 nF+0.5]); ylabel('cross-prep CV of the ratio (%)','FontSize',13);
legend({'raw','+ rundown correction'},'Location','northwest','FontSize',12);
title('Does the correction make the preps agree? (lower is better)','FontSize',14,'FontWeight','bold');
for bb = nb, xline(bb+0.5,'-','Color',[0.75 0.75 0.75],'LineWidth',1.2,'HandleVisibility','off'); end

exportgraphics(fig, fullfile(resDir,'atp_full_feature_panel.png'), 'Resolution', 170);
save(fullfile(resDir,'atp_full_feature_panel.mat'),'FEAT','raw','cor','shp','frac','P');

%% ---- console table -----------------------------------------------------
fprintf('\n%-20s %8s %8s %8s %8s %8s\n','feature','raw','+rundown','shape','CV raw','CV cor');
for j = 1:nF
    k = ordF(j);
    fprintf('%-20s %8.3f %8.3f %8.3f %7.1f%% %7.1f%%\n', FEAT{k,1}, ...
        mean(raw(:,k),'omitnan'), mean(cor(:,k),'omitnan'), mean(shp(:,k),'omitnan'), ...
        cvR(j), cvC(j));
end
disp('FULL FEATURE PANEL DONE');
