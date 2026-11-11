% RunPoolingStrategies.m
% Six ways to pool 8 mM and 2 mM records into ONE fitting target, compared on
% bias, precision, and what each has to assume.
%
% THE GOAL
% Identification is expensive, so we want a single 8 mM target and a single 2 mM
% target, pooled over preparations, such that the DIFFERENCE between them is the
% ATP effect and not a rundown artefact.
%
% THE PROPOSAL THAT MOTIVATED THIS (P2 below)
% "Average three 8 mM records from 8->2 with three from 2->8, and the same for
%  2 mM. Both averages are then equally damaged, so the damage cancels."
% Right in spirit -- but it does not cancel, because the two orders do NOT
% inflict equal damage:
%     8->2 : the 8 mM run damages the 2 mM run     f ~ 0.05
%     2->8 : the 2 mM run damages the 8 mM run     f ~ 0.15   (x3)
% Averaging a 5 % deficit against a 15 % deficit leaves ~5 %.
%
% THE FIX (P6) : both orders together let you SOLVE for the damage instead of
% assuming it. With f_i = phi * g_i and g_i fully measured per prep,
%     8->2 :  R_i = a / (1 + phi*g_i)
%     2->8 :  R_j = a * (1 + phi*g_j)
% is 2 unknowns (a, phi) in n equations. phi no longer has to be transferred
% from 03/27's bracket -- the design itself identifies it.
%
% Outputs -> results/pooling_strategies.png, .mat

clear; close all; clc;
cd(fileparts(which('RunPoolingStrategies')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end

%% ── Measured inputs, per preparation ───────────────────────────────────────
% R      : own-peak-referenced 2 mM / 8 mM force ratio, as measured
% g      : |slope_earlier| * T_dec_earlier / F_later   (phi factored out)
% is82   : true if the prep ran 8 mM -> 2 mM
prep = struct( ...
 'day',  {'03/27 M','04/03 F','04/10 M'}, ...
 'R',    {1.269, 1.233, 1.649}, ...
 'g',    {0.4592*9.50/83.59, 0.7570*9.58/63.54, 1.0169*10.67/43.03}, ...
 'is82', {true, true, false});
R = [prep.R];  g = [prep.g];  is82 = [prep.is82];
sdPrep = 0.053;                      % residual between-prep SD after correction
PHI_BR = 0.581;                      % phi from the 03/27 bracket (own-peak footing)

fprintf('%-9s %-7s %8s %8s %10s\n','prep','order','R meas','g','phi*g');
for i=1:3
    if is82(i); o='8->2'; else; o='2->8'; end
    fprintf('%-9s %-7s %8.3f %8.4f %10.3f\n', prep(i).day, o, R(i), g(i), PHI_BR*g(i));
end

%% ── P1 naive: ignore order, average all ratios ─────────────────────────────
P(1).name = 'P1 naive (ignore order)';
P(1).a    = mean(R);
P(1).ass  = 'none';

%% ── P2 balanced-count pooling (the proposal) ───────────────────────────────
% Average the ratios within each order group, then combine the two groups.
r82 = mean(R(is82));  r28 = mean(R(~is82));
P(2).name = 'P2 balanced-count (arith)';
P(2).a    = (r82 + r28)/2;
P(2).ass  = 'equal n both orders';
P(3).name = 'P3 balanced-count (geom)';
P(3).a    = sqrt(r82*r28);
P(3).ass  = 'equal n both orders';

%% ── P4 first-activations only (unpaired, unbiased) ─────────────────────────
% Every prep's FIRST activation is undamaged. Take undamaged 8 mM from the 8->2
% preps and undamaged 2 mM from the 2->8 preps. Unbiased -- but UNPAIRED, so the
% ~20 % between-prep force-scale variation no longer cancels.
P(4).name = 'P4 first activations only';
P(4).a    = NaN;                      % value needs absolute forces; bias is 0 by construction
P(4).ass  = 'unbiased, but unpaired';

%% ── P5 per-prep correction with a TRANSFERRED phi (current practice) ───────
corr5 = R;
corr5(is82)  = R(is82) .*(1 + PHI_BR*g(is82));
corr5(~is82) = R(~is82)./(1 + PHI_BR*g(~is82));
P(5).name = 'P5 correct with phi from bracket';
P(5).a    = mean(corr5);
P(5).ass  = 'phi transfers between preps';

%% ── P6 two-order solve: fit a AND phi from the design itself ───────────────
% Least squares on logs:  ln R_i = ln a -+ ln(1 + phi*g_i)
obj = @(v) sum(( log(R(:)) - ( v(1) - (2*is82(:)-1).*log(1+v(2)*g(:)) ) ).^2);
v0  = [log(mean(R)), PHI_BR];
opt = optimset('Display','off','TolX',1e-10,'TolFun',1e-12);
vh  = fminsearch(obj, v0, opt);
aHat = exp(vh(1));  phiHat = vh(2);
P(6).name = 'P6 two-order solve (a AND phi)';
P(6).a    = aHat;
P(6).ass  = 'both orders present; phi NOT assumed';

fprintf('\n=== P6: solved from the design ===\n');
fprintf('  a   = %.3f\n  phi = %.3f   (bracket gave %.3f)\n', aHat, phiHat, PHI_BR);
pred = aHat * (1 + phiHat*g).^(1-2*is82);
fprintf('  predicted vs measured R: %s  vs  %s\n', mat2str(round(pred,3)), mat2str(round(R,3)));
fprintf('  residuals: %s   (these are BIOLOGICAL scatter, n=3 so it is fragile)\n', ...
        mat2str(round(pred-R,3)));

%% ── Bias of each strategy, evaluated against a known truth ─────────────────
% Simulate: true a, per-prep g drawn near the measured ones, biological noise.
aTrue = 1.344;  phiTrue = 0.581;
rng(7);  nRep = 4000;  nEach = 3;
g82pop = [0.052 0.114];  g28pop = 0.252;      % measured g values by order
bias = nan(1,6);  sem = nan(1,6);
for s = [1 2 3 5 6]
    est = nan(1,nRep);
    for m = 1:nRep
        g8 = g82pop(randi(2,1,nEach)).*(1+0.15*randn(1,nEach));
        g2 = g28pop      .*(1+0.15*randn(1,nEach));
        R8 = aTrue./(1+phiTrue*g8) .* (1+sdPrep/aTrue*randn(1,nEach));
        R2 = aTrue.*(1+phiTrue*g2) .* (1+sdPrep/aTrue*randn(1,nEach));
        switch s
            case 1, est(m) = mean([R8 R2]);
            case 2, est(m) = (mean(R8)+mean(R2))/2;
            case 3, est(m) = sqrt(mean(R8)*mean(R2));
            case 5, est(m) = mean([R8.*(1+PHI_BR*g8), R2./(1+PHI_BR*g2)]);
            case 6
                RR=[R8 R2]; gg=[g8 g2]; ii=[true(1,nEach) false(1,nEach)];
                o2 = @(v) sum(( log(RR(:)) - ( v(1) - (2*ii(:)-1).*log(1+v(2)*gg(:)) ) ).^2);
                vv = fminsearch(o2,[log(mean(RR)) 0.5],opt);  est(m)=exp(vv(1));
        end
    end
    bias(s) = 100*(mean(est)/aTrue - 1);  sem(s) = std(est);
end
bias(4) = 0;  sem(4) = sqrt(2)*0.20*aTrue/sqrt(nEach);   % unpaired: 20 % prep scale

fprintf('\n=== STRATEGY COMPARISON (simulated, 3 preps per order, truth a=%.3f) ===\n', aTrue);
fprintf('%-34s %9s %9s %9s   %s\n','strategy','bias','SEM','RMSE','assumes');
for s = 1:6
    rmse = sqrt((bias(s)/100*aTrue)^2 + sem(s)^2);
    fprintf('%-34s %8.1f%% %9.3f %9.3f   %s\n', P(s).name, bias(s), sem(s), rmse, P(s).ass);
end

fprintf('\n  P2/P3 (the balanced-count proposal) leave ~%.0f%% because the two orders do\n', bias(2));
fprintf('  NOT inflict equal damage: f ~ 0.05 (8->2) vs ~0.15 (2->8).\n');
fprintf('  P4 is unbiased but UNPAIRED -- the 20%% prep-to-prep force scale stops\n');
fprintf('  cancelling, and its SEM (%.3f) is %.0fx worse than the paired strategies.\n', sem(4), sem(4)/sem(6));
fprintf('  P6 is unbiased AND paired: it uses both orders to SOLVE for phi.\n');

%% ── Figure ─────────────────────────────────────────────────────────────────
fig = figure(1030); clf; set(fig,'Position',[20 20 1450 700]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nm = {P.name};

nexttile([1 2]); hold on; box on;
b = bar([bias; 100*sem/aTrue]'); b(1).FaceColor=[.85 .33 .1]; b(2).FaceColor=[.3 .5 .8];
set(gca,'XTick',1:6,'XTickLabel',{'P1 naive','P2 bal-arith','P3 bal-geom','P4 first-act','P5 phi-transfer','P6 two-order'});
xtickangle(20); ylabel('% of the true effect'); yline(0,'k-');
title('Bias and random error by pooling strategy');
legend({'bias','SEM'},'Location','northwest','FontSize',8);

nexttile; hold on; box on;
rmse = sqrt((bias/100*aTrue).^2 + sem.^2);
b = bar(rmse); b.FaceColor='flat';
b.CData = repmat([.6 .6 .6],6,1); [~,ib]=min(rmse); b.CData(ib,:)=[.2 .6 .3];
set(gca,'XTick',1:6,'XTickLabel',{'P1','P2','P3','P4','P5','P6'});
ylabel('RMSE'); title('Total error (lower is better)');

nexttile; hold on; box on;
xq = linspace(0,0.35,100);
plot(100*xq, 100*(1./(1+PHI_BR*xq)-1),'-','Color',[0 .45 .74],'LineWidth',2);
plot(100*xq, 100*((1+PHI_BR*xq)-1),'-','Color',[.85 .33 .1],'LineWidth',2);
for i=1:3
    c = [0 .45 .74]; if ~is82(i); c=[.85 .33 .1]; end
    plot(100*g(i), 100*((1+PHI_BR*g(i))^(1-2*is82(i))-1),'o','Color',c,'MarkerFaceColor',c,'MarkerSize',10);
    text(100*g(i)+0.8, 100*((1+PHI_BR*g(i))^(1-2*is82(i))-1), prep(i).day,'FontSize',7);
end
yline(0,'k:'); xlabel('g = |slope|\cdotT / F_{later}   (%)'); ylabel('bias (%)');
title('Why the orders do not cancel');
legend({'8\rightarrow2','2\rightarrow8'},'Location','northwest','FontSize',7);

nexttile; hold on; box on;
nn = 1:6;
semA = nan(size(nn));
for k = nn
    est = nan(1,1500);
    for m=1:1500
        g8 = g82pop(randi(2,1,k)).*(1+0.15*randn(1,k)); g2 = g28pop.*(1+0.15*randn(1,k));
        R8 = aTrue./(1+phiTrue*g8).*(1+sdPrep/aTrue*randn(1,k));
        R2 = aTrue.*(1+phiTrue*g2).*(1+sdPrep/aTrue*randn(1,k));
        RR=[R8 R2]; gg=[g8 g2]; ii=[true(1,k) false(1,k)];
        o2=@(v) sum((log(RR(:))-(v(1)-(2*ii(:)-1).*log(1+v(2)*gg(:)))).^2);
        vv=fminsearch(o2,[log(mean(RR)) 0.5],opt); est(m)=exp(vv(1));
    end
    semA(k)=std(est);
end
plot(nn, semA,'ko-','LineWidth',2,'MarkerFaceColor','k');
xlabel('preps PER ORDER'); ylabel('SEM of a (P6)');
title('P6 precision vs design size');
for k=[1 3 6]; text(k, semA(k)+0.004, sprintf('%.3f',semA(k)),'FontSize',8); end

nexttile; hold on; box on;
txt = {'\bfRECOMMENDATION\rm','', ...
       '\bfP6\rm: solve a and \phi jointly from', ...
       'both orders. Needs no transferred \phi.','', ...
       'Build ONE 8 mM and ONE 2 mM target:', ...
       ' 1. per prep, ratio (cancels scale)', ...
       ' 2. correct the SECOND record by', ...
       '    the fitted \phi\cdotg', ...
       ' 3. normalise, then average','', ...
       '\bfDesign\rm: 3 preps per order.', ...
       'A sandwich (8-2-8) makes one prep', ...
       'self-correcting and pins \phi directly.'};
text(0.02,0.97,txt,'VerticalAlignment','top','FontSize',9.5); axis off;

sgtitle('Pooling strategies — how to build one 8 mM and one 2 mM target');
exportgraphics(fig, fullfile(resDir,'pooling_strategies.png'), 'Resolution', 150);
save(fullfile(resDir,'pooling_strategies.mat'),'P','bias','sem','aHat','phiHat','R','g','is82');
fprintf('\nSaved %s\n', fullfile(resDir,'pooling_strategies.png'));
