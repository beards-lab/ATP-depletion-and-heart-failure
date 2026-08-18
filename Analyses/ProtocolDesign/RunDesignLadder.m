% RunDesignLadder.m — the acquisition ladder under two user constraints:
%   (1) a bracketing repeat in EVERY session (so no phi is ever transferred), and
%   (2) 4 mM squeezed into an existing session rather than bought separately.
%
% Answers three things RunDesignPower.m did not:
%   A. Where does 4 mM go inside a 4-activation session? (8-2-4-8 vs 8-4-2-8)
%   B. Does balancing the orders let you skip the rundown correction entirely?
%      strategy.md §2.4 says no (+5.1 % residual) for plain 8-2 / 2-8. Does that
%      still hold once the sessions are stacked and bracketed?
%   C. Cumulative statistical power, session by session, for 2/8 and for 4/8.
%
% Reproduce: RunDesignLadder -> results/design_ladder.png + design_ladder.mat

clear; close all;
here   = fileparts(mfilename('fullpath'));
resdir = fullfile(here,'results');
if ~exist(resdir,'dir'); mkdir(resdir); end

Phi = @(x) 0.5*(1+erf(x/sqrt(2)));          % normal CDF, no toolbox needed

%% ------------------------------------------------------------------ inputs
ATP   = [8 4 2];
slope = [0.560 0.782 1.092];                % kPa/s; 4 mM log-interpolated
Tact  = [10.4  11.0  11.6 ];                % s
F8    = 50;  R2 = 1.34;  R4 = sqrt(R2);     % 4 mM = geometric mean of 2 and 8
Fc    = F8 * [1 R4 R2];
phi_set = [0.45 0.581 0.889];  phi_c = 0.581;
CV_set  = [0.04 0.08 0.11];    CV_plan = 0.08;

D_of = @(a,phi) phi * slope(ATP==a) * Tact(ATP==a);

%% ---------------------------------------- A. where does 4 mM go in the session
cand = {[8 2 4 8], [8 4 2 8], [8 4 8], [8 2 8], [2 8 2]};
cname = {'8-2-4-8 (proposed)','8-4-2-8 (counter)','8-4-8','8-2-8','2-8-2'};
sig   = [NaN R4-1 R2-1];                    % true signal size per condition

fprintf('\n=== A. Placement of 4 mM inside a bracketed 4-run session (phi=%.3f) ===\n', phi_c);
fprintf('%-20s %-28s %14s %14s\n','sequence','damage at start of each run','4mM corr/sig','2mM corr/sig');
A = struct([]);
for k = 1:numel(cand)
    [f,b] = runSeq(cand{k}, phi_c, ATP, Fc, D_of);
    A(k).seq=cname{k}; A(k).o=cand{k}; A(k).f=f; A(k).bias=b;
    s = '';
    for j=1:numel(cand{k}); s=[s sprintf('%d:%4.1f%% ',cand{k}(j),100*f(j))]; end %#ok
    c4 = abs(b(2)-1)/sig(2);  c2 = abs(b(3)-1)/sig(3);
    A(k).c4=c4; A(k).c2=c2;
    fprintf('%-20s %-28s %14s %14s\n', cname{k}, s, ...
        tern(isnan(c4),'-',sprintf('%.2f',c4)), tern(isnan(c2),'-',sprintf('%.2f',c2)));
end

fprintf('\n  What the uncorrected 4 mM ratio would READ:\n');
for k = 1:2
    fprintf('    %-20s  R4_measured = %.3f   (truth %.3f)%s\n', A(k).seq, ...
        R4*A(k).bias(2)/1, R4, tern(R4*A(k).bias(2)<1,'   <-- SIGN FLIPS',''));
end

%% ------------------------------- B. does balancing remove the need to correct?
% Naive average = unweighted mean of log ratios across sessions, no correction.
pairs = { {[8 2],[2 8]},          '8-2  vs  2-8      (strategy.md §2.4 case)'; ...
          {[8 2 8],[2 8 2]},      '8-2-8 vs 2-8-2    (bracketed, 3 runs)'; ...
          {[8 4 2 8],[2 8 2]},    '8-4-2-8 vs 2-8-2  (YOUR stacked design)'; ...
          {[8 2 4 8],[2 8 2]},    '8-2-4-8 vs 2-8-2  (4 mM late)' };

fprintf('\n=== B. Residual bias of a NAIVE 1:1 average, no rundown correction ===\n');
fprintf('%-38s %10s %10s %10s\n','design pair','phi=0.45','phi=0.58','phi=0.89');
B = zeros(size(pairs,1), numel(phi_set));
for p = 1:size(pairs,1)
    for q = 1:numel(phi_set)
        [~,b1] = runSeq(pairs{p,1}{1}, phi_set(q), ATP, Fc, D_of);
        [~,b2] = runSeq(pairs{p,1}{2}, phi_set(q), ATP, Fc, D_of);
        B(p,q) = exp(mean(log([b1(3) b2(3)]))) - 1;      % 2 mM/8 mM ratio bias
    end
    fprintf('%-38s %9.1f%% %9.1f%% %9.1f%%\n', pairs{p,2}, 100*B(p,:));
end

%% ------------------------------------------------------ C. the session ladder
% Alternating orders, bracket in every session, sex balanced at every rung.
L = struct( ...
 'id',   {'S1','S2','S3','S4','S5','S6'}, ...
 'seq',  {[8 4 2 8],[2 8 2],[8 4 2 8],[2 8 2],[8 4 2 8],[2 8 2]}, ...
 'name', {'8-4-2-8','2-8-2','8-4-2-8','2-8-2','8-4-2-8','2-8-2'}, ...
 'sex',  {'F','F','M','M','M','F'});

% Existing holdings: 2/8 from three preps; no 4 mM.
n28_0 = 3; n48_0 = 0; sexM_0 = 2; sexF_0 = 1; desc_0 = 2; rev_0 = 1;

fprintf('\n=== C. The ladder (cumulative, CV = %.0f %%) ===\n', 100*CV_plan);
fprintf('%-4s %-9s %-4s %4s | %s\n','','session','sex','#act', ...
        'n(2/8) SEM   pwr_sex20   | n(4/8) SEM   pwr(R4~=1)  pwr(15% err) | M/F   desc/rev');
fprintf('%-4s %-9s %-4s %4s | %s\n','--','existing','-','-', ...
        sprintf('%4d %5.1f%% %8s   | %4d %5s %11s %13s | %d/%d   %d/%d', ...
        n28_0, 100*CV_plan/sqrt(n28_0), '-', n48_0, '-','-','-', sexM_0,sexF_0,desc_0,rev_0));

n28=n28_0; n48=n48_0; sM=sexM_0; sF=sexF_0; nd=desc_0; nr=rev_0; act=0;
Lad = struct([]);
for k = 1:numel(L)
    o = L(k).seq;
    n28 = n28 + 1;                                   % every session gives 2/8
    if any(o==4); n48 = n48 + 1; end
    if o(1)==8; nd = nd+1; else; nr = nr+1; end
    if L(k).sex=='M'; sM=sM+1; else; sF=sF+1; end
    act = act + numel(o);

    sem28 = CV_plan/sqrt(n28);  sem48 = CV_plan/sqrt(n48);
    % power to detect a 20 % sex difference in the 2/8 ratio (two groups)
    g = floor(min(sM,sF));
    pw_sex = tern(g<1, 0, Phi(log(1.20)/(CV_plan*sqrt(2/max(g,1))) - 1.96));
    % power that 4 mM differs from 8 mM at all, and that a 15 % model error shows
    pw_R4  = tern(n48<1, 0, Phi(log(R4)/sem48 - 1.96));
    pw_err = tern(n48<1, 0, Phi(log(1.15)/(sem48*sqrt(1.25)) - 1.96));

    Lad(k).id=L(k).id; Lad(k).n28=n28; Lad(k).n48=n48; Lad(k).sem28=sem28;
    Lad(k).sem48=sem48; Lad(k).pw_R4=pw_R4; Lad(k).pw_err=pw_err; Lad(k).act=act;

    fprintf('%-4s %-9s %-4s %4d | %4d %5.1f%% %8.0f%%   | %4d %5s %11s %13s | %d/%d   %d/%d\n', ...
        L(k).id, L(k).name, L(k).sex, numel(o), ...
        n28, 100*sem28, 100*pw_sex, ...
        n48, tern(n48<1,'-',sprintf('%.1f%%',100*sem48)), ...
        tern(n48<1,'-',sprintf('%.0f%%',100*pw_R4)), ...
        tern(n48<1,'-',sprintf('%.0f%%',100*pw_err)), sM,sF,nd,nr);
end
fprintf('\n  total new activations after S6: %d across %d preparations\n', act, numel(L));

% What curvature would need, for contrast.
d_curv = abs(log(1+(R2-1)*(1/4-1/8)/(1/2-1/8)) - log(R4));
fprintf('  curvature test (vs hyperbolic): needs n = %d at CV %.0f %% -- not on this ladder\n', ...
        ceil(1.25*(2.802*CV_plan/d_curv)^2), 100*CV_plan);

%% -------------------------------------------------------------- the figure
fig = figure('Position',[70 50 1480 820],'Color','w');
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl,'Bracketed, order-alternating ladder: where 4 mM goes and what each session buys', ...
      'FontWeight','bold','FontSize',13);
cols = lines(6);

% 1. placement of 4 mM
nexttile; hold on; grid on; box on;
cm = containers.Map({8,4,2},{cols(1,:),cols(2,:),cols(3,:)});
for k=1:2
    for j=1:numel(A(k).o)
        bar(k+0.22*(j-(numel(A(k).o)+1)/2), 100*A(k).f(j), 0.19, ...
            'FaceColor',cm(A(k).o(j)),'EdgeColor','none');
    end
end
set(gca,'XTick',1:2,'XTickLabel',{'8-2-4-8','8-4-2-8'});
ylabel('force lost when the run starts (%)'); title('Where the damage lands');
h=[patch(nan,nan,cols(1,:)),patch(nan,nan,cols(2,:)),patch(nan,nan,cols(3,:))];
legend(h,{'8 mM','4 mM','2 mM'},'Location','northwest','Box','off');

% 2. correction relative to signal
nexttile; hold on; grid on; box on;
bar([A(1).c4 A(2).c4; A(1).c2 A(2).c2]');
yline(1,'r--','correction = signal','LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'8-2-4-8','8-4-2-8'});
ylabel('|correction| / true signal'); title('Correction as a fraction of the effect');
legend({'4 mM','2 mM'},'Location','northwest','Box','off');

% 3. naive-average residual
nexttile; hold on; grid on; box on;
bar(100*B); set(gca,'XTick',1:size(B,1),'XTickLabel',{'8-2 / 2-8','8-2-8 / 2-8-2','8-4-2-8 / 2-8-2','8-2-4-8 / 2-8-2'}, ...
    'XTickLabelRotation',22);
yline(0,'k-'); ylabel('residual bias, uncorrected (%)');
title('Can you skip the correction?'); legend({'\phi=0.45','\phi=0.58','\phi=0.89'},'Box','off','Location','northeast');

% 4. SEM ladder
nexttile; hold on; grid on; box on;
plot(0:numel(Lad), [CV_plan/sqrt(n28_0) [Lad.sem28]]*100,'o-','LineWidth',1.8,'Color',cols(1,:));
plot(find([Lad.n48]>0), [Lad([Lad.n48]>0).sem48]*100,'s-','LineWidth',1.8,'Color',cols(2,:));
set(gca,'XTick',0:numel(Lad),'XTickLabel',[{'now'} {Lad.id}]);
ylabel('SEM of the ratio (%)'); title('Precision, session by session');
legend({'2 mM / 8 mM','4 mM / 8 mM'},'Box','off');

% 5. power ladder for 4 mM
nexttile; hold on; grid on; box on;
ii = find([Lad.n48]>0);
plot(ii, 100*[Lad(ii).pw_R4],'o-','LineWidth',1.8,'Color',cols(2,:));
plot(ii, 100*[Lad(ii).pw_err],'s-','LineWidth',1.8,'Color',cols(4,:));
yline(80,'k--','80 %'); ylim([0 100]);
set(gca,'XTick',1:numel(Lad),'XTickLabel',{Lad.id});
ylabel('power (%)'); title('4 mM: what each session buys');
legend({'4 mM differs from 8 mM','catches a 15 % model error'},'Location','southeast','Box','off');

% 6. cost
nexttile; hold on; grid on; box on;
bar(1:numel(Lad), [Lad.act], 'FaceColor',[.6 .6 .75],'EdgeColor','none');
set(gca,'XTick',1:numel(Lad),'XTickLabel',{Lad.id});
ylabel('cumulative activations'); xlabel('stop after ...');
title('Cost in activations');

exportgraphics(fig, fullfile(resdir,'design_ladder.png'),'Resolution',150);
save(fullfile(resdir,'design_ladder.mat'),'A','B','Lad','L','pairs','phi_set','CV_set');
fprintf('\nWrote %s\n', fullfile(resdir,'design_ladder.png'));

function y = tern(c,a,b); if c; y=a; else; y=b; end; end

% Damage suffered by each run of a sequence, and the bias on each ratio taken
% against the session's FIRST 8 mM run.
function [f, bias] = runSeq(o, phi, ATP, Fc, D_of)
    Dc = 0;  f = zeros(1,numel(o));
    for j = 1:numel(o)
        f(j) = Dc / Fc(ATP==o(j));
        Dc   = Dc + D_of(o(j), phi);
    end
    i8 = find(o==8,1,'first');
    bias = nan(1,numel(ATP));
    for c = 1:numel(ATP)
        j = find(o==ATP(c) & (1:numel(o))~=i8, 1, 'first');
        if ~isempty(j) && ~isempty(i8)
            bias(c) = (1-f(j)) / (1-f(i8));
        end
    end
end
