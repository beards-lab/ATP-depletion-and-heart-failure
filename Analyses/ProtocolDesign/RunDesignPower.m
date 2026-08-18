% RunDesignPower.m — how many sessions, in what order, with which repeats, and
% what a 4 mM arm would cost.
%
% Extends Analyses/RundownCorrection/strategy.md, which settled ORDER (majority
% 8->2, keep >=1 reversed to validate) and POOLING (P6 two-order solve, fit
% ratios with a free per-prep scale) but did not address:
%     (a) SEX  — nowhere in the repo is sex treated as a variable; it exists
%                only as a day label (03/27 M, 04/03 F, 04/10 M).
%     (b) N    — "3 preps per order" is asserted as the design target without a
%                precision requirement behind it.
%     (c) 4 mM — a third dose has never been costed against rundown.
%
% Everything below is arithmetic on measured quantities. Sources:
%   slope, T_act, phi, F        RundownCorrection/strategy.md §2.3, methods.md
%   between-prep SD of ratio    RundownCorrection/strategy.md §5.5  (s = 0.053)
%   ATP effect anchor           ATPEffectReconciliation/conclusions.md (x1.34)
%
% Reproduce:  RunDesignPower  ->  results/design_power.png + design_power.mat

clear; close all;
here    = fileparts(mfilename('fullpath'));
resdir  = fullfile(here,'results');
if ~exist(resdir,'dir'); mkdir(resdir); end

%% ------------------------------------------------------------------ inputs
% Measured within-run decay (kPa/s) and activated time (s) per condition.
% 8 and 2 mM are measured; 4 mM is INTERPOLATED three ways and is the single
% untested input in this whole analysis (see 'dose law' below).
ATP     = [8 4 2];
slope   = [0.560 NaN 1.092];        % kPa/s, strategy.md §2.3
Tact    = [10.4  NaN 11.6 ];        % s,     strategy.md §2.3

% phi = fraction of each run's within-run decline that is PERMANENT.
% The repo quotes three values from three routes; the spread is ~x2 and is the
% largest single uncertainty in any rundown correction.
phi_set = [0.45 0.581 0.889];       % [bracket-consistent, strategy §2.3, P6 fit]
phi_lbl = {'0.45 (03/27 bracket)','0.58 (strategy §2.3)','0.89 (P6 joint fit)'};
phi_c   = 0.581;                    % central value used for the design tables

% Reference fresh 8 mM plateau force. Absolute force spans 37.5-62.2 kPa across
% preps (CV 20 %); f scales as 1/F, so this is a representative, not a constant.
F8      = 50;

% The measured ATP effect on force, rundown-corrected (ATPEffectReconciliation).
R2      = 1.34;

% Between-prep SD of the corrected ratio, in ratio units (strategy §5.5).
% CV = 0.053/1.34 = 4.0 %. That estimate has 2 df, so it is planned against
% three values: optimistic (as measured), realistic, pessimistic (the raw
% pooled force CV of 11 % from ATPEffectReconciliation).
CV_set  = [0.04 0.08 0.11];
CV_lbl  = {'4 % (corrected, n=3)','8 % (planning)','11 % (raw pooled)'};
CV_plan = 0.08;

%% ------------------------------------------- 4 mM: three dose-response laws
% 4 mM is the GEOMETRIC mean of 2 and 8 (sqrt(2*8) = 4), so a log-linear dose
% dependence predicts R4 = sqrt(R2) exactly. That makes 4 mM the natural test
% point for log-linearity, and it is why the three laws separate so little.
law.name = {'hyperbolic (linear in 1/ATP)','log-linear (linear in ln ATP)','linear in ATP'};
law.R4   = [ 1 + (R2-1)*(1/4-1/8)/(1/2-1/8), ...
             sqrt(R2), ...
             1 + (R2-1)*(8-4)/(8-2) ];
% Same three laws applied to the decay slope, to fill the NaN above.
law.s4   = [ interp1([1/8 1/2], slope([1 3]), 1/4), ...
             exp(interp1(log([8 2]), log(slope([1 3])), log(4))), ...
             interp1([8 2], slope([1 3]), 4) ];

slope(2) = law.s4(2);                       % central: log-linear
Tact(2)  = interp1([8 2], Tact([1 3]), 4);
R_dose   = [1 law.R4(2) R2];                % force ratio vs 8 mM, central law
F_cond   = F8 * R_dose;                     % fresh plateau force per condition

fprintf('\n=== 4 mM interpolation (the one untested input) ===\n');
for i = 1:3
    fprintf('  %-32s  slope_4 = %.3f kPa/s   R_4 = %.3f\n', law.name{i}, law.s4(i), law.R4(i));
end
fprintf('  spread in R_4 = %.3f  (%.1f %% of the ratio)\n', ...
        max(law.R4)-min(law.R4), 100*(max(law.R4)-min(law.R4))/mean(law.R4));

%% ------------------------------------------------------- candidate sequences
% Each entry: the order of ATP conditions in the session. Relax precedes and
% PNB+Mava follows every one of them; those are not activations of interest.
seq.order = { [8 2], [2 8], [8 2 8], [8 4], [4 8], [8 4 8], [8 4 2], [8 4 2 8], [8 2 4] };
seq.name  = { '8-2','2-8','8-2-8','8-4','4-8','8-4-8','8-4-2','8-4-2-8','8-2-4' };
nseq      = numel(seq.order);

% Permanent damage inflicted BY one run of each condition, at the central phi.
D_of = @(a,phi) phi * slope(ATP==a) * Tact(ATP==a);

S = struct([]);
for k = 1:nseq
    o = seq.order{k};  Dcum = 0;  f = zeros(1,numel(o));  Fmeas = zeros(1,numel(o));
    for j = 1:numel(o)
        Ftrue    = F_cond(ATP==o(j));
        f(j)     = Dcum / Ftrue;            % fraction of TRUE force already lost
        Fmeas(j) = Ftrue - Dcum;
        Dcum     = Dcum + D_of(o(j), phi_c);
    end
    S(k).name = seq.name{k};  S(k).order = o;  S(k).f = f;  S(k).Fmeas = Fmeas;
    S(k).nact = numel(o);     S(k).Dtot = Dcum;

    % Bias on each non-8 ratio measured in this session, using the FIRST 8 mM
    % run as the reference:  R_meas/R_true = (1-f_c)/(1-f_8).
    i8 = find(o==8, 1, 'first');
    S(k).bias = struct('cond',{},'val',{});
    if ~isempty(i8)
        for j = 1:numel(o)
            if o(j) ~= 8
                b = (1-f(j)) / (1-f(i8));
                S(k).bias(end+1) = struct('cond', o(j), 'val', b);
            end
        end
    else
        S(k).bias(end+1) = struct('cond', NaN, 'val', NaN);
    end
end

fprintf('\n=== Rundown budget per session (phi = %.3f) ===\n', phi_c);
fprintf('%-9s %5s %8s   %s\n','sequence','#act','tot kPa','per-run damage already suffered (%% of true force)');
for k = 1:nseq
    fprintf('%-9s %5d %8.1f   ', S(k).name, S(k).nact, S(k).Dtot);
    for j = 1:numel(S(k).order)
        fprintf('%d mM:%5.1f%%  ', S(k).order(j), 100*S(k).f(j));
    end
    fprintf('\n');
end

fprintf('\n=== Uncorrected ratio bias, first 8 mM run as reference ===\n');
for k = 1:nseq
    if isnan(S(k).bias(1).val); continue; end
    fprintf('%-9s  ', S(k).name);
    for b = S(k).bias
        fprintf('%d mM/8 mM: %+6.1f%%   ', b.cond, 100*(b.val-1));
    end
    fprintf('\n');
end

% Sensitivity of the worst case to phi.
fprintf('\n=== phi sensitivity of the 2 mM correction (how much is being assumed) ===\n');
fprintf('%-9s  %s\n','sequence', sprintf('%-26s', phi_lbl{:}));
for k = [1 7]
    fprintf('%-9s  ', S(k).name);
    for p = phi_set
        o = S(k).order;  Dc = 0;  fv = 0;
        for j = 1:numel(o)
            if o(j)==2; fv = Dc / F_cond(ATP==2); end
            Dc = Dc + D_of(o(j), p);
        end
        fprintf('%-26s', sprintf('f(2 mM) = %.1f %%', 100*fv));
    end
    fprintf('\n');
end

%% ----------------------------------------------------------- sample sizes
% Work in LOG ratio units: sd(ln R) ~= CV, so every requirement below is a
% distance in log units divided by a CV.
z = 2.802;                                  % z_{0.025} + z_{0.20}: 95 %, 80 % power

% (1) Place the 4 mM point to a stated relative precision (SEM as a fraction).
prec_target = [0.10 0.05 0.035];
n_place = zeros(numel(CV_set), numel(prec_target));
for a = 1:numel(CV_set)
    n_place(a,:) = ceil((CV_set(a)./prec_target).^2);
end

% (2) Discriminate the dose-response laws. The test is on ln R4 - 0.5*ln R2
% (zero under log-linearity). With n preps contributing to each of R4 and R2,
% var = CV^2/n + 0.25*CV^2/n = 1.25*CV^2/n.
delta_law = abs(log(law.R4) - 0.5*log(R2));         % separation from log-linear
delta_hard = min(delta_law(delta_law > 1e-9));      % the harder of the two
n_law = zeros(numel(CV_set), 2);
d2 = delta_law([1 3]);
for a = 1:numel(CV_set)
    n_law(a,:) = ceil(1.25 * (z*CV_set(a)./d2).^2);
end

% (3) Detect a gross model prediction failure at 4 mM (model fitted to 8+2,
% predicts R4, is wrong by X %).
model_err = [0.25 0.15 0.10];
n_model = zeros(numel(CV_set), numel(model_err));
for a = 1:numel(CV_set)
    n_model(a,:) = ceil(1.25 * (z*CV_set(a)./log(1+model_err)).^2);
end

% (4) Sex: two-group comparison of the ATP ratio, n per group.
sex_eff = [0.30 0.20 0.10];
n_sex = zeros(numel(CV_set), numel(sex_eff));
for a = 1:numel(CV_set)
    n_sex(a,:) = ceil(2 * (z*CV_set(a)./log(1+sex_eff)).^2);
end

pr = @(hdr, cols, M, fmtc) deal();
show = @(hdr, colhdr, M) printTable(hdr, colhdr, CV_lbl, M);

fprintf('\n=== (1) n needed to PLACE the 4 mM ratio (SEM as %% of the ratio) ===\n');
show('CV', arrayfun(@(x) sprintf('+-%.1f%%', 100*x), prec_target, 'uni', 0), n_place);

fprintf('\n=== (2) n needed to DISCRIMINATE dose-response curvature (per dose, 80%% power) ===\n');
show('CV', {sprintf('vs %s', law.name{1}), sprintf('vs %s', law.name{3})}, n_law);
fprintf('    separation from log-linear: %.4f and %.4f log units (%.1f %% and %.1f %%)\n', ...
        delta_law(1), delta_law(3), 100*delta_law(1), 100*delta_law(3));

fprintf('\n=== (3) n needed to CATCH a model prediction error at 4 mM (80%% power) ===\n');
show('CV', arrayfun(@(x) sprintf('%.0f%% off', 100*x), model_err, 'uni', 0), n_model);

fprintf('\n=== (4) n PER SEX GROUP to detect a sex difference in the ATP ratio ===\n');
show('CV', arrayfun(@(x) sprintf('%.0f%% diff', 100*x), sex_eff, 'uni', 0), n_sex);

% Precision actually delivered by the recommended n, at each CV.
n_grid = 1:12;
sem_grid = CV_set(:) * (1./sqrt(n_grid));

%% -------------------------------------------------------------- the figure
fig = figure('Position',[80 60 1500 900],'Color','w');
tl  = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(tl, 'Protocol design: order, repeats, sample size, and the cost of a 4 mM arm', ...
      'FontWeight','bold','FontSize',13);

% --- 1. damage suffered, per sequence
nexttile; hold on; grid on; box on;
cols = lines(3);  cmap = containers.Map({8,4,2}, {cols(1,:), cols(2,:), cols(3,:)});
for k = 1:nseq
    for j = 1:numel(S(k).order)
        bar(k + 0.26*(j-(numel(S(k).order)+1)/2), 100*S(k).f(j), 0.22, ...
            'FaceColor', cmap(S(k).order(j)), 'EdgeColor','none');
    end
end
set(gca,'XTick',1:nseq,'XTickLabel',seq.name,'XTickLabelRotation',35);
ylabel('force already lost when the run starts (%)');
title('Rundown suffered by each run'); xlim([0.4 nseq+0.6]);
h = [patch(nan,nan,cols(1,:)), patch(nan,nan,cols(2,:)), patch(nan,nan,cols(3,:))];
legend(h, {'8 mM','4 mM','2 mM'}, 'Location','northwest','Box','off');

% --- 2. uncorrected bias per sequence
nexttile; hold on; grid on; box on;
for k = 1:nseq
    for bi = 1:numel(S(k).bias)
        b = S(k).bias(bi);
        if isnan(b.val); continue; end
        bar(k + 0.2*(bi-(numel(S(k).bias)+1)/2), 100*(b.val-1), 0.18, ...
            'FaceColor', cmap(b.cond), 'EdgeColor','none');
    end
end
yline(0,'k-');
set(gca,'XTick',1:nseq,'XTickLabel',seq.name,'XTickLabelRotation',35);
ylabel('bias on the ratio vs 8 mM (%)');
title('Uncorrected order bias'); xlim([0.4 nseq+0.6]);

% --- 3. the three dose laws
nexttile; hold on; grid on; box on;
aa = linspace(1.6, 8.6, 200);
plot(aa, 1 + (R2-1)*(1./aa-1/8)/(1/2-1/8), '-', 'LineWidth',1.8, 'Color',cols(1,:));
plot(aa, exp(log(R2)*(log(8)-log(aa))/(log(8)-log(2))), '-','LineWidth',1.8,'Color',cols(2,:));
plot(aa, 1 + (R2-1)*(8-aa)/(8-2), '-', 'LineWidth',1.8, 'Color',cols(3,:));
plot([8 2],[1 R2],'ko','MarkerFaceColor','k','MarkerSize',8);
semA = CV_plan*law.R4(2)/sqrt(3);
errorbar(4, law.R4(2), 2*semA, 'ks','MarkerFaceColor',[.85 .85 .85],'MarkerSize',9,'LineWidth',1.5);
text(4.15, law.R4(2)-0.055, sprintf('n=3, \\pm2 SEM\n(CV %.0f%%)',100*CV_plan), 'FontSize',9);
set(gca,'XDir','reverse','XScale','log','XTick',[2 4 8],'XTickLabel',{'2','4','8'});
xlabel('[ATP] (mM)'); ylabel('force ratio vs 8 mM');
title('What 4 mM would test');
legend([law.name, {'measured anchors'}], 'Location','northwest','Box','off','FontSize',8);

% --- 4. n vs CV, the three purposes
nexttile; hold on; grid on; box on;
plot(CV_set*100, n_place(:,2), 'o-','LineWidth',1.8,'Color',cols(1,:));
plot(CV_set*100, n_model(:,2), 's-','LineWidth',1.8,'Color',cols(2,:));
plot(CV_set*100, n_law(:,1),   'd-','LineWidth',1.8,'Color',cols(3,:));
set(gca,'YScale','log'); yline(3,'k--','recommended n = 3');
xlabel('between-prep CV of the ratio (%)'); ylabel('preparations required');
title('Required n by purpose');
legend({'place 4 mM to \pm5 %','catch a 15 % model error','resolve curvature'}, ...
       'Location','northwest','Box','off','FontSize',8);

% --- 5. SEM vs n
nexttile; hold on; grid on; box on;
for a = 1:numel(CV_set)
    plot(n_grid, 100*sem_grid(a,:), 'o-','LineWidth',1.6,'Color',cols(a,:));
end
xlabel('preparations per condition'); ylabel('SEM of the ratio (%)');
title('Precision delivered'); legend(CV_lbl,'Location','northeast','Box','off','FontSize',8);
xline(3,'k--');

% --- 6. phi sensitivity: how much of the answer is assumed
nexttile; hold on; grid on; box on;
seqs_show = [1 4 7 8];
for si = 1:numel(seqs_show)
    k = seqs_show(si);  fv = zeros(size(phi_set));
    for pi_ = 1:numel(phi_set)
        o = S(k).order; Dc = 0; last = 0;
        for j = 1:numel(o)
            if o(j) ~= 8; last = Dc / F_cond(ATP==o(j)); end
            Dc = Dc + D_of(o(j), phi_set(pi_));
        end
        fv(pi_) = 100*last;
    end
    plot(phi_set, fv, 'o-','LineWidth',1.8,'DisplayName',S(k).name);
end
xlabel('\phi (permanent fraction of within-run decline)');
ylabel('correction on the last non-8 run (%)');
title('\phi is a \times2 unknown — and it multiplies the correction');
legend('Location','northwest','Box','off');

exportgraphics(fig, fullfile(resdir,'design_power.png'), 'Resolution',150);
save(fullfile(resdir,'design_power.mat'), 'S','n_place','n_law','n_model','n_sex', ...
     'law','CV_set','phi_set','delta_law','R_dose','slope','Tact');
fprintf('\nWrote %s\n', fullfile(resdir,'design_power.png'));

%% ------------------------------------------------------------------ helper
function printTable(rowhdr, colhdr, rowlbl, M)
    fprintf('%-24s', rowhdr);
    for c = 1:numel(colhdr); fprintf('%16s', colhdr{c}); end
    fprintf('\n');
    for r = 1:size(M,1)
        fprintf('%-24s', rowlbl{r});
        for c = 1:size(M,2); fprintf('%16d', M(r,c)); end
        fprintf('\n');
    end
end
