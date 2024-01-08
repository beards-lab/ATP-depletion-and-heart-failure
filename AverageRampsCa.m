% processes the Ca data, expects peaks from normal
% clearvars -except dataset;

% 
% dataset{1} = load('DataStruct20230518_renamed.mat');
dataset{1} = load('DataStruct20230919.mat');
dataset{2} = load('DataStruct20230927.mat');
dataset{3} = load('DataStruct20230928.mat');
dataset{4} = load('DataStruct20231027.mat');
dataset{5} = load('DataStruct20231102.mat');
dataset{6} = load('DataStruct20231107.mat');
%% Get fmaxes - for each dataset there is a Fmax
% index: dataset, 
% measurement, sequence;
figure(1);clf;hold on;
fmax_pointers = [8, 1;2 1;2, 1;2, 1;2, 1;2,1];
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_maxF(i) = NaN;
        continue;
    end
    dtst = dataset{i}.dsc{fmax_pointers(i, 1), fmax_pointers(i, 2)};
    rmp = dtst.datatable;
    % subplot(122);hold on;
    % plot(rmp.t, rmp.L);
    % subplot(121);hold on;
    pfm(i) = plot(rmp.t, rmp.F);

    %  the pCa starts about halfway the dataset
    i_start = find(rmp.t > 100, 1, 'first');
    i_firstSteps = find(rmp.L(i_start:end) > 1.001, 1, 'first') - 5 + i_start;
    % rmp.t([i_start, i_firstSteps])

    % plot(rmp.t(i_start:i_firstSteps), rmp.F(i_start:i_firstSteps), 'x')
    [maxF i_max] = max(rmp.F(i_start:i_firstSteps));
    dataset_maxF(i) = maxF;
    plot(rmp.t(i_max+i_start), maxF, 'rx', LineWidth=4, MarkerSize=12)
end
legend(pfm, '20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F');
%% Get the ss F of the slowest active ramp
ss_pointers = [2, 2;3,7;3,7;3,7;3,7;3,7];
dataset_ssF = [];
clf;hold on;
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_ssF(i) = NaN;
        continue;
    end
    dttbl = dataset{i}.dsc{ss_pointers(i, 1), ss_pointers(i, 2)}.datatable;
    i_ss = find(dttbl.L > max(dttbl.L)*0.999, 40, 'last') - 10;
    dataset_ssF(i) = mean(dttbl.F(i_ss));
    pl(i) = plot(dttbl.t, dttbl.F); plot(dttbl.t(i_ss), dttbl.F(i_ss), '|');
end
legend(pl, '20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F');
%% cell for each ramp
% ramp durations
rds = [100, 10, 1, 0.1];
% rds = [1, 0.1];
% dataset, logtrace, ramp for each ramp duration
% only final dataset
pCa4_4Data{1} = [1, 2, 2;2, 3, 7;3, 3, 7;4,3, 7;5, 3, 7;6, 3, 7];
pCa4_4Data{2} = [1, 2, 3;2, 3, 8;3, 3, 8;4,3, 8;5, 3, 8;6, 3, 8];
pCa4_4Data{3} = [1, 2, 4;2, 3, 9;3, 3, 9;4,3, 9;5, 3, 9;6, 3, 9];
pCa4_4Data{4} = [1, 2, 5;2, 3,10;3, 3,10;4,3,10;5, 3,10;6, 3,10];
dsName = 'AvgpCa4.4';

%% Get corresponding relaxed for scaling options
% Only final dataset with the new protocol
relaxed{1} = [1, 1, 9;2, 1, 6;3, 1, 6;4,1,9;5, 1, 9;6, 1, 9];
relaxed{2} = [1, 1, 8;2, 1, 7;3, 1, 7;4,1,8;5, 1, 8;6, 1, 8];
relaxed{3} = [1, 1, 7;2, 1, 8;3, 1, 8;4,1,7;5, 1, 7;6, 1, 7];
relaxed{4} = [1, 1, 6;2, 1, 9;3, 1, 9;4,1,6;5, 1, 6;6, 1, 6];

%% Adjustment of remaining active force
figure(6); clf; 
clf;

% linear
FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*FremMax(1).*(rd-t)./rd) + ... % 
            (t >= rd).*FremMax(2).*(t-rd)/(t(end)*0 + 300 -rd);
% first order
k = .015;
s = 1;

% exponential
FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*((1-s).*FremMax(1) + s.*FremMax(1).*(rd-t)./rd)) +... % 
            (t >= rd).*(...
            s*FremMax(2).*...
            (1-exp(-k*(t-rd))) +...
            (1-s).*FremMax(2)); % exponential approximation
% s.*FremMax(2)/(1-exp(-k*(t(end) - rd))).*... % scaling to hit the FremMax at the end

% indexes of active PNB 300s hold datasets
i_pca4_4_100ms_hold300 = [2 5; 3 10; 3 10;3 11; 3 11; 3 11];
% get the pCa 11 long hold for decay fitting 
% row = dataset#, col 1 measurement set, col 2 ramp index
i_Relax_100ms_hold300 = [1 10;1 9;1 9;1 6;1 6;1 6];
% indexes of validation datasets
% same or one step slower
i_validation = [2 4;3 9;3 9;3 10;3 10;3 10];
% 10s
i_validation = [2 3;3 8;3 8;3 8;3 8;3 8];
tiledlayout(2,3)
for i_dtst = 1:length(i_pca4_4_100ms_hold300)
    nexttile();
    dtst = dataset{i_dtst}.dsc{i_pca4_4_100ms_hold300(i_dtst, 1), i_pca4_4_100ms_hold300(i_dtst, 2)};
    dtst_relaxed = dataset{i_dtst}.dsc{i_Relax_100ms_hold300(i_dtst, 1), i_Relax_100ms_hold300(i_dtst, 2)};
    rmp = dtst.datatableZDCorr;
    rmpRel = dtst_relaxed.datatableZDCorr;
   
    
    i_0 = find(rmp.t >= 10, 1);
    % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
    dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
    % i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
    % int_avg
    rmp.F = movmean(rmp.F, [8 8]);
    rmp.t = rmp.t - 10;
    rmpRes = rmp(1:10:end, :);
    rmpRelRes = rmpRel(1:10:end, :);
    rmpRelRes.t = rmpRelRes.t - 10;
    fitRelRng = rmpRelRes.t > dtst_relaxed.rd & rmpRelRes.t < 300;

    % fit the relaxed decay first
    f_fitRelax = @(a, b, c, x) a*(x-dtst_relaxed.rd).^(-b) + c;
    [rae rgood] = fit(rmpRelRes.t(fitRelRng), rmpRelRes.F(fitRelRng), f_fitRelax, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0]);    

    % plot(rmpRel.t-10, rmpRel.F, rmpRelRes.t(fitRelRng), rmpRelRes.F(fitRelRng),...
    %     rmpRelRes.t(fitRelRng), f_fitRelax(rae.a, rae.b, rae.c, rmpRelRes.t(fitRelRng)))
        
    %% remaining force at the beginning and at the end
    FremMax = [mean(rmp.F(1:i_0-1/dt)), mean(rmp.F(length(rmp.F)-1/dt:end))];

    % identify the exponential coefficient by fitting
    t_fit_start = 30;    
    if t_fit_start > rmp.t(end) - 10
        % fit first and last only
        fitrng = rmpRes.t < 0 | rmpRes.t > rmpRes.t(end) - 1;
        % plot(rmpRes.t, rmpRes.F, rmpRes.t, fitrng)
    else
        % we have hold, we can identify the remaining force rise
        fitrng = rmpRes.t > t_fit_start & rmpRes.t < 300;
    end
    
    f_k = @(a, b, c, x) a*(1-exp(-b*x)) + c;
    f_fit = @(a, b, c, x) f_k(a, b, c, x) + f_fitRelax(rae.a, rae.b, rae.c, x);
    [ae good] = fit(rmpRes.t(fitrng), rmpRes.F(fitrng), f_fit, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0], 'Upper', [100, 0.1, 100]);    
        
    % Match that by the fit function
    FremFitShift = FremMax(2) - f_k(ae.a, ae.b, 0, rmpRes.t(end));
    % minRemForce = f_fit(ae.a, ae.b, FremFitShift, dtst.rd);
    FremFun = @(t, rd) (t <= 0)*FremMax(1) + ...
            (t > 0 & t < rd).*max(-10, min(FremMax(1), (f_k(ae.a, ae.b, FremFitShift, rd) - FremMax(1))/rd*t + FremMax(1))) +... % 
            (t >= rd).*(...
            f_k(ae.a, ae.b, FremFitShift, t)); % exponential approximation

    % plot the raw force
    pltF = plot(rmp.t, rmp.F, ':');
    hold on;
    % show the remaining force approximation
    pltApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd), '--', LineWidth=3);
    % show the corrected dataset
    pltCorr = plot(rmpRes.t, rmpRes.F - FremFun(rmpRes.t, dtst.rd), '-');    
    % show the fit
    x_ax = 0.1:0.1:300;
    plot(x_ax, f_fit(ae.a, ae.b, ae.c, x_ax), '--', LineWidth=2);
    plot(x_ax, f_k(ae.a, ae.b, ae.c, x_ax), '--', LineWidth=2);
    plot(x_ax, f_fitRelax(rae.a, rae.b, rae.c, x_ax), '--', LineWidth=2);
    
    % approximation zones
    pltB = plot(rmp.t(1:i_0-1/dt), rmp.F(1:i_0-1/dt), rmp.t(length(rmp.F)-1/dt:end), rmp.F(length(rmp.F)-1/dt:end), LineWidth=2, Color=[0.8500    0.3250    0.0980]);
    xlim([rmp.t(1), rmp.t(end)]);
    legend([pltF, pltB(1), pltApprx, pltCorr], {'Raw tension data', 'Zones used for averaging', 'Linear approximation of F_{rem}', 'Corrected data'}, 'AutoUpdate','off')
    title(['Correction for dataset ' num2str(i_dtst)])        
    FremFunArr{i_dtst} = FremFun;
    %% repeat for the validation ramp, using the same parametrization of the FremFun
    dtst = dataset{i_dtst}.dsc{i_validation(i_dtst, 1), i_validation(i_dtst, 2)};
    rmp = dtst.datatableZDCorr;
    set(gca,'ColorOrderIndex',1)
    F = movmean(rmp.F, [16 16]);
    rmp.t = rmp.t - 10;
       
    % plot the raw force
    pltF = plot(rmp.t, F, ':');
    % show the remaining force approximation
    pltApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd), '--');
    % show the corrected dataset
    pltCorr = plot(rmp.t, F - FremFun(rmp.t, dtst.rd), '-');    
end

%% Average all data
% figure(5);clf;
% shift of peaks to have the same tail - just guessed
rampShift = [-94, -7.2, -0.35, 0];

figure(4);clf;hold on;

clear peaks peaks_norm leg;

for i_rds = 1:length(rds)
    % sp =subplot(4, 4, (i_rds-1)*4 +  (1:2));cla;
    % excluding 100s
    if i_rds ~= 1 % skip plotting of 100s
        sp =subplot(3, 4, (i_rds-2)*4 +  (1:2));cla;
    end
    outF = [];sum_squared_diff = []; n = 1;clear leg;
    rampSet = pCa4_4Data{i_rds};
    clin = lines(size(rampSet, 1)+1);

    for i_logtrace = 1:size(rampSet, 1)
        % if rampSet(i_logtrace, 1) == 0
        %     % placeholder to match the relaxed peaks, skipping            
        %     continue;
        % end
        % experiment dataset - whole measurement session
        eds = dataset{rampSet(i_logtrace, 1)}.dsc;
        % dataset - particular conditions
        dtst = eds{rampSet(i_logtrace, 2), rampSet(i_logtrace, 3)};
        rmp = dtst.datatableZDCorr;
        i_0 = find(rmp.t >= 10, 1);
        % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
        dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
        i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
        
        % trim the ramp - start to end
        rmp = rmp(1:i_end, :);
        % rmp = rmp(i_0:i_end, :);
        rmp.t = rmp.t - 10;
        % i_0 = 1;

        
        
        %% filter out the remaining force
        % saving the original 
        rmp.Fraw = rmp.F;
        rmp.F = rmp.F - FremFunArr{i_logtrace}(rmp.t, dtst.rd);            
        
        % base on nothing - just absolute peak values
        base_rel = 1;
        % base on absolute peak of each ramp
        % base_rel = max(rmp.F);
        % base on absolute peak of fastest ramp
        % base_rel = ...
        %     dataset{pCa4_4Data{4}(i_logtrace, 1)}.dsc{pCa4_4Data{4}(i_logtrace, 2), pCa4_4Data{4}(i_logtrace, 3)}.peak;
        % base on steady state
        % base_rel = rmp.F(end);
        % base on relaxed peaks for each ramp
        % base_rel = eds{relaxed{i_rds}(i_logtrace, 2),...
                        % relaxed{i_rds}(i_logtrace, 3)}.peak;
        % base on relaxed peaks for fastest ramp only
        % base_rel = eds{relaxed{4}(i_logtrace, 2),...
                        % relaxed{4}(i_logtrace, 3)}.peak;
        % base on Fmax
        base_rel = dataset_maxF(rampSet(i_logtrace, 1));
        % base on staeady state of the slowest ramp (in AverageRampsCa)
        % base_rel = dataset_ssF(i_logtrace);

        if isnan(base_rel) || base_rel == 0
            continue;
        end

        % adjusted for Fremaining, normalized
        F = rmp.F/base_rel;
        % just raw force, normalized
        % F = rmp.Fraw/base_rel;
        % just raw force
        % F = rmp.F;

        % prepare legend
        leg{i_logtrace} = [dtst.folder ':' dtst.datasetTitle];        
        % semilogx(rmp.t, movmean(F, [16 16]), ':', Color=clin(i_logtrace, :));hold on;
        semilogx(rmp.t, F, '-', Color=[clin(i_logtrace, :), 0.05]);hold on;
        % semilogx(rmp.t, FremFunArr{i_logtrace}(rmp.t, dtst.rd)/base_rel, '--', Color=clin(i_logtrace, :));
        % semilogx(dtst.datatableZDCorr.t - 10, dtst.datatableZDCorr.F, '-', Color=clin(i_logtrace, :));hold on;
        % semilogx(rmp.t, FremFun(rmp.t, dtst.rd), '--', Color=clin(i_logtrace, :));
        set(gca, 'FontSize', 14);
        % base on absolute peak
        peaks(i_rds, i_logtrace) = max(F*base_rel);
        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        
        % normalization        
        peaks_norm(i_rds, i_logtrace) = max(F);

        if isempty(outF)
            outF = F;
            outFr = rmp.Fraw;
            Fmax = base_rel;
            n = 1;
            sum_squared_diff = outF*0;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF  = outF(1:L)  + (F(1:L)     - outF(1:L))/n;
            outFr = outFr(1:L) + (rmp.Fraw(1:L) - outFr(1:L))/n;
            outT = rmp.t(1:L);
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;

            % estimate the standard error out of sum squared
            sum_squared_diff = sum_squared_diff(1:L) + (F(1:L) - outF(1:L)).^2;            
        end    

        %% getting stiffnes based on linear fit
        % fitrg = rmp.L > 1.1 & rmp.t < dtst.rd;
        % rmpFitrg = rmp(fitrg, :);
        % % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        % fitfun = @(a, b, x) a.*(x + b);
        % 
        % [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1]);
        % as(i_rds, i_logtrace) = ae.a;
        figure(5);
        semilogx(rmp.t + rampShift(i_rds), rmp.Fraw, '-', Color=[clin(i_logtrace, :), 0.1]);hold on;
        figure(4);
        
        
    end
%% resample and save
    t_ignore = 0.01; 
    % resample up the peak and for the tail separately
    % t_s = [linspace(0, 1, 20)*(rmp.t(i_0) + rds(i_rds) - t_ignore) ...
    %     logspace(log10(rmp.t(i_0) + rds(i_rds) + t_ignore), log10(rmp.t(i_0) + rds(i_rds) + 10), 20) ...
    %     logspace(log10(rmp.t(i_0) + rds(i_rds) + 10), log10(outT(end)), 5)];

    % resample up the peak and for the tail separately
    t_s = [linspace(0, 1, 20)*(rmp.t(i_0) + rds(i_rds) - t_ignore) ...
        logspace(log10(rmp.t(i_0) + rds(i_rds) + t_ignore), log10(outT(end)), 40)];

    % resample log equally
    % t_s = [logspace(log10(1e-3), log10(outT(end)), 40)];
    % remove consequent duplicates at joints
    t_s = t_s(~[false t_s(2:end) == t_s(1:end-1)]);
    SD = sqrt(sum_squared_diff / (n * (n - 1)));
    % force and length interpolation
    FLSDint = interp1(outT, [outF, rmp.L(1:L), SD], t_s, "pchip", 'extrap');
    
    tab_rmpAvg = table(t_s' + 2, FLSDint(:, 2), FLSDint(:, 1)*Fmax);
    tab_rmpAvg.Properties.VariableNames = {'Time', 'L', 'F'};
    % writetable(tab_rmpAvg, ['data/' dsName '_' num2str(rds(i_rds)) 's.csv']);
%% plot the AVG

    % Fill the area between upper and lower bounds to show the confidence interval
    % fill([outT', fliplr(outT')], [(outF*Fmax - SD*1.96)', fliplr((outF*Fmax + SD*1.96)')], 'b');

    % semilogx(t_s, FLSDint(:, 1) + FLSDint(:, 3), '-', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1) - FLSDint(:, 3), '-', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1), '-|', LineWidth=2, Color=clin(end, :))
    
    errorbar(t_s,FLSDint(:, 1),FLSDint(:, 3), '-', LineWidth=2, Color=clin(end, :))
    xlim([1e-2 2e2])
    yl = ylim;
    ylabel('\itT_{m,rel}');
    yyaxis right; ylim(yl*Fmax);
    ylabel('\itT_m (kPa)');
    yyaxis left;
    % plot(rmp.t(1:L), outF(1:L) + (rmp.F(1:L) - outF(1:L))/n)
    leg{length(leg) +1} = 'Averaged';
    title(['Ramp ' num2str(rds(i_rds)) 's'])
    if i_rds == 1
        % x = 0.1;y = 0.03;
        % legend(leg, 'Interpreter','none', 'Position', ...
        %     [sp.Position(1) + sp.Position(3) + x, sp.Position(2) - y, 1 - sp.Position(1) - sp.Position(3) - x, sp.Position(4)])
    elseif i_rds == 4
        xlabel('Time (s)')
    end    
    
    % figure(5);
    % rampShift = [-92, -8.6, -0.58, 0];
    % semilogx(outT + rampShift(i_rds), outFr, 'k-', LineWidth=1);hold on;
    % xlim([1e-2 1e2])
    % figure(4);
    % Tarr{i_rds} = t_s;
    % Farr{i_rds} = interp1(outT, outFr, t_s, "pchip", 'extrap');
    Tarr{i_rds} = outT;
    Farr{i_rds} = outFr;
    

end
% peak sum up
% normtype = 'Norm to highest peak'
% normtype = 'Norm to steady state'
% normtype = 'Norm to relaxed peaks for each ramp'
% normtype = 'Norm to highest relaxed peak'
normtype = 'Norm to Fmax'
subplot(4, 4, [3 16]);cla;

% figure(11);
% subplot(1, 5, 5);cla reset;
% remove zero peaks as NaNs to fix the average - only when something was missing
peaks_norm(peaks_norm == 0) = NaN;

% for x axis we go from fastest to slowest
x_ax = length(rds):-1:1;


boxplot(fliplr(peaks_norm'), PlotStyle="traditional", Notch="off");hold on;

% coefficient of variance
cv = nanstd(peaks_norm')./nanmean(peaks_norm');
fprintf('Coefficient of variance: %1.3f \n', sum(cv));

clin = lines(size(peaks_norm, 2)+1);
for i_pk = 1:size(peaks_norm, 2)
    if all(isnan(peaks_norm(:, i_pk)))
        % it was just a placeholder
        leg{i_pk} = '';
        continue;
    end
    % set(gca, 'colororderindex', 1);
    
    plot(x_ax, peaks_norm(:, i_pk), '.:', 'MarkerSize',12, LineWidth=1.5, Color=clin(i_pk, :));hold on;    
    % semilogx(rds, peaks_relaxed(:, i_pk), 'x-', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
    % semilogx(rds, as(:, i_pk), 'x:', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
end

% just for the legend
plot(NaN, NaN, '-|',LineWidth=3, Color=clin(end, :))
%# find non-empty legend
% validLeg = ~cellfun(@isempty,leg);
% %# remove empty cells
% leg_cleared = leg(validLeg);
leg_cleared = {'20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};

legend(leg_cleared, 'Interpreter','none', 'AutoUpdate','off', 'Location','northeast')

plot(x_ax, nanmean(peaks_norm, 2)', '-', LineWidth=3, Color=clin(end, :))
plot(x_ax, nanmean(peaks_norm, 2)', '|', LineWidth=5, Color=clin(end, :), MarkerSize=12)

% plot(x_ax, nanmedian(peaks_norm, 2)', '--', LineWidth=3, Color=clin(end, :))
% plot(x_ax, nanmedian(peaks_norm, 2)', '_', LineWidth=5, Color=clin(end, :), MarkerSize=12)

xlabel('Ramp duration (s)');
% xlim([0.08, 150])
% ylabel('Peak tension (kPa)')
% semilogx(rds, mean(peaks, 2)', '_', LineWidth=3, MarkerSize=12)
% set(gca, 'XTick', fliplr(rds));
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'});

set(gca, 'FontSize', 14);
title(sprintf('Averaging peaks (%1.3f)\n %s', sum(cv), normtype))

yl = ylim();
ylim([0, yl(2)]);

% boxplot(peaks', 'Positions',rds)
% end
%% 
% F_corrw = F_corr(win);
% Fdet = detrend(F_corrw);
% Fremdiff = [0 ;diff(movmean(F_corr(win), 4e5))];
% Fremdiff(Fremdiff<0) = 0;
% FremIncrS = cumsum(Fremdiff);
% subplot(133);hold on;
% plot(rmp.t(win), [0 ;diff(movmean(F_corr(win), 4e5))], '+-', LineWidth=2);
% plot(rmp.t(win), FremIncrS)

%% Compare the shift

% rampShift = [-94, -7.2, -0.35, 0];
figure(5);
% clf;
for i_rds = 4:-1:1
    semilogx(Tarr{i_rds} + rampShift(i_rds), Farr{i_rds}, 'k', LineWidth=5-i_rds);hold on;
end
xlim([1e-1, 50])

return;
%% compare pca and relaxed

figure(5);clf;
peaks(peaks==0) = NaN;
boxplot(fliplr(peaks'), PlotStyle="traditional", Notch="off");hold on;
plot(4:-1:1, peaks, '.:', 'LineWidth',1);
plot(4:-1:1, nanmean(peaks, 2)', '-', LineWidth=3, Color=clin(end, :))
plot(4:-1:1, nanmean(peaks, 2)', '|', LineWidth=5, Color=clin(end, :), MarkerSize=12)
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'});
xlabel('Ramp duration (s)')
ylabel('Tension peak (kPa)');
title('Averaging the peaks for pCa 4.4')

%%
% semilogx(rds, peaks_relaxed, 'x-');
semilogx(rds, peaks./peaks_relaxed(4, :), 'x-');
% invalidLeg = cellfun(@isempty,leg);
%# remove empty cells
leg(3) = {'skip'};

legend(leg, 'Interpreter','none', 'AutoUpdate','off', 'Location','northeast')
%%
f = @(t)8/90*(t-10);
plot(rmp.t, rmp.F, rmp.t, f(rmp.t), rmp.t, rmp.F - f(rmp.t))