% processes the Ca data, expects peaks from normal
clearvars -except dataset;

% 
% dataset{1} = load('DataStruct20230518_renamed.mat');
% dataset{2} = load('DataStruct20230919.mat');
% dataset{3} = load('DataStruct20230927.mat');
% dataset{4} = load('DataStruct20230928.mat');
% dataset{5} = load('DataStruct20231027.mat');
% dataset{6} = load('DataStruct20231102.mat');
% dataset{7} = load('DataStruct20231107.mat');
%% Get fmaxes - for each dataset there is a Fmax
% index: dataset, 
% measurement, sequence;
fmax_pointers = [0, 0;8, 1;2 1;2, 1;2, 1;2, 1;2,1];
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_maxF(i) = NaN;
        continue;
    end
    dataset_maxF(i) = max(...
        dataset{i}.dsc{fmax_pointers(i, 1), fmax_pointers(i, 2)}.datatable.F);
end

%% cell for each ramp
% ramp durations
rds = [100, 10, 1, 0.1];
% rds = [1, 0.1];
% dataset, logtrace, ramp for each ramp duration
% all
pCa4_4Data{1} = [1, 4, 2;2, 2, 2;3, 3, 7;4, 3, 7;5,3, 7;6, 3, 7;7, 3, 7];
pCa4_4Data{2} = [1, 4, 3;2, 2, 3;3, 3, 8;4, 3, 8;5,3, 8;6, 3, 8;7, 3, 8];
pCa4_4Data{3} = [1, 4, 4;2, 2, 4;3, 3, 9;4, 3, 9;5,3, 9;6, 3, 9;7, 3, 9];
pCa4_4Data{4} = [1, 4, 5;2, 2, 5;3, 3,10;4, 3,10;5,3,10;6, 3,10;7, 3,10];
dsName = 'AvgpCa4.4';

%% Get corresponding relaxed for scaling options
relaxed{1} =    [1, 1, 2;2, 1, 2;3, 1, 6;4, 1, 6;5,1, 9;6, 1, 9;7, 1, 9];
relaxed{2} =    [1, 1, 3;2, 1, 3;3, 1, 7;4, 1, 7;5,1, 8;6, 1, 8;7, 1, 8];
relaxed{3} =    [1, 1, 4;2, 1, 4;3, 1, 8;4, 1, 8;5,1, 7;6, 1, 7;7, 1, 7];
relaxed{4} =    [1, 1, 5;2, 1, 5;3, 1, 9;4, 1, 9;5,1, 6;6, 1, 6;7, 1, 6];


% %% 
% dtst = dataset{2}.dsc;
% dtst{1, 1}.datasetTitle
% rmp = dtst{1, 1}.datatable;

%
figure(4);clf;hold on;
clear peaks peaks_norm leg;

for i_rds = 1:length(rds)
    sp =subplot(4, 4, (i_rds-1)*4 +  (1:2));%cla;
    outF = [];n = 1;clear leg;
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
        
        % remaining force at the beginning and at the end
        FremMax = [mean(rmp.F(i_0 -100:i_0+100)), mean(rmp.F(length(rmp.F)-1/dt:end))];

        % trim the ramp - start to end
        rmp = rmp(1:i_end, :);
        % rmp = rmp(i_0:i_end, :);
        % rmp.t = rmp.t - 10;
        % i_0 = 1;


        %% filter out the remaining force
        % linear
        % FremFun = @(t, rd) (t < rd).*FremMax0.*(rd-t)./rd + ... 
        %    (t >= rd).*FremMaxEnd.*(t-rd)/(t(end)-rd);
        
        % first order exponential
        k = .015; s = 1;
        FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*((1-s).*FremMax(1) + s.*FremMax(1).*(rd-t)./rd)) +... % 
            (t >= rd).*(...
            s*FremMax(2).*...
            (1-exp(-k*(t-rd))) +...
            (1-s).*FremMax(2)); % exponential approximation
        
        rmp.Fraw = rmp.F;
        rmp.F = rmp.F - FremFun(rmp.t - rmp.t(i_0), dtst.rd, FremMax);        
        
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
        % base_rel = dataset_maxF(rampSet(i_logtrace, 1));
        if isnan(base_rel) || base_rel == 0
            continue;
        end

        F = rmp.F/base_rel;
        % F = rmp.F;

        % prepare legend
        leg{i_logtrace} = [dtst.folder ':' dtst.datasetTitle];
        plot(rmp.t, F, ':', Color=clin(i_logtrace, :));hold on;
        % semilogx(dtst.datatableZDCorr.t - 10, dtst.datatableZDCorr.F, '-', Color=clin(i_logtrace, :));hold on;
        % semilogx(rmp.t, FremFun(rmp.t, dtst.rd), '--', Color=clin(i_logtrace, :));
        set(gca, 'FontSize', 14);
        % base on absolute peak
        peaks(i_rds, i_logtrace) = max(rmp.F);
        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        
        % normalization        
        peaks_norm(i_rds, i_logtrace) = max(F);

        if isempty(outF)
            outF = F;
            outFr = rmp.Fraw;
            Fmax = base_rel;
            n = 1;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF = outF(1:L) + (F(1:L) - outF(1:L))/n;
            outFr = outFr(1:L) + (rmp.F(1:L) - outFr(1:L))/n;
            outT = rmp.t(1:L);
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;
            
        end    

        %% getting stiffnes based on linear fit
        % fitrg = rmp.L > 1.1 & rmp.t < dtst.rd;
        % rmpFitrg = rmp(fitrg, :);
        % % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        % fitfun = @(a, b, x) a.*(x + b);
        % 
        % [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1]);
        % as(i_rds, i_logtrace) = ae.a;

        
    end
%% resample and save
    
    % resample up the peak and for the tail separately
    t_s = [linspace(0, 1, 20)*(rmp.t(i_0) + rds(i_rds)) ...
        logspace(log10(rmp.t(i_0) + rds(i_rds)), log10(outT(end)), 40)];
    % resample log equally
    % t_s = [logspace(log10(1e-3), log10(outT(end)), 40)];
    % remove consequent duplicates at joints
    t_s = t_s(~[false t_s(2:end) == t_s(1:end-1)]);
    % force and length interpolation
    FLint = interp1(outT, [outF*Fmax, rmp.L(1:L)], t_s, "pchip", 'extrap');
    
    tab_rmpAvg = table(t_s' + 2, FLint(:, 2), FLint(:, 1));
    tab_rmpAvg.Properties.VariableNames = {'Time', 'L', 'F'};
    writetable(tab_rmpAvg, ['data/' dsName '_' num2str(rds(i_rds)) 's.csv']);
%% plot the AVG
    semilogx(t_s, FLint(:, 1)/Fmax, '-|', LineWidth=2, Color=clin(end, :))
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
end
%% peak sum up
normtype = 'Norm to highest peak'
normtype = 'Norm to steady state'
normtype = 'Norm to relaxed peaks for each ramp'
normtype = 'Norm to highest relaxed peak'
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
% set(gca, 'YAxisLocation', 'right');
% just for the legend
plot(NaN, NaN, '-|',LineWidth=3, Color=clin(end, :))
%# find non-empty legend
validLeg = ~cellfun(@isempty,leg);
%# remove empty cells
leg_cleared = leg(validLeg);
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

%% representative case of adjustment of active force
figure(6); clf; 
subplot(1, 3, [1 2])
i_dtst = 6;
% linear
FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*FremMax(1).*(rd-t)./rd) + ... % 
            (t >= rd).*FremMax(2).*(t-rd)/(t(end)*0 + 300 -rd);
% first order
k = .015;
s = 1;

FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*((1-s).*FremMax(1) + s.*FremMax(1).*(rd-t)./rd)) +... % 
            (t >= rd).*(...
            s*FremMax(2).*...
            (1-exp(-k*(t-rd))) +...
            (1-s).*FremMax(2)); % exponential approximation
% s.*FremMax(2)/(1-exp(-k*(t(end) - rd))).*... % scaling to hit the FremMax at the end
            


dtst = dataset{i_dtst}.dsc{3,11}
rmp = dtst.datatableZDCorr;

i_0 = find(rmp.t >= 10, 1);
% i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
% int_avg
rmp.t = rmp.t - 10;

% identify the exponential coefficient
fitrng = rmp.t > 100 & rmp.t < 300;
fit_F = rmp.F(fitrng);
fit_t = rmp.t(fitrng);

f_fit = @(a, b, c, x) a*(1-exp(-b*x)) + c;
[ae good] = fit(fit_t, fit_F, f_fit, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0], 'Upper', [100, 0.1, 100]);
figure(21);plot(fit_t, fit_F, ':', 1:300, f_fit(ae.a, ae.b, ae.c, 1:300), '-|');
FremFun
ae.b
%%

% remaining force at the beginning and at the end
FremMax = [mean(rmp.F(1:i_0-1/dt)), mean(rmp.F(length(rmp.F)-1/dt:end))]
pltF = plot(rmp.t, rmp.F, ':', LineWidth=0.5);hold on;
pltB = plot(rmp.t(1:i_0-1/dt), rmp.F(1:i_0-1/dt), rmp.t(length(rmp.F)-1/dt:end), rmp.F(length(rmp.F)-1/dt:end), LineWidth=2, Color=[0.8500    0.3250    0.0980]);
pltCorr = plot(rmp.t, rmp.F - FremFun(rmp.t, dtst.rd, FremMax), '-', LineWidth=2);
pltLinApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd, FremMax), '--', LineWidth=3);
legend([pltF, pltB(1), pltLinApprx, pltCorr], {'Raw tension data', 'Zones used for averaging', 'Linear approximation of F_{rem}', 'Corrected data'})

dtst = dataset{i_dtst}.dsc{3,11}
rmp = dtst.datatableZDCorr;
set(gca,'ColorOrderIndex',1)

i_0 = find(rmp.t >= 10, 1);
% i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
% int_avg
rmp.t = rmp.t - 10;

% the non-increasing window
win = rmp.t > dtst.rd + 20 & rmp.t < rmp.t(end) - 40;
i_FremDecEnd = find(rmp.t > rmp.t(end)-41, 1);

% remaining force at the beginning and at the end

FremMax = [mean(rmp.F(1:i_0-1/dt)), mean(rmp.F(length(rmp.F)-1/dt:end))]
F_corr = rmp.F - FremFun(rmp.t, dtst.rd, FremMax);

pltF = plot(rmp.t, rmp.F, ':', LineWidth=0.5);hold on;
pltB = plot(rmp.t(1:i_0-1/dt), rmp.F(1:i_0-1/dt), rmp.t(length(rmp.F)-1/dt:end), rmp.F(length(rmp.F)-1/dt:end), LineWidth=2, Color=[0.8500    0.3250    0.0980]);
pltCorr = plot(rmp.t, rmp.F - FremFun(rmp.t, dtst.rd, FremMax), '-', LineWidth=2);
% plot(rmp.t, movmean(rmp.F, 1000), '-', LineWidth=3);
% plot(rmp.t, movmean(F_corr, 10000), '--', LineWidth=3);
% plot(rmp.t(win), movmean(F_corr(win), 4e5), '--', LineWidth=3);
plot(rmp.t(win), movmean(detrend(F_corr(win), 1), 1e4), '-', LineWidth=2);

plot([0 rmp.t(end)], repmat(F_corr(i_FremDecEnd), [2 1]), 'k-')
pltLinApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd, FremMax), '--', LineWidth=3);

legend([pltF, pltB(1), pltLinApprx, pltCorr], {'Raw tension data', 'Zones used for averaging', 'Approximation of F_{rem}', 'Corrected data'})

F_corrw = F_corr(win);
Fdet = detrend(F_corrw);
Fremdiff = [0 ;diff(movmean(F_corr(win), 4e5))];
Fremdiff(Fremdiff<0) = 0;
FremIncrS = cumsum(Fremdiff);
subplot(133);hold on;
plot(rmp.t(win), [0 ;diff(movmean(F_corr(win), 4e5))], '+-', LineWidth=2);
plot(rmp.t(win), FremIncrS)
return
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