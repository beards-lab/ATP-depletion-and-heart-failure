% processes the relaxed data
% clear -except dataset;

% using only the final dataset
dataset{1} = load('DataStruct20230919.mat');
dataset{2} = load('DataStruct20230927.mat');
dataset{3} = load('DataStruct20230928.mat');
dataset{4} = load('DataStruct20231027.mat');
dataset{5} = load('DataStruct20231102.mat');
dataset{6} = load('DataStruct20231107.mat');
isMale = [1, 1, 0, 0, 1, 0];
%% cell for each ramp
% ramp durations
rds = [100, 10, 1, 0.1];
dsName = 'AvgRelaxed';

% Only final dataset with the new protocol
relaxed{1} = [1, 1, 9;2, 1, 6;3, 1, 6;4,1,9;5, 1, 9;6, 1, 9];
relaxed{2} = [1, 1, 8;2, 1, 7;3, 1, 7;4,1,8;5, 1, 8;6, 1, 8];
relaxed{3} = [1, 1, 7;2, 1, 8;3, 1, 8;4,1,7;5, 1, 7;6, 1, 7];
relaxed{4} = [1, 1, 6;2, 1, 9;3, 1, 9;4,1,6;5, 1, 6;6, 1, 6];

% peaks = relaxed;
% % Only with decay 60s+
% relaxed{1} = [3, 1, 6;4, 1, 6;5,1,9];
% relaxed{2} = [3, 1, 7;4, 1, 7;5,1,8];
% relaxed{3} = [3, 1, 8;4, 1, 8;5,1,7];
% relaxed{4} = [3, 1, 9;4, 1, 9;5,1,6];

% % with PNB
% relaxed{1} = [1, 1, 2;2, 1, 2;2, 1, 9;3, 1, 6;4, 1, 6;5,1,9];
% relaxed{2} = [1, 1, 3;2, 1, 3;2, 1, 8;3, 1, 7;4, 1, 7;5,1,8];
% relaxed{3} = [1, 1, 4;2, 1, 4;2, 1, 7;3, 1, 8;4, 1, 8;5,1,7];
% relaxed{4} = [1, 1, 5;2, 1, 5;2, 1, 6;3, 1, 9;4, 1, 9;5,1,6];


% %%
% dtst = dataset{2}.dsc;
% dtst{1, 1}.datasetTitle
% rmp = dtst{1, 1}.datatable;

%
figure(3);clf;hold on;
clear peaks peaks_norm;

for i_rds = 1:length(rds)
    % sp =subplot(4, 4, (i_rds-1)*4 +  (1:2));cla;
    % excluding 100s
    if i_rds ~= 1 % skip plotting of 100s
        sp =subplot(3, 4, (i_rds-2)*4 +  (1:2));cla;
    end
    outF = [];sum_squared_diff = []; n = 1;clear leg;
    rampSet = relaxed{i_rds};
    clin = lines(size(rampSet, 1)+1);

    for i_logtrace = 1:size(rampSet, 1)
        % experiment dataset - whole measurement session
        eds = dataset{rampSet(i_logtrace, 1)}.dsc;
        % dataset - particular conditions
        dtst = eds{rampSet(i_logtrace, 2), rampSet(i_logtrace, 3)};
        leg{i_logtrace} = [dtst.folder ':' dtst.datasetTitle];
        rmp = dtst.datatableZDCorr;
        i_0 = find(rmp.t >= 10, 1);
        % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
        dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
        i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
        % cut out
        rmp = rmp(i_0:i_end, :);
        rmp.t = rmp.t - 10;
        
        % base rebase not base
        base_rel = 1;
        % base on absolute peak
        % base_rel = max(rmp.F);
        % base on steady state
        % base_rel = rmp.F(end);
        % base on max relaxed peak
        % base_rel = peaks_relaxed(4, i_logtrace);
        % base on Fmax, if available (in AverageRampsCa)
        base_rel = dataset_maxF(rampSet(i_logtrace, 1));
        % base on staeady state of the slowest ramp (in AverageRampsCa)
        % base_rel = dataset_ssF(i_logtrace);

        if isnan(base_rel) || base_rel == 0
            continue;
        end

        F = rmp.F/base_rel;
        % F = rmp.F;
        
        % if isMale(i_logtrace)
        if i_rds ~= 1 % skip 100s
            semilogx(rmp.t, F, '-', LineWidth=0.5, Color=[clin(i_logtrace, :), 0.1]);hold on;
        end
        % else
        %     semilogx(rmp.t, F, ':', LineWidth=1.5);hold on;
        % end
        set(gca, 'FontSize', 14);
        
        % base on obsolute peak
        peaks(i_rds, i_logtrace) = max(rmp.F);
        % normalization        
        peaks_norm(i_rds, i_logtrace) = max(F);

        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        if isempty(outF)
            outF = F;
            Fmax = base_rel;
            n = 1;
            sum_squared_diff = outF*0;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF = outF(1:L) + (F(1:L) - outF(1:L))/n;
            outT = rmp.t(1:L);
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;
             
			% estimate the standard error out of sum squared
            sum_squared_diff = sum_squared_diff(1:L) + (F(1:L) - outF(1:L)).^2;
        end    

        %%
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
    t_s = [linspace(0, 1, 20)*rds(i_rds) ...
        logspace(log10(rds(i_rds)), log10(outT(end)), 40)];
    % resample log equally
    % t_s = [logspace(log10(1e-3), log10(outT(end)), 40)];
    % remove consequent duplicates at joints
    t_s = t_s(~[false t_s(2:end) == t_s(1:end-1)]);
    SD = sqrt(sum_squared_diff / (n * (n - 1)));
    % Force, Length and SD interpolation (relative to Fmax)
    FLSDint = interp1(outT, [outF, rmp.L(1:L), SD], t_s, "pchip", 'extrap');
    
    tab_rmpAvg = table(t_s' + 2, FLSDint(:, 2), FLSDint(:, 1)*Fmax);
    tab_rmpAvg.Properties.VariableNames = {'Time', 'L', 'F'};
    writetable(tab_rmpAvg, ['data/' dsName '_' num2str(rds(i_rds)) 's.csv']);
%% plot the AVG
    % semilogx(outT, outF, 'k-');hold on;
    % semilogx(t_s, FLSDint(:, 1) + FLSDint(:, 3), '--', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1) - FLSDint(:, 3), '--', Color=[clin(end, :), 1], LineWidth=2);
    % plot(t_s, FLSDint(:, 1), '-|', Color=[clin(end, :), 0.5], LineWidth=2)
    
    % Fill the area between upper and lower bounds to show the confidence interval
    % fill([outT', fliplr(outT')], [(outF - SE*1.96)', fliplr((outF + SE*1.96)')], 'b');

    errorbar(t_s,FLSDint(:, 1),FLSDint(:, 3), '-', LineWidth=2, Color=clin(end, :))
    
    xlim([1e-2 2e2])
    yl = ylim;
    ylabel('\itT_{m,rel}');
    yyaxis right; ylim(yl*Fmax);
    ylabel('\itT_m (kPa)');
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
% peak sum up
subplot(4, 4, [3 16]);cla;
% remove zero peaks as NaNs to fix the average - only when something was missing
peaks(peaks == 0) = NaN;


% for x axis we go from fastest to slowest
x_ax = length(rds):-1:1;

boxplot(fliplr(peaks_norm'), PlotStyle="traditional", Notch="off");hold on;

clin = lines(size(peaks, 2)+1);
for i_pk = 1:size(peaks, 2)
    if all(isnan(peaks(:, i_pk))) || all (peaks(:, i_pk) == 0)
        % it was just a placeholder
        leg{i_pk} = '';
        continue;
    end    
    % set(gca, 'colororderindex', 1);
    plot(x_ax, peaks_norm(:, i_pk)', '.:', 'MarkerSize',12, LineWidth=1.5, Color=clin(i_pk, :));hold on;    
    % semilogx(rds, as(:, i_pk), 'x:', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
end
% set(gca, 'YAxisLocation', 'right');
% just for the legend
plot(NaN, NaN, '-|',LineWidth=3, Color=clin(end, :))
% validLeg = ~cellfun(@isempty,leg);
% %# remove empty cells
% leg_cleared = leg(validLeg);

leg_cleared = {'20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};
legend(leg_cleared, 'Interpreter','none', 'AutoUpdate','off', 'Location','northeast')

plot(x_ax, nanmean(peaks_norm, 2)', '-', LineWidth=3, Color=clin(end, :))
plot(x_ax, nanmean(peaks_norm, 2)', '|', LineWidth=5, Color=clin(end, :), MarkerSize=12)


xlabel('Ramp duration (s)');
% xlim([0.08, 150])
ylabel('Peak tension (kPa)')
% semilogx(rds, mean(peaks, 2)', '_', LineWidth=3, MarkerSize=12)
% set(gca, 'XTick', fliplr(rds));
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'});

set(gca, 'FontSize', 14);
title('Absolute peak height')

% boxplot(peaks', 'Positions',rds)
% end
%% the same, but for pCa
% ds = readtable(['Data\AvgpCa4.4_0.1s.csv']);
  