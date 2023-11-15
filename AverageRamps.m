% processes the data
clear;

dataset{1} = load('DataStruct20230518_renamed.mat');
dataset{2} = load('DataStruct20230919.mat');
dataset{3} = load('DataStruct20230927.mat');
dataset{4} = load('DataStruct20230928.mat');
dataset{5} = load('DataStruct20231027.mat');
dataset{6} = load('DataStruct20231102.mat');
%% cell for each ramp
% ramp durations
rds = [100, 10, 1, 0.1];
% dataset, logtrace, ramp for each ramp duration
% all
relaxed{1} = [1, 1, 2;2, 1, 2;2, 1, 9;3, 1, 6;4, 1, 6;5,1,9;6, 1, 9];
relaxed{2} = [1, 1, 3;2, 1, 3;2, 1, 8;3, 1, 7;4, 1, 7;5,1,8;6, 1, 8];
relaxed{3} = [1, 1, 4;2, 1, 4;2, 1, 7;3, 1, 8;4, 1, 8;5,1,7;6, 1, 7];
relaxed{4} = [1, 1, 5;2, 1, 5;2, 1, 6;3, 1, 9;4, 1, 9;5,1,6;6, 1, 6];
dsName = 'AvgRelaxed';
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
figure(2);clf;hold on;
clear peaks;

for i_rds = 1:length(rds)
    sp =subplot(4, 4, (i_rds-1)*4 +  (1:2));cla;
    outF = [];n = 1;clear leg;
    rampSet = relaxed{i_rds};

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
        
        % base on absolute peak
        base_rel = max(rmp.F);
        % base on steady state
        % base_rel = rmp.F(end);
        F = rmp.F/base_rel;
        % F = rmp.F;

        plot(rmp.t, F, ':');hold on;
        set(gca, 'FontSize', 14);
        % base on obsolute peak
        peaks(i_rds, i_logtrace) = max(rmp.F);
        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        if isempty(outF)
            outF = F;
            Fmax = base_rel;
            n = 1;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF = outF(1:L) + (F(1:L) - outF(1:L))/n;
            outT = rmp.t(1:L);
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;
            
        end    

        %%
        fitrg = rmp.L > 1.1 & rmp.t < dtst.rd;
        rmpFitrg = rmp(fitrg, :);
        % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        fitfun = @(a, b, x) a.*(x + b);

        [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1]);
        as(i_rds, i_logtrace) = ae.a;

        
    end
%% resample and save
    
    % resample up the peak and for the tail separately
    t_s = [linspace(0, 1, 20)*rds(i_rds) ...
        logspace(log10(rds(i_rds)), log10(outT(end)), 40)];
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
    plot(t_s, FLint(:, 1)/Fmax, '-|', LineWidth=2)
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
%% peak sum up
subplot(4, 4, [3 16]);cla;
clin = lines(size(peaks, 2)+1);
for i_pk = 1:size(peaks, 2)
    % set(gca, 'colororderindex', 1);
    semilogx(rds, peaks(:, i_pk), '.:', 'MarkerSize',12, LineWidth=1.5, Color=clin(i_pk, :));hold on;    
    % semilogx(rds, as(:, i_pk), 'x:', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
end
% set(gca, 'YAxisLocation', 'right');
% just for the legend
semilogx(NaN, NaN, '-|',LineWidth=3, Color=clin(end, :))
legend(leg, 'Interpreter','none', 'AutoUpdate','off', 'Location','northeast')

semilogx(rds, mean(peaks, 2)', '-', LineWidth=3, Color=clin(end, :))
semilogx(rds, mean(peaks, 2)', '|', LineWidth=5, Color=clin(end, :), MarkerSize=12)


xlabel('Ramp duration (s)');
xlim([0.08, 150])
ylabel('Peak tension (kPa)')
% semilogx(rds, mean(peaks, 2)', '_', LineWidth=3, MarkerSize=12)
set(gca, 'XTick', fliplr(rds));
set(gca, 'FontSize', 14);
title('Absolute peak height')

% boxplot(peaks', 'Positions',rds)
% end
%% the same, but for pCa
