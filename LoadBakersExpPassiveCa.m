%% read new Anthony's data
close all;

rds = [0.1 1 10 100];
% ramp height 0.95 - 1.175, i.e. 1.9 - 2.35um
S1 = dir('data/PassiveCaSrc/20230518');
S1 = S1(~[S1.isdir]);
[~,idx] = sort([S1.datenum]);
S1 = S1(idx);

% ramp height 0.95 - 1.175, i.e. 1.9 - 2.35um
S1 = dir('data/PassiveCaSrc2/20230518');
S1 = S1(~[S1.isdir]);
[~,idx] = sort({S1.name});
S1 = S1(idx);
rds = fliplr([0.1 1 10 100]);

%{
% ramp height 0.95 - 1.175, i.e. 1.9 - 2.35um
% to be upscaled by 1.5x
S2 = dir('data/PassiveCaSrc/20230608');
S2 = S2(~[S2.isdir]);
[~,idx] = sort([S2.datenum]);
S2 = S2(idx);

% ramp height 1.0 to 1.1 ML
S3 = dir('data/PassiveCaSrc/20230616');
S3 = S3(~[S3.isdir]);
S3 = S3(~contains({S3.name}, '0.02'));
[~,idx] = sort([S3.datenum]);
S3 = S3(idx);
%}

% mergedTables = [struct2table(S1);struct2table(S2)];
% S = table2struct(mergedTables);
S = S1;
%%


close all;
clear peaks ss;
clear legnames;
clear timecourses;


% deifning the sequence
% seq = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
seq = [1 1 1 1 1]'*[1:12];
seq = seq(:)';

skipPlots = false;
for i = 1:length(S)
    datatable = readtable([S(i).folder '/' S(i).name], 'filetype', 'text', 'NumHeaderLines',4);
    datatable.Properties.VariableNames = {'t', 'L','F', 'SL'};

    if strcmp(S(i).name, '30_extracted_Log.txt')
        % some strange end messes with zero drift correction
        datatable = datatable(datatable.t < 1000, :);
    % elseif strcmp(S(i).name, '15_pCa6_PNB_Log.txt')
    end

    % subtrack zero drift
    % i_fzero = datatable.t > datatable.t(end) - 35 & datatable.t < datatable.t(end) - 25;
    % Fzero = mean(datatable.F(i_fzero));
    % datatable.F = (datatable.F - Fzero);
    datatable.F = movmean(datatable.F, [8 8]);
    i_ss = datatable.t > datatable.t(end) - 46.0 & datatable.t < datatable.t(end) - 45.5;% index of steady states
    ss_cur = mean(datatable.F(i_ss));

    seq_cur = i - seq(i)*5 + 5;
    % timecourses{seq(i), seq_cur} = [datatable.t, datatable.F];
    timecourses{seq(i), seq_cur} = datatable;
    

    if seq_cur == 5
        % log of the whole experiment
        % continue;
    end

    % cut out the title
    tit = split(S(i).name, {'_', '.txt'});
    legnames{seq(i)} = [num2str(seq(i)) ':' strjoin(tit(3), '_')];    


    peaks(seq(i), seq_cur) = max(datatable.F);        
    ss(seq(i), seq_cur) = ss_cur;

    if skipPlots
        continue
    end
    figure(seq(i));
    set(gcf, 'Position',  [769.8000   41.8000  766.4000 740.8000]);

    subplot(211);hold on;plot(datatable.t, datatable.L);
    subplot(212);hold on;
    plot(datatable.t, datatable.F);
    % plot(datatable.t, F);
    plot(datatable.t(i_ss), repmat(ss_cur, [sum(i_ss), 1])+2, 'LineWidth',2);
    title(legnames{seq(i)}, 'Interpreter','None');
    % ylim([0 10])
    % pause
end

%% Filter out zeros in a separate pass
for i_logtrace = 1:size(timecourses, 1)
    %%
    % i_logtrace = 6
    % is slack? ML below 0.85 definitely is
    is_slack = timecourses{i_logtrace, 5}.L < 0.85; 
    % allow first 2s to be considered as slack too. Usually over zero
    % though
    % is_slack = is_slack | timecourses{i_logtrace, 5}.t < 2; % is slack?

    zdt = timecourses{i_logtrace, 5}.t(is_slack); % zero drift time
    zdF = timecourses{i_logtrace, 5}.F(is_slack); % zero drift Force

    % Force drift as a function of time
    f_Fdt = @(a, b, c, x)0*a.*x.^2 +b.*x + c;        
    [ae be] = fit(zdt, zdF, f_Fdt, 'StartPoint', [1e-6, 1e-4, 1]);        
    % plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
    zd = f_Fdt(ae.a, ae.b, ae.c, timecourses{i_logtrace, 5}.t); % zero drift

    figure(i_logtrace); clf; 
    subplot(222)
    plot(timecourses{i_logtrace, 5}.t, timecourses{i_logtrace, 5}.F, ...
        timecourses{i_logtrace, 5}.t(is_slack), timecourses{i_logtrace, 5}.F(is_slack), ...
        timecourses{i_logtrace, 5}.t, zd);
    legend('Raw Force reading', 'Zero regions', 'Zero drift');

    %% Identify position of individual ramps
    % We use start of the slack, which beggins 4s from the end of indi ramp
    % of the individual cut-out
    i_slackStart = find(diff(is_slack) > 0);
    t_slackStart = (timecourses{i_logtrace, 5}.t(i_slackStart));    
    % take last 4, some traces have some weird beggining
    t_slackStart = t_slackStart([end-3:end]);
    % compare to the 100s, 10s, 1s and 0.1s ramp-up cutouts
    figure(i_logtrace);subplot(221);cla;hold on;
    zd = cell(1, 5);
    for i_ramp = 1:4
        % i_ramp = 3;
        dt(i_ramp) = t_slackStart(i_ramp) - (timecourses{i_logtrace, i_ramp}.t(end) - 36.4);
        t_slack = timecourses{i_logtrace, i_ramp}.t(end) - 36.4 + [2, 10];
        i_slack = timecourses{i_logtrace, i_ramp}.t > t_slack(1) & timecourses{i_logtrace, i_ramp}.t < t_slack(2);
        avg_Fslack(i_ramp) = mean(timecourses{i_logtrace, i_ramp}.F(i_slack));
        plot(timecourses{i_logtrace, i_ramp}.t + dt(i_ramp), timecourses{i_logtrace, i_ramp}.L)
        plot(timecourses{i_logtrace, i_ramp}.t(i_slack) + dt(i_ramp), timecourses{i_logtrace, i_ramp}.L(i_slack), '*--')
    end
    plot(timecourses{i_logtrace, 5}.t, timecourses{i_logtrace, 5}.L, ':')
    legend('Muscle length', 'Zero regions')
    
    figure(i_logtrace);subplot(212);cla;hold on;
    for i_ramp = 1:4
        zd{i_ramp} = f_Fdt(ae.a, ae.b, ae.c, timecourses{i_logtrace, i_ramp}.t + dt(i_ramp)); % zero drift
        % zero drift fitted
        plot(timecourses{i_logtrace, i_ramp}.t + dt(i_ramp), timecourses{i_logtrace, i_ramp}.F - zd{i_ramp})
        % compare with simple cut-out
        plot(timecourses{i_logtrace, i_ramp}.t + dt(i_ramp), timecourses{i_logtrace, i_ramp}.F - avg_Fslack(i_ramp), ':')        

        % save the force corrected for the zero drift
        % timecourses{i_logtrace, i_ramp}.F = timecourses{i_logtrace, i_ramp}.F - zd;
    end
    % plot(timecourses{i_logtrace, 5}.t(is_slack), timecourses{i_logtrace, 5}.F(is_slack), '.', 'LineWidth',2)
    % plot(timecourses{i_logtrace, 5}.t, f_Fdt(ae.a, ae.b, ae.c, timecourses{i_logtrace, 5}.t))
    % plot(timecourses{i_logtrace, 5}.t, timecourses{i_logtrace, 5}.F, ':')
    legend('Zero drift removal', 'zero-value removal')

    axis tight;


end

%% plot All the peaks and ss

% rds = [0.1 1 10 100];
figure(21);clf; 
markers = 'sd<^>vox.ph*';
colors = jet(size(peaks, 1));

for i = 1:size(peaks, 1)
    % subplot(211);
    semilogx(rds, (peaks(i, :)), [markers(i) '-'], 'LineWidth',2, 'Color', colors(i, :));
    hold on;    
end
title('All peaks and SS')
legend(legnames, 'AutoUpdate',false)

% separate loop to separate the legends
for i = 1:size(peaks, 1)
% subplot(212);
    semilogx(rds, (ss(i, :)), [markers(i) '--'], 'LineWidth',1, 'Color', colors(i, :));
    hold on;
end

%% compare just those of interest: relaxed vs extracted
intrs = [1 4 5]
figure(22);clf;

for i = intrs
    % subplot(211);
    semilogx(rds, fliplr(peaks(i, :)), [markers(i) '-'], 'LineWidth',2, 'Color', colors(i, :));
    hold on;    
end
title('Relaxed vs extracted')
legend(legnames{intrs}, 'AutoUpdate',false, 'Interpreter', 'None')

% separate loop to separate the legends
for i = intrs
% subplot(212);
    semilogx(rds, fliplr(ss(i, :)), [markers(i) '--'], 'LineWidth',1, 'Color', colors(i, :));
    hold on;
end

%% compare just those of interest: scale-up the second experiment
intrs = [1 3 6 7 11]
figure(23);
subplot(211);cla;subplot(212);cla;
mult = 1.5; % multiplier for the second data set
for i = intrs
    subplot(211);
    if i < 6 % first experiment series
        semilogx(rds, fliplr(peaks(i, :)), [markers(i) '-'], 'LineWidth',2, 'Color', colors(i, :));
    else 
        semilogx(rds, fliplr(peaks(i, :))*mult, [markers(i) '-'], 'LineWidth',2, 'Color', colors(i, :));
    end
    hold on;    
end
title(sprintf('Second set of experiments: scaling up by %0.0f %%', mult*100));
legend(legnames{intrs}, 'AutoUpdate',false, 'Interpreter', 'None')

% separate loop to separate the legends
for i = intrs
subplot(212);
    if i < 6 % first experiment series
        semilogx(rds, fliplr(ss(i, :)), [markers(i) '--'], 'LineWidth',1, 'Color', colors(i, :));
    else
        semilogx(rds, fliplr(ss(i, :))*mult, [markers(i) '--'], 'LineWidth',1, 'Color', colors(i, :));
    end
    hold on;
end

%% Compare extracted time-courses
figure(24);clf; hold on;
mult = 1.5;
for i = [1 3 6 7]
    for j = 1:4
        subplot(2,2,j);hold on;
        % Scale equally
        % plot(timecourses{i, j}(:, 1), movmean(timecourses{i, j}(:, 2), [16*j^1.5*2]))
        
        % compare scaling the other experiment set
        if i < 6 % first experiment series
            plot(timecourses{i, j}(:, 1), movmean(timecourses{i, j}(:, 2), [8]))
        else
            plot(timecourses{i, j}(:, 1), mult*movmean(timecourses{i, j}(:, 2), [8]))
        end
        title(sprintf('Ramp %0.1fs', rds(5-j)));
        xlim([10 40 + rds(5-j)])
        ylim([0 40])
    end        
end
subplot(2,2,1);legend(legnames{[1 3 6 7]}, 'Interpreter', 'None');

%% stiffness at Ca levels
figure(28);clf;hold on;
colors = lines(5);
pCa = [4.4, 6, 11, 12]

indx = [4, 3, 5, 1]
% pCa = [4.4, 5.5, 5.75, 6, 6.25, 11, 12];
% indx = [2 3 4 5 6 7 1];
plot(-repmat(pCa, [4,1])', peaks(indx, :), 's', 'Linewidth', 1); 

% TODO fix me
ramp01 = peaks(indx, 4);

yf = @(L, c, k, x0, x)L./(1+exp(-k*(-x-x0)))+c;
t_f = 3:0.2:13;% time fit 

[ae, be] = fit(pCa', ramp01, yf, 'StartPoint', [6.0000, 3.0000, 5.7500, -5.0000]);
plot(-t_f, yf(ae.L, ae.c, ae.k, ae.x0, t_f), '--', 'LineWidth',1, 'Color', colors(4, :));

[ae, be] = fit(pCa', ramp02, yf, 'StartPoint', [6.0000, 3.0000, 5.7500, -5.0000]);
plot(-t_f, yf(ae.L, ae.c, ae.k, ae.x0, t_f), '--', 'LineWidth',1, 'Color', colors(3, :))

[ae, be] = fit(pCa', ramp03, yf, 'StartPoint', [6.0000, 3.0000, 5.7500, -5.0000]);
plot(-t_f, yf(ae.L, ae.c, ae.k, ae.x0, t_f), '--', 'LineWidth',1, 'Color', colors(2, :))

[ae, be] = fit(pCa', ramp04, yf, 'StartPoint', [6.0000, 3.0000, 5.7500, -5.0000]);
plot(-t_f, yf(ae.L, ae.c, ae.k, ae.x0, t_f), '--', 'LineWidth',1, 'Color', colors(1, :))

% t_i = 4:0.1:12;
% F_i = interp1(pCa,peaks(indx, end),t_i,'makima');
% plot(-t_i, F_i);
xlabel('-pCa');
ylabel('Tension (kPa)')
legend('100s', '10s', '1s', '0.1s');

%% Overlap the decays
for i = 1:size(timecourses, 1)
    %%
    i = 1;
    figure(i);clf;hold on;
    for i_rd = 1:length(rds)
      plot(timecourses{i, i_rd}.t - rds(i_rd), (timecourses{i, i_rd}.F  - ss(i, i_rd))/(peaks(i, i_rd) - ss(i, i_rd)));
    end
end



%% Resample and save

intrs = [2 3 4]
pCas = [11, 6, 4];
% rds = [100, 10, 1, 0.1];
% ts_d = 
clf; hold on;
for i = 1:length(intrs)
    for i_rd = 1:length(rds)
        %%
        clf; hold on;
        rd = rds(i_rd);
        datatable_cur = timecourses{intrs(i), i_rd};

        % correct for zero drift
        datatable_cur.F = datatable_cur.F - zd{i_rd};
        
        % cut out the ramp and decay
        % validIds = datatable_cur(:, 1) >= 8 & datatable_cur(:, 1) < rd + 10 + 30;
        validIds = datatable_cur.t >= 8 & datatable_cur.t < rd + 10 + 30;
        datatable_cur = datatable_cur(validIds, :);
        % shift the time base
        datatable_cur(:, 1) =  datatable_cur(:, 1) - 8;


        plot(datatable_cur.t,datatable_cur.F, '--');
        % sample times for datas; baseline, ramp-up, peak decay, long tail
        dwnsmpl = [0:0.5:2, 2 + (0:max(2,rd/30):rd), 2 + rd + (0:max(2,rd/30):min(30,rd)), (2 + 2*rd):0.5:(30 + rd + 2 - 1.0)];

        datatable_cur.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
        i_data = interp1(datatable_cur.Time, 1:length(datatable_cur.Time), dwnsmpl, 'nearest', 'extrap');
        
        dsf = 8;scaleF = 1;
        datatable_cur.F = movmean(datatable_cur.F*scaleF,[dsf/2 dsf/2]); % force filtered
        % datatable_cur = 
        
        plot(dwnsmpl, datatable_cur.F(i_data), 'x-', 'Linewidth', 2);
        % 
        % plot(dwnsmpl, zeros(1, length(dwnsmpl)), '|-', 'Linewidth', 2);
        % DownSampleAndSplit(datatable, dwnsmpl, [], 2.0, 1, 1.5, filename, 0, 1);

        saveAs = ['bakers_passiveStretch_pCa' num2str(pCas(i)) '_' num2str(rds(i_rd)*1000) 'ms'];
        writetable(datatable_cur, ['data/PassiveCa_3/' saveAs '.csv']);
        disp(['Saved as ' saveAs ' and csv'])
    end
end


%% fft?
ms = movstd(datatable.F, 5);
figure(33); plot(datatable.t, ms); 
Y = fft(datatable.F);
L = length(datatable.F);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = 1/mean(diff(datatable.t));
f = Fs*(0:(L/2))/L;
plot(f,P1)