function [dSLpc, ktr, df, del, E, SL, x0lin] = fitRecovery(datatable, zones, zeroTreshold, fixed_df, plotData)


    if nargin < 4
        fixed_df = [];
    end

    if nargin < 5
        plotData = true;
    end

    % subplot(2, 1, [1]);hold on;
    % plot(datatable(:, 1), datatable(:, 2));
    % xlim([datatable(1, 1), datatable(end, 1)]);
    % subplot(2, 1, [2]);
    if plotData
        hold on;
        plot(datatable(:, 1), datatable(:, 3))
    end
%     plot(datatable([1, end], 1), [zeroTreshold zeroTreshold], 'k-');
    
    ML0 = 2.0; % sarcomere length at nominal muscle length(um)

    % zones = [1180, 1209;1485 1519;1839 1889;2290 2355;2794 2900];
    for z = 1:size(zones, 1)
        
        %% exp fit
        bt = 5; % buffer time (ms)
        z1 = zones(z, :);
        i1 = find(datatable(:, 1) > (z1(1) + bt)/1000, 1);
        i2 = find(datatable(:, 1) > z1(2)/1000, 1);
        z1 = [i1:i2];

        to = datatable(z1(1), 1); % time offset
        timebase_exp = datatable(z1, 1)-to;
        
        if isempty(fixed_df)
            y_exp = @(df, ktr, s, x)df*(1-exp(-(x-s)*ktr));
        else
            y_exp = @(df, ktr, s, x)fixed_df(z)*(1-exp(-(x-s)*ktr)) + 0*df;
        end

        [ae be] = fit(timebase_exp, datatable(z1, 3)- zeroTreshold, y_exp, 'StartPoint', [50, 2, bt/1000]);
        ktr(z) = ae.ktr; 
        df(z) = ae.df; % difference in force
        del(z) = ae.s; % delay time [s]
        E(z) = be.rmse;
        SL(z) = datatable(i1, 2);
        
        timebase_exp = (-(bt+15)/1000:0.01:0.3);

        
        %% linear approx
        ilin1 = find(datatable(:, 1) > zones(z, 1)/1000, 1); % start at the zone
        zlin = ilin1:ilin1+30; % linear zone length
        timebase_lin = datatable(zlin, 1)-datatable(zlin(1), 1);
        y_line = @(k, x0, x)k.*(x-x0);
        [al bl] = fit(timebase_lin, datatable(zlin, 3)-zeroTreshold, y_line, 'StartPoint', [1000, -0.001]);
        x0lin(z) = al.x0 + datatable(ilin1);
        timebase_lin = (-10/1000:0.01:0.05); % extending the timebase

        if plotData
            plot(datatable(z1, 1), datatable(z1, 3), 'x', timebase_exp + to, y_exp(ae.df, ae.ktr, ae.s, timebase_exp), '--', 'Linewidth', 2);
            ci = get(gca,'ColorOrderIndex');
            set(gca,'ColorOrderIndex', max(ci-2, 1));
            plot(datatable(zlin, 1), datatable(zlin, 3), 'o', timebase_lin + datatable(zlin(1), 1), y_line(al.k, al.x0, timebase_lin), ':', 'Linewidth', 2);
        end
%%
        % cursors - exp zone start and end

        % find max slack velocity from data
        si = 100;% search indices
        if i1-si < 0
            si = i1 - 1;
        end
        [~, imin] = min(diff(datatable(i1-si:i1, 2))); % SL reached the bottom;
%         datatable(imin + i1 - si, 1);
        sc = find(...
            datatable(i1:-1:i1-si, 2) > datatable(i1-si, 2)*0.99 & ... % one percent drop
            (si:-1:0 < imin)', ... % search until the bottom is reached
            1, 'first') -1;
        iss = i1- sc; % index of slack start (from ML)
        % get the index of fitted delay
        i_del = find(datatable(iss:i2, 1) >= datatable(i1, 1) + del(z), 1);
        ise = iss + i_del; % index of slack end (the extrapolation crosses zero)
        
        % cursors - start slack (g) and end slack (r)
%         plot([datatable(iss, 1) datatable(iss, 1)], [0, 100], 'g')
%         plot([datatable(i1- si + i_del2, 1) datatable(i1-si + i_del2, 1)], [0, 100], ':r', 'Linewidth', 3)
        dSL =  datatable(1, 2) - datatable(iss, 2);
        
        % surprisingly, it is percent of nominal ML, not actual elongated
        % one
        % dSLpc(z) = dSL/datatable(iss, 2)*100; % dSL in percent of stretcehd ML
        dSLpc(z) = dSL/ML0*100; % dSL in percent of nominal ML
        
        dT(z) =  datatable(ise, 1) - datatable(iss, 1);
        vel(z) = dSL/dT(z);% um/s
        
%         i_del_lin = find(datatable(iss:i2, 1) >= ix0lin(z), 1);
%         dTlin(z) = x0lin(z) - datatable(iss, 1);
        vel_lin(z) = 0;%sdSL/dTlin(z);
        if plotData
            plot([datatable(ise, 1) datatable(ise, 1)], [0, 100], 'r')
            
            plot([x0lin(z) x0lin(z)], [0, 100], 'm--');
            xlabel('Time (s)');
            ylabel('Force (kPa)')
            xlim([datatable(1, 1), datatable(end, 1)]);
        end
        
%         fprintf('At ML %1.2f, and time %1.2f: ' , datatable(i1- sc, 2), datatable(i1- sc, 1));
%         fprintf(' dSL = %1.3f (%1.1f pct), dT = %1.1f ms and vel is %1.1f (%1.1f) \n', dSL, dSLpc(z), dT(z)*1000, -vel(z), -vel_lin(z));
        
%         text(to + 0.001, datatable(z1(1), 3), sprintf('ktr = %1.1f, fm=%0.1f,\n with rmse %0.3f \n, vel = %0.1f', ae.ktr, ae.df, be.rmse, -vel(z)), 'Color', [1 0 0])
    end
% %% Summary graph    
% % 8mM	
% % X dT (ms),	Y dML (% ML)
% bp8 = [2.735152222	8;
% 4.380311224	10;
% 6.466291958	12;
% 9.056710864	14;
% 12.05655714	16];
% 
% % 2mM	
% % X dT (ms),	Y dML (% ML)
% bp2 = [1.191034575	7.954971857;
% 3.025150961	9.981238274;
% 5.552426411	11.96998124;
% 7.953651442	13.9587242;
% 11.48933025	15.94746717];
% 
% % 02mM	
% % X dT (ms),	Y dML (% ML)
% bp02 = [-9.491620288	7.99249531
% -6.018966686	9.981238274
% -3.302497359	12.00750469
% -1.090347959	13.99624765
% 0.2394485	15.98499062];
% 
% bp = bp2; % show just one Baker's Plot
% 
% dSL_B = bp(:, 2)/100*ML0;
% v_B = dSL_B./bp(:, 1)*1000;
% 
%     subplot(3, 4, [4 8]);hold on;
%     plot(-dSLpc, ktr, 'o-', 'Linewidth', 2);plot(-dSLpc, df, 'x-', 'Linewidth', 2);
%     ylabel('Ktr and max force');
%     yyaxis right;
%     plot(bp(:,2), v_B/ML0, 'x-', 'Linewidth', 2);
%     plot(-dSLpc, -vel_lin/ML0, 's-', 'Linewidth', 2)
%     plot(-dSLpc, -vel/ML0, 'o--', 'Linewidth', 2);
%     legend('ktr (s^{-1})', 'max force (kPa)', 'slack velocity - Bakers fit(um/s)*', 'slack velocity - lin fit(um/s)*', 'slack velocity - exp fit (um/s)*');
%     ylabel('Velocity (ML/s)');
%     xlabel('Experiment nr')
% % plot(df, ktr, 'o-');
% 
% subplot(3, 4, 12);hold on;title('% ML to delay time')
% plot(bp(:, 1),bp(:, 2),'x-',  dTlin*1000, -dSLpc, 's-', dT*1000, -dSLpc, 'o--');
% xlabel('dt (ms)')
% ylabel('Step down (% of ML)')
% ylim([0, 20])
% legend('Bakers fit', 'Linear fit', 'Exponential fit');

% vel = (dSLpc/100 + 1)./dT

end