function [ktr, df, st] = fitRecovery(datatable, zones, zeroTreshold)
    
    clf; subplot(1, 4, [1 3]);hold on;
    plot(datatable(:, 1), datatable(:, 3))
    plot(datatable([1, end], 1), [zeroTreshold zeroTreshold], 'k-');

    % zones = [1180, 1209;1485 1519;1839 1889;2290 2355;2794 2900];
    for z = 1:size(zones, 1)
        
        %% exp fit
        bt = 15; % buffer time (ms)
        z1 = zones(z, :);
        i1 = find(datatable(:, 1) > (z1(1) + bt)/1000, 1);
        i2 = find(datatable(:, 1) > z1(2)/1000, 1);
        z1 = [i1:i2];

        to = datatable(z1(1), 1); % time offset
        timebase_exp = datatable(z1, 1)-to;
        
        y_exp = @(df, ktr, s, x)df*(1-exp(-(x-s)*ktr));
        [ae be] = fit(timebase_exp, datatable(z1, 3)- zeroTreshold, y_exp, 'StartPoint', [50, 20, bt/1000]);
        ktr(z) = ae.ktr; 
        df(z) = ae.df; % difference in force
        del(z) = ae.s; % delay time [s]
        
        timebase_exp = (-(bt+5)/1000:0.01:0.3);
        plot(datatable(z1, 1), datatable(z1, 3), 'x', timebase_exp + to, y_exp(ae.df, ae.ktr, ae.s, timebase_exp), '--', 'Linewidth', 2);
        
        %% linear approx
        ilin1 = find(datatable(:, 1) > zones(z, 1)/1000, 1); % start at the zone
        zlin = ilin1:ilin1+20; % linear zone length
        timebase_lin = datatable(zlin, 1)-datatable(zlin(1), 1);
        y_line = @(k, x0, x)k.*(x-x0);
        [al bl] = fit(timebase_lin, datatable(zlin, 3)-zeroTreshold, y_line, 'StartPoint', [1000, -0.001]);
        x0lin(z) = al.x0 + datatable(ilin1);
        timebase_lin = (-10/1000:0.01:0.05); % extending the timebase
        
        ci = get(gca,'ColorOrderIndex');
        set(gca,'ColorOrderIndex', max(ci-2, 1));
        plot(datatable(zlin, 1), datatable(zlin, 3), 'o', timebase_lin + datatable(zlin(1), 1), y_line(al.k, al.x0, timebase_lin), ':', 'Linewidth', 2);
%%

        
        
        
        % cursors - exp zone start and end
%         plot([datatable(i1, 1) datatable(i1, 1);datatable(i2, 1) datatable(i2, 1)]', [0, 50; 0 100]')

        % find max slack velocity from data
        si = 100;% search indices
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
        plot([datatable(iss, 1) datatable(iss, 1)], [0, 100], 'g')
        plot([datatable(ise, 1) datatable(ise, 1)], [0, 100], 'r')
%         plot([datatable(i1- si + i_del2, 1) datatable(i1-si + i_del2, 1)], [0, 100], ':r', 'Linewidth', 3)
        dSL =  datatable(i1, 2) - datatable(iss, 2);
        dSLpc(z) = dSL/datatable(iss, 2)*100; % dSL in percent of ML
        dT(z) =  datatable(ise, 1) - datatable(iss, 1);
        vel(z) = dSL/dT(z);% um/s
        
%         i_del_lin = find(datatable(iss:i2, 1) >= ix0lin(z), 1);
        dTlin(z) = x0lin(z) - datatable(iss, 1);
        vel_lin(z) = dSL/dTlin(z);
        plot([x0lin(z) x0lin(z)], [0, 100], 'm--');
        
        fprintf('At ML %1.2f, and time %1.2f: ' , datatable(i1- sc, 2), datatable(i1- sc, 1));
        fprintf(' dSL = %1.3f (%1.1f pct), dT = %1.1f ms and vel is %1.1f (%1.1f) \n', dSL, dSLpc(z), dT(z)*1000, -vel(z), -vel_lin(z));
        
        text(to + 0.001, datatable(z1(1), 3), sprintf('ktr = %1.1f, fm=%0.1f,\n with rmse %0.3f \n, vel = %0.1f', ae.ktr, ae.df, be.rmse, -vel(z)), 'Color', [1 0 0])
        xlabel('Time (s)');
        ylabel('Force (kPa)')
    end
    yyaxis right;
    plot(datatable(:, 1), datatable(:, 2));
    yyaxis left;
    title('1-exp(-x) exponential recovery identification')
    subplot(1, 4, 4);hold on;
    plot(-dSLpc, ktr, 'o-', 'Linewidth', 2);plot(-dSLpc, df, 'x-', 'Linewidth', 2);
    ylabel('Ktr and max force');
    yyaxis right;
    plot(-dSLpc, -vel/2.0, '|-', 'Linewidth', 2);
    plot(-dSLpc, -vel_lin/2.0, 'x--', 'Linewidth', 2)
    legend('ktr (s^{-1})', 'max force (kPa)', 'slack velocity - exp fit (um/s)*', 'slack velocity - lin fit(um/s)*');
    ylabel('Velocity (ML/s)');
    xlabel('Experiment nr')
% plot(df, ktr, 'o-');

figure;
plot(dT*1000, -dSLpc, 'o-', dTlin*1000, -dSLpc, 'x-');
xlabel('dt (ms)')
ylabel('dsl (%)')
ylim([0, 20])

vel = (dSLpc/100 + 1)./dT

end