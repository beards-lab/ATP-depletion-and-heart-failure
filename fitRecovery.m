function [ktr, df, st] = fitRecovery(datatable, zones)

    clf; subplot(1, 4, [1 3]);hold on;
    plot(datatable(:, 1), datatable(:, 3))

    % zones = [1180, 1209;1485 1519;1839 1889;2290 2355;2794 2900];
    for z = 1:size(zones, 1)
        
        %% exp fit
        bt = 15; % buffer time (ms)
        z1 = zones(z, :);
        i1 = find(datatable(:, 1) > (z1(1) + bt)/1000, 1);
        i2 = find(datatable(:, 1) > z1(2)/1000, 1);
        z1 = [i1:i2];

        to = datatable(z1(1), 1); % time offset
        timebase = datatable(z1, 1)-to;
        
        y_exp = @(df, ktr, s, x)df*(1-exp(-(x-s)*ktr));
        [ae be] = fit(timebase, datatable(z1, 3), y_exp, 'StartPoint', [50, 20, bt/1000]);
        ktr(z) = ae.ktr; 
        df(z) = ae.df; % difference in force
        del(z) = ae.s; % delay time [s]
        
        
        %% linear approx
%         ilin1 = (find(datatable(:, 1) > (z1(1))/1000, 1));
%         zlin = ilin1:ilin1+20; % linear zone length
%         timebase_lin = datatable(zlin, 1)-datatable(zlin(1), 1);
%         y_line = @(k, x0, y0, x)k.*(x-x0) + y0;
%         [al bl] = fit(timebase_lin, datatable(zlin, 3), y_line, 'StartPoint', [1000, 0.05, 1]);
%         fprintf('Delay %1.1f ms (exp) vs %1.1f ms (lin)', ae.s*1000, al.x0*1000);
%         del2(z) = al.y0/al.k +al.x0 + datatable(z1(1));


        timebase_fit = (-bt/1000:0.01:0.3);
%         timebase2_fit = (-30/1000:0.01:0.1);
        plot(timebase+to, datatable(z1, 3), 'x', timebase_fit + to, y_exp(ae.df, ae.ktr, ae.s, timebase_fit), '--', 'Linewidth', 2);
        ci = get(gca,'ColorOrderIndex');
        set(gca,'ColorOrderIndex', max(ci-1, 1));
%         plot(timebase_lin+to, datatable(zlin, 3), 'o', timebase2_fit + to, y_line(al.k, al.x0, al.y0, timebase2_fit), ':', 'Linewidth', 2);
        
        
        plot([datatable(i1, 1) datatable(i1, 1);datatable(i2, 1) datatable(i2, 1)]', [0, 2.2; 0 100]')
        % find max slack velocity from data
        si = 100;% search indices
        [~, imin] = min(diff(datatable(i1-si:i1, 2))); % SL reached the bottom;
        datatable(imin + i1 - si, 1);
        sc = find(...
            datatable(i1:-1:i1-si, 2) > datatable(i1-si, 2)*0.99 & ... % one percent drop
            (si:-1:0 < imin)', ... % search until the bottom is reached
            1, 'first') -1;
%         datatable(i1 - si, 1)
%         plot([datatable(i1- si + imin, 1) datatable(i1-si + imin, 1)], [0, 100])
        plot([datatable(i1- sc, 1) datatable(i1-sc, 1)], [0, 100], 'r')
        % get the velocity
%         i_del = find(datatable(i1-si:i1+si, 1) >= datatable(i1) + del(z), 1);
        
        i_del = find(datatable(i1-si:i2+si, 1) >= datatable(i1) + del(z), 1);
%         i_del2 = find(datatable(i1-si:i2+si, 1) >= datatable(i1) + del2(z), 1);
        plot([datatable(i1- si + i_del, 1) datatable(i1-si + i_del, 1)], [0, 100], 'r')
%         plot([datatable(i1- si + i_del2, 1) datatable(i1-si + i_del2, 1)], [0, 100], ':r', 'Linewidth', 3)
        dSL = -datatable(i1- sc, 2) + datatable(i1 - si + i_del, 2);
        dT = (datatable(i1- sc, 1) - datatable(i1 - si + i_del, 1));
        vel = dSL/dT;
        
        fprintf('At ML %1.2f, and time %1.2f, dSL = %1.3f (%1.1f %%), dT = %1.1f ms and vel is %1.1f \n', ...
            datatable(i1- sc, 2), datatable(i1- sc, 1), dSL, dSL/2.0*100, dT*1000, vel);
        
        sv(z) = vel; % slack velocity
        text(to + 0.001, datatable(z1(1), 3), sprintf('ktr = %1.1f, fm=%0.1f,\n with rmse %0.3f \n, vel = %0.1f', ae.ktr, ae.df, be.rmse, vel), 'Color', [1 0 0])
        xlabel('Time (s)');
        ylabel('Force (kPa)')
    end
    title('1-exp(-x) exponential recovery identification')
    subplot(1, 4, 4);hold on;
    plot(ktr, 'o-', 'Linewidth', 2);plot(df, 'x-', 'Linewidth', 2);
    ylabel('Ktr and max force');
    yyaxis right;plot(sv, '|-', 'Linewidth', 2)
    legend('ktr (s^{-1})', 'max force (kPa)', 'slack velocity (um/s)*');
    ylabel('Velocity (um/s)');
    xlabel('Experiment nr')
% plot(df, ktr, 'o-');
end