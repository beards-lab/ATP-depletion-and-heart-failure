function [cost_dt cost_ktr] = fitSlackForceOnset(datatable, velocitytable,t, SL, Force, plotData)

% modeldatatable = [out.t; out.SL; out.Force]';
plotData_fit = plotData;
% plotData_fit = false;

modeldatatable = [t; SL; Force]';

dropstart = velocitytable([3, 7, 11, 15, 19], 1);

zones = [1162, 1209;1464 1519;1816 1889;2268.5 2359.5;2766.5 2900];
zones1 = zones(:, 1);
for z1i = 1:length(zones1)
    z1u(z1i) = find(t > dropstart(z1i) + 0.002 & t < zones(z1i, 2)/1000 & Force > 1e-3, 1, 'first');    
end
zones(:, 1) = t(z1u)*1000;
%
if plotData_fit
    nexttile;hold on;
end
% zones = [1162, 1209;1464 1519;1816 1889;2269 2359.5];
[dSLpc, ktr, df, del, E, SL, x0lin]  = fitRecovery(modeldatatable, zones, 0, [], plotData_fit);
    
    % zones_d = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];

if plotData_fit
    plot([1 3], [0 0], 'k-')    
end
% zones_d = [1162, 1209;1464 1519;1816 1889;2269 2359.5;2774 2900];
% [dSLpc_d, ktr_d, df_d, del_d, E_d, SL_d, x0lin_d]  = fitRecovery(datatable, zones_d, 0, [], plotData);
ktr_d = [39.7759   34.0831   31.0308   27.6087   26.2495];

SL_d = [2.0400    2.0000    1.9600    1.9200    1.8800];
x0lin_d = [1.1616    1.4637    1.8162    2.2691    2.7724];

dt_d = x0lin_d' - dropstart; 
dL_d = 2.2 - SL_d;    
v_d = dL_d'./dt_d;


% times of start of the SL drop
% dropstart = velocitytable([3, 7, 11, 15], 1);

dt = x0lin' - dropstart; 
dL = 2.2 - SL;

v = dL'./dt;
%

if plotData    
    nexttile;hold on;
    plot(modeldatatable(:, 1)-dropstart', modeldatatable(:, 2));
    leg_m = plot([3e-4 dt'], [2.2 SL], '*-', LineWidth=2);
    xlim([-0.05, 0.25])

   leg_d =  plot([3e-4 dt_d'], [2.2 SL_d], '*-', LineWidth=2);
   legend([leg_m, leg_d], 'model', 'Data');

    nexttile;
    plot(datatable(:, 1)-dropstart', datatable(:, 3), '-', LineWidth=1);
    hold on; set(gca,'ColorOrderIndex',1);

    plot(modeldatatable(:, 1)-dropstart', modeldatatable(:, 3), ':', LineWidth=2);
    xlim([-0.05, 0.25])
end

slack_x = [3e-4 dt'] - 3e-4;
slack_y = [2.2 SL];

cost_dt = sum((dt_d - dt).^2);
cost_ktr = sum((ktr_d - ktr).^2);