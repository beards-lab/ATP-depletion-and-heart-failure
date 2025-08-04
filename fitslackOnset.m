clf;
dropstart = velocitytable([3, 7, 11, 15, 19], 1);

ktr_d = [39.7759   34.0831   31.0308   27.6087   26.2495];

SL_d = [2.0400    2.0000    1.9600    1.9200    1.8800];
x0lin_d = [1.1616    1.4637    1.8162    2.2691    2.7724];

dt_d = x0lin_d' - dropstart; 
dL_d = 2.2 - SL_d;    
v_d = dL_d'./dt_d;

dt_set = 0:0.001:0.015;
FOfun = @(t) 402*t.^2 - 20.8*t + 2.085;


leg_d =  plot([1e-6 dt_d'], [2.2 SL_d], '*-', dt_set, FOfun(dt_set), '--', LineWidth=2);

