load pca4data.mat
load pca11data.mat
% load pca4dataAdj.mat
% load pca4dataAdj60s.mat

% dsc = load('DataStruct20240705.mat').dsc;
i = 1;rd = {100, 10, 1, 0.1};
figure(80085);
clf
vz = Tarr{i} > rd{i}*0.1 & Tarr{i} < rd{i}*1.2;
ML = 0.95 + (1.175 - 0.95)*max(0, min(rd{i}, Tarr{i}))/rd{i};



nexttile;
plot(Tarr{i}, Farr{i}, Tarr{i}(vz), Farr{i}(vz))

% nexttile;
% plot(Tarr{i}, ML);

nexttile;
plot(ML, Farr{i}, ML(vz), Farr{i}(vz), LineWidth=2);hold on;

y = @(k_pas, x0, gamma, z, x) k_pas.*(x-x0).^gamma + z*0 + z*1e9.*(x-0.95).^13;
% y = @(k_pas, x0, gamma, x) k_pas.*exp((x-x0)*gamma);

[a, b, c] = fit(ML(vz), Farr{i}(vz), y, 'StartPoint', [50, 0.5, 4, 0.5]);

plot(ML, y(a.k_pas, a.x0, a.gamma, a.z, ML), ML, y(0, 0, 0, a.z, ML),LineWidth=2)

x_ax = 0.8:0.01:1.2;
% plot(x_ax, y(a.k_pas, a.x0, a.gamma, x_ax))
Force_ML_BestFit = @(x_ax) y(0.4, -0.4, 7.9, 0.5, x_ax);
plot(x_ax, Force_ML_BestFit(x_ax), '--', LineWidth=2)

SL = ML*2;

% plot(ML, params.k_pas*max(SL - params.Lsc0, 0).^params.gamma, '--'); 

nexttile(1);
hold on;plot(Tarr{i}, y(a.k_pas, a.x0, a.gamma, a.z, ML), '--', LineWidth=2)
plot(Tarr{i}, Force_ML_BestFit(ML), '--', LineWidth=2)
xlim([1e-1, 200])
a
b.rmse
% plot(Tarr{i}(vz), y(a.k_pas, a.x0, a.gamma, ML(vz)))

%% high ca for comparison

load pca4data.mat
vz = Tarr{i} > rd{i}*0.1 & Tarr{i} < rd{i}*1.2;
ML = 0.95 + (1.175 - 0.95)*max(0, min(rd{i}, Tarr{i}))/rd{i};

[a_4, b_4, c_4] = fit(ML(vz), Farr{i}(vz), y, 'StartPoint', [50, 0.5, 4, 0.5]);

nexttile(2)
plot(ML, Farr{i}, ML(vz), Farr{i}(vz), LineWidth=2);hold on;
plot(ML, y(a_4.k_pas, a_4.x0, a_4.gamma, a_4.z, ML),LineWidth=2)
a
b.rmse
legend('pCa 11 - data', 'pCa 11 - fit area', ...
    sprintf('Fit %1.1f(ML + %1.1f)^{%1.1f} + %1.1f*F_{coll}', a.k_pas, -a.x0, a.gamma, a.z), 'F_{coll}',...
    sprintf('Simplified function %1.1f(ML + %1.1f)^{%1.1f} + %1.1f*F_{coll}', 0.4, 0.4, 7.9, 0.5), 'pCa 4 - data', 'pCa 4 - fit area', ...
    sprintf('pCa 4 fit %1.1f(ML + %1.1f)^{%1.1f} + %1.1f*F_{coll}', a_4.k_pas, -a_4.x0, a_4.gamma, a_4.z))


fontsize(14, "points")
xlabel('ML (L/L0)');ylabel('\Phi (kPa)')