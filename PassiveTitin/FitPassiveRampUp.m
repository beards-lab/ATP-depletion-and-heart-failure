load pca4data.mat
load pca11data.mat
% load pca4dataAdj.mat
% load pca4dataAdj60s.mat

% dsc = load('DataStruct20240705.mat').dsc;
i = 2;rd = {100, 10, 1, 0.1};
figure(80085);clf;
vz = Tarr{i} > rd{i}*0.1 & Tarr{i} < rd{i}*0.9;
ML = 0.95 + (1.175 - 0.95)*min(rd{i}, Tarr{i})/rd{i};



nexttile;
plot(Tarr{i}, Farr{i}, Tarr{i}(vz), Farr{i}(vz))

nexttile;
plot(Tarr{i}, ML);

nexttile;
plot(ML, Farr{i}, ML(vz), Farr{i}(vz), LineWidth=2);hold on;

y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;
% y = @(k_pas, x0, gamma, x) k_pas.*exp((x-x0)*gamma);


[a, b, c] = fit(ML(vz), Farr{i}(vz), y, 'StartPoint', [50, 0.7, 3]);

plot(ML, y(a.k_pas, a.x0, a.gamma, ML), ML, y(0, 0, 0, ML),LineWidth=2)

x_ax = 0.8:0.01:1.2;
% plot(x_ax, y(a.k_pas, a.x0, a.gamma, x_ax))
plot(x_ax, y(0.4, -0.4, 7.9, x_ax), '--', LineWidth=2)

SL = ML*2;

% plot(ML, params.k_pas*max(SL - params.Lsc0, 0).^params.gamma, '--'); 

nexttile(1);
hold on;plot(Tarr{i}, y(a.k_pas, a.x0, a.gamma, ML), '--', LineWidth=2)
a

% plot(Tarr{i}(vz), y(a.k_pas, a.x0, a.gamma, ML(vz)))