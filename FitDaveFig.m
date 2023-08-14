% fit dave's output

figure(1);clf;

% % check Dave's peaks
% semilogx(0.1, 1, 'r+');
% hold on;
% semilogx(1, 0.76, 'r+');
% semilogx(10, 0.5, 'r+');
% semilogx(100, 0.35, 'r+');
% semilogx(5000, 0.15, 'r+');

% hold on;



rds = [0.02, 0.1, 1, 10, 100];
colors = colormap(lines(length(rds)));
normalizer = 10.6611; % peak value for 100ms ramp - this is our normalized max
% normalizer = 25.6611; % peak value for 100ms ramp - this is our normalized max
for rd_i = 2:length(rds)
    rd = rds(rd_i);
    datatable = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']).datatable;
    semilogx(datatable(:, 1) - 2, datatable(:, 3)/normalizer, 'color',colors(rd_i, :), LineWidth=2);
hold on;
end

% Draw the Ca experiments
data = readtable('data/PNB_dataset.csv');
semilogx(data.RampDuration, data.Relax/max(data.Relax), 'ro', LineWidth=2);
semilogx(data.RampDuration, data.pCa6/max(data.pCa6), 'gs', LineWidth=2);
semilogx(data.RampDuration, data.pCa4/max(data.pCa4), 'bd', LineWidth=2);

% draw steady state of PNBs
semilogx(500, 4.5/max(data.Relax), 'r_', LineWidth=2);
semilogx(500, 4.5/max(data.pCa6), 'g_', LineWidth=2);
semilogx(500, 4.5/max(data.pCa4), 'b_', LineWidth=2);

legend('Data 0.1s', 'Data 1s', 'Data 10s', 'Data 100s',...
    'AutoUpdate','off')


xlim([5e-3 1e4]);
ylim([0 1]);

% draw Dave's plot
I = imread('data\FitDave.png'); 
h = image(xlim,flip(ylim),I); 
uistack(h,'bottom');

title('Ramp-up (0.8 - 1.2 ML = 1.6 - 2.4 um SL) force response');
xlabel('Time (s)');ylabel('Normalize tension');
