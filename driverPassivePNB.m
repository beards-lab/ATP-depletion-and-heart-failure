%% DriverPassivePNB
% Follows experiments from Baker using paranitroblebistation at various Ca
% levels

data = readtable('data/PNB_dataset.csv');

figure(1); clf;
semilogx(data.RampDuration, data.Relax, 'o-', 'LineWidth',2);
hold on;
semilogx(data.RampDuration, data.pCa6, '^-', 'LineWidth',2);
semilogx(data.RampDuration, data.pCa4, 's-', 'LineWidth',2);

xlabel('Ramp duration');
ylabel('Peak force');

% compare with normal
semilogx(rds, peaks_data, 'x-', 'LineWidth',2);
semilogx(rds, ss_data, 'x-', 'LineWidth',2);


%% do a SA on parameters - which increase the peak?
figure(2);
clf;hold on;
% Documentation use
mods = {'r_a', 'r_d', 'mu', 'ks', 'k1', 'c', 'gamma', 'alpha1', 'e'};
% sign of correlation with the Ftot peak - positive value means increasing
% the param value results in increased peak in Ftot
dir =  [  1      -1     1     1    1     1     -1          1      -1];
opt_mods = [1.1568    0.7497    .20208    2.414/5  0.5852    1.0600    1.1421    1.6024    1.0790];

rd = 1;
plotEach = false;
evaluatePassive;

% Ftot_base = Ftot;
opt_mods_base = opt_mods;
plotEach = false;

plot(Tsim, Ftot, '--', 'LineWidth',3);
%%
for i = 1:numel(opt_mods)
    opt_mods = opt_mods_base;
    opt_mods(i) = opt_mods_base(i)*(1 + dir(i)/2);
    evaluatePassive;
    plot(Tsim, Ftot, 'LineWidth',0.5);
end

legend(cat(2, 'Baseline', mods))
cmap = colormap(colorcube(10));hold on;
set(gca, 'ColorOrder', cmap)

%% plot the r_d to peak characteristics
% select only fastest ramp - the first one
i_rd = 1;
rd = data.RampDuration(i_rd);
opt_mods = opt_mods_base;

mrds = 0.2:0.2:1; % set of detachment modifiers
passivePeaks = zeros(1, length(mrds));
for iterDetach = 1:length(mrds)
    opt_mods(2) = opt_mods_base(2)*mrds(iterDetach);
    evaluatePassive;
    passivePeaks(iterDetach) = max(Ftot);
end
clf;
plot(mrds*100, passivePeaks, 's-');
xlabel('r_d modifier (%)');
ylabel('Max peak [kPa]');
title('Detachment rate characteristics on peak response');
%% Plot the data characteristics - this is tough because of the resting zero
figure(4);
clf; 
i_rd = 1;
rd = data.RampDuration(i_rd);

semilogx([1e-9 1e-6 1e-4], [data.Relax(i_rd) data.pCa6(i_rd) data.pCa4(i_rd)], 'o-', 'LineWidth',2);
hold on;
i_rd = i_rd +1;
semilogx([1e-9 1e-6 1e-4], [data.Relax(i_rd) data.pCa6(i_rd) data.pCa4(i_rd)], 'o-', 'LineWidth',2);
i_rd = i_rd +1;
semilogx([1e-9 1e-6 1e-4], [data.Relax(i_rd) data.pCa6(i_rd) data.pCa4(i_rd)], 'o-', 'LineWidth',2);


%% Attempt on fit
% lets start with fastest pCa4 peak 
i_rd = 1;
rd = data.RampDuration(i_rd);

opt_mods = opt_mods_base;

% detachment
% mrd = 0.006; % modifier
% opt_mods(2) = opt_mods_base(2)*mrd;

% e
% mrd = 0.41; % modifier
% opt_mods(9) = opt_mods_base(9)*mrd;

% C
% mrd = 11; % modifier
% opt_mods(6) = opt_mods_base(6)*mrd;

% % attach only is not converging
% mrd = 10000; % modifier
% opt_mods(1) = opt_mods_base(1)*mrd;

% % ks not converging
% mrd = 0.1000; % modifier
% opt_mods(4) = opt_mods_base(4)*mrd;

% k1 - good fit?
mrd = 9; % modifier for pca4
opt_mods(5) = opt_mods_base(5)*mrd;

% % gamma - Does not konverge
% mrd = 0.00001; % modifier
% opt_mods(7) = opt_mods_base(7)*mrd;
% 
% % alpha 1
% mrd = 1.32; % modifier
% opt_mods(8) = opt_mods_base(8)*mrd;

plotEach = false;
evaluatePassive;

E = max(Ftot) - data.pCa4(i_rd)
%%
i_rd = 1;
rd = data.RampDuration(i_rd);
opt_mods = opt_mods_base;

mrd = 9;
% mrd = mrd*4/6*0.7; % log-relation to pca6
opt_mods(5) = opt_mods_base(5)*mrd;

plotEach = false;
evaluatePassive;

E = max(Ftot) - data.pCa6(i_rd)

% Reevaluate
rds = [0.02, 0.1, 1, 10, 100];
peaksFtot = zeros(1, 4)
for iterPeaks = 1:5
    rd = rds(iterPeaks);
    figure(100 + iterPeaks);
    plotEach = true;

    evaluatePassive;
    
    peaksFtot(iterPeaks) = max(Ftot);
end
%
figure(1);
semilogx(rds, peaksFtot, 'x--', 'LineWidth',2);
% legend('Peaks - Relax (Data)', 'Peaks - pCa6', 'Peaks - pCa4', 'Peaks - previous data set', 'Steady state - previous data set', 'Peaks - model')


