%% Control room for evaluating the simulations
% init
load gopt;
LoadData;

% Set up environment
MgADP = 0; 
Pi    = 0;
F_active_0 = zeros(length(vel), length(MgATP));
t_sl0 = [0 0.11]; % time at which the SL = 2.2 - time to stop the experiment
t_ss = [0 0.2]; %% steady state time

fcn = @dPUdTCa;
opts = struct('N', 40, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 1);
%% Evaluate state distribution and different speeds
figure(1);clf;
params = struct('Pi', 0, 'MgATP', 2, 'MgADP', 0, 'Ca', 100,'Velocity', 0,'UseCa', true,'UseOverlap', false);
[F, outs] = evaluateModel(fcn, t_ss, params,g, opts);
%% Force velocity profile

figure(4);clf;
subplot(211);
plot(outs.t, outs.F);
title(['Force at velocity ' num2str(params.Velocity)]);
xlabel('time')

subplot(223);
plot(outs.t, outs.NR, outs.t,outs.NP, outs.t, outs.p1_0,outs.t, outs.p2_0, outs.t, outs.p3_0);
title('Zero moment to time')
legend('NR', 'NP', 'P1', 'P2', 'P3');

subplot(224);
plot(outs.t, outs.p1_1,outs.t, outs.p2_1, outs.t, outs.p3_1);
title('First moment to time');


if isfield(outs, 'ps0_t')

    figure(5);clf;
    plot(outs.t, outs.ps0_t);
    title('S(1) maximal value (checking the boundary)')
end

%% draw the force-pCa plot
% 1e-3:100
% x axis: 10 to 2
pRange = 4.5:0.1:7.5;
cRange = 10.^(-pRange+6);% our [Ca] in nM
% ATP = [ MgATP_iso(1:end-1) 2 4 8];
ATP = [0.0100  0.0150  0.0200    0.0500    0.1000 0.5000    2.0000    4.0000    8.0000];
F = zeros(length(ATP), length(cRange));
xbtor = zeros(length(ATP), length(cRange));
params = struct('Pi', 0, 'MgATP', 2, 'MgADP', 0, 'Ca', 100,'Velocity', 0,'UseCa', true,'UseOverlap', false);
for a = 1:length(ATP)
    params.MgATP = ATP(a);
    for i = 1:length(cRange)
        params.Ca = cRange(i);
        [f out] = evaluateModel(fcn, t_ss, params, g, opts); 
        F(a, i) = f(end);
        xbtor(a, i) = out.XB_TOR(end);
    end
end
%%
figure(5);clf;hold on;
subplot(231);cla;hold on;
plotstyle = {'o-', '-', '-', '-', '-', '-', '-', '-', 'x-'};
hf = zeros(size(F, 1), 4);
for i = 1:size(F, 1)
    if i > 5 
        lw = 2;
    else
        lw = 1;
    end
    plot(pRange, F(i,:), plotstyle{i}, 'Linewidth', lw, 'MarkerSize', 4);
    leg{i} = "ATP " + num2str(ATP(i) + " mM");
    c50(i) = interp1(F(i,:), pRange,  max(F(i, :))/2);
    f_max(i) = max(F(i,:));
    f_min(i) = min(F(i, :));
    
    % MAPPING hill coefficient
%     hill_fit = @(b,x)  (x.^b(1))./(b(2)+x.^b(1));

    % Hill func from https://www.sciencedirect.com/science/article/pii/S1056871914002470
%     hill_fit = @(p, x) p(1) + p(2)./(1 + (p(3)./x).^p(4));
    
    % Hill func recommended by Kenneth Campbell
    hill_fit = @(p, Ca) p(1) + p(2) * Ca.^p(4) ./ (p(3).^p(4) + Ca.^p(4));
    b0 = [0; 1; 6;10];                                  % Initial Parameter Estimates
    scaled = F(i, :)/max(F(i, :));
%     B = lsqcurvefit(hill_fit, b0, pRange, scaled);
    B = lsqcurvefit(hill_fit, b0, cRange, scaled);
%     HC(i) = B(1);
    hf(i, :) = B;
    
    
end
set(gca, 'XDir','reverse');
legend(leg, 'Location', 'Best')
title("Force-pCa relation");
xlabel('pCa');
ylabel('Force')

subplot(234);cla;hold on;
for a = 1:size(xbtor, 1)
    if a > 5 
        lw = 2;
    else
        lw = 1;
    end
    plot(pRange, xbtor(a, :), plotstyle{a}, 'LineWidth', lw);title('XB TOR for ATPs');
end
title('XB turnover rate drops with ATP');
xlabel('pCa')
ylabel('XB TOR')
set(gca, 'XDir','reverse');
legend(leg, 'Location', 'Best')

subplot(232);cla;
plot(F(:, end), 'o--');
ylabel('Force');
title('Force at minimal Ca (pCa = 7.5)')
xticks(1:length(ATP))
xticklabels(leg);
xtickangle(45);

subplot(235);cla;
plot(F(:, 1), 'o--');
ylabel('Force');
title('Force at maximal Ca (pCa = 4.5)')
xticks(1:length(ATP))
xticklabels(leg);
xtickangle(45);

subplot(233);cla;
plot(hf(:,3), 'o--');
ylabel('Ca_{50} (\muM)');
title('Ca curve shifts with ATP concentration');
xticks(1:length(ATP))
xticklabels(leg);
xtickangle(45);
yyaxis right;
plot(-log10(hf(:,3))+6, 'o--');
ylabel('pCa_{50}');
% - log10(cRange) + 6 = pRange;
% log10(10^2) = 2

subplot(236);cla;
plot(hf(:,4), 'o--');
ylabel('Hill coefficient (-)');
xlabel('ATP');
xticks(1:length(ATP))
xticklabels(leg);
xtickangle(45);
title('pCa steepness changes with ATP');

%% test the hill fit
figure(6);clf;hold on;
% hf =  @(p, x)(1)./(1 + abs(p(1)/x)^p(2))
c = jet(size(hf, 1));
for i = 1:size(hf, 1)
    subplot(211);%hold on;
    semilogx(cRange, hill_fit(hf(i, :), cRange) - F(i,:)/max(F(i, :)), 'x--', 'Color', c(i, :));
    hold on;
    title('Hill Fit Error');
    subplot(212);hold on;
    semilogx(cRange, F(i,:)/max(F(i, :)), 'o:', 'Color', c(i, :));    
    hold on;
    semilogx(cRange, hill_fit(hf(i, :), cRange), 'x-', 'linewidth', 1, 'Color', c(i, :));
    title('Hill Fit (Abs)');
end

%% Test stepping up the SL
% cant use zero velocity
zv = 1e-3; % zero velocity
sv = 10; % step velocity
st = 1/sv; % step time
durVel = [1, zv;0.1/sv,sv;1, zv]; % duration x velocity
timesVel = [[0; cumsum(durVel(:, 1))], [durVel(:, 2);0]]
ML = 2.2;
dML = durVel(1:end, 1).*durVel(1:end, 2); % change in ML per each step
MLt = [0 ; cumsum(dML)]

plot(timesVel(:, 1), MLt*ML)

% timesVel(2:end,1) = durVel(1:end, 1) + durVel(2:end-1, 1) 