%% Runs an vmax to atp char
clear;
% original dans data
g = [0.7009    1.2243    1.0965    1.8390    0.4718 2.3357    0.3960    0.2372    0.1465    0.9817 0.8737    1.8333    0.2916    0.9513    1.0085]

% "optimized" set
g = [1.0695    3.4444    1.1394    0.8159    0.0215 1.0595    0.1061    1.3011    0.0920    2.6525 1.1631    2.2622    0.0997    1.0826    3.7106];


%% Set up environment
MgATP = [0.01 0.05 0.1 0.15 0.2 0.3 1 2 4 8];
MgADP = 0; 
Pi    = 0; 
fcn = @dPUdT;
vel = 0:0.5:8;
for v = 1:length(vel)
    for a = 1:length(MgATP)
%         F_active(v, a) = vel(v) + MgATP(a);
        F_active(v, a) = evaluateModel(fcn, vel(v), 1, MgATP(a),Pi,MgADP,g);
        P_active(v, a) = F_active(v, a)*v;
    end
end
disp('Done with the shit')
%%
% plot the shit out
set_atp = [1 7:length(MgATP)]

clf;

subplot(131);cla;hold on;title('Force to velocity');
for a = 1:length(set_atp)
    if MgATP(set_atp(a)) >= 2
        plot(F_active(:, set_atp(a)), vel, '-*', 'linewidth', 2);
    else
        plot(F_active(:, set_atp(a)), vel, '-*', 'linewidth',0.5);
    end
end
ylabel('velocity')
xlabel('force')
legend(strcat("ATP: ",string(num2cell(MgATP(set_atp))), " mM"))


subplot(132);cla;hold on;title('Power to velocity');
for a = 1:length(set_atp)
    if MgATP(set_atp(a)) >= 2
    plot(vel,P_active(:, set_atp(a)), '-*', 'linewidth', 2);
    else
            plot(vel,P_active(:, set_atp(a)), '-*', 'linewidth', 0.5);
    end
end
xlabel('velocity')
ylabel('Power')
legend(strcat("ATP: ",string(num2cell(MgATP(set_atp))), " mM"))

p_max = max(P_active(2:end, :));
subplot(133);cla;hold on;title('maximal P with ATP');
set(gca,'xscale','log');
plot(MgATP, p_max, '.-');
plot(MgATP(set_atp), p_max(set_atp), 'r.');

xlabel('MgATP (log)');
ylabel('Maximal power for velocity > 0');
