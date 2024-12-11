% plot state transitions
myl = 500;

% figure(22); clf;
nexttile; title('R1D'); hold on;
plot(s, R1D, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('R12 and R21'); hold on;
plot(s, R12, 'x-', s, R21, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

% nexttile; title('R21'); hold on;
% plot(s, R21, 'x-'); % p2 to p1
% ylim([0, myl])
% xlim([s(1) s(end)])
% xlabel('s (\mum)');ylabel('Transition rate (1/s)');
legend('R12', 'R21')

nexttile; title('R2T'); hold on;
plot(s, R2T, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('T2SR and SR2T'); hold on;
plot(F_SR, RPT2SR, 'x-', F_SR, RSR2PT, '+-');
plot(F_SR, RPD2SRD, 'o-', F_SR, RSRD2PD, 'o-');

legend('T to SR', 'SR to T','D to SRD', 'SRD to D');
% ylim([0, myl])
xlim([F_SR(1) F_SR(end)])
xlabel('Force (kPa)');ylabel('Transition rate (1/s)');

% s - Distance from attachment point