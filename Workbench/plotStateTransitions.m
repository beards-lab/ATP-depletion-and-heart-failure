% plot state transitions
myl = 500;

% figure(22); clf;
% ATTACHMENT - DETACHMENT
nexttile; title('Hydrolyses and ATTACHMENT'); cla;hold on;
plot([s(1), 0, s(end)], [RD1, RD1, RD1], '-', LineWidth=1.5);
plot(s, R1D,'x-.', LineWidth=1.5);
plot([s(1), 0, s(end)], [RTD,RTD, RTD], 'x-', LineWidth=1, MarkerSize=10);
plot([s(1), 0, s(end)], [RDT,RDT, RDT], 'x-.', LineWidth=1, MarkerSize=10);

plot([s(1), 0, s(end)], [RSR2SRD,RSR2SRD, RSR2SRD], 's-', LineWidth=0.5, MarkerSize=14);
plot([s(1), 0, s(end)], [RSRD2SR,RSRD2SR, RSRD2SR], 's-.', LineWidth=0.5, MarkerSize=14);

legend('R1D - Attachment', 'RD2 - detachment', 'RTD - hydrolysis', 'RDT - rev hydrolysis', ...
    'SRX RTD - hydrolysis', 'SRX - rev hydrolysis');
% ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('Attached states'); hold on;
plot(s, R12, 'x-', s, R21, 'x-.');
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
plot(s, R2, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('T2SR and SR2T'); hold on;
plot(F_SR, RPT2SR, 'x-', F_SR, RSR2PT, '+-.', LineWidth=1.5, MarkerSize=10);
plot(F_SR, RPD2SRD, 'o-', F_SR, RSRD2PD, 's-.', LineWidth=1, MarkerSize=14);

legend('T to SR', 'SR to T','D to SRD', 'SRD to D');
% ylim([0, myl])
xlim([F_SR(1) F_SR(end)])
xlabel('Force (kPa)');ylabel('Transition rate (1/s)');
return

%% test Hill SRD transition
% F_SR, RPD2SRD,
figure(24);clf;hold on;
plot(F_SR, RPD2SRD, 'o-', F_SR, RSRD2PD, 'o-');

D2SRD_hill = @(K, n, x) K./(1+exp(n*(x)));

plot(F_SR, D2SRD_hill(200, 0.4, F_SR), 'x-');

% s - Distance from attachment point