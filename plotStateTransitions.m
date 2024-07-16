% plot state transitions
myl = 500;

figure(22); clf;
nexttile; title('R1D'); hold on;
plot(s, R1D, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('R12'); hold on;
plot(s, R12, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('R21'); hold on;
plot(s, R21, 'x-'); % p2 to p1
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

nexttile; title('R2T'); hold on;
plot(s, R2T, 'x-');
ylim([0, myl])
xlim([s(1) s(end)])
xlabel('s (\mum)');ylabel('Transition rate (1/s)');

% s - Distance from attachment point