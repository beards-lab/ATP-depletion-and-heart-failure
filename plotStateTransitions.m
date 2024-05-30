% plot state transitions
myl = 500;

figure(22); clf;
nexttile; title('R1D'); hold on;
plot(s, params.kd*strainDep(params.alpha0, params.dr0), 'x-');
ylim([0, myl])
xlim([-.1 0.1])

nexttile; title('R12'); hold on;
plot(s, params.k1*exp(-params.alpha1*s), 'x-');
ylim([0, myl])
xlim([-.1 0.1])


nexttile; title('R21'); hold on;
plot(s, f1*params.k_1*strainDep(params.alpha_1, params.dr_1), 'x-'); % p2 to p1
ylim([0, myl])
xlim([-.1 0.1])

nexttile; title('R2T'); hold on;
plot(s, (params.k2 + kL + kR), 'x-');
ylim([0, myl])
xlim([-.1 0.1])
