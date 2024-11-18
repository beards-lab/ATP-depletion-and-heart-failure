function plotRates(out)

t = out.t;
%% Define forward fluxes between states
ST2UT = out.RSR2PT;
SD2UD = out.RSRD2PD;
UT2UD = out.RTD;  % Previously d2t, then s12s2
UD2A1 = out.RD1;  % Previously t2a1, then s22a1
A12A2 = out.R12;
A22UT = out.R2T;  % Previously a22d, then a22s1

% Define backward fluxes between states
UT2ST = out.RPT2SR;
UD2SD = out.RPD2SRD;
UD2UT = 0;  % Previously t2d, then s22s1
A12UD = out.R1D;  % Previously a12t, then a12s2
A22A1 = out.R21;

% Define flux between SD and ST
SD2ST = out.RSRD2SR;

% Assemble forward and backward fluxes into arrays
forward_fluxes = [UT2UD, UD2A1, A12A2, A22UT, ST2UT, SD2UD];
backward_fluxes = [UD2UT, A12UD, A22A1, UT2ST, UD2SD, SD2ST];

plot(t, UD2A1, t, A12A2, t, A22UT, t, UT2UD, LineWidth=2); hold on;
plot(t, UT2ST, 'g--', t, ST2UT, 'g:',t, UD2SD, 'c--', t, SD2UD,'c:', LineWidth=2);
legend('Attachment', 'Ratcheting', 'Unattachment', 'Hydrolysis', 'SR_T+', 'SR_T-', 'SR_D', 'SR_D-')
hold off;
