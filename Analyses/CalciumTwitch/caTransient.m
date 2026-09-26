function Ca = caTransient(kind, t, opt)
%CATRANSIENT  Candidate prescribed cytosolic [Ca2+](t) driving signals, in uM.
%
%   CA = CATRANSIENT(KIND, T) evaluates transient KIND at times T (s). Every
%   transient is periodic with period OPT.T and a beat starts at each multiple
%   of OPT.T. CA = CATRANSIENT(KIND, T, OPT) overrides the defaults below.
%
%   KIND
%     'rice2008'     bi-exponential; rat, 24 C (Rice et al. 2008, CellML defaults):
%                    rise tau1 = 20 ms, decay tau2 = 110 ms, 0.09 -> 1.45 uM, onset 5 ms.
%     'marzban2020'  log-normal in cycle phase phi = mod(t,T)/T; rat, 2-10 Hz table,
%                    pchip in rate (Marzban 2020 PhysiomeS000007; identical table in
%                    Multiscale-Hx-P21-Model/model/compare_ca_curve.m).
%     'tewari_tanh'  tanh - exponential; human, Piacentino 1 Hz shape time-scaled by
%                    the Kang restitution Td(T) (cross-bridge-model XB.mo, block Ca2).
%     'fibersim2c'   FiberSim two-compartment SR/myofilament pulse, closed form
%                    (defaults of demo_files/electrical_stimulation/twitch).
%     'pCa'          constant 10^(6-OPT.pCa) uM.
%
%   OPT fields: T (period, s; default 1), pCa, and per-kind overrides.

if nargin < 3, opt = struct(); end
T  = getOpt(opt, 'T', 1);
tt = mod(t, T);                                   % time since beat onset

switch kind
    case 'rice2008'
        tau1 = getOpt(opt, 'tau1', 0.020);  tau2 = getOpt(opt, 'tau2', 0.110);
        Cad  = getOpt(opt, 'Cadia', 0.09);  Cam  = getOpt(opt, 'Camax', 1.45);
        t0   = getOpt(opt, 't0', 0.005);
        r    = tau1/tau2;
        beta = r^(tau1/(tau2-tau1)) - r^(tau2/(tau2-tau1));   % normalises the peak to Camax
        x    = max(0, tt - t0);
        Ca   = Cad + (Cam - Cad)/beta*(exp(-x/tau2) - exp(-x/tau1));

    case 'marzban2020'
        tab = [2      3      4      5      6      7      8      9      10;      % Hz
               0.0838 0.1306 0.1802 0.2557 0.3099 0.3613 0.4080 0.4539 0.4602; % a
               0.7513 0.8756 1.0274 1.4988 1.6107 1.6741 1.7902 2.1398 1.9832; % b
               2.2370 2.0486 1.9480 1.8520 1.6932 1.6773 1.5988 1.4952 1.4524; % c
               0.1865 0.1815 0.1709 0.1693 0.1610 0.1661 0.1425 0.1458 0.1222]; % Ca0
        c4  = pchip(tab(1,:), tab(2:5,:), 1/T);
        phi = max(tt/T, 1e-9);
        Ca  = c4(1)./phi.*exp(-c4(2)*(log(phi) + c4(3)).^2) + c4(4);

    case 'tewari_tanh'
        Td1 = 0.36/(0.39^2.2 + 1);                    % Kang 80 % relaxation time at 1 Hz
        r   = (0.36*T^2.2/(0.39^2.2 + T^2.2))/Td1;    % duration relative to 1 Hz
        a   = 0.56 + 0.1655*(1/T - 1);
        Ca  = max(0, a*(0.5*(1 - tanh(7.66/r*(tt - 0.289*r))) - exp(-162.2/r*tt))) ...
              + getOpt(opt, 'Ca0', 0.040);

    case 'fibersim2c'
        C  = getOpt(opt, 'Ca_content', 1e-3);  ks = getOpt(opt, 'k_serca', 20);
        kl = getOpt(opt, 'k_leak', 6e-4);      ka = getOpt(opt, 'k_act', 8.2e-2);
        D  = getOpt(opt, 'stim', 0.010);
        myo = @(x0, ko, s) C*ko/(ks+ko) + (x0 - C*ko/(ks+ko)).*exp(-(ks+ko)*s);  % Ca_myofil, M
        m0  = C*kl/(ks+kl);                                                     % resting
        mD  = myo(m0, kl+ka, D);                                                % end of stimulus
        Ca  = 1e6*((tt <= D).*myo(m0, kl+ka, tt) + (tt > D).*myo(mD, kl, tt - D));

    case 'pCa'
        Ca = 10^(6 - opt.pCa)*ones(size(t));

    otherwise
        error('caTransient:kind', 'Unknown transient ''%s''.', kind);
end
end

function v = getOpt(opt, name, default)
if isfield(opt, name), v = opt.(name); else, v = default; end
end
