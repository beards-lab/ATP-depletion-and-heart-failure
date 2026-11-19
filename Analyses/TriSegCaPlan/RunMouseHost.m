%% RunMouseHost  The lab's mouse TriSeg host (beards-lab/Rat_Triseg, private) as an alternative to the P21 host (#82, Q7).
%
% Part 1 runs the host's own two fitted cases (the literature-target mouse and one
% echo-fitted mouse) with the host's scripts and keeps the sarcomere kinematics.
% Part 2 drives one wall of the cross-bridge model open loop on those beats, as
% RunHostProbe does on the P21 beat. The mouse host has no Ca transient (its
% activation is a prescribed function of time), so Ca(t) is the P21 host's analytic
% transient compressed to the mouse beat period: an assumption.
%
% Needs the host repository next to this one (../Rat_Triseg) or MOUSEROOT, and the
% P21 host for the Ca transient. Results -> results/mouseHost_*.mat, mouseHostProbe.mat.
%   matlab -batch "run('Analyses/TriSegCaPlan/RunMouseHost.m')"      (about 5 min)

here = fileparts(mfilename('fullpath'));  root = fileparts(fileparts(here));
if ~exist('MOUSEROOT', 'var'), MOUSEROOT = fullfile(fileparts(root), 'Rat_Triseg'); end
if ~exist('HOSTROOT', 'var'),  HOSTROOT  = fullfile(fileparts(root), 'Multiscale-Hx-P21-Model'); end
SP = fullfile(here, 'results');  if ~exist(SP, 'dir'), mkdir(SP); end
D  = MOUSEROOT;

%% Part 1: the host's own beats
cd(D); addpath(D);
cases = {'mouse', 'm22146'};
for ic = 1:2
    clearvars -except D SP cases ic here root MOUSEROOT HOSTROOT
    cs = cases{ic};
    try
        if strcmp(cs, 'mouse')
            [targets, inputs, mods] = targetVals_mouse();
            modifiers = readmatrix('modifiers_mouse.csv', 'range', '2:2');
        else
            [targets, inputs, mods] = targetVals_22146();
            modifiers = readModifiers('modifiers_22146.csv', mods);
            kill_flag = [1 1 1 1 0 0 0];
        end
        MRI_flag = 1;
        [INIparams, INIinit] = estiminiParams(targets, inputs);
        [params, init] = optParams(INIparams, INIinit, mods, modifiers);
        print_sim = true;
        tic; runSimOnGL; el = toc;
        S = struct('t', t(:), 'T', params.T, 'HR', params.HR, 'Y', Y(:), 'P_LV', P_LV(:), 'V_LV', V_LV(:), 'el', el);
        for w = {'LV', 'SEP', 'RV'}
            S.(['Ls_' w{1}])  = 2*exp(eval(['eps_' w{1}]));  S.(['Lsc_' w{1}]) = eval(['Lsc_' w{1}]);
            S.(['spas_' w{1}]) = eval(['sigma_pas_' w{1}]);  S.(['sact_' w{1}]) = eval(['sigma_act_' w{1}]);  S.(['s_' w{1}]) = eval(['sigma_' w{1}]);
        end
        S.params = params;
        save(fullfile(SP, ['mouseHost_' cs '.mat']), 'S');
        k = S.t >= S.t(end) - S.T - 1e-9;  tt = S.t(k) - S.t(find(k, 1));
        fprintf('CASE %s: HR %.0f bpm, T %.1f ms, %d samples in the last beat, sim %.0f s\n', cs, S.HR, 1e3*S.T, sum(k), el);
        fprintf('   LV: EDV %.1f uL, ESV %.1f uL, P %.1f..%.1f mmHg; tau_TS %.3f tau_TR %.3f; k_act_LV %.4g, k_pas_LV %.4g\n', 1e3*max(S.V_LV(k)), 1e3*min(S.V_LV(k)), min(S.P_LV(k)), max(S.P_LV(k)), params.tau_TS, params.tau_TR, params.k_act_LV, params.k_pas_LV);
        for w = {'LV', 'SEP', 'RV'}
            Ls = S.(['Ls_' w{1}])(k);  v = gradient(Ls, tt);
            fprintf('   %-3s: Ls %.3f..%.3f um, peak shortening %.1f um/s (%.2f lengths/s) at %.0f %% of the beat, peak lengthening %.1f um/s; stress act max %.1f, pas max %.2f, total max %.1f (host units)\n', ...
                w{1}, min(Ls), max(Ls), -min(v), -min(v)/max(Ls), 100*tt(v == min(v))/S.T, max(v), max(S.(['sact_' w{1}])(k)), max(S.(['spas_' w{1}])(k)), max(S.(['s_' w{1}])(k)));
        end
    catch ME
        fprintf('CASE %s FAILED: %s\n', cs, ME.message);
        if ~isempty(ME.stack), fprintf('   at %s line %d\n', ME.stack(1).name, ME.stack(1).line); end
    end
end

%% Part 2: the cross-bridge wall on those beats
rmpath(D);
cd(root); addpath(genpath(fullfile(root, 'Model'))); addpath(genpath(fullfile(root, 'Auxiliary'))); addpath(fullfile(root, 'params'));
addpath(fullfile(root, 'Analyses', 'TriSegCaPlan'));
Hp = hostSignals('Nx51', HOSTROOT);
res = struct();
for cs = {'mouse', 'm22146'}
    L = load(fullfile(SP, ['mouseHost_' cs{1} '.mat'])); S = L.S; T = S.T;
    k = S.t >= S.t(end) - T - 1e-9; tb = S.t(k) - S.t(find(k, 1)); [tb, iu] = unique(tb);
    tu = (0:0.25e-3:T)';
    H = struct('T', T, 'HR', S.HR);
    for w = {'LV', 'SEP', 'RV'}
        Ls = S.(['Ls_' w{1}])(k); Ls = Ls(iu); Lu = interp1(tb, Ls, min(tu, tb(end)), 'pchip');
        pp = spline(tu, Lu); [br, c, ~, ord] = unmkpp(pp); dpp = mkpp(br, c(:, 1:ord-1) .* (ord-1:-1:1));
        H.(w{1}).Ls = @(tt) ppval(pp, mod(tt, T)); H.(w{1}).v = @(tt) ppval(dpp, mod(tt, T));
        sa = S.(['sact_' w{1}])(k); sp = S.(['spas_' w{1}])(k);
        H.raw.(w{1}).sigma_XB = interp1(tb, sa(iu), min(tu, tb(end)))/7.50062; H.raw.(w{1}).sigma_pas = interp1(tb, sp(iu), min(tu, tb(end)))/7.50062;
    end
    H.Ca = @(tt) Hp.Ca(mod(tt, T)*Hp.T/T);  H.raw.t = tu;
    arms = {'kon3', struct('ss', 60, 'velMode', 'imposed', 'konMult', 3); 'kon10', struct('ss', 60, 'velMode', 'imposed', 'konMult', 10); ...
            'kon3_x5', struct('ss', 60, 'velMode', 'imposed', 'konMult', 3, 'xrate', 5)};
    for a = 1:size(arms, 1)
        nm = [cs{1} '_' arms{a, 1}];
        try
            P = setupWall(H, 'LV', arms{a, 2}); R = runHostBeats(P, 4, struct('tPre', 0.5));
            kk = R.t >= 3*T - 1e-9; tk = R.t(kk) - 3*T; dtk = [diff(tk); 0];
            F = R.F(kk); Fa = R.Factive(kk); Fp = R.Fpassive(kk);
            [pk, ip] = max(Fa); hx = H.raw.LV.sigma_XB;
            ej = tk >= 0.25*T & tk <= 0.45*T;          % the ejection part of the host's beat (by its Ls)
            fprintf('%-16s: %.1f s/beat, slides %d | active peak %.1f kPa at %.0f ms, min %.2f | total %.2f..%.1f | passive at ED %.2f (host %.2f) | slack (<0.05 kPa) %.1f ms | gate max %.2f | host active peak %.1f at %.0f ms\n', ...
                nm, R.wallPerBeat, R.nShift, pk, 1e3*tk(ip), min(Fa), min(F), max(F), Fp(1), H.raw.LV.sigma_pas(1), 1e3*sum(dtk(F < 0.05)), max(R.gate(kk)), max(hx), 1e3*tu(hx == max(hx)));
            R = rmfield(R, 'Y'); res.(nm) = R; res.(nm).tk = tk; res.(nm).hostT = T;
        catch ME
            fprintf('%-16s FAILED: %s\n', nm, ME.message);
        end
    end
    res.(['H_' cs{1}]) = H.raw; res.(['T_' cs{1}]) = T;
    tt = linspace(0, T, 600)'; res.(['Ls_' cs{1}]) = [tt, H.LV.Ls(tt), H.LV.v(tt), H.Ca(tt)];
end
save(fullfile(SP, 'mouseHostProbe.mat'), 'res');
