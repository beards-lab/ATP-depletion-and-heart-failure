% RunCharacterize.m — the two model constants the protocol port depends on.
%
% (1) Vmax, the unloaded shortening velocity (ML/s). van der Zee et al. 2026
%     (PLOS Comput Biol 10.1371/journal.pcbi.1014748) state their stretch rate
%     of 0.45 L0/s is "~35% of the expected unloaded shortening velocity" of
%     rat soleus (~1.3 L0/s). v/Vmax = 0.35 is therefore the paper's OWN
%     normalisation, and is what "kinetically matched" means here.
%
% (2) ktr, the force-redevelopment rate (1/s). Their fibres regenerate tension
%     at < 2 1/s; cardiac is ~an order of magnitude faster. This sets how many
%     cross-bridge turnovers elapse during the ramp, and in particular whether
%     the 10 ms SRS regression window is still "instantaneous elasticity".
%
% Writes results/characterize.mat

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));

params0 = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
params0.PlotEachSeparately = false;
modelFcn = str2func(params0.modelFcn);

%% [1] Force-velocity curve, fine near the zero-force end ------------------
vs = [0 -0.25 -0.5 -0.75 -1 -1.5 -2 -3 -4 -5 -5.5 -6 -6.5 -7 -7.5 -8 -9];
F_tot = nan(size(vs)); F_act = nan(size(vs)); F_pas = nan(size(vs));

for j = 1:numel(vs)
    p = params0;
    p.SL0 = 2.0;
    dsl = 80 / p.kSE;
    p.Slim_r = 2.0 + p.A1AttachmentWidth + p.dS;
    p.Slim_l = 2.0 - 0.25 - dsl - p.A1AttachmentWidth - p.dS;
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = vs(j);
    tspan = [0 min(0.15, 0.2 / max(abs(vs(j)), 0.5))];
    [F, out] = evaluateModel(modelFcn, tspan, p);
    F_tot(j) = F; F_act(j) = out.FXB(end); F_pas(j) = out.FXBPassive(end);
    fprintf('v=%6.2f  F_tot=%8.3f  F_act=%8.3f  F_pas=%7.3f kPa\n', ...
            vs(j), F_tot(j), F_act(j), F_pas(j));
end

F0 = F_tot(1);

% Unloaded shortening velocity = |v| at which TOTAL force reaches zero (that is
% what "unloaded" means: zero external load). Active force alone only asymptotes
% to a small positive floor, so a zero-crossing of F_act does not exist.
iz = find(F_tot <= 0, 1);
if isempty(iz)
    Vmax = NaN;
else
    Vmax = abs(interp1(F_tot(iz-1:iz), vs(iz-1:iz), 0, 'linear'));
end

% Robustness: also report the (definition-insensitive) velocity at 5% F0.
i5 = find(F_tot <= 0.05*F0, 1);
V05 = abs(interp1(F_tot(i5-1:i5), vs(i5-1:i5), 0.05*F0, 'linear'));

fprintf('\nF0 (isometric, SL 2.0) = %.2f kPa\n', F0);
fprintf('Vmax (F_tot -> 0)      = %.3f ML/s  (%.2f um/s SL)\n', Vmax, Vmax*params0.ML);
fprintf('V at 5%% F0             = %.3f ML/s\n', V05);

%% [2] Force redevelopment rate -------------------------------------------
ktr = ktrProbe(params0);
fprintf('ktr = %.2f 1/s   (tau = %.1f ms)\n', ktr, 1000/ktr);

%% [3] Derive the two protocol cases --------------------------------------
% Both start at SL 2.0 um (L0 = 1.0 ML, since params0.ML = 2) and stretch by
% 4% of L0, so the AMPLITUDE is identical; only the ramp rate differs.
fracAmp = 0.04;
dL      = fracAmp * 1.0;                  % ML

v_verb  = fracAmp / 0.100;                % 4% L0 in 100 ms  -> ML/s
v_kin   = 0.35 * Vmax;                    % same v/Vmax as the paper
dur_kin = dL / v_kin;

% Secondary, for the discussion only: matching TURNOVERS during the ramp.
% Their ramp spans 0.100 s * ktr_sol turnovers; ours must span the same.
ktr_soleus = 2.0;                         % upper bound from the paper (<2 1/s)
dur_tor    = 0.100 * ktr_soleus / ktr;
v_tor      = dL / dur_tor;

fprintf('\n--- protocol cases (both: SL0 = 2.0 um, +4%% L0 = +0.08 um SL = 40 nm/HS) ---\n');
fprintf('A verbatim          : v = %.3f ML/s   ramp = %6.1f ms   v/Vmax = %.3f\n', ...
        v_verb, 1000*dL/v_verb, v_verb/Vmax);
fprintf('B kinetically match : v = %.3f ML/s   ramp = %6.1f ms   v/Vmax = %.3f\n', ...
        v_kin, 1000*dur_kin, 0.35);
fprintf('  (turnover-matched : v = %.3f ML/s   ramp = %6.1f ms   v/Vmax = %.3f)\n', ...
        v_tor, 1000*dur_tor, v_tor/Vmax);
fprintf('\nSRS window (first 10 ms of the ramp) spans:\n');
fprintf('  A: %.4f ML = %.1f nm/HS, and %.2f ktr time constants\n', ...
        v_verb*0.010, v_verb*0.010*params0.ML*1000/2, 0.010*ktr);
fprintf('  B: %.4f ML = %.1f nm/HS  (ramp is only %.1f ms, window truncates)\n', ...
        v_kin*0.010, v_kin*0.010*params0.ML*1000/2, 1000*dur_kin);

save(fullfile(here, 'results', 'characterize.mat'), ...
     'vs', 'F_tot', 'F_act', 'F_pas', 'F0', 'Vmax', 'V05', 'ktr', ...
     'v_verb', 'v_kin', 'dur_kin', 'v_tor', 'dur_tor', 'ktr_soleus', 'fracAmp', 'dL');
