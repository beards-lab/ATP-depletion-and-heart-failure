function H = hostP21Load(fitId, hostRoot)
%HOSTP21LOAD  One saved fit of the P21 host: parameters, the state at the start of
% its last beat, and that beat's traces (the reference the tests compare against).
%
% PROVENANCE: HOST DATA. Reads beards-lab/Multiscale-Hx-P21-Model,
% outputs/outputs_<fitId>.mat (written by the host's model_sol.m). Nothing here is a
% model equation; the index maps below follow model/model.m @ 80f86d3:
%   pars (stored as logs, model_sol.m l.20): 1-4 compliances C_SA C_SV C_PA C_PV,
%     5-6 resistances R_SA R_PA, 7-10 valve resistances m a t p, 11-13 Amref,
%     14-16 Vw (LV SEP RV), 17-18 k_TS k_TR                       (model.m l.44-72)
%   fixpars: 13 Lsref, 14 Lsc0, 18-19 k_passive LV/RV, 20 gamma  (model.m l.75-84)
%            75 k_on, 76-77 k_off LV/RV, 78 K_coop               (model.m l.110-113)
% Stored units: pressures mmHg, Lsc in 1e-4 um, everything else as in the model
% (kPa, cm, cm^3, um, uM, s).
%
%   H = hostP21Load('Nx51')            host repository at ../Multiscale-Hx-P21-Model
%   H = hostP21Load('Nx51', hostRoot)

if nargin < 2 || isempty(hostRoot)
    hostRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Multiscale-Hx-P21-Model');
end
S  = load(fullfile(hostRoot, 'outputs', ['outputs_' fitId '.mat']));
p  = exp(S.pars(:)');  fx = S.data.fixpars(:)';  o = S.outputs;
MMHG = 7.50062;                                   % mmHg per kPa

H.fitId = fitId;  H.hostRoot = hostRoot;  H.hostCommit = '80f86d3';
H.HR = S.data.HR;  H.T = 60/S.data.HR;
H.C  = struct('SA', p(1), 'SV', p(2), 'PA', p(3), 'PV', p(4));
H.R  = struct('SA', p(5), 'PA', p(6), 'm', p(7), 'a', p(8), 't', p(9), 'p', p(10));
H.Amref = p(11:13);  H.Vw = p(14:16);  H.k_TS = p(17);  H.k_TR = p(18);
H.Lsref = fx(13);  H.Lsc0 = fx(14);  H.k_passive = [fx(18), fx(18), fx(19)];  H.gamma = fx(20);
H.tf = struct('k_on', fx(75), 'k_off', [fx(76), fx(76), fx(77)], 'K_coop', fx(78));
Ca = o.XB.Ca_i(:);
H.Ca = struct('amplitude', max(Ca) - min(Ca), 'diastole', min(Ca));

% ---- the last full beat of the stored run
t = o.time(:);  k = find(t >= t(end) - H.T - 1e-9);  k1 = k(1);
walls = {'LV', 'SEP', 'RV'};
r.t = t(k) - t(k1);
H.tPhase = mod(t(k1), H.T);     % the stored beat does not start on a beat boundary: host time = r.t + tPhase
cv = @(x) reshape(x(k), [], 1);                    % stored traces are rows or columns
r.g = [cv(o.displacements.xm_LV), cv(o.displacements.xm_SEP), cv(o.displacements.xm_RV), cv(o.displacements.ym)];
for i = 1:3
    w = walls{i};
    r.eps(:, i)       = cv(o.strains.(['eps_' w]));
    r.Lsc(:, i)       = cv(o.lengths.(['Lsc_' w]))/1e4;
    r.sigma(:, i)     = cv(o.stresses.total.(['sigma_' w]));
    r.sigma_XB(:, i)  = cv(o.stresses.active.(['sigma_XB_' w]));
    r.sigma_pas(:, i) = cv(o.stresses.passive.(['sigma_pas_' w]));
end
r.Ls = H.Lsref*exp(r.eps);
P = @(n) cv(o.pressures.(n))/MMHG;                 % kPa
r.P = [P('P_LV'), P('P_SA'), P('P_SV'), P('P_RV'), P('P_PA'), P('P_PV')];
% volumes of the lumped compartments from their pressures (P = V/C, model.m l.324-327)
r.V = [cv(o.volumes.V_LV), cv(o.volumes.V_RV), r.P(:, 3)*H.C.SV, r.P(:, 6)*H.C.PV, r.P(:, 2)*H.C.SA, r.P(:, 5)*H.C.PA];
r.Ca = cv(Ca);
H.ref = r;
H.g0 = r.g(1, :)';  H.V0 = r.V(1, :)';  H.Ls0 = r.Ls(1, :);
H.data = struct('EDV_LV', S.data.EDV_LV, 'ESV_LV', S.data.ESV_LV, 'EDP_LV', S.data.EDP_LV, 'ESP_LV', S.data.ESP_LV, ...
                'EDV_RV', S.data.EDV_RV, 'ESV_RV', S.data.ESV_RV, 'EDP_RV', S.data.EDP_RV, 'ESP_RV', S.data.ESP_RV);
end
