function H = hostSignals(fitId, hostRoot, interp)
%HOSTSIGNALS  One periodic beat of a P21 TriSeg fit: Ls(t), dLs/dt and Ca(t) per wall.
%
%   H = HOSTSIGNALS(FITID) loads Multiscale-Hx-P21-Model/outputs/outputs_<FITID>.mat
%   (e.g. 'Nx51') and returns periodic piecewise-polynomial signals for ONE beat:
%     H.T              beat period (s) = 60/HR
%     H.<w>.Ls(t)      sarcomere length (um), w in {LV, SEP, RV}; Ls = Lsref*exp(eps)
%     H.<w>.v(t)       dLs/dt (um/s), the derivative of the same pchip
%     H.Ca(t)          cytosolic Ca (uM) the host drove its own cross-bridges with
%     H.raw            the stored traces (t, Ls, Lsc, sigma_pas, P_LV, ...)
%   Every signal is evaluated at mod(t, T), so it repeats forever.
%   The stored output spans 2 beats; the LAST full beat is used (the first may
%   still carry the host's own start-up transient).

if nargin < 2 || isempty(hostRoot), hostRoot = fullfile('..', '..', '..', 'Multiscale-Hx-P21-Model'); end
if nargin < 3, interp = 'pchip'; end       % 'pchip' (C1, kinks in dv/dt at every stored sample) or 'spline' (C2)
if strcmp(interp, 'spline'), fit = @spline; else, fit = @pchip; end
S = load(fullfile(hostRoot, 'outputs', ['outputs_' fitId '.mat']));
o = S.outputs;  t = o.time(:);
T = 60/S.data.HR;
Lsref = S.data.fixpars(13);

k  = t >= t(end) - T - 1e-9;                % last full beat
tb = t(k) - t(find(k, 1));
H.T = T;  H.HR = S.data.HR;  H.fitId = fitId;
walls = {'LV', 'SEP', 'RV'};
for i = 1:3
    w  = walls{i};
    Ls = Lsref*exp(o.strains.(['eps_' w])(:));
    Ls = Ls(k);
    pp = fit(tb, Ls);
    dpp = fnder_pp(pp);
    H.(w).Ls = @(tt) ppval(pp, mod(tt, T));
    H.(w).v  = @(tt) ppval(dpp, mod(tt, T));
    H.raw.(w).Ls  = Ls;
    H.raw.(w).Lsc = o.lengths.(['Lsc_' w])(:);  H.raw.(w).Lsc = H.raw.(w).Lsc(k)/1e4;  % stored in 1e-4 um
    H.raw.(w).sigma_pas = o.stresses.passive.(['sigma_pas_' w])(:);  H.raw.(w).sigma_pas = H.raw.(w).sigma_pas(k);
    H.raw.(w).sigma_XB  = o.stresses.active.(['sigma_XB_' w])(:);    H.raw.(w).sigma_XB  = H.raw.(w).sigma_XB(k);
    H.raw.(w).sigma     = o.stresses.total.(['sigma_' w])(:);        H.raw.(w).sigma     = H.raw.(w).sigma(k);
end
Ca  = o.XB.Ca_i(:);  Ca = Ca(k);
ppC = fit(tb, Ca);
H.Ca = @(tt) ppval(ppC, mod(tt, T));
H.raw.t = tb;  H.raw.Ca = Ca;
H.raw.P_LV = o.pressures.P_LV(:);  H.raw.P_LV = H.raw.P_LV(k);   % mmHg in the stored output
H.raw.V_LV = o.volumes.V_LV(:);    H.raw.V_LV = H.raw.V_LV(k);
H.data = S.data;
end

function d = fnder_pp(pp)
% derivative of a piecewise cubic (base MATLAB, no Curve Fitting Toolbox)
[br, c, ~, ord] = unmkpp(pp);
d = mkpp(br, c(:, 1:ord-1) .* (ord-1:-1:1));
end
