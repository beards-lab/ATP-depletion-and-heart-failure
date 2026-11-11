function K = xbKernels(s, params, vel)
% XBKERNELS  Per-head strain kernels of the active cross-bridge mechanism set.
%
%   K = XBKERNELS(S, PARAMS, VEL) evaluates, at strains S (column), every
%   strain-dependent quantity that dPUdT_CombinedTransitions integrates against
%   p1 or p2. Each field is a PER-HEAD rate (1/s) or per-head force (kPa per unit
%   probability), so a flux is dS*sum(K.x .* p) on the grid, or sum(w .* K.x) on
%   quadrature nodes.
%
%   This is the ONE thing a moments backend and the grid backend must share to
%   stay in sync. It replicates exactly the branches active in fastgate_best:
%   UsePieceWiseStrainDep + UseFastPPval, UseNegativeKstiff, SRXFromR2HighStrain,
%   UseA2Reattaching (hop). checkKernelParity() verifies it against the dxdt.
%
%   Fields:
%     F1, F2   force kernels on p1, p2
%     k1D      p1 -> PD (detachment of the weak state)
%     k12      p1 -> p2 (power stroke)
%     k21      p2 -> p1 (reverse stroke)
%     k2       p2 -> detached (PT; part routed to SRX by fSRX)
%     fSRX     fraction of k2 flux routed to SRX (dimensionless)
%     kHop     p2 high-strain hop back to strain dA2re (only when VEL > 0)

if ~isfield(params, 'kernelsChecked'), assertSupported(params); end   % hot path skips it
s = s(:);
MAX_RATE = 1e4;

% force kernels (UseNegativeKstiff branch)
s1 = s + params.dr1;
e2 = sign(s + params.dr) .* abs(s + params.dr).^params.estiff;
K.F1 = params.kstiff1 * s1 .* (s1 >= 0) + params.kstiff1_n * s1 .* (s1 < 0);
K.F2 = params.kstiff2 * e2 .* (s >= -params.dr) + params.kstiff2_n * e2 .* (s < -params.dr);

% transition kernels (UsePieceWiseStrainDep, fast path)
c  = params.cached;
cl = params.ClampStrainDepExtrapolation;
K.k1D = params.kd  * pwEval(c.pwsdR1D_brk, c.pwsdR1D_cof, s, cl);
K.k12 = params.k1  * pwEval(c.pwsd_brk,    c.pwsd_cof,    s, cl);
K.k21 = params.k_1 * pwEval(c.pwsdR21_brk, c.pwsdR21_cof, s, cl);
K.k2  = params.k2  * pwEval(c.pwsd2_brk,   c.pwsd2_cof,   s + params.dr2 - params.dr, cl);

% detachment routed to SRX (high-strain rupture branch)
f = params.SRXFromR2 + params.SRXFromR2HighStrain * ...
    min(1, max(0, (s - params.sSRXrip) / max(params.dSRXrip, 1e-9)));
K.fSRX = min(1, f);

% A2 hop (lengthening only)
if params.UseA2Reattaching && vel > 0 && params.kA2hop > 0
    K.kHop = min(MAX_RATE, params.kA2hop * max(0, s - params.sA2hop));
else
    K.kHop = zeros(size(s));
end
end

function assertSupported(p)
% Fail loudly if the parameter set uses a branch this kernel set does not model.
off = {'UseYieldGate','UseCatchBond','UseR2Ceiling','UseAtpDetachR2','UseD0State', ...
       'UseA2AttachmentShift','UseA2MechanicalRecocking','UseA2Popping', ...
       'UseTargetZoneSaturation','UseVernierVelocity','UseA1AttachmentKernel', ...
       'UseDirectSRXTransition','UseUniformTransitionFunc','UseVelGaussAttachment', ...
       'UseStretchActivation','LegacyStrainFlipping'};
for i = 1:numel(off)
    if isfield(p, off{i}) && any(p.(off{i}))
        error('xbKernels:unsupported', '%s is on but not modelled by the moments backend', off{i});
    end
end
on = {'UsePieceWiseStrainDep','UseFastPPval','UseNegativeKstiff'};
for i = 1:numel(on)
    if ~p.(on{i}), error('xbKernels:unsupported', '%s must be on', on{i}); end
end
if p.k2rip ~= 0 || p.UseStrictDetachmentAt > 0
    error('xbKernels:unsupported', 'XB_Ripped / strict detachment not modelled');
end
end

function y = pwEval(b, c, x, clampToKnots)
% verbatim copy of the local function in dPUdT_CombinedTransitions
if clampToKnots
    x = min(max(x, b(1)), b(end));
end
si = max(1, min(size(c,1), sum(x >= b(1:end-1)', 2)));
dx = x - b(si);
y = ((c(si,1).*dx + c(si,2)).*dx + c(si,3)).*dx + c(si,4);
end
