function [eff, dCa, diag] = calciumActivation(t, xCa, bound, params, Nov)
%CALCIUMACTIVATION  Ca source + thin-filament activation, as effectors on the XB cycle.
%
%   [EFF, DCA, DIAG] = CALCIUMACTIVATION(T, XCA, BOUND, PARAMS, NOV) is called once
%   per RHS evaluation by dPUdT_CombinedTransitions when PARAMS.ca.model ~= 0.
%   It never runs for CaModel = 'none', which is the calibrated saturating-Ca
%   model and stays byte-identical.
%
%   Inputs
%     T       time (s)
%     XCA     the Ca-layer states, PU(params.ca.idx) (empty for 'clamp')
%     BOUND   attached-head fraction p1_0 + p2_0 + p3_0. Proxy for the bound-site
%             fraction (one binding site per head); the overlap factor stays
%             in the RHS's own RD1 = ka_eff*PD*N_overlap*f_lattice.
%     NOV     geometric filament overlap (fraction of sites in the overlap zone),
%             before the RHS turns it into its overlap factor; read by 'myosim' only
%     PARAMS  read only. PARAMS.ca is resolved by getParams (buildCaLayer):
%             numeric ids, precomputed constants and transient coefficients.
%             Writing into PARAMS here would copy its ~300-field header on
%             every call (measured 10.6 us, ~7 % of an RHS call); reading is free.
%
%   Outputs
%     EFF     effector multipliers on the cross-bridge cycle. EFF(1) multiplies
%             ka_eff (attachment availability). Further effectors (SRX exit,
%             attachment kinetics) take the next indices when implemented.
%     DCA     d(XCA)/dt, same length as XCA
%     DIAG    [Ca (uM, NaN for 'clamp'), EFF(1)], appended to the RHS outputs
%
%   Layers (see Analyses/CalciumTwitch/conclusions.md):
%     source  params.ca.input  1 const (pCa) | 2 analytic (caTransient kind,
%                              fired every CaPeriod from CaFirstBeat) | 3 signal
%     model   params.ca.model  1 clamp    : rectangular availability pulse per
%                                            beat, no state (bandwidth probe)
%                              2 twoState : Tewari 2016 / Campbell 2018 /
%                                            Marzban 2020, one state N (non-
%                                            permissive fraction of sites):
%         Jon  = k_on*Ca*N*(1 + K_coop*(1-N))
%         Joff = k_off*Pu*(1 + K_coop*N),   Pu = max(0, 1 - N - bound)
%         gate = Pu/(1 - bound)
%       Sites held by bound heads cannot switch off (cross-bridge feedback);
%       at N = 0 the gate is exactly 1.
%                              3 myosim   : the same kinetics written as in MyoSim
%                                            (Campbell lab MATMyoSim, update_*.m), one
%                                            state N_on (permissive fraction of ALL
%                                            sites), switchable only inside the overlap:
%         Jon  = k_on*Ca*(Nov - N_on)*(1 + K_coop*N_on/Nov)
%         Joff = k_off*(N_on - bound)*(1 + K_coop*(Nov - N_on)/Nov)
%       MyoSim attaches in proportion to the free ON sites (N_on - bound). The RHS
%       already applies its overlap factor (Nov - bound)/unattached (UseOverlapFactor),
%       so the gate is (N_on - bound)/(Nov - bound): the product is MyoSim's site
%       term, and the gate is 1 when every overlap site is on. Without
%       UseOverlapFactor the RHS applies Nov, and the gate is (N_on - bound)/Nov.
%       At Nov = 1 the model reduces exactly to twoState with N = 1 - N_on, except
%       for the gate normalisation (1 - bound) vs (Nov - bound) = the same at Nov = 1.
%
%   See also: getParams, caTransient, dPUdT_CombinedTransitions

ca = params.ca;
switch ca.model
    case 1  % clamp
        tt = t - ca.t1;
        on = tt >= 0 && mod(tt, ca.T) < ca.clampD;
        eff = double(on);
        dCa = zeros(0, 1);
        diag = [NaN, eff];

    case 2  % twoState
        Ca = caSource(t, ca);
        N  = min(1, max(0, xCa(1)));
        Pu = max(0, 1 - N - bound);
        Jon  = ca.k_on*Ca*N*(1 + ca.K_coop*(1 - N));
        Joff = ca.k_off*Pu*(1 + ca.K_coop*N);
        dCa  = Joff - Jon;
        eff  = Pu/max(1 - bound, eps);
        diag = [Ca, eff];

    case 3  % myosim
        Ca  = caSource(t, ca);
        Nov = max(Nov, eps);
        Non = min(Nov, max(0, xCa(1)));
        Pu  = max(0, Non - bound);                      % free ON sites
        Jon  = ca.k_on*Ca*(Nov - Non)*(1 + ca.K_coop*Non/Nov);
        Joff = ca.k_off*Pu*(1 + ca.K_coop*(Nov - Non)/Nov);
        dCa  = Jon - Joff;
        if ca.useOverlapFactor
            eff = Pu/max(Nov - bound, eps);
        else
            eff = Pu/Nov;
        end
        diag = [Ca, eff];
end
end

function Ca = caSource(t, ca)
% Cytosolic [Ca2+] (uM) at time t.
switch ca.input
    case 1  % const
        Ca = ca.Ca_const;
    case 2  % analytic, one fixed frequency
        tt = t - ca.t1;
        if tt < 0
            Ca = ca.Cadia;
        else
            Ca = caTransient(ca.kind, tt, ca.opt);
        end
    case 3  % measured
        Ca = sigLookup(ca.sig, t);
end
end
