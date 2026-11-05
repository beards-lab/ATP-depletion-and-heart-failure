function [D, fPCr] = adpScenario(A, pH, path)
%ADPSCENARIO  Free cytosolic ADP (mM) at ATP A (mM) and pH from the creatine-kinase equilibrium
%       PCr + ADP + H+ <-> Cr + ATP,     ADP = A * Cr / (K' * PCr),   K' = 166 * 10^-(pH - 7)
%   with the human pools of TWorld@lowAtp getCKEquilibriumADP (Weiss 2005: PCr = 10.1/5.7 x 8 mM,
%   total Cr = 27.6/5.7 x 8 mM; K_CK = 1.66e9 /M). H+ is a substrate of the ATP-forming direction, so
%   at a fixed creatine phosphorylation state ADP scales with 10^(pH - 7): acidosis LOWERS free ADP.
%   PCr is spent before ATP, so the PCr fraction fPCr(A) left at reduced ATP is the open variable;
%   PATH sets it, piecewise linear in A through (8, 1):
%     'control'  fPCr = 1 (fully CK-buffered; ADP falls in proportion to ATP)
%     'mild'     anchored to the lowest published ischaemic free ADP, guinea pig (Terkildsen 2007 /
%                Ferrero 2023): 182 uM at 4 mM, 109 uM at 2 mM
%     'severe'   anchored to rat (Bak & Ingwall 1994): 284 uM at 4 mM, 767 uM at 2 mM
%   The anchors are read as CK-derived at pH 7.05 (TWorld@lowAtp's reference pH) and inverted for fPCr.
Pn = 10.1/5.7*8;  Crt = 27.6/5.7*8;
K = @(p) 166*10.^(-(p - 7));
inv = @(a, d) (Crt/(1 + d*K(7.05)/a))/Pn;            % fPCr giving ADP d at ATP a, pH 7.05
switch path
    case 'control', f4 = 1;  f2 = 1;
    case 'mild',    f4 = inv(4, 0.182);  f2 = inv(2, 0.109);
    case 'severe',  f4 = inv(4, 0.284);  f2 = inv(2, 0.767);
    otherwise, error('adpScenario:path', 'unknown path %s', path);
end
fPCr = interp1([2 4 8], [f2 f4 1], A, 'linear', 'extrap');
PCr = fPCr*Pn;  Cr = Crt - PCr;
D = A .* Cr ./ (K(pH) .* PCr);
end
