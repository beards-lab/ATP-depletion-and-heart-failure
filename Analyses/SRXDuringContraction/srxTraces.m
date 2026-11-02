% srxTraces  isometric state (SL 2.2) plus two traced protocols for params q:
%   TRC  slack-restretch cycle 1 of the 03/27 8 mM protocol (release at t = 0)
%   TRK  the 03/27 ktr protocol (release at t = 0)
% Trace columns: [t, Force, SRX (SR+SRD), SL, attached]. Each run starts from the
% variant's own isometric SRX. See RunSrxRateKicks.m / RunSrxCompensation.m.
function [s, TRC, TRK] = srxTraces(q, root)
    s = srxIsoState(q);
    dP = load(fullfile(root, 'data', 'protocol_03_27_2026_8mM_slack.mat'));
    vt0 = dP.velocitytable; iRel = find(vt0(:, 2) < -1, 2); tR = vt0(iRel(1), 1);
    vtC = [-2, 0; vt0(iRel(1):iRel(2)-1, 1:2) - [tR 0]; vt0(iRel(2), 1) - tR - 0.05, 0];
    vtC(end-1:end, :) = [vtC(end-1, 1), 0; vtC(end, 1), 0];
    qc = q; qc.SRXT_0 = s.SR; qc.SRXD_0 = s.SRD; qc.SL0 = 2 * vt0(1, 4);
    qc.Velocity = vtC(:, 2); qc = getParams(qc, qc.g, true, false);
    [~, o] = evaluateModel(str2func(qc.modelFcn), vtC(:, 1), qc);
    TRC = iThin(o, (-0.05:5e-4:vtC(end, 1))');
    dK = load(fullfile(root, 'data', 'protocol_03_27_2026_velocitytable_ktr.mat'));
    vtK = dK.velocitytable(:, 1:2); vtK(:, 1) = vtK(:, 1) - vtK(2, 1);
    qk = q; qk.SRXT_0 = s.SR; qk.SRXD_0 = s.SRD; qk.SL0 = 2 * dK.velocitytable(1, 4);
    qk.Velocity = vtK(:, 2); qk = getParams(qk, qk.g, true, false);
    [~, o] = evaluateModel(str2func(qk.modelFcn), vtK(:, 1), qk);
    TRK = iThin(o, (-0.01:2.5e-4:vtK(end, 1))');
end

function T = iThin(o, tq)
    k = makeMonotonous(o.t(:));
    T = [tq, interp1(o.t(k), o.Force(k), tq), interp1(o.t(k), o.SR(k) + o.SRD(k), tq), ...
         interp1(o.t(k), o.SL(k), tq), interp1(o.t(k), o.p1_0(k) + o.p2_0(k), tq)];
end
