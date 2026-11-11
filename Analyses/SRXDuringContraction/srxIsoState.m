% srxIsoState  isometric steady state at SL 2.2 after a 3 s warm start: force, SRX pool,
% attached fraction and the net ATPase R2T (strong-detachment flux = net turnover; the
% objective's XTOR is the gross DRX flux RTD). See RunSrxRateKicks.m.
function [s, o] = srxIsoState(p)
    q = p; q.Velocity = [0; 0]; q.SL0 = 2.2; q = getParams(q, q.g, true, false);
    [~, o] = evaluateModel(str2func(q.modelFcn), [-3 0], q);
    s = struct('F', o.Force(end), 'SR', o.SR(end), 'SRD', o.SRD(end), 'SRX', o.SR(end) + o.SRD(end), ...
               'att', o.p1_0(end) + o.p2_0(end), 'R2T', o.R2T(end), 'RTD', o.RTD(end));
end
