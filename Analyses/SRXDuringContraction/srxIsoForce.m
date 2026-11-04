% srxIsoForce  isometric force and [SR SRD] at SL 2.2 after a Tw-second warm start.
% See RunSrxFamily.m.
function [F, ss] = srxIsoForce(p, ka, ic, Tw)
    % isometric force at SL 2.2 after a Tw-second warm start; ss = [SR SRD] at the end
    q = p; q.ka = ka; q.Velocity = [0; 0]; q.SL0 = 2.2;
    if ~isempty(ic); q.SRXT_0 = ic(1); q.SRXD_0 = ic(2); end
    q = getParams(q, q.g, true, false);
    [~, o] = evaluateModel(str2func(q.modelFcn), [-Tw 0], q);
    F = o.Force(end); ss = [o.SR(end), o.SRD(end)];
end
