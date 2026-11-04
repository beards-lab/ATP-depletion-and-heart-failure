function err = checkKernelParity(ref, params, idx)
% CHECKKERNELPARITY  xbKernels vs the real dxdt, on states from a PDE run.
%
%   For each sample idx of REF, rebuild PU, call dPUdT_CombinedTransitions and
%   compare its integrated rates/force with dS*sum(kernel.*p). This is the sync
%   test any moments backend must pass after a dxdt edit.
ss = params.ss; dS = params.dS;
err = zeros(numel(idx), 6);
for j = 1:numel(idx)
    i = idx(j);
    pr = params; pr.LXBpivot = ref.LXBPivot(i);
    pr.Vums = 0;                                   % rates are velocity-free here (hop is vel>0 only)
    PU = [ref.p1(i,:), ref.p2(i,:), ref.tail(i,:)]';
    [~, outp, rates] = dPUdT_CombinedTransitions(ref.t(i), PU, pr);
    SL = PU(2*ss+3); LSE = PU(2*ss+4);
    s = pr.s(:) - (-(SL - LSE) + pr.LXBpivot)/2;
    p1 = max(0, ref.p1(i,:)'); p2 = max(0, ref.p2(i,:)');
    K = xbKernels(s, pr, 0);
    mine = dS*[sum(K.F1.*p1 + K.F2.*p2), sum(K.k1D.*p1), sum(K.k12.*p1), sum(K.k21.*p2), sum(K.k2.*p2), sum(K.fSRX.*K.k2.*p2)];
    theirs = [outp(2), rates(4), rates(5), rates(6), rates(7), rates(17)];
    err(j,:) = abs(mine - theirs) ./ max(abs(theirs), 1e-9);
end
end
