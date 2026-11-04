% srxApplyTheta  apply log multipliers th to the knob list KN ({name, {fields}} rows).
% See RunSrxCompensation.m.
function q = srxApplyTheta(p, KN, th)
    q = p;
    for j = 1:size(KN, 1)
        if th(j) ~= 0; q = srxScaleRates(q, KN{j, 2}, exp(th(j))); end
    end
end
