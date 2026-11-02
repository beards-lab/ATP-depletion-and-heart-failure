% srxApplyVariant  apply one SRX family variant (fam L/S/M, factor x) to a params struct.
% See RunSrxFamily.m.
function q = srxApplyVariant(p, v)
    q = p;
    switch v.fam
        case 'L'; q.ksr0 = q.ksr0 * v.x; q.ksrd = q.ksrd * v.x;
        case 'S'; q.sigma2 = q.sigma2 * v.x; q.sigma_srd1 = q.sigma_srd1 * v.x;
        case 'M'; q.ksr0 = q.ksr0 * v.x; q.ksrd = q.ksrd * v.x; q.kmsrd = q.kmsrd * v.x;
        % MA: every SRX rate incl. the hydrolysis inside SRX, so the pool's steady
        % composition is kept exactly. M alone shrinks it: SRD is still set by
        % influx/kmsrd (mu cancels), but SR -> SRD hydrolysis stays at 130 /s while the
        % influx that keeps SR stocked is slowed, so SR empties (0.023 -> ~0) and the
        % pool falls 0.030 -> 0.008.
        case 'MA'; q.ksr0 = q.ksr0 * v.x; q.ksrd = q.ksrd * v.x; q.kmsrd = q.kmsrd * v.x;
                   q.ksr2srd = q.ksr2srd * v.x; q.ksrd2sr = q.ksrd2sr * v.x;
    end
end
