% srxScaleRates  multiply the named params fields by x (the SRX rate kicks).
% See RunSrxRateKicks.m.
function q = srxScaleRates(p, names, x)
    q = p;
    for i = 1:numel(names); q.(names{i}) = q.(names{i}) * x; end
end
