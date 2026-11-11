function params = updateRates(params)
% UPDATERATES  Scale all cross-bridge turnover rates by a single multiplier.
%
%   PARAMS = UPDATERATES(PARAMS) checks for the field PARAMS.XRATE. If
%   present, it multiplies all kinetic rate constants (kah, kamh, ka, kd,
%   k1, k_1, k2) by PARAMS.XRATE, then removes the field so the scaling
%   is not applied again on subsequent calls.
%
%   If PARAMS.XRATE is absent, the function returns PARAMS unchanged.
%
%   This is called automatically by getParams (Step 2) and should not
%   normally be called directly.
%
%   Input / Output:
%     PARAMS - parameter struct; modified in-place if xrate is set
%
%   See also: getParams

if ~isfield(params, 'xrate')
    return;
end

xr = params.xrate;
params = rmfield(params, 'xrate');


params.kah = params.kah*xr;
params.kamh = params.kamh*xr;

params.ka = params.ka*xr;
params.kd = params.kd*xr;

params.k1 = params.k1*xr;
params.k_1 = params.k_1*xr;

params.k2 = params.k2*xr;

% params.ksr0 = params.ka*params.xr;
% params.kmsr = params.ka*params.xr;
