function params = updateRates(params)
% update turnover rates at once

if ~isfield(params, 'xrate')
    return;
end

xr = params.xrate;
params = rmfield(params, 'xrate');


params.kah = params.kah*xr;
params.ka = params.ka*xr;
% params.kd

params.k1 = params.k1*xr;
% params.k_1

params.k2 = params.k2*xr;

% params.ksr0 = params.ka*params.xr;
% params.kmsr = params.ka*params.xr;
