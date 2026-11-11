function [f, outputs, rates] = dPUdT_CombinedTransitions_mex(t, PU, params)
% DPUDT_COMBINEDTRANSITIONS_MEX  Thin wrapper calling the compiled C MEX RHS.
%
%   Same signature as dPUdT_CombinedTransitions. Set as model function via:
%     params0.modelFcn = '@dPUdT_CombinedTransitions_mex';
%
%   Requires:
%     - dPUdT_core_mex.mexw64 compiled (run compileMex.m once)
%     - params.mex_vals set (UseCompiledMex=true + UseFastPPval=true in getParams)
%
%   params.Vums and params.LXBpivot are read each call (change between
%   ode15s segments; only these two scalars need live lookup).
%
%   See also: compileMex, packMexVals, dPUdT_CombinedTransitions_mex_flat

[f, outputs, rates] = dPUdT_core_mex(t, PU, params.mex_vals, params.Vums, params.LXBpivot);
end
