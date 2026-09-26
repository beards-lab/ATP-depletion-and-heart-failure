function varargout = dPUdT_CountedTyped(t, PU, params)
%DPUDT_COUNTEDTYPED Diagnostic wrapper: counts RHS calls by caller type.
%   ode15s asks for one output; evaluateModel's storeOutputs asks for three.
%   Serial runs only (globals). See RunCostCount.m.
global RHS_T RHS_K RHS_N %#ok<GVMIS>
RHS_N = RHS_N + 1;
if RHS_N > numel(RHS_T), RHS_T(2*numel(RHS_T) + 1e5) = 0; RHS_K(numel(RHS_T)) = 0; end
RHS_T(RHS_N) = t; RHS_K(RHS_N) = nargout > 1;
[varargout{1:nargout}] = dPUdT_CombinedTransitions(t, PU, params);
end
