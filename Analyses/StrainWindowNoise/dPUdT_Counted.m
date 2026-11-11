function varargout = dPUdT_Counted(t, PU, params)
%DPUDT_COUNTED  Diagnostic wrapper: records the time of every RHS call.
%   Serial runs only. See ProbeStiffness.m.
global RHS_T RHS_N
RHS_N = RHS_N + 1;
if RHS_N > numel(RHS_T), RHS_T(2*numel(RHS_T) + 1e5) = 0; end
RHS_T(RHS_N) = t;
[varargout{1:nargout}] = dPUdT_CombinedTransitions(t, PU, params);
end
