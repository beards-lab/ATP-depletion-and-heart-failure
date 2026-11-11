function stop = myoutput(gr, optimvalues, state)
% MYOUTPUT  Output function for fminsearch/fmincon: plots fit every 10 iterations.
%
%   STOP = MYOUTPUT(GR, OPTIMVALUES, STATE) is called by MATLAB optimizers
%   after each iteration. Every 10 iterations it re-runs evaluateProblem
%   with the current modifiers GR and displays the results.
%
%   Inputs:
%     GR          - current modifier vector (passed directly to evaluateProblem)
%     OPTIMVALUES - struct provided by the optimizer (fields: iteration, fval, etc.)
%     STATE       - 'init', 'iter', or 'done'
%
%   Output:
%     STOP - always false (does not interrupt the optimization)
%
%   See also: evaluateProblem, fminsearch

    stop = false;
    if ~isequal(state,'iter') || mod(optimvalues.iteration, 10) > 0
        return;
    end
%     g_all = [gr(1) 3 gr(2) 0.8 gr(3:end)];
    evaluateProblem(@dPUdTCa, gr, true, [1 1 1 1 0 1]);
%     evaluateProblem(@dPUdT, [gr(1) 0.0555 gr(2) 3.6484 gr(3:end)], true);
end
