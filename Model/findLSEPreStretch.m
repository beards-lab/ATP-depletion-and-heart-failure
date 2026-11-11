function [rsl] = findLSEPreStretch(SL0, ForceAtSL0, params)

params.SL0 = SL0;
tol = 1e-6;        % tolerance for stopping
maxIter = 50;      % safety cap
sL = 0;            % lower bound
sR = SL0;          % upper bound

% Evaluate forces at bounds
fL = simulateForce(SL0, sL, params) - ForceAtSL0;
fR = simulateForce(SL0, sR, params) - ForceAtSL0;

if fL * fR > 0
    error('Bisection requires a sign change. Adjust [sL, sR].');
end

for iter = 1:maxIter
    sMid = (sL + sR)/2;
    fMid = simulateForce(SL0, sMid, params) - ForceAtSL0;

    if abs(fMid) < tol
        break; % solution found
    end

    if fL * fMid < 0
        sR = sMid; fR = fMid;
    else
        sL = sMid; fL = fMid;
    end
end

% prestretch ratio
rsl = (SL0 + sMid)/SL0;
% fprintf('Solution: s = %.6f, force ≈ %.4f\n', s_solution, simulateForce(s_solution, params));


function f = simulateForce(SL0, LSE, params)
    params.SL0 = SL0 + LSE;
    params.rsl0 = 1;
    params = getParams(params, [], true);
    params.Vums = 0;
    params.PU0(params.NumberOfStates*params.ss + 4) = LSE;

    modelFcn = str2func(params.modelFcn);
    [~, outs] = modelFcn(0, params.PU0', params);
    f = outs(1);
end

end