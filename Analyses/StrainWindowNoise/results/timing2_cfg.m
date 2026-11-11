fx = @(varargin) struct('SmoothRectifierEps', 0.01, 'ClampStrainDepExtrapolation', 1, varargin{:});
VAR = { 'old numerics',                         struct('SmoothRectifierEps', 0, 'odeRelTol', 1e-2, 'odeAbsTol', 1e-4, 'odeOutputDt', 0, 'ClampStrainDepExtrapolation', 0)
        'fix 1e-3, dense 0.5 ms everywhere',    fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4)
        'fix 1e-3, dense only segments <= 5 s', fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5)
        'fix 1e-3, no dense',                   fx('odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 0)
        'fix 1e-4, dense 0.5 ms everywhere',    fx('odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4)
        'fix 1e-4, dense only segments <= 5 s', fx('odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4, 'odeOutputMaxSegment', 5) };
OUTTAG = 'refit9_opt_clean'; NEV = 3;
