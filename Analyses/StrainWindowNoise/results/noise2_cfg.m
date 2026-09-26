CFG = { 'SiLU + dense (RelTol 1e-2)',              struct('SmoothRectifierEps', 0.01, 'odeOutputDt', 5e-4)
        'SiLU + RelTol 1e-3 + dense',              struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4)
        'SiLU + RelTol 1e-4 + dense',              struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-4, 'odeAbsTol', 1e-7, 'odeOutputDt', 5e-4)
        'SiLU + 1e-3 + dense + wide 128 + clamp',  struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6, 'odeOutputDt', 5e-4, 'MaxStrainArraySize', 128, 'ClampStrainDepExtrapolation', true) };
OUTTAG = 'refit8_opt_round2';
