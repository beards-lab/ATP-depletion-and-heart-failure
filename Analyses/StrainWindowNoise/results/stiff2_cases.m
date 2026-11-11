CASES = { 'SiLU, RelTol 1e-2',           struct('SmoothRectifierEps', 0.01)
          'SiLU, RelTol 1e-3',           struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-3, 'odeAbsTol', 1e-6)
          'SiLU, RelTol 1e-4',           struct('SmoothRectifierEps', 0.01, 'odeRelTol', 1e-4, 'odeAbsTol', 1e-7) };
BUDGET = 300;
