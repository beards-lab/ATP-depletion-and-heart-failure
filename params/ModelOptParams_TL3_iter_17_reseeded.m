% ModelOptParams_TL3_iter_17_reseeded.m
% =========================================================================
% Rate-reseeding overlay — load ON TOP of ModelOptParams_TL3_iter_17.m
%
% Purpose: bring kah, kd, k2, ka from their below-floor iter_17 values up
% to physiologically motivated starting points for mouse/rat alpha-MHC at
% ~20-25 deg C before any optimizer run. Forces coupled mechanics to shift
% (especially attached fraction and force amplitude); kstiff1/kstiff2 may
% need manual adjustment once the new steady state is evaluated.
%
% Sources (see Model/parameterBounds.m and Docs/sources.md for details):
%   kah = 80   Mijailovich/Walklate 2017 alpha-cardiac kH = 77/s at 20C;
%              matches getParams default
%   kd  = 100  Wang 2013 mouse alpha skinned ~111/s at 17C;
%              Little 2012 rat alpha ~159/s at 0-ADP; matches getParams default
%   k2  = 350  Lowey 2013 mouse alpha-S1 ADP release ~350/s at 20C
%   ka  = 375  Beard 2022 parent-model fit to rat cardiac FV/ktr at 25C;
%              matches getParams default
%   kamh       NOT set here - resolved via formula '=0.1*kah' in getParams;
%              iter_17 hard-codes 2.82418 which is overridden by the getParams
%              formula only if the field is removed or set to the formula string.
%              Reseeding restores it to 0.1*80 = 8/s (within Tier B [1,30]/s).
%
% Parameter violations fixed:
%   kah: 20.1 -> 80  (was below Tier A lb=40, weight=100)
%   kd:  10.0 -> 100 (was below Tier A lb=20, weight=100)
%   k2:  69.5 -> 350 (was below Tier B lb=100, weight=10)
%   ka:  37.7 -> 375 (was below Tier B lb=50, weight=10)
% =========================================================================

% --- Kinetic rates (all base rates, xrate=1 in iter_17) ------------------
params0.kah  = 80;          % ATP hydrolysis:  was 20.1   floor=40   s^-1
params0.kamh = '=0.1*kah';  % reverse hydrolysis: restore linked formula
                             % (iter_17 hard-coded 2.82418; formula gives 8/s)
params0.ka   = 375;         % attachment:      was 37.7   floor=50   s^-1
params0.kd   = 100;         % P1 detachment:   was 10.0   floor=20   s^-1
params0.k2   = 350;         % P2 ADP release:  was 69.5   floor=100  s^-1

% NOTE: k1, k_1, k_2 are NOT changed here. They may need rebalancing after
% evaluating the new steady state (especially k1 which sets the power-stroke
% flux and couples to attached fraction via P1->P2 forward rate).

% --- Mechanical: placeholder for rebalancing ----------------------------
% Raising kah and ka will increase attached fraction; force amplitude may
% overshoot. Adjust kstiff1/kstiff2 here once DriverBoundedFit is evaluated.
% (Uncomment and set as needed after first scoring run.)
% params0.kstiff1 = 60383.9;   % current iter_17 value -- adjust if needed
% params0.kstiff2 = 6.9794e+04*0.5;  % already halved in DriverBoundedFit

% --- No mods/g changes: this is a manual reseed, not an optimizer run ---
% params0.mods and params0.g are left as-is (empty from iter_17 tail).
