function features = extractForceVelocityAttributes(velocities, F_active, features, evaluateSpeeds, plotResults)

if nargin < 5 
    plotResults = false;
end
if nargin < 4 
    evaluateSpeeds   = -[0.5 1 2 3 4 5];  
end

scaleVal  = -0.5;

win = find(ismember(velocities, evaluateSpeeds));

if isempty(win)
    error('Velocities do not include any from the window');
end

scaleTo = find(velocities == scaleVal, 1);

% win = 2:8;
% win = [2 3 4 5 6];
% win = 2:length(velocities);
% scaleTo = 2;

FV_x = -velocities(win)';

if isempty(scaleTo)
    FV_y = F_active(win)';
else
    FV_y = F_active(win)'/F_active(scaleTo);
end

FVfit_fun = @(a, b, c, x) a +b*exp(-(x)*c);
[fvfit fvgood] = fit(FV_x, FV_y, FVfit_fun, 'StartPoint', [1, 100, 0.5], 'Lower', [0 0 0]);
FV_t = (0:0.2:7);
if plotResults
    nexttile;
    plot(FV_y, FV_x, 'x', fvfit(FV_t), FV_t, '-', LineWidth=1.5, MarkerSize=12);
end

features.FV_expFit = fvfit;
features.FV_expFit_v = FV_t;
features.FV_expFit_f = fvfit(features.FV_expFit_v);
features.FV_v = FV_x;
features.FV_f = F_active(win)';
% force normalized
features.FV_fnorm = features.FV_f./features.FV_f(1);

% ---- POWER --------------------------------------------------------------
% Power = force x shortening velocity. It is the observable the ATP effect
% actually shows up in (peak power and the velocity it occurs at both move
% between 8 and 2 mM), and unlike FV_fnorm it is not dominated by the
% isometric point, which carries no information about cross-bridge kinetics.
%
% The MODEL has a single FV curve, so the per-dataset fields are aliases of
% that one curve. They exist so the same fn entry can score the model against
% the two-dataset average built on the data side (see runFVExperiment), which
% is where FV_power1/2 genuinely differ. FV_normpowerVar is a property of the
% DATA spread only; on the model side it is zero by construction and is
% emitted purely so a stray fn entry cannot score MISSING_FEATURE_COST.
% The v = 0 point is EXCLUDED: power there is identically zero in every dataset
% and in every model, by definition. Including it would contribute exactly zero
% error while soaking up weight mass, making the inverse-variance weights on the
% informative points smaller than they appear. FV_powerv records the velocities
% these arrays are on, since they are shorter than FV_v.
ipw = features.FV_v > 0;
features.FV_powerv    = features.FV_v(ipw);
features.FV_power     = features.FV_f(ipw)     .* features.FV_powerv;
features.FV_normpower = features.FV_fnorm(ipw) .* features.FV_powerv;

features.FV_power1       = features.FV_power;
features.FV_power2       = features.FV_power;
features.FV_normpower1   = features.FV_normpower;
features.FV_normpower2   = features.FV_normpower;
features.FV_normpowerAvg = features.FV_normpower;
features.FV_normpowerVar = zeros(size(features.FV_normpower));
