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
features.FV_fnorm = FV_y;
