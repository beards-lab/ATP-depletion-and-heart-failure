% RunFullBatteryCheck.m -- does the new gate form leave the 8 mM fit alone?
%
% RunHillGateConversion checks the slack features only. This runs the FULL high-ATP
% battery (force-velocity, slack, passive) that the 8 mM fit was scored on, under
% three configurations that must be indistinguishable:
%
%   (a) rskR2_w025_opt as-is        -- the pure 8 mM fit, no ATP gates at all
%   (b) + power-law gates @ 8 mM    -- the incumbent low-ATP fit, evaluated at 8 mM
%   (c) + saturable Hill gates @ 8 mM -- this work, evaluated at 8 mM
%
% Both gate forms are pinned to 1 at MgATP_ref, and evaluateModel short-circuits the
% whole block when MgATP == MgATP_ref, so identity is expected BY CONSTRUCTION; this
% is the empirical proof that the implementation actually honours it -- per feature,
% not just in the total.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunFullBatteryCheck.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

HI  = 'rskR2_w025_opt';          % the 8 mM fit the low-ATP work was built on
PWR = 'lowATP_powergate2_opt';   % power-law gates
HIL = 'lowATP_hillgate_opt';     % saturable gates (this work)

lab = {'(a) 8 mM fit, no gates', '(b) + power gates @ 8 mM', '(c) + hill gates @ 8 mM'};
snp = {HI, PWR, HIL};
Et = nan(1,3); FC = cell(1,3); FNM = cell(1,3);

for i = 1:3
    q = getParams(loadParams(snp{i}), [], true, false);
    % every configuration is scored as the 8 mM battery, whatever the snapshot's own
    % condition was -- the gates are what differ, not the protocol
    q.MgATP = 8; q.MgADP = 0; q.Pi = 0;
    q.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
    q.EvalFeatures = true; q.OptimizeOn = 'Feats'; q.BreakOnODEUnstable = false;
    q.RunSlack = 1; q.RunSlackSegments = 'AllPar'; q.RunForceVelocity = 1;
    q.RunSlackPassive = 1; q.RunKtr = 0; q.RunStairs = 0;
    q.PlotEachSeparately = 0; q.PlotFeatureFitting = 0; q.MaxRunTime = 300;
    if i == 1, q.UseAtpPowerGates = 0; end
    t = tic;
    [Et(i), ~, FC{i}, FNM{i}] = evaluateBakersExp([], q, false, false);
    fprintf('%-28s  Et = %.8f   (%.0f s)\n', lab{i}, Et(i), toc(t));
end

fprintf('\n=== 8 mM total cost ===\n');
fprintf('  (b) - (a) = %.3g\n  (c) - (a) = %.3g\n  (c) - (b) = %.3g\n', ...
    Et(2)-Et(1), Et(3)-Et(1), Et(3)-Et(2));

fprintf('\n=== per-feature cost, 8 mM battery ===\n');
fprintf('%-34s %12s %12s %12s\n','feature','(a) no gate','(b) power','(c) hill');
n = min(cellfun(@numel, FC));
worst = 0;
for k = 1:n
    d = max(abs(FC{2}(k)-FC{1}(k)), abs(FC{3}(k)-FC{1}(k)));
    if isfinite(d), worst = max(worst, d); end
    if FC{1}(k) > 1e-6 || d > 1e-9
        fprintf('%-34s %12.6f %12.6f %12.6f\n', FNM{1}{k}, FC{1}(k), FC{2}(k), FC{3}(k));
    end
end
fprintf('\nlargest per-feature deviation from the ungated 8 mM fit: %.3g\n', worst);
if worst == 0
    fprintf('=> the 8 mM battery is BIT-IDENTICAL under both gate forms.\n');
end

save(fullfile(resdir,'fullbattery8mM.mat'),'Et','FC','FNM','lab','snp');
fprintf('\nFULL BATTERY CHECK DONE\n');
