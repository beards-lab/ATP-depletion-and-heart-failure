% RunBest.m -- the current best fit, run twice: once at 8 mM ATP, once at 2 mM.
%
%   matlab -batch "run('RunBest.m')"
%
% ONE parameter file gives BOTH conditions. The low-ATP gate is normalised to 1
% at MgATP_ref = 8 mM, so at 8 mM it is exactly inert: the two runs differ by
% [MgATP] and the protocol file, not by a second fit.
%
% Scored against Model/pooledTarget8mM.m / pooledTarget2mM.m -- three
% preparations (03/27, 04/03, 04/10), each normalised by its own pre-slack
% isometric and rundown-corrected. runSlackExperiment picks the target that
% matches params0.MgATP.
%
% CAVEAT before quoting the 2 mM number: its single gate was fitted with the
% weights and floors DERIVED at 8 mM and rescaled to the 2 mM spread, not
% re-derived. C0, C4 and P1 all converged on this same one gate, every extra
% gate landing at zero -- so k2 is currently the whole low-ATP mechanism here.
% See Docs/experiments/status-2026-09-16.md sec. 7.

cd(fileparts(mfilename('fullpath'))); addpath(genpath(pwd));
if isempty(gcp('nocreate')); parpool(5); end   % slack and FV are parfor

%% ---- 8 mM ATP ----------------------------------------------------------
params0 = getParams(loadParams('lean5_a_opt'), [], true, false);
params0.PlotEachSeparately = true;      % model-vs-data traces, one fig per protocol

RunBakersExp;                           % a SCRIPT: reads params0, writes features_*
fm8 = features_model;  fd8 = features_data;

plotFeatures(fd8, fm8, [], params0.fn, params0);
fprintf('\n8 mM objective %.4f\n', sum(evalFeatureCost(fd8, fm8, params0.fn)));

%% ---- 2 mM ATP ----------------------------------------------------------
% One gate, one protocol file. Everything else is the 8 mM fit untouched.
params0.MgATP = 2;
params0.UseAtpPowerGates = true;
params0.pAtp_k2 = 0.44;                 % k2 -> k2*(MgATP/8)^0.44
params0.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
params0.PlotEachSeparately = true;
figure;
RunBakersExp;
figure;
fm2 = features_model;  fd2 = features_data;

plotFeatures(fd2, fm2, [], params0.fn, params0);
fprintf('\n2 mM objective %.4f\n', sum(evalFeatureCost(fd2, fm2, params0.fn)));

%% ---- the ATP effect, which is why both runs exist ----------------------
fprintf('\n%-8s %8s %8s\n', '2/8 mM', 'model', 'data');
for nm = {'steady','ktr','rsK80'}
    fprintf('%-8s %8.3f %8.3f\n', nm{1}, ...
        mean([fm2.(nm{1})],'omitnan') / mean([fm8.(nm{1})],'omitnan'), ...
        mean([fd2.(nm{1})],'omitnan') / mean([fd8.(nm{1})],'omitnan'));
end
