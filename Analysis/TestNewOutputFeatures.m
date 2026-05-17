% TestNewOutputFeatures.m
% One-shot test of new output-bounding features extracted from extractSlackAttributes.
% Runs a single param file through the slack experiment and displays
% the new feature values + the feature cost figure.
%
% Verifies:
%   - SRX_ss, attached_ss, PT_ss extracted at pre-slack steady state
%   - XTOR_vmax extracted at maximum shortening speed during slack
%   - t0_crossing extracted as SL-based slack time
%   - Grouped bar tile renders correctly in figure 80085

clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(root));

p = parpool('Threads', 5);

%% Load param file
params0 = getParams();
% param_file = fullfile(root, 'params', 'ModelParamsWPassive_PreOpt.m');
% run(param_file);
ModelParamsWPassive_PreOpt
% ModelParamsFeats_FVSlackUpdateValley2

%% Configure for slack experiment only
params0.RunForceVelocity             = true;
params0.RunKtr                       = false;
params0.RunSlack                     = true;
params0.RunStairs                    = false;
params0.RunForceLengthEstim          = false;
params0.RunForceVelocityTime         = false;
params0.RunMinStretch                = false;
params0.justPlotStateTransitionsFlag = false;
params0.BreakOnODEUnstable           = false;
params0.ghostSave                    = '';
params0.ghostLoad                    = '';
params0.PlotEachSeparately           = true;
params0.ShowStatePlots               = false;
params0.velocitytableonfile          = 'protocol_03_27_2026_8mM_slack.mat';
params0.RunSlackSegments             = 'All';
params0.MaxRunTime                   = 120;
params0.OptimizeOn                   = 'Feats';

params0.fn = {
    'FV_f|FV_v', 'ktr|2', 'A|50', ...
    'ktr_rmse|0-0.1', ...
    'XTOR|2-8', ...
    'XTOR_vmax|4-30', ...
    't0_crossing|SLdiff|2', ...
    'SRX_ss|0.25-0.75,attached_ss|0.02-0.45,PT_ss|0.02-0.70', ...
    'restretchSlopeStart|0.1', 'restretchSlopeEnd|0.1', ...
    'peak1_y|10', 'peak1_dSL|0.2', ...
    'vall_y|10', 'vall_t|0.2', 'peak2|5', ...
    'steady|50', 'vall2_dy|0.1', 'ovrsht_dy|0.1'
};

%% Run
features_data  = struct();
features_model = struct();
ModelParams_Passive;
params0.xrate = 2;

params0.SRXD_0 = .15;

params0.UseSuperRelaxed = false;
params0.UseSuperRelaxedADP = false;

figure(1); clf;
params0.RunSlackSegments = 'AllPar';
tic
RunBakersExp;
toc

%% Print new features
fprintf('\n=== New output features (model) ===\n');
if isfield(features_model, 'SRX_ss')
    fprintf('SRX_ss:       %s\n', mat2str(features_model.SRX_ss, 3));
end
if isfield(features_model, 'attached_ss')
    fprintf('attached_ss:  %s\n', mat2str(features_model.attached_ss, 3));
end
if isfield(features_model, 'PT_ss')
    fprintf('PT_ss:        %s\n', mat2str(features_model.PT_ss, 3));
end
if isfield(features_model, 'XTOR_vmax')
    fprintf('XTOR_vmax:    %s\n', mat2str(features_model.XTOR_vmax, 3));
end
if isfield(features_model, 'XTOR')
    fprintf('XTOR (ss):    %s\n', mat2str(features_model.XTOR, 3));
end
if isfield(features_model, 't0_crossing')
    fprintf('t0_crossing:  %s\n', mat2str(features_model.t0_crossing, 3));
end
if isfield(features_model, 't0')
    fprintf('t0 (force):   %s\n', mat2str(features_model.t0, 3));
end
fprintf('\n=== Conservation check (should sum to ~1 per segment) ===\n');
if isfield(features_model, 'SRX_ss') && isfield(features_model, 'attached_ss') && isfield(features_model, 'PT_ss')
    row_sum = features_model.SRX_ss + features_model.attached_ss + features_model.PT_ss;
    fprintf('SRX+attached+PT: %s (remaining = PuR + NP)\n', mat2str(row_sum, 3));
end

%% Feature cost figure
if ~isempty(fieldnames(features_data)) && ~isempty(fieldnames(features_model))
    figure(80085); clf;
    cost_vec = plotFeatures(features_data, features_model, [], params0.fn, params0);
    annotation('textbox', [0 0.96 1 0.04], ...
        'String', sprintf('Total cost: %.3f', sum(cost_vec, 'omitnan')), ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 11);
    fprintf('\nTotal cost: %.3f\n', sum(cost_vec, 'omitnan'));
end
