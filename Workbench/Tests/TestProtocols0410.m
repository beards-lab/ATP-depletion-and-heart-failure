% TestProtocols0410.m
% Validation pass for the 04/10/2026 (Male 2-8) protocol files.
%
% Two things are checked for every new .mat:
%
%   1. VELOCITYTABLE / DATATABLE CONSISTENCY (no simulation involved).
%      The velocitytable drives the model through its VELOCITY column, while
%      the cost is scored against the datatable's measured length/force. If the
%      two disagree, the model is silently driven along a different length
%      trajectory than the data was recorded at, and every feature is off. So
%      integrate the velocity column into a commanded length and compare it,
%      on the datatable's own time grid, against the measured length.
%
%   2. THAT IT SIMULATES. Each file is driven through its runXExperiment with
%      the ThursdayNightFever baseline parameters, and the run is timed and
%      error-trapped.
%
% Output: three figures (slack / ktr / stairs+FV2), each with a length row
% (commanded + measured + simulated SL) and a force row (measured vs simulated),
% exported to Docs/figures/.
%
% The FV2 files are driven through runKtrExperiment: their datatable force is
% normalised by the pre-ramp isometric, exactly the convention that runner
% expects, so no dedicated experiment function is needed.

clear; close all; clc;
cd(fileparts(mfilename('fullpath')));
addpath(genpath('../..'));

figDir = '../../Docs/figures/';
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% ── Baseline parameters ────────────────────────────────────────────────────
params0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
params0.PlotEachSeparately = 0;      % this script does its own plotting
params0.PlotFullscreen     = 0;
params0.PlotFeatureFitting = 0;
params0.ghostSave          = '';
params0.ghostLoad          = '';
params0.EvalFeatures       = 1;
params0.RunSlackSegments   = 'All';  % 'AllPar' needs a parpool; 'All' is equivalent

%% ── What to validate ───────────────────────────────────────────────────────
% {label, kind, file, figure group}
specs = {
  '8mM slack',       'slack',   'protocol_04_10_2026_8mM_slack.mat',            1
  '2mM slack',       'slack',   'protocol_04_10_2026_2mM_slack.mat',            1
  'PNB+Mava slack',  'passive', 'protocol_04_10_2026_ActivePNBMava_slack.mat',  1
  '8mM ktr',         'ktr',    'protocol_04_10_2026_velocitytable_ktr.mat',    2
  '2mM ktr',         'ktr',    'protocol_04_10_2026_velocitytable_ktr_2mM.mat',2
  'PNB+Mava ktr',    'ktr',    'protocol_04_10_2026_velocitytable_ktr_PNBMava.mat', 2
  '8mM stairs',      'stairs', 'protocol_04_10_2026_velocitytable_stairs.mat', 3
  '8mM FV2',         'ktr',    'protocol_04_10_2026_velocitytable_FV2.mat',    3
  '2mM FV2',         'ktr',    'protocol_04_10_2026_velocitytable_FV2_2mM.mat',3
  'PNB+Mava FV2',    'ktr',    'protocol_04_10_2026_velocitytable_FV2_PNBMava.mat', 3
};
nSpec = size(specs, 1);

groupNames = {'Slack', 'Ktr', 'Stairs + FV2'};
groupCols  = arrayfun(@(g) sum([specs{:,4}] == g), 1:3);
for g = 1:3
    figure(300+g); clf;
    set(gcf, 'Position', [40 40 420*groupCols(g) 640], 'Name', groupNames{g});
    tiledlayout(2, groupCols(g), 'TileSpacing', 'compact', 'Padding', 'compact');
end
groupSlot = zeros(1, 3);

R = struct('label', {}, 'kind', {}, 'ok', {}, 'E', {}, 'secs', {}, ...
           'vtMismatch', {}, 'vtHold', {}, 'slMod', {}, 'nVT', {}, 'nDT', {}, 'msg', {});

%% ── Run ────────────────────────────────────────────────────────────────────
for i = 1:nSpec
    label = specs{i,1};  kind = specs{i,2};  fname = specs{i,3};  grp = specs{i,4};
    fprintf('\n===== %s  (%s) =====\n', label, fname);

    D  = load(['data/' fname]);
    vt = D.velocitytable;  dt = D.datatable;

    % ── check 1: does the velocitytable's velocity column reproduce the data
    %             length trace it will be scored against? ────────────────────
    Lcmd_dt = iCommandedL(vt, dt(:,1));          % ML
    Lmeas   = dt(:,2) / 2;                        % datatable length is um -> ML
    inWin   = dt(:,1) >= vt(1,1) & dt(:,1) <= vt(end,1);
    vtMismatch = max(abs(Lcmd_dt(inWin) - Lmeas(inWin)));

    % Split out the HOLD segments. The overall max is dominated by corner
    % rounding: the velocitytable is piecewise-linear while the servo has a
    % finite rise time, so every transition contributes ~1e-2 ML that means
    % nothing. A mismatch on the holds is the one that matters -- it is a real
    % length offset between what drives the model and what the data was at.
    isHold = false(size(dt,1), 1);
    for k = 1:size(vt,1)-1
        if abs(vt(k,2)) < 0.05
            isHold = isHold | (dt(:,1) > vt(k,1) + 0.002 & dt(:,1) < vt(k+1,1) - 0.002);
        end
    end
    holdSel = isHold & inWin;
    if any(holdSel)
        vtHoldMismatch = max(abs(Lcmd_dt(holdSel) - Lmeas(holdSel)));
    else
        vtHoldMismatch = NaN;
    end

    % ── check 2: does it simulate? ────────────────────────────────────────
    p = params0;
    ok = true; msg = ''; E = NaN; out = [];
    tic;
    try
        switch kind
            case 'slack'
                p.velocitytableonfile = fname;
                [E, out] = runSlackExperiment(p);
                E = E(1);
            case 'passive'
                % The PNB/Mava file's consumer is runPassiveExperiment (ka=0),
                % not the active slack runner. It returns features rather than
                % a cost, so score the trace RMSE here for display.
                p.passive_velocitytableonfile = fname;
                [~, ~, out] = runPassiveExperiment(p);
                mono = makeMonotonous(out.t);
                Fi   = interp1(out.t(mono), out.Force(mono), dt(:,1), 'linear', NaN);
                E    = sqrt(mean((dt(:,3) - Fi).^2, 'omitnan'));
            case 'ktr'
                p.ktr_velocitytableonfile = fname;
                p.Slim_l = 1.4; p.Slim_r = 2.3;   % protocol reaches SL 1.60-2.20 um
                [E, out] = runKtrExperiment(p, []);
            case 'stairs'
                p.stairs_velocitytableonfile = fname;
                [E, out] = runStairsExperiment(p);
        end
    catch ME
        ok = false; msg = ME.message;
        fprintf('  *** FAILED: %s\n', msg);
    end
    secs = toc;

    % Did the model actually follow the commanded length?
    if ok
        slMod = max(abs(interp1(out.t, out.SL, dt(holdSel,1), 'linear', NaN) ...
                        - 2*Lcmd_dt(holdSel)), [], 'omitnan');
    else
        slMod = NaN;
    end

    R(end+1) = struct('label', label, 'kind', kind, 'ok', ok, 'E', E, ...
        'secs', secs, 'vtMismatch', vtMismatch, 'vtHold', vtHoldMismatch, ...
        'slMod', slMod, 'nVT', size(vt,1), 'nDT', size(dt,1), 'msg', msg); %#ok<SAGROW>

    fprintf('  VT %dx4, DT %dx3 | commanded-vs-measured L: max=%.2e, on holds=%.2e ML | model SL vs commanded (holds)=%.2e um\n', ...
            size(vt,1), size(dt,1), vtMismatch, vtHoldMismatch, slMod);
    if ok
        fprintf('  simulated OK in %.1f s, E = %.4g\n', secs, E);
    end

    % ── plot ────────────────────────────────────────────────────────────────
    groupSlot(grp) = groupSlot(grp) + 1;
    col = groupSlot(grp);
    figure(300+grp);

    % model force normalisation matches what each runner scores against
    if ok
        if any(strcmp(kind, {'slack', 'passive'}))
            Fmod = out.Force;  ylab = 'Force (kPa)';
        else
            Fss_m = interp1(out.t, out.Force, vt(2,1), 'linear', 'extrap');
            if Fss_m == 0; Fss_m = 1; end
            Fmod = out.Force / Fss_m;  ylab = 'Force (norm.)';
        end
    else
        ylab = 'Force';
    end

    xl = [vt(2,1) - 0.02*(vt(end,1)-vt(2,1)), vt(end,1)];

    % top row: lengths
    nexttile(col); hold on; box on;
    tq = linspace(vt(1,1), vt(end,1), 4000);
    plot(tq, 2*iCommandedL(vt, tq), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5);
    plot(dt(:,1), dt(:,2), 'k-', 'LineWidth', 0.75);
    if ok; plot(out.t, out.SL, 'b--', 'LineWidth', 1.25); end
    ylabel('SL (\mum)');
    title(sprintf('%s  (\\Delta L hold = %.1e ML)', label, vtHoldMismatch));
    xlim(xl);
    if col == 1; legend({'VT commanded', 'data', 'model SL'}, 'Location', 'best', 'FontSize', 7); end

    % bottom row: forces
    nexttile(groupCols(grp) + col); hold on; box on;
    plot(dt(:,1), dt(:,3), 'k-', 'LineWidth', 0.75);
    if ok; plot(out.t, Fmod, 'b-', 'LineWidth', 1.25); end
    ylabel(ylab); xlabel('Time (s)');
    if ok
        title(sprintf('E = %.4g   (%.1f s)', E, secs));
    else
        title('SIMULATION FAILED', 'Color', 'r');
    end
    xlim(xl);
    if col == 1; legend({'data', 'model'}, 'Location', 'best', 'FontSize', 7); end
end

%% ── Summary ────────────────────────────────────────────────────────────────
fprintf('\n\n================ 04/10/2026 PROTOCOL VALIDATION ================\n');
fprintf('%-16s %-8s %-5s %5s %7s %10s %10s %10s %11s\n', ...
    'protocol', 'kind', 'runs', 'VT', 'DT', 'dL max', 'dL hold', 'SLmod err', 'E');
for i = 1:numel(R)
    if R(i).ok; oks = 'yes'; else; oks = 'NO'; end
    fprintf('%-16s %-8s %-5s %5d %7d %10.2e %10.2e %10.2e %11.4g\n', ...
        R(i).label, R(i).kind, oks, R(i).nVT, R(i).nDT, ...
        R(i).vtMismatch, R(i).vtHold, R(i).slMod, R(i).E);
end
nFail = sum(~[R.ok]);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('%d/%d simulate.\n', numel(R)-nFail, numel(R));
fprintf('dL max  : commanded vs measured length, worst point (dominated by servo corner rounding)\n');
fprintf('dL hold : same, restricted to held segments -- this is the one that must be ~0\n');
fprintf('SLmod   : simulated SL vs commanded SL on holds (um) -- catches a wrong SL0\n');
fprintf('worst dL on holds = %.2e ML,  worst model SL error = %.2e um\n', ...
        max([R.vtHold]), max([R.slMod]));

for g = 1:3
    figure(300+g);
    sgtitle(sprintf('04/10/2026 Male 2-8 — %s protocols (ThursdayNightFever params)', groupNames{g}));
    exportgraphics(gcf, fullfile(figDir, sprintf('validate_0410_%d_%s.png', g, ...
        matlab.lang.makeValidName(groupNames{g}))), 'Resolution', 150);
end
fprintf('Figures exported to %s\n', figDir);

% ---------------------------------------------------------------------------
function Lq = iCommandedL(vt, tq)
%ICOMMANDEDL Length (ML) the velocitytable commands, by integrating column 2.
%   The model is driven by the VELOCITY column, so this is the trajectory it
%   actually follows -- not column 4, which is only a bookkeeping copy and can
%   drift from the integral if a row was edited.
    tk = vt(:,1);  vk = vt(:,2);
    Lk = zeros(numel(tk), 1);  Lk(1) = vt(1,4);
    for k = 2:numel(tk)
        Lk(k) = Lk(k-1) + vk(k-1) * (tk(k) - tk(k-1));
    end
    Lq = interp1(tk, Lk, tq(:), 'linear', 'extrap');
end
