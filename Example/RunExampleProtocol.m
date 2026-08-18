% RunExampleProtocol.m
% =========================================================================
% Run the synthetic HOLD - STRETCH - HOLD protocol, two ways.
%
%   WAY 1 — DIRECT.   Build a velocity table in memory and hand it straight to
%                     evaluateModel. No files, no cost function, no feature
%                     extraction. This is what you want while exploring: it is
%                     the shortest path from an idea to a force trace.
%
%   WAY 2 — THROUGH THE MACHINERY.  Save the same protocol as a .mat under
%                     data/ and let runStairsExperiment drive it. Now you get
%                     the standard plotting, the data-vs-model overlay, and a
%                     cost term — everything the optimizer would see.
%
% Because way 2 is scored against a "measurement" that way 1 produced, the
% cost must come out at essentially ZERO. That is the point: it proves the
% build -> save -> load -> run round trip is lossless. Once you replace the
% synthetic datatable with real measurements, any cost you then see is
% physics, not plumbing.
%
%   WHY runStairsExperiment AS THE HOST?
%   ------------------------------------
%   Of the four protocol runners, it is the one with no structural
%   expectations: plain trace MSE, no feature extraction, no assumptions about
%   how many rows the velocity table has. runSlackExperiment cannot host this
%   protocol — it segments the table on velocities below -1 ML/s (it is looking
%   for the fast release of a slack cycle) and its internal switch indexes rows
%   literally, so a table of the wrong shape is silently mis-read.
%
% See also: Example/BuildExampleProtocol.m, Model/experiments/runStairsExperiment.m
% =========================================================================

clear; clc; close all;
root = fullfile(fileparts(mfilename('fullpath')), '..');
cd(root);
addpath(genpath(root));


%% [1] BUILD + RUN DIRECTLY ------------------------------------------------
% BuildExampleProtocol is a script, so it runs in THIS workspace and leaves
% behind: vt, out_direct, params_direct, protocolFile.
BuildExampleProtocol;


%% [2] LOOK AT THE DIRECT RESULT ------------------------------------------
% `out` carries far more than force — the strain-resolved state probabilities,
% the serial-elastic extension, the super-relaxed fraction, and so on. A few
% of the most useful traces:
%
%   out.t       time (s)
%   out.Force   total force (kPa)
%   out.SL      sarcomere length (um)
%   out.LXB     mean cross-bridge extension
%   out.SR      super-relaxed (SRX) fraction
figure(10); clf;
set(gcf, 'Position', [80 80 1100 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Window covering the stretch and the hold, i.e. everything except the long
% warm-start. Most of what is worth looking at happens in here.
tZoom = [-0.1, vt(end, 1)];

nexttile;
plot(out_direct.t, out_direct.SL, 'k-', 'LineWidth', 1.5);
ylabel('SL (\mum)'); xlabel('t (s)'); grid on;
title('Commanded length (full run)');
% Everything at t < 0 is the warm-start: the model is activated from rest and
% settles to its isometric steady state before the protocol proper begins.
xline(0, 'r--', 'stretch', 'LabelVerticalAlignment', 'bottom');

nexttile;
plot(out_direct.t, out_direct.Force, 'b-', 'LineWidth', 1.5);
ylabel('Force (kPa)'); xlabel('t (s)'); grid on;
title('Force (full run) — note the activation transient at t < 0');
xline(0, 'r--');

nexttile;
plot(out_direct.t, out_direct.Force, 'b-', 'LineWidth', 1.5);
xlim(tZoom); ylabel('Force (kPa)'); xlabel('t (s)'); grid on;
% WHAT YOU ARE LOOKING AT
% -----------------------
% Force rises almost vertically at the onset of the ramp — that is the elastic
% loading of the cross-bridges that were already attached — and then SATURATES
% rather than continuing to climb. The saturation level is the lengthening
% YIELD point: strained bridges detach as fast as the stretch loads them.
%
% Worth checking for yourself: that peak is set by the ramp VELOCITY, not by
% how far you stretch. Re-run with dL = 0.01 instead of 0.05 in
% BuildExampleProtocol.m and the peak lands in the same place (~72.9 kPa here);
% change V instead and it moves.
%
% The oscillation during the ramp is real model dynamics, not numerical noise
% — repeated detach/re-attach of the strained population. Its amplitude is
% unchanged when the strain-grid spacing params.dS is halved or quartered,
% which is the test that distinguishes dynamics from discretisation error. If
% you push into more violent protocols, run that check again.
title('Force, zoomed: elastic rise, yield, and settling');
xline(0, 'r--');

nexttile;
yyaxis left;
plot(out_direct.t, out_direct.SR, '-', 'LineWidth', 1.5);
ylabel('SRX fraction');
yyaxis right;
plot(out_direct.t, out_direct.LXB, '-', 'LineWidth', 1.5);
ylabel('mean XB extension (\mum)');
xlim(tZoom); xlabel('t (s)'); grid on;
title('Internal states: super-relaxed pool and mean XB extension');

sgtitle('Synthetic hold - stretch - hold (direct evaluateModel call)');


%% [3] RUN THE SAME PROTOCOL THROUGH THE MACHINERY -------------------------
% Point the stairs runner at the file we just wrote. Everything else about
% params0 is unchanged, which is the whole idea: a custom protocol is just
% another file name.
p = params0;
p.stairs_velocitytableonfile = protocolFile;
p.PlotEachSeparately = true;
p.PlotFullscreen     = false;   % false -> the runner calls nexttile
p.ghostLoad = ''; p.ghostSave = '';

figure(11); clf;
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
[E_stairs, out_hosted] = runStairsExperiment(p);
sgtitle('Same protocol, driven by runStairsExperiment');


%% [4] THE ROUND-TRIP SELF-TEST -------------------------------------------
% E_stairs is mean((F_data - F_model)^2) * 10 over the sampled trace. Since
% F_data was generated from an identical simulation, it should be at solver
% round-off, not merely "small".
TOL = 1e-8;
fprintf('\n===== ROUND-TRIP SELF-TEST =====\n');
fprintf('  runStairsExperiment cost E = %.3e   (tolerance %.0e)\n', E_stairs, TOL);
if E_stairs < TOL
    fprintf('  PASS — build -> save -> load -> run is lossless.\n');
else
    fprintf(['  FAIL — the hosted run does not reproduce the direct run.\n' ...
             '  The usual cause is a mismatch in SL0 / Slim_l / Slim_r between\n' ...
             '  BuildExampleProtocol.m cell [2] and runStairsExperiment.m:27-29;\n' ...
             '  a different strain grid gives a slightly different solution.\n']);
end

fprintf(['\nNow make it yours: change L0 / dL / V / T_HOLD at the top of\n' ...
         'Example/BuildExampleProtocol.m and re-run. Replace the synthetic\n' ...
         'datatable with real [t, SL, F] measurements and the same cost term\n' ...
         'becomes a real fit.\n']);
