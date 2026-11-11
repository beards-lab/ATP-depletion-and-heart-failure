% BuildExampleProtocol.m
% =========================================================================
% Build a protocol of your own: HOLD - STRETCH - HOLD.
%
% This script shows the whole path from "I want the muscle to do X" to a
% simulated force trace, and then packages the result as a protocol file that
% the standard machinery can consume. Run it directly, or let
% Example/RunExampleProtocol.m call it.
%
%   THE ONE THING TO INTERNALISE
%   ----------------------------
%   The model is driven by a VELOCITY TABLE, and only by its second column:
%
%       velocitytable = [ t(s), v(ML/s), v_um(bookkeeping), L(ML) ]
%
%   Row i's velocity is held constant over the interval [t_i, t_i+1]. Length
%   is whatever integrating that piecewise-constant velocity produces.
%   Columns 3 and 4 are BOOKKEEPING ONLY — nothing reads them during the
%   simulation, so they can silently disagree with column 2 if you edit a row
%   and forget to update them. When you need the commanded length, integrate
%   column 2 (see iCommandedL at the bottom of this file).
%
%   Units: velocity and length are in muscle lengths (ML); the model's
%   sarcomere length SL is in um, with SL = 2 * L for this preparation
%   (params0.ML = 2).
%
% Leaves in the workspace:  vt, out_direct, protocolFile, params_direct
%
% See also: Example/RunExampleProtocol.m, Model/evaluateModel.m,
%           Analyses/RestretchRecoveryFit/RunAmplitudeTest.m (a real use of
%           this same pattern, with a release/re-stretch instead of a stretch)
% =========================================================================

root = fullfile(fileparts(mfilename('fullpath')), '..');
cd(root);
addpath(genpath(root));

% Reuse params0 if a driver already built one; otherwise load the frozen set.
if ~exist('params0', 'var') || ~isstruct(params0)
    params0 = getParams(loadParams('params/ExampleHighATP.m'), [], true, false);
    params0.MaxRunTime = 320;
    params0.ghostLoad = ''; params0.ghostSave = ''; params0.SaveBest = false;
end


%% [1] DESIGN THE PROTOCOL -------------------------------------------------
L0     = 1.00;    % starting length (ML)   -> SL0 = 2.0 um
dL     = 0.05;    % stretch amplitude (ML) -> +5 % of muscle length
V      = 0.50;    % stretch velocity (ML/s)
T_WARM = 2.00;    % warm-start hold before the stretch (s)
T_HOLD = 1.00;    % hold at the stretched length (s)

% Ramp duration follows from amplitude and velocity — deriving it (rather than
% typing a third number) is what keeps column 4 consistent with column 2.
dt     = dL / V;

t_str  = 0;                % stretch starts at t = 0 by convention
t_hold = t_str + dt;       % stretch ends / hold begins
t_end  = t_hold + T_HOLD;

% The warm-start row matters: it lets the cross-bridge population reach its
% isometric steady state before the perturbation, so what you see afterwards
% is the response to the stretch and not to the initial condition.
vt = [ -T_WARM, 0, 0, L0        % warm-start hold
        t_str,  V, 0, L0        % STRETCH ramp   (v applies over [t_str, t_hold])
        t_hold, 0, 0, L0 + dL   % HOLD at length (v applies over [t_hold, t_end])
        t_end,  0, 0, L0 + dL ];% terminal row: its velocity is never used

fprintf('Protocol: hold %.2g s -> stretch %+.3g ML in %.3g s (%.3g ML/s) -> hold %.2g s\n', ...
        T_WARM, dL, dt, V, T_HOLD);


%% [2] RUN IT --------------------------------------------------------------
% The strain grid has to span the excursion, and the initial state vector has
% to be rebuilt for the new SL0 — that is what the `true` (updateInit) argument
% to getParams does. Forgetting it is the classic silent failure: the model
% runs, but from the wrong initial condition.
%
% These three numbers are copied deliberately from runStairsExperiment.m
% (lines 27-29). Cell [3] hands this protocol back to that runner, and the
% round-trip check there is only exact if both runs use the same grid.
params_direct        = params0;
params_direct.SL0    = 2.0;    % um  (= 2 * L0)
params_direct.Slim_l = 1.4;    % um  left edge of the strain grid
params_direct.Slim_r = 2.3;    % um  right edge
params_direct = getParams(params_direct, params_direct.g, true);

% Velocity is assigned AFTER getParams (getParams does not set it, and
% evaluateModel re-runs getParams internally with updateInit = false, so it
% survives).
params_direct.Velocity = vt(:, 2);

modelFcn = str2func(params_direct.modelFcn);
fprintf('Simulating... ');
tSim = tic;
[~, out_direct] = evaluateModel(modelFcn, vt(:, 1), params_direct);
fprintf('%.1f s\n', toc(tSim));


%% [3] PACKAGE IT AS A PROTOCOL FILE ---------------------------------------
% A protocol file is a .mat holding `velocitytable` and `datatable`:
%
%   datatable = [ t(s), SL(um), F ]
%
% For the stairs/ktr runners F is NORMALISED by the isometric force at the
% instant the first length change begins. We reproduce that normalisation
% exactly as runStairsExperiment.m does it (lines 39-43), so that re-running
% this file through that runner reproduces this very trace.
%
% Here the "measurement" IS the model's own output. That is obviously circular
% as science, but it is exactly what you want as a SELF-TEST: it proves the
% build -> save -> load -> run round trip is lossless, so when you later swap
% in real measurements you know any residual is physics, not plumbing.
i_move = find(vt(:, 2) ~= 0, 1);
t_move = vt(i_move, 1);
Fm_ss  = interp1(out_direct.t, out_direct.Force, t_move, 'linear', 'extrap');
if Fm_ss == 0; Fm_ss = 1; end
F_norm = out_direct.Force / Fm_ss;

tq = (vt(1, 1) : 1e-3 : vt(end, 1))';                  % sampling grid (1 kHz)
datatable = [ tq, ...
              2 * iCommandedL(vt, tq), ...             % SL (um) = 2 * L (ML)
              interp1(out_direct.t, F_norm, tq, 'linear', 'extrap') ];
velocitytable = vt;
features_data = struct();                              % none for a synthetic run

% data/ is gitignored, so this generated file stays out of version control.
protocolFile = 'protocol_example_holdstretchhold.mat';
save(fullfile('data', protocolFile), 'velocitytable', 'datatable', 'features_data');
fprintf('Wrote data/%s  (%d samples)\n', protocolFile, size(datatable, 1));

fprintf(['Isometric force at stretch onset: %.3g kPa\n' ...
         'Peak (yield) force during ramp:   %.3g kPa (%.2fx isometric)\n' ...
         'Force at end of hold:             %.3g kPa (%.2fx isometric)\n'], ...
        Fm_ss, max(out_direct.Force), max(F_norm), out_direct.Force(end), F_norm(end));


% -------------------------------------------------------------------------
function Lq = iCommandedL(vt, tq)
%ICOMMANDEDL Length (ML) the velocitytable commands, by integrating column 2.
%   The model is driven by the VELOCITY column, so this is the trajectory it
%   actually follows -- not column 4, which is only a bookkeeping copy and can
%   drift from the integral if a row was edited.
%   (Same helper as Workbench/Tests/TestProtocols0410.m:228.)
    tk = vt(:, 1);  vk = vt(:, 2);
    Lk = zeros(numel(tk), 1);  Lk(1) = vt(1, 4);
    for k = 2:numel(tk)
        Lk(k) = Lk(k-1) + vk(k-1) * (tk(k) - tk(k-1));
    end
    Lq = interp1(tk, Lk, tq(:), 'linear', 'extrap');
end
