function tw = tworldCaLayer(params, ca, prev)
%TWORLDCALAYER  Resolve the T-World myocyte (ThirdParty/TWorld) as the Ca source of the Ca layer.
%
%   TW = TWORLDCALAYER(PARAMS, CA, PREV) is called once by getParams/buildCaLayer
%   when PARAMS.CaInput = 'tworld'. It returns the struct CA.tw read by
%   calciumActivation.m. PREV is the previous CA.tw (or []); when PARAMS.CaTworld
%   and the beat period are unchanged it is returned as is, so the pre-pacing is
%   not repeated on every getParams call (evaluateModel calls getParams again).
%
%   Options, all in the struct PARAMS.CaTworld (a struct is never written into a
%   snapshot by writeParamsToMFile, so none of this leaks into params/*.m):
%     coupling  'signal'  T-World runs standalone for one beat here, at the period
%                         CaPeriod; its cytosolic [Ca2+] of that beat is stored as a
%                         periodic table and read in the ODE (no extra state).
%                         Troponin Ca buffering stays Land's (T-World's own).
%               'land'    the 93 T-World states are integrated INSIDE the cross-
%                         bridge ODE (appended after the thin-filament state).
%                         [Ca2+]_i drives the thin filament; T-World keeps its own
%                         (Land 2017) troponin Ca buffering. Open loop: the same
%                         Ca as 'signal' up to integrator tolerance.
%               'xb'      as 'land', but the troponin Ca buffering flux of T-World
%                         (Bmax_TnClow * dCa_TRPN/dt) is REPLACED by the flux into
%                         the permissive fraction of this model's thin filament
%                         (Bmax_TnClow * dA/dt, A = 1 - N for 'twoState', N_on for
%                         'myosim'). Closed loop: the cross-bridge model's
%                         thin-filament kinetics, including the bound-head feedback
%                         on the off rate, shape the Ca transient.
%               Default 'land'.
%     prebeats  number of standalone beats run here at CaPeriod before the
%               starting state is taken (0 = use X0 / the stored state as is).
%               Default 0. The stored states are 1 Hz steady states, so set
%               prebeats > 0 when CaPeriod ~= 1 (about 0.8 s per beat).
%     X0        explicit T-World starting state (1x93). Default: getStartingState(state).
%     state     getStartingState code, default 'TW_endo' ('TW_epi', 'TW_mid').
%     param     T-World parameter struct overrides (cellType, *_Multiplier, sex,
%               stimAmp, stimDur, nao/cao/ko, ...; see getDefaultParameters.m).
%               bcl is set from CaPeriod; runSignallingPathway is not supported.
%     model     T-World model function handle, default @model_TWorld. Any function with
%               model_TWorld's argument list works, e.g. @model_TWorld_shock_ATP (TW_ATP/,
%               pH/ATP/ADP via param.extraParams = [isHypoxic, T, pH, ATP, ADP, GKATP]);
%               put that model's own helpers (modelRunner, getStartingState, ...) on the path.
%     absTol    AbsTol for the embedded T-World states (mM, mV, fractions);
%               default 1e-6, T-World's own ode15s default. Applied per state by
%               evaluateModel (the cross-bridge states keep their scalar AbsTol).
%
%   Timing inside the cross-bridge ODE (t in s): beats start at
%   CaFirstBeat + k*CaPeriod; T-World sees its own time in ms since the beat onset,
%   so the 1 ms stimulus (stimDur) starts every beat. Before CaFirstBeat the
%   embedded states are frozen (derivative 0) at the pre-paced state and Ca is its
%   diastolic value, so beat 1 starts from exactly the standalone steady state.
%
%   Units: T-World keeps mM and ms internally; the layer converts to uM and 1/s.
%
%   See also: calciumActivation, getParams, ThirdParty/TWorld/README.md

opt = struct();
if isfield(params, 'CaTworld') && isstruct(params.CaTworld), opt = params.CaTworld; end
coupName = getOpt(opt, 'coupling', 'land');
couplings = {'signal', 'land', 'xb'};
coupling = find(strcmp(coupName, couplings));
if isempty(coupling)
    error('tworldCaLayer:coupling', 'Unknown CaTworld.coupling ''%s''. Use one of: %s.', ...
        coupName, strjoin(couplings, ', '));
end
prebeats = getOpt(opt, 'prebeats', 0);
stateCode = getOpt(opt, 'state', 'TW_endo');
X0in = getOpt(opt, 'X0', []);
twp = getOpt(opt, 'param', struct());
absTol = getOpt(opt, 'absTol', 1e-6);
model = getOpt(opt, 'model', @model_TWorld);

% reuse the previous resolution when nothing that shaped it changed
key = struct('T', ca.T, 'coupling', coupling, 'prebeats', prebeats, 'state', stateCode, ...
    'X0', X0in, 'param', twp, 'absTol', absTol, 'model', func2str(model));
if ~isempty(prev) && isstruct(prev) && isfield(prev, 'key') && isequaln(prev.key, key)
    tw = prev;
    return;
end

if exist(func2str(model), 'file') ~= 2
    error('tworldCaLayer:path', ['%s.m not found. Add ThirdParty/TWorld to the path ' ...
        '(addpath(genpath(''.'')) from the repo root).'], func2str(model));
end

twp.model = model;
twp.bcl = 1e3*ca.T;                       % ms
twp = getDefaultParameters(twp);
if twp.runSignallingPathway
    error('tworldCaLayer:signalling', 'runSignallingPathway = 1 (57 PKA states) is not supported in the coupled model.');
end

X0 = X0in;
if isempty(X0), X0 = getStartingState(stateCode); end
X0 = X0(:)';
n = numel(X0);
if n ~= 93
    error('tworldCaLayer:state', 'T-World state has %d entries, expected 93 (no signalling block).', n);
end
if prebeats > 0
    [~, X] = modelRunner(X0, [], twp, prebeats, prebeats - 1);
    X0 = X{1}(end, :);
end

% the argument list of model_TWorld after (t, y, flag_ode); same order as modelRunner.m
args = {twp.cellType, twp.ICaL_Multiplier, twp.INa_Multiplier, twp.Itos_Multiplier, twp.Itof_Multiplier, ...
    twp.INaL_Multiplier, twp.IKr_Multiplier, twp.IKs_Multiplier, twp.IK1_Multiplier, twp.IKb_Multiplier, ...
    twp.INaCa_Multiplier, twp.INaK_Multiplier, twp.INab_Multiplier, twp.ICab_Multiplier, twp.IpCa_Multiplier, ...
    twp.IClCa_Multiplier, twp.IClb_Multiplier, twp.Jrel_Multiplier, twp.Jup_Multiplier, ...
    twp.nao, twp.cao, twp.ko, twp.ICaL_fraction_junc, twp.INaCa_fraction_junc, twp.PKA_P, twp.Whole_cell_PP1, ...
    twp.runSignallingPathway, [], twp.sex, twp.stimAmp, twp.stimDur, twp.vcParameters, twp.apClamp, twp.extraParams};

tw = struct('key', key, 'coupling', coupling, 'model', model, 'n', n, 'x0', X0, 'args', {args}, ...
    'stimDur', twp.stimDur/1e3, ...          % s
    'Cadia', 1e3*X0(8), ...                  % uM, cytosolic Ca of the starting state
    'Bmax_TnClow', 70e-3, ...                % mM, troponin C low-affinity sites (model_TWorld.m, Buffering)
    'iCa', 8, 'iTnC', 13, 'iCaTRPN', 71, ... % state indices used by the 'xb' coupling
    'absTol', absTol*ones(1, n), 'param', twp, 'sig', [], 'ref', []);

% one standalone beat from the starting state: the periodic Ca table for 'signal',
% and the reference traces (Ca, Land tension) every coupling can be compared with
[tm, X] = modelRunner(X0, [], twp, 1, 0);
c = getCurrentsStructure(tm, X, twp, 0);
t_s = c.time/1e3;  Ca_uM = 1e3*c.Ca_i;
tw.ref = struct('t', t_s, 'Ca', Ca_uM, 'Ta', c.Land_Ta, 'V', c.V, 'X', X{1});
Fs = 4000;
tg = (0:1/Fs:ca.T)';
[tu, iu] = unique(t_s);
y = interp1(tu, Ca_uM(iu), tg, 'pchip');
tw.sig = struct('t0', 0, 'Fs', Fs, 'N', numel(tg), 'y', y(:), 'm', gradient(y(:))*Fs, 'hermite', true);
end

function v = getOpt(opt, name, default)
if isfield(opt, name) && ~isempty(opt.(name)), v = opt.(name); else, v = default; end
end
