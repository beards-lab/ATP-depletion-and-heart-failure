function cfg = tierConfigs(name)
%TIERCONFIGS  Fitting tiers for the human-iPSC EHT twitch with a FIXED literature Ca transient.
%
%   Stage 2 ("humanizing trajectories only"): the mouse incumbent (robust1_seed) is moved along
%   the literature directions of Analyses/HumanizeXB/conclusions.md:
%     r_iso  isoform (alpha -> beta MHC) ratio on k2, k_2, kah, kamh, ka, kd   (literature 0.1)
%     r_Pi   Pi-release / power-stroke ratio on k1, k_1                         (literature 0.3, weak)
%     r_T    temperature (rig 20-25 C -> 37 C, Q10 3-3.5) on ALL cycle rates    (literature x3.7-8.4)
%     r_SRX  SRX/OFF exchange rates (no species evidence; Q10 would apply)
%     r_kSE  series compliance (the EHT pulls on compliant posts; unknown stiffness)
%     r_pas  passive titin (fetal N2BA in hiPSC-CM: softer)
%     SL0    sarcomere length (unknown in the EHT)
%   The thin filament (tf_k_on, tf_k_off, tf_K_coop) was never calibrated on any data, so it is
%   free in every tier within the bounds used by the EHT report (K_coop <= 20), as is the beat
%   offset (stimulus record not synchronised). sigma (uN per kPa) is solved in closed form.
%   Stage 3 ("anything but the Ca transient"): individual multipliers on every rate, stiffness,
%   SRX and gate parameter, wide thin-filament bounds.
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
fullName = name;  prof = regexp(name, '^(?<base>.+?)__(?<par>\w+)__(?<k>\d+)$', 'names');
if ~isempty(prof), name = prof.base; end        % profile point: build the base tier, then fix one parameter
cfg.name = fullName;  cfg.seed = 'robust1_seed';  cfg.pre = 2;  cfg.nb = 2;  cfg.fdStep = 2e-3;  cfg.maxIter = 40;
if endsWith(name, '_ex')            % variant on the EHT report's core (ExampleHighATP = spreadw_opt) instead of the HumanizeXB incumbent
    cfg.seed = fullfile(here, '..', '..', 'params', 'ExampleHighATP.m');  name = extractBefore(name, '_ex');
end
b = loadParams(cfg.seed);  b.mods = {};  b.g = [];  b.SL0 = 2.0;
cfg.base = getParams(b, [], true, false);
cfg.ca = struct('kind', 'risedecay', 'opt', struct('Cadia', 0.15, 'Camax', 0.60, 'tau_r', 0.025, 'm', 1, 't50', 0.265, 'hdec', 3, 't0', 0.005, 'nsup', 3));
cfg.CaModel = 'twoState';
cyc = {'k2', 'k_2', 'kah', 'kamh', 'ka', 'kd', 'k1', 'k_1'};
iso = cyc(1:6);  pi_ = cyc(7:8);
srx = {'ksr0', 'kmsrd', 'ksrd', 'ksr2srd', 'ksrd2sr'};
P = @(nm, type, fields, lb, ub, x0, islog) struct('name', nm, 'type', type, 'fields', {fields}, 'lb', lb, 'ub', ub, 'x0', x0, 'log', islog);
% thin filament + offset (shared by every tier); x0 from the EHT fit C rescaled to the new Ca gauge
tf = [P('tf_k_on', 'tf', {}, log10(0.1), log10(3000), log10(7), true), ...
      P('tf_k_off', 'tf', {}, log10(0.5), log10(2000), log10(5), true), ...
      P('tf_K_coop', 'tf', {}, 1, 20, 12, false), ...
      P('shift', 'shift', {}, -0.150, 0.060, 0, false)];
hum = [P('r_iso', 'scale', iso, log10(0.03), log10(0.5), log10(0.1), true), ...
       P('r_T', 'scale', cyc, log10(1), log10(15), log10(5), true)];
switch name
    case 't0_v0'      % literature core, no temperature factor: human_v0 ratios fixed
        cfg.fixedScale = fixed(iso, 0.1, pi_, 0.3);
        cfg.par = tf;
    case {'t1_isoT', 't1_caFast', 't1_caSlow', 't1_caHiLvl'}   % tier 1; the ca* variants test the Ca-transient assumption
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, tf];  cfg.maxIter = 50;
        switch name
            case 't1_caFast',  cfg.ca.opt.tau_r = 0.018;  cfg.ca.opt.t50 = 0.200;   % TTP ~70 ms, peak->50 % ~140 ms (fast edge of the range)
            case 't1_caSlow',  cfg.ca.opt.tau_r = 0.035;  cfg.ca.opt.t50 = 0.330;   % TTP ~120 ms, peak->50 % ~230 ms (slow edge)
            case 't1_caHiLvl', cfg.ca.opt.Cadia = 0.25;   cfg.ca.opt.Camax = 1.0;   % high levels (adult-clamp-like), diastolic/peak ratio 0.25
        end
    case 't2a_Pi'
        cfg.par = [hum, P('r_Pi', 'scale', pi_, log10(0.1), log10(1), log10(0.3), true), tf];
    case 't2b_SRX'
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, P('r_SRX', 'scale', srx, log10(0.1), log10(15), 0, true), tf];
    case 't2b_SRXc'   % SRX direction alone: cycle rates held at the tier-1 optimum (r_iso 0.0655, r_T 6.93)
        T1 = load(fullfile(res, 't1_isoT_fit.mat'));  x1 = T1.R.x;
        cfg.fixedScale = fixed(iso, 10^(x1(1) + x1(2)), pi_, 0.3*10^x1(2));
        cfg.par = [P('r_SRX', 'scale', srx, log10(0.1), log10(15), 0, true), tf];
    case 't2c_SE'
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, P('r_kSE', 'scale', {'kSE'}, log10(0.01), log10(3), 0, true), tf];
    case 't2d_pas'
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, P('r_pas', 'scale', {'k_pas'}, log10(0.1), log10(3), 0, true), tf];
    case 't2e_SL'
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, P('SL0', 'abs', {'SL0'}, 1.9, 2.3, 2.0, false), tf];
    case 't2_all'
        cfg.par = [hum, P('r_Pi', 'scale', pi_, log10(0.1), log10(1), log10(0.3), true), ...
                   P('r_SRX', 'scale', srx, log10(0.1), log10(15), 0, true), ...
                   P('r_kSE', 'scale', {'kSE'}, log10(0.01), log10(3), 0, true), ...
                   P('r_pas', 'scale', {'k_pas'}, log10(0.1), log10(3), 0, true), ...
                   P('SL0', 'abs', {'SL0'}, 1.9, 2.3, 2.0, false), tf];  cfg.maxIter = 50;
    case {'tB', 'tB_caE1', 'tB_caE2', 'tB_caE3', 'tB_caE4'}
        % case-B identifiability set (Ca transient treated as MEASURED): net factor on the six isoform
        % rates, net factor on k1,k_1, SRX and series factors, thin filament, offset. caE1-4: the Ca shape
        % perturbed by a realistic indicator error (TTP -/+ ~5 ms, peak->50 % -/+ ~10 ms).
        cfg.par = [P('r_net', 'scale', iso, log10(0.05), log10(20), log10(0.45), true), ...
                   P('r_Pi', 'scale', pi_, log10(0.05), log10(20), log10(2.08), true), ...
                   P('r_SRX', 'scale', srx, log10(0.1), log10(15), 0, true), ...
                   P('r_kSE', 'scale', {'kSE'}, log10(0.01), log10(3), 0, true), tf];
        cfg.maxIter = 12;      % profile points start warm; 12 iterations suffice (wave 1 ran with 25)
        switch name
            case 'tB_caE1', cfg.ca.opt.tau_r = 0.022;  cfg.ca.opt.t50 = 0.253;
            case 'tB_caE2', cfg.ca.opt.tau_r = 0.028;  cfg.ca.opt.t50 = 0.277;
            case 'tB_caE3', cfg.ca.opt.tau_r = 0.022;  cfg.ca.opt.t50 = 0.277;
            case 'tB_caE4', cfg.ca.opt.tau_r = 0.028;  cfg.ca.opt.t50 = 0.253;
        end
        cfg.profGrid = struct('r_net', log10([0.15 0.25 0.45 0.8 1.5 3 6 10]), ...
                              'r_SRX', log10([0.15 0.3 1 3 8]), 'r_kSE', log10([0.1 0.3 0.6 1 2]));
    case 't2_allc'    % the combined tier INSIDE the literature corridor: r_iso 0.05-0.15 (isoform), r_T 3.7-8.4 (Q10 3-3.5 over 12-17 C)
        cfg.par = [P('r_iso', 'scale', iso, log10(0.05), log10(0.15), log10(0.1), true), ...
                   P('r_T', 'scale', cyc, log10(3.7), log10(8.4), log10(4.5), true), ...
                   P('r_Pi', 'scale', pi_, log10(0.1), log10(1), log10(0.3), true), ...
                   P('r_SRX', 'scale', srx, log10(0.1), log10(15), 0, true), ...
                   P('r_kSE', 'scale', {'kSE'}, log10(0.01), log10(3), 0, true), ...
                   P('r_pas', 'scale', {'k_pas'}, log10(0.1), log10(3), 0, true), ...
                   P('SL0', 'abs', {'SL0'}, 1.9, 2.3, 2.0, false), tf];  cfg.maxIter = 50;
    case {'s3_free', 's3_free_b', 's3_free_c'}
        rate = @(f) P(f, 'scale', {f}, log10(0.01), log10(20), 0, true);
        cfg.par = [cellfun(rate, cyc), ...
                   P('ksr0', 'scale', {'ksr0'}, log10(0.05), log10(20), 0, true), ...
                   P('kmsrd', 'scale', {'kmsrd'}, log10(0.05), log10(20), 0, true), ...
                   P('ksrd', 'scale', {'ksrd'}, log10(0.05), log10(20), 0, true), ...
                   P('sigma2', 'scale', {'sigma2'}, log10(0.3), log10(3), 0, true), ...
                   P('sigma_srd1', 'scale', {'sigma_srd1'}, log10(0.3), log10(3), 0, true), ...
                   P('kstiff1', 'scale', {'kstiff1'}, log10(0.3), log10(3), 0, true), ...
                   P('kstiff2', 'scale', {'kstiff2'}, log10(0.3), log10(3), 0, true), ...
                   P('kSE', 'scale', {'kSE'}, log10(0.01), log10(3), 0, true), ...
                   P('k_pas', 'scale', {'k_pas'}, log10(0.1), log10(3), 0, true), ...
                   P('A0', 'scale', {'A0'}, log10(0.5), log10(4.5), 0, true), ...
                   P('dr', 'scale', {'dr'}, log10(0.5), log10(2), 0, true), ...
                   P('SL0', 'abs', {'SL0'}, 1.9, 2.3, 2.0, false), ...
                   P('tf_k_on', 'tf', {}, log10(0.01), log10(1e4), log10(7), true), ...
                   P('tf_k_off', 'tf', {}, log10(0.1), log10(1e4), log10(5), true), ...
                   P('tf_K_coop', 'tf', {}, 1, 100, 12, false), ...
                   P('shift', 'shift', {}, -0.150, 0.060, 0, false)];
        cfg.maxIter = 60;
    case {'tw1_isoT', 'tw1_isoT_xb', 'tw2_caH'}
        % T-World (ThirdParty/TWorld, CaInput 'tworld') as the Ca source instead of the de Lange
        % transient: the adult human endocardial 1 Hz beat (0.097 -> 0.42 uM, TTP 54 ms) taken as
        % granted ('signal' = its standalone-beat Ca; '_xb' = the 93 states embedded, the XB thin
        % filament buffering the troponin Ca). Same free set as t1_isoT: iso x T + thin filament + shift.
        cfg.CaInput = 'tworld';
        cfg.CaTworld = struct('coupling', 'signal', 'prebeats', 0, 'state', 'TW_endo');
        if strcmp(name, 'tw1_isoT_xb'), cfg.CaTworld.coupling = 'xb'; end
        cfg.fixedScale = fixed(pi_, 0.3);
        cfg.par = [hum, tf];  cfg.maxIter = 50;
        if strcmp(name, 'tw2_caH')
            % tier 2: T-World's Ca-handling levers freed (SERCA uptake = decay, ICaL = release trigger/amplitude,
            % NCX = extrusion), each pre-paced 20 beats towards its own steady state (~16 s per evaluation).
            cfg.CaTworld.prebeats = 20;
            cfg.par = [cfg.par, P('Jup_Multiplier',   'tw', {}, log10(0.1), log10(3), 0, true), ...
                                P('ICaL_Multiplier',  'tw', {}, log10(0.3), log10(3), 0, true), ...
                                P('INaCa_Multiplier', 'tw', {}, log10(0.3), log10(3), 0, true)];
            cfg.maxIter = 30;
        end
        % resolve the T-World layer once here, so every evaluation reuses it (tworldCaLayer caches on
        % params.ca.tw; a changed 'tw' parameter changes the key and re-paces, as it must)
        b2 = cfg.base;  b2.CaModel = 'twoState';  b2.CaInput = 'tworld';  b2.CaTworld = cfg.CaTworld;
        b2.CaPeriod = 1;  b2.CaFirstBeat = cfg.pre;
        b2 = getParams(b2, [], true, false);  cfg.base.ca = b2.ca;
    otherwise
        error('tierConfigs:name', 'unknown tier %s', name);
end
% warm start from an earlier tier's optimum where available
switch name
    case {'t2a_Pi', 't2b_SRX', 't2b_SRXc', 't2c_SE', 't2d_pas', 't2e_SL', 't2_all', 't2_allc', 't1_caFast', 't1_caSlow', 't1_caHiLvl'}
        cfg = warm(cfg, fullfile(res, 't1_isoT_fit.mat'));
    case {'tB', 'tB_caE1', 'tB_caE2', 'tB_caE3', 'tB_caE4'}
        cfg = warm(cfg, fullfile(res, 't1_isoT_fit.mat'));
    case {'tw1_isoT', 'tw1_isoT_xb'}
        cfg = warm(cfg, fullfile(res, 't1_isoT_fit.mat'));
    case 'tw2_caH'
        f = fullfile(res, 'tw1_isoT_fit.mat');  if ~exist(f, 'file'), f = fullfile(res, 't1_isoT_fit.mat'); end
        cfg = warm(cfg, f);
    case 's3_free'          % from the best stage-2 tier (t2_all if present, else t1), expanded to individual rates
        f = fullfile(res, 't2_all_fit.mat');  if ~exist(f, 'file'), f = fullfile(res, 't1_isoT_fit.mat'); end
        cfg = warmS3(cfg, f, iso, pi_, srx);
    case 's3_free_c'        % third start: the fast-cycle / slow-SRX regime found by t2_all
        cfg = warmS3(cfg, fullfile(res, 't2_all_fit.mat'), iso, pi_, srx);
    case 's3_free_b'        % second start: literature human_v0 ratios at Q10 factor 5
        cfg = warmS3(cfg, '', iso, pi_, srx);
end
if ~isempty(prof)                                % fix the profiled parameter at its grid value
    j = find(strcmp({cfg.par.name}, prof.par), 1);
    if isempty(j), error('tierConfigs:prof', 'no parameter %s in %s', prof.par, prof.base); end
    v = cfg.profGrid.(prof.par)(str2double(prof.k));  p = cfg.par(j);  cfg.par(j) = [];
    if p.log, vv = 10^v; else, vv = v; end
    switch p.type
        case 'scale'
            if ~isfield(cfg, 'fixedScale'), cfg.fixedScale = struct(); end
            for f = string(p.fields), if isfield(cfg.fixedScale, f), cfg.fixedScale.(f) = cfg.fixedScale.(f)*vv; else, cfg.fixedScale.(f) = vv; end, end
        case 'tf'
            cfg.tfFixed.(p.name) = vv;
        otherwise
            error('tierConfigs:prof', 'cannot profile %s', p.type);
    end
    cfg.profFixed = struct('par', p.name, 'value', v, 'k', str2double(prof.k));
end
end

function s = fixed(varargin)
s = struct();
for k = 1:2:numel(varargin), for f = varargin{k}, s.(f{1}) = varargin{k+1}; end, end
end

function cfg = warm(cfg, f)
if ~exist(f, 'file'), return; end
R0 = load(f);  R0 = R0.R;
for i = 1:numel(cfg.par)
    j = find(strcmp({R0.par.name}, cfg.par(i).name), 1);
    if ~isempty(j), cfg.par(i).x0 = min(max(R0.x(j), cfg.par(i).lb), cfg.par(i).ub); end
end
end

function cfg = warmS3(cfg, f, iso, pi_, srx)
% expand group multipliers (r_iso, r_T, r_Pi, r_SRX, r_kSE, r_pas) into the individual ones
g = struct('r_iso', log10(0.1), 'r_T', log10(5), 'r_Pi', log10(0.3), 'r_SRX', 0, 'r_kSE', 0, 'r_pas', 0, 'SL0', 2.0, ...
           'tf_k_on', log10(7), 'tf_k_off', log10(5), 'tf_K_coop', 12, 'shift', 0);
if ~isempty(f) && exist(f, 'file')
    R0 = load(f);  R0 = R0.R;
    for j = 1:numel(R0.par), g.(R0.par(j).name) = R0.x(j); end
end
for i = 1:numel(cfg.par)
    n = cfg.par(i).name;  v = [];
    if any(strcmp(n, iso)), v = g.r_iso + g.r_T;
    elseif any(strcmp(n, pi_)), v = g.r_Pi + g.r_T;
    elseif any(strcmp(n, srx)), v = g.r_SRX;
    elseif strcmp(n, 'kSE'), v = g.r_kSE;
    elseif strcmp(n, 'k_pas'), v = g.r_pas;
    elseif isfield(g, n), v = g.(n);
    end
    if ~isempty(v), cfg.par(i).x0 = min(max(v, cfg.par(i).lb), cfg.par(i).ub); end
end
end
