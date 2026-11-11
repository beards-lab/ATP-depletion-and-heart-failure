% BuildLeanSeed.m -- the ablation campaign's verdict, turned into a seed.
%
% Three changes, all of them things the campaign or the degeneracy analysis
% earned. Written with %.17g, NOT writeParamsToMFile's %g (see §5 below).
%
% 1. MECHANISMS DELETED (free, per Analyses/MechanismAblation/conclusions.md)
%      UseA2Reattaching   -> 0   : 3 params for +0.0001
%      UseMaxwellTensionOnly -> 0: 1 param  for +0.037, rsK80 0.06 BETTER
%      UseDynamicPassive  -> 0   : exactly inert (c_titin_visc = 0)
%    UseSerialStiffness is left alone: it is dead code, so setting it changes
%    nothing. Deleting the FLAG is a source edit, not a seed edit.
%
% 2. OBJECTIVE
%      coolDownLS REMOVED. It is an unbounded raw MSE on the cool-down trace
%      (extractSlackAttributes.m:456), not a chi-2 term, and it dominated the
%      movers of nearly every ablation arm -- a guard term masquerading as a
%      data term. rsK80 already scores post-re-stretch recovery on a proper
%      floored chi-2, so it inherits the role, at a raised weight.
%      FV share 0.40 -> FV_SHARE (default 0.25). The pin has been arbitrary
%      since it was introduced; the campaign gave no reason to keep it high.
%
% 3. POOL (written into params0.leanPool for RunLeanOptim)
%      Dropped as INERT   : sigma1, sigma_srd2, LXBpivot
%      Dropped as DELETED : kA2re, kA2hop, sA2hop, x_M_slack, c_titin_visc
%      Dropped as REDUNDANT (degeneracy analysis, conclusions §7):
%        d10_ref, SL_ref_lattice, R_thick, R_thin  -- these four and d_optimal
%          enter f_lattice only through two ratios. Fixed at literature values;
%          d_optimal and sigma_lattice carry the whole mechanism.
%        L_hbare -- in the operating regime it appears only in the denominator
%          (L_thick - L_hbare), a pure gain that ka already spans.
%
% Cost is NOT comparable to 0.9081 afterwards: the objective changed. The
% printout reports the seed under both objectives so the change is auditable.
%
%   matlab -batch "run('Analyses/MechanismAblation/BuildLeanSeed.m')"

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('BASE','var');      BASE      = 'floored_seed'; end
if ~exist('SEEDTAG','var');   SEEDTAG   = 'lean_seed';    end
if ~exist('FV_SHARE','var');  FV_SHARE  = 0.25;           end
if ~exist('RSK_BOOST','var'); RSK_BOOST = 2.0;            end  % rsK80 inherits coolDownLS's role
if ~exist('NWORKERS','var');  NWORKERS  = 5;              end

p = getParams(); s = loadParams(BASE);
f = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunKtr=false; p.RunSlack=true; p.RunStairs=false;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.mods = {}; p.g = ones(1,0);

if isempty(gcp('nocreate')); parpool('local', NWORKERS); end

fprintf('=== seed under the OLD objective ===\n');
[c_old, ~, ct_old, nm_old] = evaluateBakersExp(1, p, true, false);
fprintf('  %s = %.4f\n', BASE, c_old);
iShow(ct_old, nm_old, c_old, 8);

%% ---- 1. mechanism deletions --------------------------------------------
p.UseA2Reattaching = 0;      p.kA2re = 0; p.kA2hop = 0;
p.UseMaxwellTensionOnly = 0; p.x_M_slack = 0;
p.UseDynamicPassive = 0;     p.c_titin_visc = 0;

%% ---- 2. objective surgery ----------------------------------------------
fn = p.fn;

% (a) drop coolDownLS
keep = ~strcmp(iPrim(fn), 'coolDownLS');
if any(~keep); fprintf('\ndropped fn entry: %s\n', strjoin(fn(~keep), ', ')); end
fn = fn(keep);

% (b) rsK80 inherits the role
iR = find(strcmp(iPrim(fn), 'rsK80'), 1);
if ~isempty(iR)
    fn{iR} = iSetW(fn{iR}, iGetW(fn{iR}) * RSK_BOOST);
    fprintf('rsK80 weight x%.3g -> %s\n', RSK_BOOST, fn{iR});
end

% (c) FV share. Contribution is linear in the weight, so to move FV from its
% current share to FV_SHARE, scale its weight by the ratio the two shares imply.
iF = find(strcmp(iPrim(fn), 'FV_normpowerAvg'), 1);
ctFV = 0; iOld = find(strcmp(nm_old,'FV_normpowerAvg'),1);
if ~isempty(iOld); ctFV = ct_old(iOld); end
if ~isempty(iF) && ctFV > 0
    rest  = c_old - ctFV;                       % everything that is not FV
    ctNew = FV_SHARE/(1-FV_SHARE) * rest;       % what FV must contribute
    scale = ctNew / ctFV;
    fprintf('FV share %.3f -> %.3f  (weight x%.4g)\n', ctFV/c_old, FV_SHARE, scale);
    fn{iF} = iSetW(fn{iF}, iGetW(fn{iF}) * scale);
end
p.fn = fn; p.mods = {};

%% ---- 3. the lean pool ---------------------------------------------------
% Literature values for the four lattice constants that are observations, not
% fit parameters. Rat/mouse cardiac, X-ray diffraction and EM:
%   d10 ~ 37-42 nm at SL 2.0-2.2; thick filament radius ~7-8 nm; thin ~4-5 nm.
% The incumbent's values are already inside those ranges, so FIXING them here
% changes nothing numerically -- it only stops the optimiser spending draws on
% four numbers that enter the model through two ratios.
p.d10_ref        = 0.0365153;   % um, X-ray d10 at the reference SL   [FIXED]
p.SL_ref_lattice = 2.00579;     % um, the SL that d10 was measured at [FIXED]
p.R_thick        = 0.00825173;  % um, thick filament radius           [FIXED]
p.R_thin         = 0.00389595;  % um, thin filament radius            [FIXED]

p.leanPool = { ...
    'ka','kd','k1','k_1','k2','k_2','kah','kamh', ...
    'kstiff1','kstiff2','kstiff1_n','kstiff2_n','kSE','estiff','ekSE', ...
    'kSE_M','eta_M','mu_neg','mu', ...
    'dr','dr1','dr2', ...
    'ksr0','kmsrd','ksr2srd','ksrd2sr','sigma2','sigma_srd1', ...
    'A0','v_ref_reg','tau_reg','P_bound_max', ...
    'k_cr','tau_cr','F_cr','SRXFromR2HighStrain','sSRXrip', ...
    'L_thick','L_thin', ...
    'd_optimal','sigma_lattice', ...
    'Lsc0','gamma','k_pas', ...
    'PieceWiseStrainDep2X_logOffset','PieceWiseStrainDepR1DX_logOffset', ...
    'PieceWiseStrainDepR21X_logOffset','PieceWiseStrainDepX_logOffset', ...
    'PieceWiseStrainDepX__2','PieceWiseStrainDepX__3','PieceWiseStrainDepX__4', ...
    'PieceWiseStrainDep2X__2','PieceWiseStrainDep2X__3','PieceWiseStrainDep2X__4','PieceWiseStrainDep2X__5', ...
    'PieceWiseStrainDepR1DX__1','PieceWiseStrainDepR1DX__2','PieceWiseStrainDepR1DX__3','PieceWiseStrainDepR1DX__4', ...
    'PieceWiseStrainDepR21X__2','PieceWiseStrainDepR21X__3','PieceWiseStrainDepR21X__4', ...
    'PieceWiseStrainDepParams__2','PieceWiseStrainDepParams__3','PieceWiseStrainDepParams__4', ...
    'PieceWiseStrainDep2Params__2','PieceWiseStrainDep2Params__3','PieceWiseStrainDep2Params__4', ...
    'PieceWiseStrainDep2Params__5','PieceWiseStrainDep2Params__6','PieceWiseStrainDep2Params__7', ...
    'PieceWiseStrainDepR1DParams__1','PieceWiseStrainDepR1DParams__2', ...
    'PieceWiseStrainDepR1DParams__3','PieceWiseStrainDepR1DParams__4', ...
    'PieceWiseStrainDepR21Params__2','PieceWiseStrainDepR21Params__3','PieceWiseStrainDepR21Params__4'};

%% ---- 4. the seed under the new objective -------------------------------
fprintf('\n=== seed under the NEW objective ===\n');
[c_new, ~, ct_new, nm_new] = evaluateBakersExp(1, p, true, false);
fprintf('  %s = %.4f   (was %.4f under the old one -- NOT comparable)\n', SEEDTAG, c_new, c_old);
iShow(ct_new, nm_new, c_new, 8);

%% ---- 5. write, at full precision ---------------------------------------
% writeParamsToMFile serialises scalars with '%g' -- SIX significant digits. A
% snapshot round-trip therefore perturbs every parameter by up to ~1e-6
% relative. That is not why the optimiser's own guard fires (it re-evaluates in
% memory), but it does mean a snapshot is not the point that was fitted. Write
% %.17g, which round-trips a double exactly.
out = fullfile(root, 'params', [SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['%s: A2Reattaching/MaxwellTensionOnly/DynamicPassive deleted; ' ...
    'coolDownLS dropped; rsK80 x%.3g; FV share %.2f; lattice constants fixed. ' ...
    'cost %.6f (old objective %.6f)'], BASE, RSK_BOOST, FV_SHARE, c_new, c_old));
fprintf('\nwritten: %s\n', out);

% Prove the round-trip is now exact.
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i = 1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
nm = fieldnames(p); worst = 0; wname = '';
for i = 1:numel(nm)
    a = p.(nm{i});
    if ~isnumeric(a) || ~isscalar(a) || a == 0 || ~isfinite(a); continue; end
    if ~isfield(q, nm{i}); continue; end
    r = abs(q.(nm{i}) - a)/abs(a);
    if r > worst; worst = r; wname = nm{i}; end
end
fprintf('round-trip worst relative error: %.3g (%s)\n', worst, wname);

%% ---- helpers ------------------------------------------------------------
function pr = iPrim(fn)
pr = regexprep(fn, '[|\[].*$', '');
end

function w = iGetW(entry)
tk = strsplit(entry, '|'); w = str2double(tk{end});
end

function e = iSetW(entry, w)
tk = strsplit(entry, '|'); tk{end} = sprintf('%.6g', w); e = strjoin(tk, '|');
end

function iShow(ct, nm, c, n)
[~, o] = sort(ct, 'descend');
for i = o(1:min(n,numel(o)))
    fprintf('      %-22s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c);
end
end

function iWrite17g(fname, params, comment)
fid = fopen(fname, 'w');
fprintf(fid, '%% Generated file to manipulate the simulation parameters directly\r\n');
fprintf(fid, '%% %s\r\n', comment);
nm = fieldnames(params);
for i = 1:numel(nm)
    v = params.(nm{i});
    if isa(v,'function_handle') || isstruct(v); continue; end
    if isnumeric(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %.17g;\r\n', nm{i}, v);
    elseif islogical(v) && isscalar(v)
        fprintf(fid, 'params0.%s = %d;\r\n', nm{i}, v);
    elseif ischar(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = ''%s'';\r\n', nm{i}, v);
    elseif iscell(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = {''%s''};\r\n', nm{i}, strjoin(v, ''', '''));
    elseif isnumeric(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = [%s];\r\n', nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g',x), v(:)', 'uni', 0), ' '));
    end
end
fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n');
fclose(fid);
end
