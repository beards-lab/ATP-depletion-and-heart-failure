% BuildRefit5Seed.m -- give the objective back the two things it cannot
% currently see: the shape of the force block, and the slope of ktr across the
% slack ladder.
%
% WHY (measured, Analyses/KtrSegmentSlope/RunSteadySpread.m):
%   The pooled target's between-preparation spread for every force-class
%   feature is dominated by how strong each HEART is, not by how well the
%   feature is determined:
%
%       feature        s used   s level   s shape
%       steady          0.172     0.233     0.021
%       A               0.179     0.264     0.021
%       peak1_y         0.175     0.237     0.016
%       vall_y          0.180     0.244     0.026
%       peak2           0.175     0.211     0.024
%
%   The three hearts' plateau forces are 80.6 / 68.9 / 50.6 kPa, but once each
%   is divided by its own mean they agree to 1-2 %. Because the floor rule is
%   fl = (1/3)*nEl*s^2, using the level spread makes every one of these
%   deadbands ~100x too wide, and refit4_opt sits inside all of them: steady by
%   3.7x while running a uniform +5.0 % high, peak2 by 2.0x with an +11 %
%   segment-4 error, vall_y by 1.0x at -9.7 %.
%
%   `ktr` is the same story with a different cause. Its ladder falls in all
%   three preparations (-4.9 %, -23.3 %, -15.2 %, pooled -14.8 %) and the model
%   RISES (+0.34 per segment against -1.85). That costs nothing: raw 0.021
%   against a 0.054 floor. Here the level spread is genuine biology -- hearts
%   really do differ in rate -- so the fix is not a narrower deadband on `ktr`
%   but a separate ratio feature, ktr_shape = ktr / ktr(1), exactly as
%   A_over_steady exists alongside A and steady.
%
% THE RULE APPLIED HERE (uniform, not per-feature tuning):
%   for a LEVEL-dominated feature, s := 2 * s_shape
%   then fl = (1/3)*nEl*s^2 and w = lambda/s^2 as before.
%   The factor 2 is the concession to n = 3: a shape spread from three hearts
%   is itself only known to about a factor of two.
%
% NOT touched: restretchSlopeStart and vall_t (their spread really is shape),
% rsK80, FV, the PS_* passive guards, and the `steady[1]` boundary term, which
% is the LEVEL guard and correctly keeps the wide 52.5-87.7 kPa heart range.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/BuildRefit5Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~exist('SEEDTAG','var'); SEEDTAG = 'refit5_seed'; end
LAMBDA  = 0.0106294;
SHAPE_K = 2;            % s := SHAPE_K * s_shape

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r5_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

p = getParams(loadParams('refit4_opt'), [], true, false);
p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);

fprintf('=== refit4_opt on its OWN objective, for the record ===\n');
c_before = evaluateBakersExp(1, p, true, false);
fprintf('  %.4f\n', c_before);

%% ---- measure the shape spreads straight from the three preparations -----
DAYS = {'03_27','04_03','04_10'};
FD = cell(1,3);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    FD{d} = S.features_data;
end
% shape spread per feature: divide each preparation by its own mean, then take
% the between-preparation CV element by element and average. No local function
% -- in a script those have to sit at the very end, and this needs FD.
LEVELDOM = {'steady','A','peak1_y','vall_y','peak2','vall2_dy'};
SS = containers.Map('KeyType','char','ValueType','double');
for q1 = 1:numel(LEVELDOM)
    nmq = LEVELDOM{q1};
    X = zeros(3, numel(FD{1}.(nmq)));
    for dd = 1:3; X(dd,:) = FD{dd}.(nmq)(:)'; end
    Xs = X ./ mean(X,2,'omitnan');
    SS(nmq) = mean(std(Xs,0,1,'omitnan') ./ abs(mean(Xs,1,'omitnan')), 'omitnan');
end
% ktr_shape is already level-free, so no per-preparation rescaling, and
% element 1 is identically 1 so it is excluded.
Xk = zeros(3,5);
for dd = 1:3; Xk(dd,:) = FD{dd}.ktr(:)'/FD{dd}.ktr(1); end
Xk = Xk(:,2:5);
SS('ktr_shape') = mean(std(Xk,0,1,'omitnan') ./ abs(mean(Xk,1,'omitnan')), 'omitnan');
fprintf('\n%-22s %5s %9s %9s %11s %11s %9s %9s\n', ...
        'term','nEl','s_old','s_new','fl_old','fl_new','w_old','w_new');
fn = p.fn;
for i = 1:numel(fn)
    tk = split(fn{i},'|');
    if contains(fn{i}, ','); continue; end            % grouped boundary entry
    base = regexprep(tk{1}, '\[.*$', '');
    if ~ismember(base, LEVELDOM); continue; end
    if any(startsWith(tk,'0-')) || any(~cellfun(@isempty, regexp(tk,'^[\d.]+-[\d.]+$'))); continue; end % boundary term
    nEl = 5;
    if contains(tk{1},'[')
        inner = extractBetween(tk{1},'[',']'); nEl = numel(eval(['[' inner{1} ']']));
    end
    jf = find(startsWith(tk,'fl:'),1);
    flold = NaN; if ~isempty(jf); flold = str2double(extractAfter(tk{jf},'fl:')); end
    wold  = str2double(tk{end});
    snew  = SHAPE_K * SS(base);
    sold  = sqrt(LAMBDA/wold);
    flnew = (1/3)*nEl*snew^2;
    wnew  = LAMBDA/snew^2;
    if isempty(jf)
        tk = [tk(1:end-1); {sprintf('fl:%.6g',flnew)}; tk(end)];
    else
        tk{jf} = sprintf('fl:%.6g', flnew);
    end
    tk{end} = sprintf('%.6g', wnew);
    fn{i} = strjoin(tk,'|');
    fprintf('%-22s %5d %9.4f %9.4f %11.5f %11.5f %9.4f %9.4f\n', ...
            base, nEl, sold, snew, flold, flnew, wold, wnew);
end

%% ---- the new ktr_shape term --------------------------------------------
s_ks  = SHAPE_K * SS('ktr_shape');
fl_ks = (1/3)*4*s_ks^2;
w_ks  = LAMBDA/s_ks^2;
fn{end+1} = sprintf('ktr_shape[2:5]|fl:%.6g|%.6g', fl_ks, w_ks);
fprintf('%-22s %5d %9s %9.4f %11s %11.5f %9s %9.4f   <- NEW\n', ...
        'ktr_shape[2:5]', 4, '-', s_ks, '-', fl_ks, '-', w_ks);
p.fn = fn;

fprintf('\nper-preparation ktr ladder (each / its own segment 1):\n');
for d = 1:3
    v = FD{d}.ktr(:)'/FD{d}.ktr(1);
    fprintf('  %-8s %7.4f %7.4f %7.4f %7.4f %7.4f\n', strrep(DAYS{d},'_','/'), v);
end

%% ---- score --------------------------------------------------------------
rep = checkFeatureCoverage(p);
if ~rep.ok
    fprintf(2,'UNSCORABLE:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildRefit5Seed:coverage','objective has unscorable terms');
end
fprintf('\n=== refit4_opt on the NEW objective ===\n');
[c_new,~,ct,nm] = evaluateBakersExp(1, p, true, false);
fprintf('  %s = %.4f   (was %.4f on the old one)\n', SEEDTAG, c_new, c_before);
[~,o] = sort(ct,'descend');
for i = o(1:min(16,numel(o)))
    if ct(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm{i}, ct(i), 100*ct(i)/c_new);
end
fprintf('\n  zero-cost terms: %s\n', strjoin(nm(ct<=1e-9), ', '));

%% ---- write --------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['refit4_opt params; force block re-floored on the SHAPE spread ' ...
    '(x%g) and ktr_shape[2:5] added. cost %.6f'], SHAPE_K, c_new));
fprintf('\nwritten: %s\n', out);
q = getParams(); s2 = loadParams(out);
f2 = fieldnames(s2); for i=1:numel(f2); q.(f2{i}) = s2.(f2{i}); end
nmf = fieldnames(p); worst = 0;
for i=1:numel(nmf)
    a = p.(nmf{i});
    if ~isnumeric(a)||~isscalar(a)||a==0||~isfinite(a); continue; end
    if ~isfield(q,nmf{i}); continue; end
    worst = max(worst, abs(q.(nmf{i})-a)/abs(a));
end
fprintf('round-trip worst relative error: %.3g\n', worst);
fprintf('pool: %d parameters\n', numel(p.leanPool));
fprintf('\nDONE\n');

function iWrite17g(fname, params, comment)
fid = fopen(fname, 'w');
fprintf(fid, '%% Generated file to manipulate the simulation parameters directly\r\n');
fprintf(fid, '%% %s\r\n', comment);
nm = fieldnames(params);
for i = 1:numel(nm)
    v = params.(nm{i});
    if isa(v,'function_handle') || isstruct(v); continue; end
    if isnumeric(v) && isscalar(v);       fprintf(fid, 'params0.%s = %.17g;\r\n', nm{i}, v);
    elseif islogical(v) && isscalar(v);   fprintf(fid, 'params0.%s = %d;\r\n', nm{i}, v);
    elseif ischar(v) && ~isempty(v);      fprintf(fid, 'params0.%s = ''%s'';\r\n', nm{i}, v);
    elseif iscell(v) && ~isempty(v);      fprintf(fid, 'params0.%s = {''%s''};\r\n', nm{i}, strjoin(v, ''', '''));
    elseif isnumeric(v) && ~isempty(v)
        fprintf(fid, 'params0.%s = [%s];\r\n', nm{i}, ...
                strjoin(arrayfun(@(x) sprintf('%.17g',x), v(:)', 'uni', 0), ' '));
    end
end
fprintf(fid, 'params0.mods = {};\r\nparams0.g = [];\r\n');
fclose(fid);
end
