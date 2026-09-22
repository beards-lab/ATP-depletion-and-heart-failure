% BuildRefit3Seed.m -- third refit, per fjezek 2026-09-21.
%
% RESTART POINT: params/ktrfence_opt.m.
%   refit2 started from ktrfence_seed and STALLED -- 9 rounds, 0 improvements,
%   4 kicks, best never left the seed value. Diagnosis: the direction that
%   earned ktrfence its -0.166 is the same one that broke A/steady, so the new
%   ratio term taxes it 0.087 and cuts the net payoff to 0.036, a 4.6x
%   reduction; a gradient that shallow is not findable with 30 simplex
%   evaluations on a 12-parameter draw, and the +-10 % kick returned 1.1681, so
%   the basin is narrow as well as flat. ktrfence_opt already scores 0.7281 on
%   that objective -- better than 2.2 h of searching reached -- AND carries an
%   explicit 0.087 debt on A_over_steady, i.e. a real gradient to descend.
%
% PASSIVE: FLOORED, NOT DOWN-WEIGHTED (fjezek: "keep it as a floored objective,
% variation within some bounds is okay").
%   The PS_* block is ~31 % of the objective and essentially INVARIANT --
%   0.2749 / 0.2374 / 0.2337 across fastgate_best / ktrfence_seed /
%   ktrfence_opt. A third of the budget that cannot move is a third wasted.
%   Cutting the weight would make it cheap to SELL; a deadband instead makes it
%   impossible to BUY and still expensive to sell, which is what it is for.
%
%   Each PS_* term gets floor = its own value at this seed x PASSIVE_HEADROOM.
%   That is frankly a freeze-where-it-is choice, not a measurement -- there is
%   no noise model for these terms, and two of them (PS_rampupRMSE,
%   PS_holdDecayRMSE) are unbounded raw RMSE against a target of zero, so no
%   chi-square reading exists to derive one from. It is defensible because the
%   passive protocol ran after four activations and carries uncorrected rundown
%   (Analyses/PassiveRundownATP: the active runs should carry ~x1.2 more passive
%   than the model has), so its absolute level is not a number to chase.
%   Record it as a judgement call.
%
%   matlab -batch "run('Analyses/RefitStrategy/BuildRefit3Seed.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
if ~exist('SEEDTAG','var');          SEEDTAG = 'refit3_seed'; end
if ~exist('PASSIVE_HEADROOM','var'); PASSIVE_HEADROOM = 1.20; end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_r3_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 5); end

% parameters from ktrfence_opt, objective from refit2_seed (A_over_steady term,
% narrowed floors, ktr fence retargeted on the pooled mean)
p  = getParams(); s = loadParams('ktrfence_opt');
f  = fieldnames(s); for i = 1:numel(f); p.(f{i}) = s.(f{i}); end
r2 = loadParams('refit2_seed');
p.fn = r2.fn; p.leanPool = r2.leanPool;
p.RunForceVelocity=true; p.RunSlack=true; p.RunStairs=false; p.RunKtr=true;
p.RunForceVelocityTime=false; p.RunSlackPassive=true; p.EvalFeatures=true;
p.BreakOnODEUnstable=false; p.PlotEachSeparately=0; p.PlotFeatureFitting=0;
p.RunSlackSegments='AllPar'; p.MaxRunTime=150; p.OptimizeOn='Feats';
p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0; p.mods={}; p.g=ones(1,0);

PASBLK = {'PS_restretchPeak','PS_rampupRMSE','PS_holdDecayRMSE','PS_steady22','PS_steady20'};

fprintf('=== ktrfence_opt on the refit2 objective, BEFORE flooring passive ===\n');
[c0,~,ct0,nm0] = evaluateBakersExp(1, p, true, false);
isP = ismember(nm0, PASBLK);
fprintf('  total %.4f   passive %.4f  (%.1f %%)\n', c0, sum(ct0(isP)), 100*sum(ct0(isP))/c0);

% raw (unweighted, unfloored) chi-2 for each PS_* term, so the floor can be set
% in the same units evalFeatureCost subtracts in
FN0 = regexprep(p.fn, '\|fl:[0-9eE.+-]+', '');
q = p; q.fn = FN0; [~,~,fcnf,nmnf] = evaluateBakersExp(1, q, true, false);
wOf = @(e) iWeight(e);
fprintf('\n  %-22s %10s %10s %12s\n','term','weighted','raw chi2','new floor');
fn = p.fn;
for i = 1:numel(fn)
    tk = split(fn{i},'|');
    base = regexprep(tk{1}, '\[.*$', '');
    if ~ismember(base, PASBLK); continue; end
    j = find(strcmp(nmnf, base), 1);
    if isempty(j); warning('no raw cost for %s', base); continue; end
    raw   = fcnf(j) / wOf(fn{i});
    newfl = raw * PASSIVE_HEADROOM;
    jf = find(startsWith(tk,'fl:'),1);
    if isempty(jf)
        tk = [tk(1:end-1); {sprintf('fl:%.6g',newfl)}; tk(end)];   % before the weight
    else
        tk{jf} = sprintf('fl:%.6g', newfl);
    end
    fn{i} = strjoin(tk,'|');
    fprintf('  %-22s %10.4f %10.4f %12.4g\n', base, fcnf(j), raw, newfl);
end
p.fn = fn;

%% ---- score --------------------------------------------------------------
rep = checkFeatureCoverage(p);
if ~rep.ok
    fprintf(2,'UNSCORABLE:\n');
    if ~isempty(rep.missingModel); fprintf(2,'  model: %s\n', strjoin(rep.missingModel,', ')); end
    if ~isempty(rep.missingData);  fprintf(2,'  data : %s\n', strjoin(rep.missingData ,', ')); end
    error('BuildRefit3Seed:coverage','objective has unscorable terms');
end
fprintf('\n=== seed on the FINAL objective ===\n');
[c1,~,ct1,nm1] = evaluateBakersExp(1, p, true, false);
isP1 = ismember(nm1, PASBLK);
fprintf('  %s = %.4f   passive now %.4f (%.1f %%), was %.4f (%.1f %%)\n', SEEDTAG, c1, ...
        sum(ct1(isP1)), 100*sum(ct1(isP1))/c1, sum(ct0(isP)), 100*sum(ct0(isP))/c0);
[~,o] = sort(ct1,'descend');
for i = o(1:min(12,numel(o)))
    if ct1(i) <= 0; continue; end
    fprintf('    %-24s %8.4f  (%4.1f %%)\n', nm1{i}, ct1(i), 100*ct1(i)/c1);
end
for q2 = {'A_over_steady','ktrP_k63','FV_normpowerAvg','A','steady'}
    i = find(strcmp(nm1,q2{1}),1);
    if ~isempty(i); fprintf('  %-22s %8.4f\n', q2{1}, ct1(i)); end
end

%% ---- write --------------------------------------------------------------
out = fullfile(root,'params',[SEEDTAG '.m']);
iWrite17g(out, p, sprintf(['ktrfence_opt params + refit2 objective (A_over_steady, narrowed floors, ' ...
    'pooled ktr fence) + PS_* FLOORED at seed x%.2f (passive %.1f%%%% -> %.1f%%%%). cost %.6f'], ...
    PASSIVE_HEADROOM, 100*sum(ct0(isP))/c0, 100*sum(ct1(isP1))/c1, c1));
fprintf('\nwritten: %s\n', out);
q3 = getParams(); s3 = loadParams(out);
f3 = fieldnames(s3); for i=1:numel(f3); q3.(f3{i}) = s3.(f3{i}); end
nmf = fieldnames(p); worst = 0;
for i=1:numel(nmf)
    a = p.(nmf{i});
    if ~isnumeric(a)||~isscalar(a)||a==0||~isfinite(a); continue; end
    if ~isfield(q3,nmf{i}); continue; end
    worst = max(worst, abs(q3.(nmf{i})-a)/abs(a));
end
fprintf('round-trip worst relative error: %.3g\n', worst);
fprintf('\nDONE\n');

function w = iWeight(entry)
% the weight is the last purely-numeric token of an fn entry (1 if none)
tk = split(entry,'|'); w = 1;
for j = 2:numel(tk)
    if startsWith(tk{j},'fl:') || startsWith(tk{j},'w:'); continue; end
    v = str2double(tk{j});
    if ~isnan(v); w = v; end
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
