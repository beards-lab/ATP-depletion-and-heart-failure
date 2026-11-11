% ReportJointRefit.m -- score the joint refit against its seed, like-for-like.
%
% Everything is scored on ONE objective (the joint 8+2 mM list carried by the
% seed), so the comparison is a comparison and not a change of yardstick. Prints:
%   1. joint cost, split into its 8 mM and 2 mM legs
%   2. per-feature deltas, both legs, sorted by what moved
%   3. the ATP effect (2 mM / 8 mM ratios) against data
%   4. which parameters actually moved, and by how much
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/ReportJointRefit.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end

if ~exist('TAG','var'); TAG = 'joint82'; end
SEEDSNAP = 'jointrefit_seed';
OPTSNAP  = [TAG '_opt'];

if ~isfile(fullfile(root,'params',[OPTSNAP '.m']))
    error('ReportJointRefit:noOpt','%s.m not found -- has the refit finished?', OPTSNAP);
end

%% ---- score both, on the SEED's objective --------------------------------
seed = getParams(loadParams(SEEDSNAP), [], true, false);
opt  = getParams(loadParams(OPTSNAP),  [], true, false);
% the objective belongs to the seed; force the optimum onto it verbatim
opt.fn = seed.fn; opt.AtpJoint_fn2 = seed.AtpJoint_fn2;
opt.UseAtpJointFit = true;
opt.AtpJoint_w8 = seed.AtpJoint_w8; opt.AtpJoint_w2 = seed.AtpJoint_w2;
opt.AtpJoint_MgATP2 = seed.AtpJoint_MgATP2; opt.AtpJoint_Pi2 = seed.AtpJoint_Pi2;
opt.AtpJoint_vt2 = seed.AtpJoint_vt2;
opt.AtpJoint_RunFV2 = false; opt.AtpJoint_RunPS2 = false;
seed.mods = {}; seed.g = []; opt.mods = {}; opt.g = [];

fprintf('\nscoring seed ...\n');
[Es, ~, fcs, nms] = evaluateBakersExp(1, seed, true, false);
fprintf('scoring refit ...\n');
[Eo, ~, fco, nmo] = evaluateBakersExp(1, opt,  true, false);
assert(isequal(nms, nmo), 'feature lists differ -- not a like-for-like comparison');

h = startsWith(nms,'hi:'); l = startsWith(nms,'lo:');
fprintf('\n================= JOINT COST =================\n');
fprintf('%-10s %10s %10s %10s\n','', 'joint', '8 mM', '2 mM');
fprintf('%-10s %10.4f %10.4f %10.4f\n','seed',  Es, sum(fcs(h)), sum(fcs(l)));
fprintf('%-10s %10.4f %10.4f %10.4f\n','refit', Eo, sum(fco(h)), sum(fco(l)));
fprintf('%-10s %+10.4f %+10.4f %+10.4f\n','delta', Eo-Es, sum(fco(h))-sum(fcs(h)), sum(fco(l))-sum(fcs(l)));

%% ---- per feature ---------------------------------------------------------
d = fco(:) - fcs(:);
[~,ix] = sort(d);
fprintf('\n============ PER-FEATURE (sorted by delta) ============\n');
fprintf('%-52s %9s %9s %9s\n','feature','seed','refit','delta');
for k = ix(:)'
    if abs(d(k)) < 5e-4 && fcs(k) < 5e-4; continue; end
    nm = nms{k}; if numel(nm) > 52; nm = nm(1:52); end
    fprintf('%-52s %9.3f %9.3f %+9.3f\n', nm, fcs(k), fco(k), d(k));
end

%% ---- the ATP effect ------------------------------------------------------
fprintf('\n============ ATP EFFECT (2 mM / 8 mM) ============\n');
R = struct();
for which = {'seed','refit'}
    p = seed; if strcmp(which{1},'refit'); p = opt; end
    q = p; q.UseAtpJointFit = false;
    q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false;
    q.RunStairs=false; q.RunSlackPassive=false; q.RunSlack=true;
    q.RunSlackSegments='All'; q.MaxRunTime=600;
    q8 = q; q8.MgATP=8; q8.Pi=0; q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
    q2 = q; q2.MgATP=2; q2.Pi=0.8; q2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
    [~,~,m8,d8] = runSlackExperiment(q8);
    [~,~,m2,d2] = runSlackExperiment(q2);
    R.(which{1}) = struct('m8',m8,'m2',m2,'d8',d8,'d2',d2);
end
keys = {'steady','A','peak1_y','vall_y','peak2','ktr','rsK'};
fprintf('%-9s %8s | %8s %7s | %8s %7s\n','feature','data', 'seedM','m/d','refitM','m/d');
for k = 1:numel(keys)
    f = keys{k};
    rd = iRat(R.seed.d2, R.seed.d8, f);
    rs = iRat(R.seed.m2,  R.seed.m8,  f);
    ro = iRat(R.refit.m2, R.refit.m8, f);
    fprintf('%-9s %8.3f | %8.3f %7.3f | %8.3f %7.3f\n', f, rd, rs, rs/rd, ro, ro/rd);
end
fprintf('\nabsolute rsK (s^-1), mean over segments:\n');
fprintf('  data   8 mM %6.1f   2 mM %6.1f\n', iM(R.seed.d8,'rsK'),  iM(R.seed.d2,'rsK'));
fprintf('  seed   8 mM %6.1f   2 mM %6.1f\n', iM(R.seed.m8,'rsK'),  iM(R.seed.m2,'rsK'));
fprintf('  refit  8 mM %6.1f   2 mM %6.1f\n', iM(R.refit.m8,'rsK'), iM(R.refit.m2,'rsK'));

%% ---- what moved ----------------------------------------------------------
fprintf('\n============ PARAMETERS THAT MOVED (>1%%) ============\n');
fn = fieldnames(seed); n = 0;
rows = {};
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(opt,f); continue; end
    a = seed.(f); b = opt.(f);
    if ~isnumeric(a) || ~isnumeric(b) || ~isequal(size(a),size(b)) || isempty(a); continue; end
    if any(~isfinite(a(:))) || any(~isfinite(b(:))); continue; end
    den = max(abs(a(:)), eps);
    rel = max(abs(b(:)-a(:))./den);
    if rel > 0.01
        n = n + 1;
        rows(end+1,:) = {f, a(1), b(1), 100*rel}; %#ok<SAGROW>
    end
end
if isempty(rows)
    fprintf('  (none)\n');
else
    [~,o] = sort(-cell2mat(rows(:,4)));
    fprintf('%-34s %12s %12s %9s\n','parameter','seed','refit','max %');
    for k = o(:)'
        fprintf('%-34s %12.5g %12.5g %8.1f%%\n', rows{k,1}, rows{k,2}, rows{k,3}, rows{k,4});
    end
end
fprintf('\n%d parameters moved by >1%%\n', n);

save(fullfile(here,'results','jointrefit_report.mat'), 'Es','Eo','fcs','fco','nms','R','rows');
fprintf('\nREPORT DONE\n');

%% ========================================================================
function r = iRat(a, b, f)
if isfield(a,f) && isfield(b,f)
    r = mean(a.(f)(:),'omitnan') / mean(b.(f)(:),'omitnan');
else
    r = NaN;
end
end

function v = iM(s, f)
if isfield(s,f); v = mean(s.(f)(:),'omitnan'); else; v = NaN; end
end
