% EvalLowATP.m -- how good is the low-ATP fit on the NEW (realign4) 8 mM baseline?
%
% Five slack simulations, each scored under BOTH feature lists, so that:
%   * old vs new baseline is like-for-like (same list), and
%   * the new number can still be tied to the incumbent's recorded 15.98 / floor 4.69.
%
% The two lists differ (the objective itself changed between the two baselines):
%   LIST_OLD  incumbent's 19 features: has ovrsht_dy, rsK weighted 0.1
%   LIST_NEW  realign4's list, slack subset: ovrsht_dy dropped (retired by
%             SlackDataAnalysis), rsK weighted 0.25
% Neither list can contain FV* or PS_* -- there is no 2 mM force-velocity or
% passive data, so those are dropped from every 2 mM score by necessity.
%
% KNOWN CONSTANT OFFSET: the composite bound term
% 'k2|100-1500|0.01,ksrd|0.-20|0.001,kmsrd|0.0-20|0.001' contributes exactly
% 1.200 to EVERY score below -- a slack-only run does not produce k2/ksrd/kmsrd,
% so each sub-entry returns MISSING_FEATURE_COST = 100 times its weight
% (100*0.01 + 100*0.001 + 100*0.001). It is not a bound violation and carries no
% information; it cancels in every comparison but inflates every absolute number.

here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root)); LoadData;
refreshPool(5);

S = load(fullfile(root,'Analyses','LowATP_AbsoluteFit','results','powergate2_state.mat'));

VT8 = 'protocol_03_27_2026_8mM_slack.mat';
VT2 = 'protocol_03_27_2026_2mM_slack.mat';

%% ---- the two feature lists -----------------------------------------------
LIST_OLD = iSlackOnly(S.FN);
q4       = getParams(loadParams('realign4_opt'), [], true, false);
LIST_NEW = iSlackOnly(q4.fn);
fprintf('LIST_OLD (%d): %s\n\n', numel(LIST_OLD), strjoin(LIST_OLD,' | '));
fprintf('LIST_NEW (%d): %s\n\n', numel(LIST_NEW), strjoin(LIST_NEW,' | '));

%% ---- the five configurations ---------------------------------------------
% NEW baseline = realign4 (Further-optim best), gates transferred verbatim.
% OLD baseline = rskR2_w025_opt, the baseline the gates were actually fitted on.
% Both OLD configs get the optimiser's full-precision exponents, so the OLD number
% here is the honest incumbent, not the %.6g snapshot.
C = struct('tag',{},'p',{});

n2 = iCfg(getParams(loadParams('lowATP_realign4_gates'), [], true, false), 2, VT2);
C(end+1) = struct('tag','NEW 2 mM  (realign4 + gates)', 'p', n2);

n8 = iCfg(getParams(loadParams('realign4_opt'), [], true, false), 8, VT8);
C(end+1) = struct('tag','NEW 8 mM  (realign4, floor)',  'p', n8);

nz = n2; nz.UseAtpPowerGates = 0;
C(end+1) = struct('tag','NEW 2 mM  (gates OFF)',        'p', nz);

o2 = iCfg(iExact(getParams(loadParams('rskR2_w025_opt'), [], true, false), S), 2, VT2);
C(end+1) = struct('tag','OLD 2 mM  (rskR2 + gates)',    'p', o2);

o8 = iCfg(getParams(loadParams('rskR2_w025_opt'), [], true, false), 8, VT8);
C(end+1) = struct('tag','OLD 8 mM  (rskR2, floor)',     'p', o8);

%% ---- run ------------------------------------------------------------------
R = struct('tag',{},'cOld',{},'cNew',{},'fm',{},'fd',{},'vOld',{},'vNew',{});
for i = 1:numel(C)
    t0 = tic;
    fprintf('--- %-30s ', C(i).tag);
    try
        [~, ~, fm, fd] = runSlackExperiment(C(i).p);
        vO = evalFeatureCost(fd, fm, LIST_OLD, 2);
        vN = evalFeatureCost(fd, fm, LIST_NEW, 2);
        R(end+1) = struct('tag', C(i).tag, 'cOld', sum(vO), 'cNew', sum(vN), ...
                          'fm', fm, 'fd', fd, 'vOld', vO, 'vNew', vN); %#ok<SAGROW>
        fprintf('OLD-list %8.4f   NEW-list %8.4f   (%.0f s)\n', sum(vO), sum(vN), toc(t0));
    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

%% ---- headline table -------------------------------------------------------
fprintf('\n\n================ 2 mM FEATURE COST ================\n');
fprintf('%-32s %12s %12s\n', 'configuration', 'OLD list', 'NEW list');
for i = 1:numel(R)
    fprintf('%-32s %12.4f %12.4f\n', R(i).tag, R(i).cOld, R(i).cNew);
end
fprintf('\nfor reference, recorded at the time of the fit:\n');
fprintf('  incumbent 2 mM best.c %.4f, 8 mM floor %.4f (OLD list)\n', S.best.c, S.FLOOR);
fprintf('  realign4 8 mM full objective (FV+slack+passive) L2 = 3.8513\n');

%% ---- per-feature, new vs old ---------------------------------------------
iPerFeature(R, LIST_NEW, 'NEW');

%% ---- the physiology: 2 mM / 8 mM ratios vs data ---------------------------
fprintf('\n\n============ ATP EFFECT (2 mM / 8 mM), model vs data ============\n');
iRatios(R, 'NEW 2 mM  (realign4 + gates)', 'NEW 8 mM  (realign4, floor)', 'NEW baseline');
iRatios(R, 'OLD 2 mM  (rskR2 + gates)',    'OLD 8 mM  (rskR2, floor)',    'OLD baseline');

save(fullfile(root,'Analyses','LowATP_Realign4Transfer','results','realign4_transfer.mat'), ...
     'R','LIST_OLD','LIST_NEW');
fprintf('\nEVAL DONE\n');

%% ===========================================================================
function fn = iSlackOnly(fn)
% Drop FV*/PS_* (no 2 mM data for either) and de-duplicate, preserving order.
[~, iu] = unique(fn, 'stable'); fn = fn(iu);
prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), fn, 'uni', 0);
fn = fn(~(startsWith(prim,'FV') | startsWith(prim,'PS_')));
end

function q = iCfg(q, atp, vt)
% Slack-only evaluation at one [ATP]. Mirrors RunIncumbentAudit's iCfg so the
% numbers stay comparable with the incumbent campaign.
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false;
q.RunStairs=false; q.RunSlackPassive=false; q.RunSlack=true;
q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=600;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0;
q.UseR2Ceiling=0; q.UseD0State=0;
q.MgATP=atp; q.MgADP=0;
if atp == 8, q.Pi = 0; else, q.Pi = 0.8; end
q.velocitytableonfile = vt;
end

function q = iExact(q, S)
% Full-precision gate exponents from the optimiser state (not the %.6g snapshot).
q.UseAtpPowerGates = 1; q.MgATP_ref = 8;
for i = 1:numel(S.P); q.(['pAtp_' S.P{i}]) = S.best.x(i); end
end

function iPerFeature(R, LIST, which)
if numel(R) < 4, return; end
iN = find(strcmp({R.tag}, 'NEW 2 mM  (realign4 + gates)'), 1);
iF = find(strcmp({R.tag}, 'NEW 8 mM  (realign4, floor)'), 1);
iO = find(strcmp({R.tag}, 'OLD 2 mM  (rskR2 + gates)'), 1);
iZ = find(strcmp({R.tag}, 'NEW 2 mM  (gates OFF)'), 1);
if isempty(iN) || isempty(iO), return; end
fprintf('\n\n======== PER-FEATURE 2 mM COST (%s list) ========\n', which);
fprintf('%-46s %10s %10s %10s %10s\n', 'feature', 'OLD 2mM', 'NEW 2mM', 'NEW nogate', 'NEW 8mM');
vN = R(iN).vNew; vO = R(iO).vNew;
vZ = []; if ~isempty(iZ), vZ = R(iZ).vNew; end
vF = []; if ~isempty(iF), vF = R(iF).vNew; end
for k = 1:numel(LIST)
    s = sprintf('%-46s %10.3f %10.3f', LIST{k}(1:min(46,end)), vO(k), vN(k));
    if ~isempty(vZ), s = [s sprintf(' %10.3f', vZ(k))]; else, s = [s blanks(11)]; end %#ok<AGROW>
    if ~isempty(vF), s = [s sprintf(' %10.3f', vF(k))]; end %#ok<AGROW>
    fprintf('%s\n', s);
end
fprintf('%-46s %10.3f %10.3f', 'TOTAL', sum(vO), sum(vN));
if ~isempty(vZ), fprintf(' %10.3f', sum(vZ)); else, fprintf('%11s',''); end
if ~isempty(vF), fprintf(' %10.3f', sum(vF)); end
fprintf('\n');
end

function iRatios(R, tag2, tag8, label)
i2 = find(strcmp({R.tag}, tag2), 1); i8 = find(strcmp({R.tag}, tag8), 1);
if isempty(i2) || isempty(i8), return; end
m2 = R(i2).fm; m8 = R(i8).fm; d2 = R(i2).fd; d8 = R(i8).fd;
keys = {'steady','ktr','rsK','peak1_y','vall_y','peak2','A'};
fprintf('\n%s\n', label);
fprintf('%-10s %10s %10s %9s | %10s %10s %9s | %8s\n', ...
        'feature','model 8mM','model 2mM','ratio','data 8mM','data 2mM','ratio','m/d');
for k = 1:numel(keys)
    f = keys{k};
    if ~isfield(m2,f) || ~isfield(m8,f), continue; end
    a8 = mean(m8.(f)(:),'omitnan'); a2 = mean(m2.(f)(:),'omitnan');
    if isfield(d2,f) && isfield(d8,f)
        b8 = mean(d8.(f)(:),'omitnan'); b2 = mean(d2.(f)(:),'omitnan');
    else
        b8 = NaN; b2 = NaN;
    end
    rm = a2/a8; rd = b2/b8;
    fprintf('%-10s %10.4g %10.4g %9.3f | %10.4g %10.4g %9.3f | %8.3f\n', ...
            f, a8, a2, rm, b8, b2, rd, rm/rd);
end
end
