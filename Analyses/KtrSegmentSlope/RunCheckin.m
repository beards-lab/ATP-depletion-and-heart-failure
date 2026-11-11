% RunCheckin.m -- cheap status read for a running campaign; writes the JSON
% the remote dashboard is built from.
%
% Reads params/<TAG>_state.mat only (plus the seed breakdown, once), so it
% costs no cores and never touches the running pool.
%
% BASELINE. Every delta is SEED -> CURRENT BEST. optimizeFeatures records a
% featCost vector only for improved iterates, so the seed is not in
% state.optIter; SaveSeedBreakdown.m evaluates it once and stores it next to
% the results. The first version of this checker diffed the first improved
% snapshot against the latest -- identical after one round -- and reported
% "*** not moving ***" on a round that had just banked -3.3 %.
%
% The two triggers that matter for the low-ATP call:
%   * RE-STRETCH block rising  -> the search is selling the features whose
%     8 mM error propagates to 2 mM; stop and re-floor
%   * ktr_shape not falling    -> the ladder slope is still not being bought
%
%   matlab -batch "TAG='refit6'; run('Analyses/KtrSegmentSlope/RunCheckin.m')"
if ~exist('TAG','var'); TAG = 'refit6'; end
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');

if ~exist('ATPC','var')
    ATPC = {'ktr80','ktr80_shape','ktrfit','ktrfit_shape','rsK80','rsT20','ktrP_k80','steady','t0_crossing'};
end
GUARD = {'A','peak1_y','vall_y','peak2','vall2_dy','restretchSlopeStart','A_over_steady','vall_t','peak1_dSL'};
RESTRETCH = ATPC;   % the block whose fall is the point of the campaign

sf = fullfile(root,'params',[TAG '_state.mat']);
if ~isfile(sf); fprintf(2,'no state file %s\n', sf); return; end
S = load(sf); st = S; if isfield(S,'state'); st = S.state; end

nm = {}; if isfield(st,'featNames'); nm = st.featNames(:)'; end
% strip both the options and the element selector, so 'rsK80[1:4]|fl:..|w'
% is recognised as rsK80 when the blocks are summed
lab = cellfun(@(s) regexprep(regexprep(s,'\|.*$',''),'\[.*$',''), nm, 'UniformOutput', false);

% improved iterates only, in order
OI = []; if isfield(st,'optIter'); OI = st.optIter; end
imp = false(1,numel(OI));
for i = 1:numel(OI)
    imp(i) = ~isfield(OI(i),'isImproved') || isempty(OI(i).isImproved) || OI(i).isImproved;
end
OI = OI(imp);
n = numel(OI);

% ---- seed baseline ------------------------------------------------------
sb = fullfile(resdir,[TAG '_seed_breakdown.mat']);
haveSeed = isfile(sb);
if haveSeed
    L = load(sb); seedC = L.seed.featCost(:)'; seedN = L.seed.featNames(:)';
    seedCost = L.seed.cost;
    % Same fn -> same order, so align by POSITION. Names are not unique: the
    % objective scores `steady` twice (the 5-element term and the steady[1]
    % level guard), and a name lookup mapped both to the first one, printing
    % a phantom -0.0365 on the guard. Fall back to names only if the lists
    % genuinely differ, and then never reuse a seed entry.
    if numel(seedN) == numel(nm) && all(strcmp(seedN, nm))
        base = seedC;
    else
        base = nan(1,numel(nm)); used = false(1,numel(seedN));
        for j = 1:numel(nm)
            k = find(strcmp(seedN, nm{j}) & ~used, 1);
            if ~isempty(k); base(j) = seedC(k); used(k) = true; end
        end
    end
else
    fprintf(2,'no seed breakdown yet -- deltas fall back to the first improved iterate\n');
    base = nan(1,numel(nm)); seedCost = NaN;
    if n >= 1 && numel(OI(1).featCost) == numel(nm); base = OI(1).featCost(:)'; end
end

cur = base; bestCost = seedCost;
if n >= 1 && numel(OI(end).featCost) == numel(nm)
    cur = OI(end).featCost(:)'; bestCost = OI(end).cost;
end
if isfield(st,'best_cost'); bestCost = st.best_cost; end

iR = ismember(lab, RESTRETCH);
iG = ismember(lab, GUARD);
iF = find(strcmp(lab,'FV_normpowerAvg'),1);
pick = @(v,j) ifelse(isempty(j), NaN, v(j));

r0 = sum(base(iR),'omitnan'); r1 = sum(cur(iR),'omitnan');
k0 = sum(base(iG),'omitnan'); k1 = sum(cur(iG),'omitnan');
fv = pick(cur,iF);

% ---- print --------------------------------------------------------------
fprintf('=== %s ===  seed %.4f -> best %.4f (%+.1f %%) | rounds %d | improvements %d\n', ...
        TAG, seedCost, bestCost, 100*(bestCost-seedCost)/seedCost, st.round, n);
fprintf('  ATP block %.4f -> %.4f (%+.1f %%)   guards %.4f -> %.4f (%+.1f %%)   FV %.4f\n', ...
        r0, r1, 100*(r1-r0)/r0, k0, k1, 100*(k1-k0)/k0, fv);
d = cur - base; [~,o] = sort(d,'ascend');
for i = o(:)'
    if ~(abs(d(i)) > 1e-6 || base(i) > 1e-6); continue; end
    fprintf('    %-22s %8.4f %8.4f %+9.4f%s\n', lab{i}, base(i), cur(i), d(i), ...
            ifelse(iR(i),'  ATP',ifelse(iG(i),'  guard','')));
end

% ---- JSON ---------------------------------------------------------------
F = struct('name',{},'seed',{},'best',{},'delta',{},'restretch',{});
for i = 1:numel(nm)
    if ~(abs(d(i)) > 1e-6 || base(i) > 1e-6 || cur(i) > 1e-6); continue; end
    F(end+1) = struct('name',lab{i},'seed',base(i),'best',cur(i), ...
                      'delta',d(i),'restretch',iR(i)); %#ok<SAGROW>
end
T = struct('round',{},'cost',{},'restretch',{},'ktr_shape',{});
for i = 1:n
    fc = OI(i).featCost(:)';
    T(end+1) = struct('round',OI(i).round,'cost',OI(i).cost, ...
                      'restretch',sum(fc(iR),'omitnan'),'ktr_shape',sum(fc(iG),'omitnan')); %#ok<SAGROW>
end
J = struct();
J.tag = TAG; J.updated = datestr(now,'yyyy-mm-dd HH:MM');
J.seed_cost = seedCost; J.best_cost = bestCost; J.rounds = st.round;
J.improvements = n; J.seed_exact = haveSeed;
J.restretch = struct('seed',r0,'now',r1,'pct',100*(r1-r0)/r0,'flag',r1 > 1.20*r0);
J.ktr_shape = struct('seed',k0,'now',k1,'pct',100*(k1-k0)/max(k0,eps),'flag',k1 > 1e-6);  % = the GUARD block; flag = a guard was breached
J.fv = fv;
J.features = F; J.trajectory = T;
fid = fopen(fullfile(resdir,['checkin_' TAG '.json']),'w');
fprintf(fid,'%s', jsonencode(J)); fclose(fid);
fprintf('wrote checkin_%s.json\n', TAG);

% one row per checkpoint for the dashboard's history table; STATE is passed
% in by the caller (the campaign log, not the state file, knows if it died)
if exist('APPEND','var') && APPEND
    if ~exist('STATE','var'); STATE = ''; end
    row = struct('time', datestr(now,'HH:MM'), 'round', st.round, 'best', bestCost, ...
                 'restretch', r1, 'ktr_shape', k1, 'state', STATE);
    fid = fopen(fullfile(resdir,['checkin_history_' TAG '.jsonl']),'a');
    fprintf(fid,'%s\n', jsonencode(row)); fclose(fid);
    fprintf('appended checkpoint row\n');
end

function out = ifelse(c,a,b)
if c; out = a; else; out = b; end
end
