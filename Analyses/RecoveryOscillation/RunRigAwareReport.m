% RunRigAwareReport.m -- consolidate the four rig-aware campaigns into one
% decision, and refuse to call any of it progress until it is verified.
%
% Four independent searches produce four snapshots and four cost histories.
% The lowest number is not automatically the answer, for two reasons this
% project has already been bitten by:
%
%   1. A SNAPSHOT HEADER IS NOT A COST. params/<tag>_opt.m carries the cost the
%      search believed it had; the objective is not bit-reproducible (the
%      MaxRunTime watchdog and the parallel pool make it rough at the ~0.076
%      level, Docs objective_noise_floor). Every snapshot here is RE-EVALUATED
%      from disk, and a header that disagrees by more than the roughness is
%      reported as such (Analyses/StairsPrediction: an iteration header was
%      fake, and a campaign was ranked on it).
%
%   2. A SHORT-BUDGET WINNER IS NOT A WINNER. Profiles differ in evaluations
%      per round, so at a fixed wall clock they have had different numbers of
%      objective evaluations. The table therefore reports rounds AND the cost
%      each profile had reached at the round count of the SHORTEST run, so the
%      ranking is budget-matched as well as wall-clock-matched.
%
% It also checks the thing that matters more than the total: whether the drop
% came from the broken features or was bought by giving back the solved ones.
% rsK80 and ktr are the ATP-carrying rates and are already inside one
% preparation-width; a campaign that lowers the total while raising those has
% traded a solved feature for a broken one.
%
% Output -> results/rigaware_report.txt
%           results/rigaware_report.png
%           results/rigaware_report.mat

clearvars -except TAGS; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('TAGS','var') || isempty(TAGS)
    TAGS = {'rigaware_a','rigaware_b','rigaware_c','rigaware_d'};
end

SEED     = 'rigaware_seed';
SEED_EXP = 2.6858;
ROUGH    = 0.25;         % peak-to-peak numerical roughness of the objective

if isempty(gcp('nocreate')); parpool('local', 5); end

%% ---- re-evaluate everything from disk ----------------------------------
C = struct([]);
C = addRow(C, SEED);
for i = 1:numel(TAGS)
    f = fullfile('params', [TAGS{i} '_opt.m']);
    if ~isfile(f)
        fprintf('  %s: no snapshot yet -- skipping\n', TAGS{i}); continue
    end
    C = addRow(C, [TAGS{i} '_opt']);
end
if isempty(C); error('RunRigAwareReport:nothing','no snapshots found'); end

%% ---- read each campaign's own history ----------------------------------
H = struct([]);
for i = 1:numel(TAGS)
    sf = fullfile('params', [TAGS{i} '_state.mat']);
    if ~isfile(sf); continue; end
    S = load(sf); st = S.state;
    h = struct('tag', TAGS{i}, 'round', st.round, 'stall', st.stall, ...
               'best', st.best_cost, ...
               'cost', arrayfun(@(x) x.cost, st.optIter));
    if isempty(H); H = h; else; H(end+1) = h; end
end

%% ---- the table ---------------------------------------------------------
fid = fopen(fullfile(resDir,'rigaware_report.txt'),'w');
p = @(varargin) tee(fid, varargin{:});

p('RIG-AWARE CAMPAIGN -- four searches, one objective\n');
p('%s\n', repmat('=',1,96));
p('seed %s = %.4f (expected %.4f)\n\n', SEED, C(1).cost, SEED_EXP);

p('%-22s %9s %9s %9s %7s %7s %7s\n', ...
  'snapshot','header','re-eval','delta','vs seed','rounds','stall');
p('%s\n', repmat('-',1,96));
for i = 1:numel(C)
    j = [];
    if ~isempty(H); j = find(strcmp({H.tag}, strrep(C(i).name,'_opt','')), 1); end
    if isempty(j); rd = NaN; stl = NaN; else; rd = H(j).round; stl = H(j).stall; end
    if isfinite(C(i).header); d = C(i).cost - C(i).header; else; d = NaN; end
    flag = '';
    if isfinite(d) && abs(d) > ROUGH; flag = '  <-- HEADER DISAGREES'; end
    p('%-22s %9s %9.4f %9s %6.1f%% %7s %7s%s\n', C(i).name, ...
      num2str(C(i).header,'%.4f'), C(i).cost, num2str(d,'%+.4f'), ...
      100*(C(i).cost/C(1).cost - 1), num2str(rd,'%d'), num2str(stl,'%d'), flag);
end

%% ---- budget-matched ranking --------------------------------------------
if ~isempty(H)
    rMin = min([H.round]);
    p('\nBUDGET-MATCHED: where each profile stood at round %d (the shortest run)\n', rMin);
    p('%-22s %10s %10s %10s\n','profile','@round','final','rounds');
    for i = 1:numel(H)
        c = H(i).cost;
        p('%-22s %10.4f %10.4f %10d\n', H(i).tag, ...
          c(min(rMin, numel(c))), H(i).best, H(i).round);
    end
    p('A profile that is behind at the matched round but ahead at the end has\n');
    p('bought its lead with evaluations, not with a better search.\n');
end

%% ---- did the drop come from the broken features? -----------------------
if numel(C) < 2
    p('\n(no campaign snapshots yet -- seed only)\n');
    iBest = 1;
else
    p('\nWHERE THE DROP CAME FROM  (per-feature cost; seed -> best)\n');
    [~, iBest] = min([C(2:end).cost]); iBest = iBest + 1;
    p('best = %s (%.4f)\n\n', C(iBest).name, C(iBest).cost);
    nm = C(1).featNames;
    d  = C(iBest).featCost - C(1).featCost;
    [~, o] = sort(d, 'ascend');
    p('%-24s %9s %9s %9s\n','feature','seed','best','change');
    for i = o
        if abs(d(i)) < 1e-4; continue; end
        warn = '';
        if any(strcmp(nm{i}, {'rsK80','ktr','FV_normpowerAvg'})) && d(i) > 0.02
            warn = '   <-- a SOLVED feature got worse';
        end
        p('%-24s %9.4f %9.4f %+9.4f%s\n', nm{i}, ...
          C(1).featCost(i), C(iBest).featCost(i), d(i), warn);
    end

    p('\nNEXT: analyzeOptIterUniqueness(''%s'') -- are the solutions around the\n', ...
      strrep(C(iBest).name,'_opt',''));
    p('best cost UNIQUE (parameters reconverge) or DEGENERATE (different\n');
    p('parameters, same cost)? A degenerate best is not a result, it is a\n');
    p('direction the data cannot see.\n');
end
fclose(fid);

save(fullfile(resDir,'rigaware_report.mat'), 'C', 'H');
plotRigAware(C, H, resDir);
fprintf('\nreport -> %s\n', fullfile(resDir,'rigaware_report.txt'));

% ------------------------------------------------------------------------
function C = addRow(C, tag)
%ADDROW  Re-evaluate one snapshot from disk, and read the cost its header
%   claims so the two can be compared.
    % The header line writeParamsToMFile emits for a gallery snapshot:
    %   % param_opt_iter #9 | cost 2.7095 | round 9 | direct
    % Matched only inside a COMMENT and only as a whole word -- a loose
    % pattern picks up `params0.drawForceOnset = 0` (that string ends in
    % "et") and invents a header of 0 for a file that has none.
    f = fullfile('params', [tag '.m']);
    hdr = NaN;
    txt = fileread(f);
    m = regexp(txt, '(?m)^\s*%[^\r\n]*?cost\s+([-\d.eE+]+)', 'tokens', 'once');
    if ~isempty(m); hdr = str2double(m{1}); end

    params0 = getParams(loadParams(tag), [], true, false);
    params0.RunForceVelocity = true; params0.RunKtr = false; params0.RunSlack = true;
    params0.RunStairs = false; params0.RunForceVelocityTime = false;
    params0.RunSlackPassive = true; params0.EvalFeatures = true;
    params0.BreakOnODEUnstable = false; params0.PlotEachSeparately = 0;
    params0.PlotFeatureFitting = 0; params0.RunSlackSegments = 'AllPar';
    params0.OptimizeOn = 'Feats'; params0.MaxRunTime = 300;
    params0.mods = {}; params0.g = ones(1,0);
    [Et, ~, featCost, featNames] = evaluateBakersExp(1, params0, true, false);
    fprintf('  %-22s %.4f\n', tag, Et);

    % {{featNames}} not {featNames}: struct() unwraps ONE cell layer, so the
    % single-brace form would build a 1 x nFeature struct array instead of one
    % row holding a cell.
    r = struct('name',tag, 'cost',Et, 'header',hdr, ...
               'featCost',featCost, 'featNames',{{featNames}});
    if isempty(C); C = r; else; C(end+1) = r; end
end

function tee(fid, varargin); fprintf(1, varargin{:}); fprintf(fid, varargin{:}); end
