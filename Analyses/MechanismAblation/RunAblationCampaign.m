% RunAblationCampaign.m -- run a list of ablation arms back to back.
%
% One MATLAB session, one parallel pool, many arms: saves ~30 s of pool startup
% per arm and keeps a failed arm from taking the rest of the list with it.
%
% SIZING. 38 runs (35 arms + 3 controls) x ~38 min = ~24 h serial. This box has
% 12 cores and one evaluation saturates at 5-6 workers, so TWO instances
% (2 x 5 = 10 cores, 2 spare for the OS) is the right parallelism: ~12 h.
% Four instances would oversubscribe and every arm would slow down.
%
%   matlab -batch "HALF=1; run('Analyses/MechanismAblation/RunAblationCampaign.m')"
%   matlab -batch "HALF=2; run('Analyses/MechanismAblation/RunAblationCampaign.m')"
%
% Or an explicit list:
%   matlab -batch "LIST={'C1','C2','B4'}; run('Analyses/MechanismAblation/RunAblationCampaign.m')"
%
% SKIP is honoured: arms the preflight found to be NO-OP or BROKEN are not worth
% a 38-minute search. Pass SKIP={} to run them anyway.
%
% See also: RunAblationArm, RunAblationPreflight, RunAblationReport

clc;
campHere = fileparts(mfilename('fullpath'));
campRoot = fullfile(campHere, '..', '..');
cd(campRoot); addpath(genpath(campRoot));

if ~exist('LIST','var') || isempty(LIST)
    A = ablationSpec();
    LIST = [{'control:1','control:2','control:3'}, {A.id}];
    if exist('HALF','var') && ~isempty(HALF)
        LIST = LIST(HALF:2:end);     % interleaved, so both halves carry controls
    end
end
if ~exist('SKIP','var'); SKIP = iPreflightSkips(campHere); end
LIST = LIST(~ismember(LIST, SKIP));

fprintf('\n=== ablation campaign: %d runs ===\n%s\n\n', numel(LIST), strjoin(LIST, ' '));
campFailed = {};
campT0 = tic;

for campI = 1:numel(LIST)
    campItem = LIST{campI};
    fprintf('\n########## [%d/%d] %s   (elapsed %.1f h) ##########\n', ...
            campI, numel(LIST), campItem, toc(campT0)/3600);
    try
        clearvars -except LIST SKIP HALF campHere campRoot campFailed campI campItem campT0
        if startsWith(campItem, 'control:')
            ARM = 'control'; REP = str2double(extractAfter(campItem, ':'));
        else
            ARM = campItem;
        end
        run(fullfile(campHere, 'RunAblationArm.m'));
    catch ME
        fprintf(2, '!!! %s FAILED: %s\n', campItem, ME.message);
        campFailed{end+1} = sprintf('%s (%s)', campItem, ME.message); %#ok<SAGROW>
    end
end

fprintf('\n=== campaign done in %.1f h ===\n', toc(campT0)/3600);
if isempty(campFailed)
    fprintf('all %d runs completed\n', numel(LIST));
else
    fprintf(2, '%d failed:\n', numel(campFailed));
    for i = 1:numel(campFailed); fprintf(2, '  %s\n', campFailed{i}); end
end

% -------------------------------------------------------------------------
function skip = iPreflightSkips(here)
% Arms the frozen-parameter preflight showed cannot move: a dead switch, an
% inert parameter, or a flip that breaks the model. Spending a search on one of
% those buys a confident null.
%
% NOT skipped: a large but FINITE frozen cost. That is an UNCOMPENSATED
% mechanism, which is precisely what stage 1 exists to test -- A1 leaves the
% frozen fit at 64.9 and a single scalar on ka may well buy all of it back.
% Reading a big frozen delta as "load-bearing" is the mistake the whole
% budget-matched design is built to avoid.
skip = {};
f = fullfile(here, 'results', 'preflight.mat');
if ~isfile(f); return; end
S = load(f, 'R');
if ~isfield(S,'R'); return; end
bad = ismember({S.R.status}, {'NO-OP','ERROR','NONFINITE'});
skip = {S.R(bad).id};
if ~isempty(skip)
    fprintf('skipping %d arm(s) on the preflight verdict: %s\n', numel(skip), strjoin(skip,', '));
end
end
