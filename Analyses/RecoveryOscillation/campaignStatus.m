function campaignStatus(tags)
%CAMPAIGNSTATUS  Where is each profile of the rig-aware campaign, really?
%
%   campaignStatus                       % rigaware_a .. rigaware_d
%   campaignStatus({'rigaware_a'})
%
% Answers the three questions you actually have when a long campaign looks
% wrong, and answers them from the STATE FILES rather than the logs, because a
% log tells you what was printed and the state tells you what was kept:
%
%   1. is it alive?      state-file mtime against the expected round time for
%                        that profile, plus a process check on Linux. A round
%                        is 25-114 min depending on profile, so a state file
%                        that is 40 minutes old is normal, not a hang.
%   2. is it descending? best_cost against the seed, and how long since the
%                        last accepted improvement.
%   3. is it stuck?      `stall` = consecutive rounds with no improvement.
%
% EXPECTED ROUND TIME. One objective evaluation is ~38 s at 5 workers, and a
% round is SIMPLEX_EVALS of them plus any surrogate sweep:
%     a  60 evals        ~38 min      c   40 evals       ~25 min
%     b  90 evals        ~57 min      d  180 evals      ~114 min
% Profile d's FIRST round is the slow one -- it spends 120 evaluations on a
% global surrogate before it ever proposes a point, so d showing nothing for
% two hours is d working as designed, not d hanging.
%
% See also RunOptimRigAware, RunRigAwareReport, analyzeOptIterUniqueness

if nargin < 1 || isempty(tags)
    tags = {'rigaware_a','rigaware_b','rigaware_c','rigaware_d'};
end
if ischar(tags) || isstring(tags); tags = cellstr(tags); end

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
SEED   = 2.6858;                       % BuildRigAwareSeed.m
MINUTES = containers.Map({'a','b','c','d'}, {38, 57, 25, 114});

fprintf('\ncampaign status   seed %.4f   %s\n', SEED, datestr(now,'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
fprintf('%s\n', repmat('=',1,104));
fprintf('%-12s %8s %8s %7s %6s %7s  %-11s %s\n', ...
        'tag','best','vs seed','round','stall','iters','state age','verdict');
fprintf('%s\n', repmat('-',1,104));

for i = 1:numel(tags)
    tag = tags{i};
    f   = fullfile(root, 'params', [tag '_state.mat']);
    p   = tag(end);
    expMin = 60; if isKey(MINUTES, p); expMin = MINUTES(p); end

    if ~isfile(f)
        fprintf('%-12s %8s %8s %7s %6s %7s  %-11s %s\n', tag, '-','-','-','-','-','-', ...
                'never started (no state file)');
        continue;
    end

    d  = dir(f);
    ageMin = (now - datenum(d.date)) * 24 * 60;                      %#ok<TNOW1,DATNM>
    S  = load(f); st = S.state;
    nIt = numel(st.optIter);

    % how long since the last ACCEPTED improvement, in rounds
    lastImpRound = NaN;
    if nIt > 0 && isfield(st.optIter, 'round'); lastImpRound = st.optIter(end).round; end

    alive = iAlive(tag);
    if ageMin > 3*expMin
        v = sprintf('STALLED/DEAD - state %.0f min old vs ~%d min/round', ageMin, expMin);
    elseif ~isempty(alive) && ~alive
        v = 'process gone - check the log tail';
    elseif st.stall >= 5
        v = sprintf('running but stuck (%d rounds no accept)', st.stall);
    else
        v = 'running normally';
    end
    if st.best_cost < SEED - 0.25
        v = ['DESCENDING - ' v];
    elseif st.best_cost > SEED + 1e-6
        v = ['ABOVE SEED?! - ' v];
    end

    fprintf('%-12s %8.4f %+7.1f%% %7d %6d %7d  %-11s %s\n', ...
        tag, st.best_cost, 100*(st.best_cost/SEED - 1), st.round, st.stall, nIt, ...
        sprintf('%.0f min', ageMin), v);
    if ~isnan(lastImpRound) && st.round > lastImpRound
        fprintf('%12s last accept was round %d of %d\n', '', lastImpRound, st.round);
    end
end

fprintf('%s\n', repmat('-',1,104));
fprintf(['Reading it:  a drop of less than 0.25 is inside the objective''s own peak-to-peak\n' ...
         'roughness (Docs objective_noise_floor), so "no decrease" only means something\n' ...
         'once several rounds have been accepted and rejected. stall counts ROUNDS, not\n' ...
         'evaluations.\n\n' ...
         'Logs:   Analyses/RecoveryOscillation/results/rigaware_<p>.log\n' ...
         'Detail: analyzeOptIterUniqueness(''<tag>'')  -- which parameters are actually moving\n\n']);
end

function tf = iAlive(tag)
%IALIVE  Is a MATLAB process running this tag? [] when it cannot be determined.
    tf = [];
    if ~isunix; return; end
    try
        [rc, out] = system(sprintf('pgrep -fa "PROFILE=''%s''" 2>/dev/null', tag(end)));
        if rc == 0 && ~isempty(strtrim(out)); tf = true; return; end
        [rc2, out2] = system('pgrep -c -u "$USER" -x MATLAB 2>/dev/null');
        if rc2 == 0 && str2double(strtrim(out2)) > 0; tf = []; else; tf = false; end
    catch
        tf = [];
    end
end
