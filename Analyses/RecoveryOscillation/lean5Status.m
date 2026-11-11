function S = lean5Status(campaign)
%LEAN5STATUS  One-screen progress read on a running campaign.
%
%   lean5Status            % the default campaign (campaignSpec)
%   S = lean5Status('lean5')
%
% Reads params/<tag>_state.mat for every profile of the campaign and reports
% where each one is, what it has bought, and -- the reason this exists -- WHICH
% FEATURES it bought it from. Two campaigns in a row put their whole gain into
% force-velocity while the features they existed to fix went backwards
% (Analyses/RecoveryOscillation sec. 7), and that is invisible in the total.
%
% A round is only counted once it is written to state, so a profile mid-round
% shows the previous round's numbers. Reading the state file while the
% optimiser is writing it is possible; a failed load is reported, not thrown.
%
% See also campaignSpec, campaignStatus, analyzeOptIterUniqueness.

if nargin < 1 || isempty(campaign); campaign = ''; end
SPEC = campaignSpec(campaign);
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

fprintf('\n=== campaign %s ===  seed %.4f  (%s)\n', ...
        SPEC.name, SPEC.seedExpect, SPEC.note);

S = struct('tag',{},'rounds',{},'best',{},'gain',{},'stall',{},'age_min',{});
for t = 1:numel(SPEC.tags)
    tag = SPEC.tags{t};
    f   = fullfile(root, 'params', [tag '_state.mat']);
    if ~isfile(f); continue; end

    d = dir(f); age = (now - datenum(d.date)) * 24 * 60;    %#ok<TNOW1,DATNM>
    try
        L = load(f, 'state');
    catch ME
        fprintf('  %-10s state unreadable (%s) -- mid-write, retry\n', tag, ME.identifier);
        continue;
    end
    if ~isfield(L,'state'); continue; end
    st = L.state;

    best   = iField(st, {'best_cost'}, NaN);
    rounds = iField(st, {'round'},     NaN);
    stall  = iField(st, {'stall'},     NaN);
    gain   = SPEC.seedExpect - best;

    fprintf(['  %-10s round %-4s  best %8.4f  (%+7.4f, %+5.1f %%)  ' ...
             'stall %-4s  state %.0f min old\n'], ...
            tag, num2str(rounds), best, -gain, -100*gain/SPEC.seedExpect, ...
            num2str(stall), age);

    % Per-feature: seed vs now. This is the column that matters.
    if isfield(st,'optIter') && ~isempty(st.optIter)
        it = st.optIter;
        nm = {}; if isfield(st,'featNames'); nm = st.featNames; end
        have = isfield(it,'featCost') && ~isempty(nm) && ...
               ~isempty(it(1).featCost) && ~isempty(it(end).featCost);
        if have
            c0 = it(1).featCost(:)';  cN = it(end).featCost(:)';
            n  = min([numel(c0), numel(cN), numel(nm)]);
            d_ = cN(1:n) - c0(1:n);
            [~, o] = sort(abs(d_), 'descend');
            fprintf('             movers (now - first round):\n');
            for i = o(1:min(6,numel(o)))
                mark = '';
                if d_(i) > 0.005; mark = '  <-- WORSE'; end
                fprintf('               %-22s %+8.4f   (%.4f -> %.4f)%s\n', ...
                        nm{i}, d_(i), c0(i), cN(i), mark);
            end
        end
    end

    S(end+1) = struct('tag',tag,'rounds',rounds,'best',best, ...
                      'gain',gain,'stall',stall,'age_min',age);  %#ok<AGROW>
end

if isempty(S)
    fprintf('  no state files yet -- the first round has not been written\n');
end
fprintf('\n');
end

function v = iField(s, names, dflt)
v = dflt;
for i = 1:numel(names)
    if isfield(s, names{i}) && ~isempty(s.(names{i}))
        v = s.(names{i})(end); return;
    end
end
end
