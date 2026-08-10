function rep = checkFeatureCoverage(params0)
%CHECKFEATURECOVERAGE Which params0.fn entries cannot actually be scored?
%
%   REP = checkFeatureCoverage(PARAMS0)
%
%   Runs the battery ONCE and checks, for every entry in PARAMS0.fn, that the
%   feature it names exists on both sides of the comparison.
%
%   WHY: evalFeatureCost returns MISSING_FEATURE_COST = 100 (a scalar
%   sentinel, applied once) when a field is absent from features_data or
%   features_model. Multiplied by the entry's weight that is a plausible-
%   looking number -- an 'rsK|0.25' whose data target is missing scores
%   exactly 25 -- so a run seeded on a machine with stale data/*.mat looks
%   like a mediocre fit rather than a broken one, and optimises a constant.
%
%   Boundary entries ('field|lb-ub') are scored against the simulation only,
%   so they are checked against features_model alone.
%
%   REP fields:
%     .missingModel  - entries whose feature is absent from features_model
%     .missingData   - entries whose feature is absent from features_data
%                      (boundary entries excluded -- they never read it)
%     .ok            - true when both lists are empty
%
%   This runs the full battery, so it costs one evaluation. Call it once at
%   the start of a long run, not inside the objective.

    ATP_c         = params0.MgATP;          %#ok<NASGU>  read by RunBakersExp
    FV_velocities = params0.FV_velocities;  %#ok<NASGU>
    params0.PlotEachSeparately = 0;
    params0.PlotFeatureFitting = 0;
    params0.EvalFeatures       = 1;

    RunBakersExp;   % script: populates features_model / features_data here

    rep = struct('missingModel', {{}}, 'missingData', {{}}, 'ok', false);

    for i_fn = 1:numel(params0.fn)
        entry = params0.fn{i_fn};
        parts = strsplit(entry, ',');        % grouped entries share one slot
        for i_p = 1:numel(parts)
            tokens = strsplit(strtrim(parts{i_p}), '|');
            tokens = tokens(~cellfun(@isempty, tokens));
            if isempty(tokens); continue; end

            % primary field name: strip any [selector] tail
            feat = regexprep(tokens{1}, '\[.*$', '');
            if isempty(feat) || strncmp(feat, '@', 1); continue; end
            if strcmp(feat, 'AssertParams'); continue; end

            isBoundary = any(cellfun(@iIsRange, tokens(2:end)));

            if ~isfield(features_model, feat)
                rep.missingModel{end+1} = sprintf('%s  (in "%s")', feat, entry); %#ok<AGROW>
            end
            if ~isBoundary && ~isfield(features_data, feat)
                rep.missingData{end+1} = sprintf('%s  (in "%s")', feat, entry); %#ok<AGROW>
            end
        end
    end

    rep.ok = isempty(rep.missingModel) && isempty(rep.missingData);
end

function tf = iIsRange(tok)
% A range token looks like '0.4-1.0' or '100-1500' (not a bare number).
    tf = ~isempty(regexp(strtrim(tok), '^-?[\d.]+\s*-\s*-?[\d.]+$', 'once'));
end
