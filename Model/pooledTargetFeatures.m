function [fd, info] = pooledTargetFeatures(fd, mode)
% POOLEDTARGETFEATURES  Swap the single-preparation slack target for the pooled one.
%
%   FD = POOLEDTARGETFEATURES(FD)          replace every pooled feature
%   FD = POOLEDTARGETFEATURES(FD, MODE)    MODE selects what gets replaced
%   [FD, INFO] = ...                       INFO reports what actually happened
%
% FD is the data-feature struct `runSlackExperiment` reads out of the curated
% protocol file — i.e. ONE preparation (03/27), in that fibre's own kPa. This
% function overwrites it with `pooledTarget8mM`: three preparations, each
% divided by its own pre-slack isometric, rundown-corrected, features pooled
% per preparation and then averaged, and multiplied back by the ONE pooled
% level Fref. See Analyses/DataPooling/conclusions.md.
%
% MODE (params0.UsePooledTarget)
%   true | 'all'    every feature the pooled target carries          (default)
%   'shape'         scale-free features + force features rescaled to keep the
%                   INCOMING level, so only the shape moves. Use this to test
%                   the shape change alone when the absolute level is being
%                   held by something else.
%   'level'         force features only, i.e. the -13 % level step alone.
%   false           no-op (the caller normally skips us entirely)
%
% WHAT IT DELIBERATELY DOES NOT TOUCH
%   SLdiff / SLslack / v_restretch — the commanded length ladder, identical
%     across the three days by protocol, so there is nothing to pool.
%   Everything outside the slack block (FV, ktr protocol, passive PS_*) — those
%     are not pooled yet (conclusions.md Sec "not ready"), so they keep coming
%     from their own sources.
%
% WHY THIS AND NOT A SEPARATE COST TERM. The objective already knows how to
% score `features_data` vs `features_model`; a pooled target is a different
% VALUE for the same quantity, not a new kind of penalty. Replacing the value
% keeps one cost path (fjezek: evaluation belongs in the feature path, not in
% disconnected side-terms) and keeps every existing fn entry, weight and
% diagnostic working unchanged.
%
% See also pooledTarget8mM, Analyses/DataPooling/WritePooledTarget,
%          runSlackExperiment, evalFeatureCost

if nargin < 2 || isempty(mode); mode = 'all'; end
% writeParamsToMFile serialises every scalar with %g, so a snapshot brings the
% switch back as the DOUBLE 1, not logical true. Accept both.
if isnumeric(mode) || islogical(mode)
    if ~any(mode); info = struct('applied',{{}}); return; end
    mode = 'all';
end
mode = validatestring(mode, {'all','shape','level'});

PT = pooledTarget8mM();
nEl = size(PT.val, 2);

info = struct('applied', {{}}, 'added', {{}}, 'kept', {{}}, 'Fref', PT.Fref, 'mode', mode);

% The pooled target is one column per slack segment of the 2026 ladder. If the
% incoming struct is a different length we are looking at a different protocol
% and quietly aligning the two would fit the wrong segment to the wrong
% amplitude — fail instead.
ref = [];
for f = {'steady','A','ktr'}
    if isfield(fd, f{1}) && isnumeric(fd.(f{1})); ref = numel(fd.(f{1})); break; end
end
if ~isempty(ref) && ref ~= nEl
    error('pooledTargetFeatures:ladderMismatch', ...
        ['incoming slack features have %d segments, the pooled target has %d — ' ...
         'the length ladder does not match, refusing to align them.'], ref, nEl);
end

isF = strcmp(PT.cls, 'F');
for k = 1:numel(PT.name)
    nm = PT.name{k};
    switch mode
        case 'shape'; if isF(k) && ~isfield(fd, nm); continue; end
        case 'level'; if ~isF(k); continue; end
    end

    v = PT.val(k, :);

    % 'shape': put the pooled shape back on the incoming preparation's own
    % level, so this call changes only the relative pattern. The incoming
    % level is read from `steady`, exactly as Fref is derived from the anchor.
    if strcmp(mode, 'shape') && isF(k)
        if ~isfield(fd, 'steady'); error('pooledTargetFeatures:noLevel', ...
                'mode ''shape'' needs `steady` in the incoming features'); end
        kS = find(strcmp(PT.name, 'steady'), 1);
        v  = v * mean(double(fd.steady), 'omitnan') / mean(PT.val(kS,:), 'omitnan');
    end

    % A NaN in the pooled target means no preparation resolved that segment.
    % Injecting it would cost a flat NAN_COST in evalFeatureCost, so keep
    % whatever was there and say so.
    if isfield(fd, nm) && isnumeric(fd.(nm)) && numel(fd.(nm)) == nEl
        old = double(fd.(nm)(:))';
        bad = ~isfinite(v);
        if any(bad)
            v(bad) = old(bad);
            info.kept{end+1} = sprintf('%s[%s]', nm, strtrim(sprintf('%d ', find(bad))));
        end
        info.applied{end+1} = nm;
    else
        if any(~isfinite(v)); continue; end   % nothing to fall back on
        info.added{end+1} = nm;
    end
    fd.(nm) = v;
end

% t0_crossing is the same measurement as t0 under a second name (the fn entry
% 't0_crossing|SLdiff' scores it against the length ladder). runSlackExperiment
% aliases it when reading the file; keep the alias true after the swap.
if isfield(fd, 't0'); fd.t0_crossing = fd.t0; end

% A_over_steady is DERIVED, so it has to be rebuilt from the pooled A and
% steady rather than left at the single preparation's ratio that
% extractSlackAttributes computed from the raw trace. Recomputing it here keeps
% the ratio target consistent with the two values it is a ratio of.
if isfield(fd,'A') && isfield(fd,'steady') && numel(fd.A) == numel(fd.steady)
    fd.A_over_steady = fd.A(:)' ./ fd.steady(:)';
end

% ktr_shape is DERIVED the same way -- rebuild it from the POOLED ktr rather
% than leave the single preparation's ladder behind. See extractSlackAttributes
% for why the ratio exists at all.
if isfield(fd,'ktr') && numel(fd.ktr) > 1 && isfinite(fd.ktr(1)) && fd.ktr(1) ~= 0
    fd.ktr_shape = fd.ktr(:)' / fd.ktr(1);
end
end
