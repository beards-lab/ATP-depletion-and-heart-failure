function [fn, info] = rescaleFnToCondition(fn, cond)
%RESCALEFNTOCONDITION  Put an 8 mM-weighted objective on another condition's scale.
%
%   FN = RESCALEFNTOCONDITION(FN, 2)
%   [FN, INFO] = ...
%
% WHY THIS IS NEEDED. Every data weight in the objective is w = lambda*R/s^2 and
% every deadband is fl = FLOOR_CHI2*nEl*s^2, with `s` the between-preparation
% spread measured at 8 mM. That makes each term a chi-2 -- "1.0 = one
% preparation-width from the data" -- but only AT 8 mM. Score the same fn at
% 2 mM and the residual comes from one condition while the scale comes from
% another, so the numbers stop being chi-2 and stop being comparable between
% terms, which was the whole point of the weighting.
%
% It is not a small effect and it does not cancel: the 2 mM spreads are
% systematically WIDER (steady 18.4 vs 15.8 %, ktr 22.8 vs 17.2 %, peak1_y 22.4
% vs 17.0 %), so the 8 mM floors are too NARROW at 2 mM and every term pays
% more than it should.
%
% WHAT IT DOES. For each data term whose feature appears in both frozen pooled
% targets, with q = s_cond / s_8mM:
%
%     w  ->  w  / q^2        (weight is lambda*R/s^2)
%     fl ->  fl * q^2        (floor is FLOOR_CHI2*nEl*s^2)
%
% Boundary entries (`lb-ub`), regularisers and anything with no counterpart in
% the targets are left ALONE -- a boundary is a hard constraint, not a chi-2,
% and rewriting one silently destroys a hand-set guardrail (the mistake
% RecoveryOscillation sec. 9 records).
%
% This rescales the SCALE, never the target VALUE. The values come from
% pooledTargetFeatures, which injects the matching frozen target.
%
% See also pooledTargetFeatures, pooledTarget8mM, pooledTarget2mM,
%          Analyses/FeatureWeighting/BuildWeightedObjective

info = struct('changed', {{}}, 'skipped', {{}}, 'q', []);
if nargin < 2 || isempty(cond) || round(cond) == 8; return; end

A = pooledTarget8mM();
switch round(cond)
    case 2; B = pooledTarget2mM();
    otherwise
        error('rescaleFnToCondition:cond', 'no pooled target for %g mM', cond);
end

for i = 1:numel(fn)
    tk = strsplit(fn{i}, '|');
    nm = regexprep(tk{1}, '\[.*$', '');

    % Boundary / grouped / custom entries are not chi-2 terms. Leave them.
    isBoundary = any(cellfun(@(t) ~isempty(regexp(t, '^-?[\d.]+-{1}-?[\d.]+$', 'once')), tk)) ...
                 || contains(fn{i}, ',') || contains(fn{i}, '@');
    ia = find(strcmp(A.name, nm), 1);
    ib = find(strcmp(B.name, nm), 1);
    if isBoundary || isempty(ia) || isempty(ib) ...
            || ~isfinite(A.s(ia)) || ~isfinite(B.s(ib)) || A.s(ia) <= 0
        info.skipped{end+1} = nm;                                   %#ok<AGROW>
        continue;
    end

    q2 = (B.s(ib) / A.s(ia))^2;

    % trailing scalar weight, if the last token is a bare number
    if numel(tk) > 1
        w = str2double(tk{end});
        if ~isnan(w) && isempty(regexp(tk{end}, '-', 'once'))
            tk{end} = sprintf('%.6g', w / q2);
        end
    end
    % floor token
    for t = 1:numel(tk)
        if startsWith(tk{t}, 'fl:')
            flv = str2double(tk{t}(4:end));
            if ~isnan(flv); tk{t} = sprintf('fl:%.6g', flv * q2); end
        end
    end

    fn{i} = strjoin(tk, '|');
    info.changed{end+1} = nm;                                       %#ok<AGROW>
    info.q(end+1)       = sqrt(q2);                                 %#ok<AGROW>
end
end
