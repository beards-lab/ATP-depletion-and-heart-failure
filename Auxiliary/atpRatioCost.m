function [cost, detail] = atpRatioCost(fm2, fm8, target)
%ATPRATIOCOST  Relative-ratio low-ATP cost (log-space, per-segment vectors).
%   [cost,detail] = atpRatioCost(fm2, fm8, target)
%   fm2/fm8  : model feature structs (5-element fields) for 2 mM and 8 mM runs.
%   target   : struct (.scalar, .vector, .weights) from the data ratios.
%
%   Scalar features  -> compare mean(fm2)/mean(fm8) to target.scalar.(f).
%   Vector features  -> compare per-segment fm2./fm8 to target.vector.(f) (5-vec).
%   Cost = sum of w * (log model_ratio - log data_ratio)^2.
    cost = 0; detail = struct();
    W = target.weights;
    sf = fieldnames(target.scalar);
    for i = 1:numel(sf)
        f = sf{i};
        if ~isfield(fm2,f) || ~isfield(fm8,f); continue; end
        r = mean(double(fm2.(f)),'omitnan') / mean(double(fm8.(f)),'omitnan');
        c = getW(W,f) * (log(r) - log(target.scalar.(f)))^2;
        cost = cost + c;
        detail.(f) = struct('ratio',r,'target',target.scalar.(f),'cost',c);
    end
    vf = fieldnames(target.vector);
    for i = 1:numel(vf)
        f = vf{i};
        if ~isfield(fm2,f) || ~isfield(fm8,f); continue; end
        r  = double(fm2.(f)) ./ double(fm8.(f));
        tv = double(target.vector.(f));
        n  = min(numel(r), numel(tv));
        c  = getW(W,f) * mean((log(r(1:n)) - log(tv(1:n))).^2, 'omitnan');
        cost = cost + c;
        detail.(f) = struct('ratio',r(1:n),'target',tv(1:n),'cost',c);
    end
end
function w = getW(W,f)
    if isfield(W,f); w = W.(f); else; w = 1; end
end
