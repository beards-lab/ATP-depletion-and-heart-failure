function out = apoPoolObjective(x, base, cfg, lo, hi)
%APOPOOLOBJECTIVE Bounded wrapper + best-so-far bookkeeping around apoPoolCost.
%
%   apoPoolObjective('reset', RESDIR, TAG)   initialise (call before fminsearch)
%   apoPoolObjective('state')                return the accumulated state struct
%   C = apoPoolObjective(X, BASE, CFG, LO, HI)   evaluate
%
%   X lives in R^n and is mapped into the box [LO, HI] by a logistic transform,
%   so fminsearch cannot walk out of the search box:
%       g = lo + (hi - lo) ./ (1 + exp(-x))
%   x = 0 corresponds to the geometric middle of the box, NOT the seed — the
%   caller computes x0 from the seed with the inverse transform.
%
%   State (best cost, best g, full history) is kept in a persistent so it
%   survives across fminsearch calls, and is saved to disk on every improvement
%   — a run that is interrupted can be resumed from the last best point.

    persistent best gbest hist resdir tag dbest

    if ischar(x) || isstring(x)
        switch lower(char(x))
            case 'reset'
                best = Inf; gbest = []; hist = []; dbest = struct();
                resdir = base; tag = cfg;      % arg slots reused for reset
                out = [];
            case 'state'
                out = struct('best',best,'gbest',gbest,'hist',hist,'dbest',dbest);
            otherwise
                error('apoPoolObjective:mode','unknown mode "%s"', char(x));
        end
        return;
    end

    g = lo(:) + (hi(:) - lo(:)) ./ (1 + exp(-x(:)));
    [c, d] = apoPoolCost(g, base, cfg);
    hist(end+1, :) = [c, g(:)'];

    if c < best
        best = c; gbest = g; dbest = d;
        state = struct('best',best,'gbest',gbest,'hist',hist,'dbest',dbest, ...
                       'mods',{cfg.mods},'lo',lo(:),'hi',hi(:));
        if ~isempty(resdir)
            save(fullfile(resdir, ['state_' tag '.mat']), 'state', 'cfg', '-v7.3');
        end
        reportRow(sprintf('BEST @%d', size(hist,1)), c, d);
        for i = 1:numel(cfg.mods)
            fprintf('        %-9s %10.4g  (x%.3f)\n', cfg.mods{i}, ...
                    cfg.seed.(cfg.mods{i})*g(i), g(i));
        end
    end
    out = c;
end

function reportRow(tag, C, D)
    if ~isfield(D, 'kR')
        fprintf('%-14s cost %.4f  (evaluation failed: %s)\n', tag, C, ...
                ternstr(isfield(D,'err'), getfieldSafe(D,'err'), 'no detail'));
        return;
    end
    fprintf(['%-14s C=%8.4f | rate %6.3f shape %6.3f trace %6.3f phys %6.3f' ...
             ' | kS %5.1f kR %5.1f kK %5.1f | peak %.3f vall %.3f' ...
             ' | E slack %7.1f ktr %6.3f fv %6.3f\n'], ...
            tag, C, D.c_rate, D.c_shape, D.c_trace, D.c_phys, ...
            D.kS, D.kR, D.kK, D.peak1, D.valley, D.E_slack, D.E_ktr, D.E_fv);
end

function s = ternstr(c, a, b); if c; s = a; else; s = b; end; end
function v = getfieldSafe(S, f); if isfield(S,f); v = S.(f); else; v = ''; end; end
