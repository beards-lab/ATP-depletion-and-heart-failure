function R = fitTwitchTier(cfg, nIterThisCall)
%FITTWITCHTIER  Resumable lsqnonlin fit of one tier (CFG from tierConfigs) to the EHT control twitch.
%   Residual (sigma*F_model - F_data)/F_peak on the 0..950 ms data grid, sigma solved in closed form.
%   Sequential (no pool), so it runs inside `matlab -batch`. Checkpoint results/<name>_ckpt.mat
%   after every iteration; every evaluation is appended to results/<name>.log. On rerun it resumes
%   from the checkpoint. Final diagnostics -> results/<name>_fit.mat, results/<name>.png.
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
if ~exist(res, 'dir'), mkdir(res); end
logf = fullfile(res, [cfg.name '.log']);  ckpt = fullfile(res, [cfg.name '_ckpt.mat']);  outf = fullfile(res, [cfg.name '_fit.mat']);
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);
lb = [cfg.par.lb];  ub = [cfg.par.ub];  x0 = [cfg.par.x0];
maxIt = cfg.maxIter;  it0 = 0;
if nargin < 2, nIterThisCall = Inf; end
if exist(ckpt, 'file')
    C = load(ckpt);  x0 = C.x;  it0 = C.iter;  maxIt = max(cfg.maxIter - C.iter, 0);
    lg(logf, sprintf('RESUME %s from iteration %d, SSR %.6f', cfg.name, C.iter, C.ssr));
end
maxIt = min(maxIt, nIterThisCall);        % a call may run only part of the budget (MCP calls must stay short)
if ~exist(ckpt, 'file')
    lg(logf, sprintf('START %s: %d free (%s); seed %s; Ca %s', cfg.name, numel(x0), strjoin({cfg.par.name}, ','), cfg.seed, cfg.ca.kind));
end
neval = 0;
    function r = resid(x)
        tic;
        try
            [base, th] = buildModel(cfg, x);
            fm = simHumanTwitch(base, th, tb);
            sig = max((fm'*Fd)/(fm'*fm), 0);
            r = (sig*fm - Fd)/max(Fd);
        catch e
            r = 10*ones(size(Fd));  sig = NaN;
            lg(logf, sprintf('EVAL ERROR %s', e.message));
        end
        neval = neval + 1;
        lg(logf, sprintf('eval %4d  SSR %.6f  rms %.3f %%  sigma %.3f  %5.1fs  x=[%s]', neval, sum(r.^2), 100*sqrt(mean(r.^2)), sig, toc, sprintf('%.4f ', x)));
    end
    function stop = ckfun(x, ov, state)
        stop = false;
        if strcmp(state, 'iter')
            iter = it0 + ov.iteration;  ssr = ov.resnorm;
            save(ckpt, 'x', 'iter', 'ssr');
            lg(logf, sprintf('ITER %d  SSR %.6f  rms %.3f %%', iter, ssr, 100*sqrt(ssr/numel(Fd))));
        end
    end
if maxIt > 0
    usePar = ~isempty(gcp('nocreate'));
    o = optimoptions('lsqnonlin', 'Display', 'off', 'UseParallel', usePar, 'FiniteDifferenceStepSize', cfg.fdStep, ...
        'StepTolerance', 1e-4, 'FunctionTolerance', 1e-6, 'MaxIterations', maxIt, 'OutputFcn', @ckfun);
    [x, ssr] = lsqnonlin(@resid, x0, lb, ub, o);
    C = load(ckpt);  x = C.x;  ssr = C.ssr;    % last checkpointed iterate (lsqnonlin's own x may differ by a rejected step)
    if isfinite(nIterThisCall) && C.iter < cfg.maxIter
        lg(logf, sprintf('PAUSE %s at iteration %d of %d', cfg.name, C.iter, cfg.maxIter));  R = [];  return
    end
else
    x = x0;  ssr = C.ssr;
end
lg(logf, sprintf('DONE %s  SSR %.6f  rms %.3f %%  x=[%s]', cfg.name, ssr, 100*sqrt(ssr/numel(Fd)), sprintf('%.4f ', x)));
R = finishTier(cfg, x, tb, Fd, D, outf, logf);
end

function lg(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), s);  fclose(fid);
end
