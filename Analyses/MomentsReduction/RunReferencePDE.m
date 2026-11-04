% RunReferencePDE.m
% Run the full strain-distributed model (fastgate_best) on every protocol in
% momentsProtocols and save the trajectories + wall clock to results/ref_<name>.mat.
% These are the ground truth for both the a-priori closure test and the
% a-posteriori moments-model comparison.

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
prefix = 'ref_';
if exist('dSoverride', 'var') && ~isempty(dSoverride)
    % grid-convergence check of the PDE itself: same strain window, finer bins
    % MaxStrainArraySize caps ss, so scale it to keep the 120 nm window width
    params0.MaxStrainArraySize = round(params0.MaxStrainArraySize * params0.dS / dSoverride);
    params0.dS = dSoverride;
    params0 = getParams(params0, [], true, false);
    prefix = sprintf('refdS%g_', dSoverride*1e3);
end
if exist('windowFactor', 'var') && ~isempty(windowFactor)
    % window-width check: same dS, window x windowFactor
    params0.MaxStrainArraySize = round(params0.MaxStrainArraySize * windowFactor);
    params0 = getParams(params0, [], true, false);
    prefix = sprintf('%sW%g_', prefix, windowFactor);
end
if exist('clampKernels', 'var') && clampKernels
    % hold every strain kernel at its outer-knot value instead of extrapolating
    % the end cubic (k2 goes NEGATIVE past +142 nm otherwise)
    params0.ClampStrainDepExtrapolation = true;
    prefix = [prefix 'clamp_'];
end
P = momentsProtocols();
if ~exist('onlyRuns', 'var'), onlyRuns = {P.name}; end
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end

for i = 1:numel(P)
    if ~ismember(P(i).name, onlyRuns), continue; end
    pr = params0;
    pr.Velocity = P(i).V(:);
    tic;
    [~, out] = evaluateModel(@dPUdT_CombinedTransitions, P(i).T, pr);
    wall = toc;
    ss = pr.ss;
    ref.name = P(i).name; ref.T = P(i).T; ref.V = P(i).V; ref.wall = wall;
    ref.t = out.t(:); ref.Force = out.Force(:); ref.FXB = out.FXB(:);
    ref.SL = out.SL(:); ref.LSE = out.LSE(:); ref.LXBPivot = out.LXBPivot(:);
    ref.p1 = out.PU(:, 1:ss); ref.p2 = out.PU(:, ss+1:2*ss);
    ref.tail = out.PU(:, 2*ss+1:end);
    ref.SR = out.SR(:); ref.PD = out.PuR(:);
    ref.p1_0 = out.p1_0(:); ref.p2_0 = out.p2_0(:);
    save(fullfile(here, 'results', [prefix P(i).name '.mat']), 'ref');
    fprintf('%-8s  %6.1f s wall  %5d samples  Fend %.2f kPa  Fmax %.2f\n', ...
        P(i).name, wall, numel(ref.t), ref.Force(end), max(ref.Force));
end
