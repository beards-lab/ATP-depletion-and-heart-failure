%% ProbeGridConvergence.m -- is the strain discretisation converged?
%
% Serial slack runs at tight integration settings (smooth rectifier, RelTol
% 1e-4, dense output, extrapolation clamp) over a centred strain grid of
% constant span (~+-60 nm), varying the bin size dS and the bin phase: an even
% bin count straddles zero strain at +-dS/2, an odd count puts a bin on zero.
% If the model is converged in dS, both phases and all dS agree.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit8_opt'; end
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');

P0 = getParams(loadParams(PARAMSET), [], true, false);
P0.RunForceVelocity = false; P0.RunSlack = true; P0.RunStairs = false; P0.RunKtr = false;
P0.RunForceVelocityTime = false; P0.RunSlackPassive = true; P0.EvalFeatures = true;
P0.BreakOnODEUnstable = false; P0.PlotEachSeparately = 0; P0.PlotFeatureFitting = 0;
P0.RunSlackSegments = 'All'; P0.MaxRunTime = 1800; P0.OptimizeOn = 'Feats';
P0.ghostSave = ''; P0.ghostLoad = ''; P0.SaveBest = 0; P0.mods = {}; P0.g = ones(1,0);
P0.SmoothRectifierEps = 0.01; P0.odeRelTol = 1e-4; P0.odeAbsTol = 1e-7; P0.odeOutputDt = 5e-4;
P0.ClampStrainDepExtrapolation = true; P0.Slim_l = 1.2; P0.Slim_r = 3.2;

% {dS (um), bins}: span ~ +-60 nm, even and odd bin counts
G = { 0.002,  60; 0.002,  61; 0.001, 120; 0.001, 121; 0.0005, 240; 0.0005, 241 };
R = struct();
for g = 1:size(G, 1)
    P = P0; P.dS = G{g,1}; P.MaxStrainArraySize = G{g,2};
    P = getParams(P, P.g, true, false);
    tic; evalc('[~, o, fm, fd] = runSlackExperiment(P);'); w = toc;
    n = makeMonotonous(o.t); F70 = interp1(o.t(n), o.Force(n), 70);
    keep = cellfun(@(e) ~contains(e, ',') && isfield(fm, regexprep(regexprep(e, '\|.*$', ''), '\[.*$', '')), P.fn);
    c = sum(evalFeatureCost(fd, fm, P.fn(keep)));
    R(g).dS = G{g,1}; R(g).ss = P.ss; R(g).F70 = F70; R(g).cost = c; R(g).fm = fm; R(g).wall = w;
    fprintf('dS %.4f  %3d bins  s0 on bin: %d  %5.0f s  F(70s) %.4f  cost %.4f  ktr80 %s  rsK80 %s  rsT20(ms) %s\n', ...
        G{g,1}, P.ss, any(abs(P.s) < 1e-12), w, F70, c, mat2str(round(fm.ktr80(:)', 2)), ...
        mat2str(round(fm.rsK80(1:4), 2)), mat2str(round(1e3*fm.rsT20(1:4), 3)));
    save(fullfile(res, ['gridconv_' PARAMSET '.mat']), 'R', 'G');
end
