%% ProbeDenseSegments.m -- does skipping the dense grid on the 94 s init hold change the chunks?
% "dense everywhere" and "dense only segments <= 5 s" differ only on the init
% hold, yet CompareOutputGrids moved rsT20's cost 0.052 -> 0.035. Compare the
% state handed to the chunks and one chunk's trace under both settings.
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(here, '..', '..')));
P = getParams(loadParams('refit9_opt'), [], true, false);
D = load(fullfile('data', P.velocitytableonfile)); vt = D.velocitytable; vt(1,1) = -20;
P = rmfield(P, 'PU0'); P.Velocity = 0;
V = {'dense everywhere', struct(); 'dense <= 5 s', struct('odeOutputMaxSegment', 5)};
for v = 1:2
    Q = P; f = fieldnames(V{v,2}); for k = 1:numel(f), Q.(f{k}) = V{v,2}.(f{k}); end
    t0 = tic; [~, oi] = evaluateModel(str2func(Q.modelFcn), [-10 vt(3,1)], Q); w = toc(t0);
    R(v).PU0 = oi.PU(end,:); R(v).piv = oi.params.LXBpivot; R(v).n = numel(oi.t); R(v).w = w; %#ok<SAGROW>
    fprintf('%-18s init: %6d stored points, %.1f s, pivot %.6f\n', V{v,1}, numel(oi.t), w, oi.params.LXBpivot);
end
fprintf('max |PU0 difference| = %.3e   (max |PU0| = %.3e)\n', max(abs(R(1).PU0 - R(2).PU0)), max(abs(R(1).PU0)));
save(fullfile(here, 'results', 'probe_dense_segments.mat'), 'R');
