function P = fv2Model(params0)
%FV2MODEL Model P/P0 on the FV2 ramp, read exactly as fv2Read reads the data.
%   Uses the 04/10 protocol table of the condition being simulated (the
%   manoeuvre is identical on both days). Same set-up as runKtrExperiment's
%   file-driven path: start length from the table, getParams, evaluateModel.
suf = ''; if params0.MgATP < 8; suf = '_2mM'; end
S = load(fullfile('data', sprintf('protocol_04_10_2026_velocitytable_FV2%s.mat', suf)));
vt = S.velocitytable;
% The protocol opens with a ~5 s isometric warm-up; 1 s reaches the same state
% (P/P0 0.1524 vs 0.152 at g = 0.54, 2 mM) at half the cost.
vt(1,1) = max(vt(1,1), -1);
params = params0;
params.SL0 = 2 * vt(1, 4);
params = getParams(params, params.g, true);
params.Velocity = vt(:, 2);
[~, out] = evaluateModel(str2func(params.modelFcn), vt(:, 1), params);
k  = makeMonotonous(out.t(:));
tm = out.t(k); Fm = out.Force(k);
tq = (0:1e-4:0.12)';                     % read on a fixed 10 kHz grid
P  = fv2Read(tq, interp1(tm, Fm, tq, 'linear'), vt);
end
