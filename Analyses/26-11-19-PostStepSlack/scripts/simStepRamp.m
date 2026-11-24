function [o, vt] = simStepRamp(pBase, step, v, opt)
% simStepRamp  One Beard2022-style cycle on the model: hold at 1.10 Lo, fast step
% of STEP (Lo), ramp at V (ML/s) down to 0.985 Lo, then hold.
%
%   [o, vt] = simStepRamp(pBase, step, v)            8 mM
%   opt.atp (mM, 8), opt.tStep (s, 0.6e-3: the measured step takes ~0.6 ms),
%   opt.tWarm (s, 2), opt.Lend (Lo, 0.985), opt.tHold (s, 0.02),
%   opt.set (struct of parameter overrides applied before getParams)
% o is evaluateModel's output with o.tc = time from the step onset (s).

if nargin < 4, opt = struct(); end
d = struct('atp', 8, 'tStep', 0.6e-3, 'tWarm', 2, 'Lend', 0.985, 'tHold', 0.02, 'set', struct());
fn = fieldnames(d);  for i = 1:numel(fn), if ~isfield(opt, fn{i}), opt.(fn{i}) = d.(fn{i}); end, end
L0 = 1.10;  tR = (L0 - step - opt.Lend)/v;
t  = [-opt.tWarm; 0; opt.tStep; opt.tStep + tR; opt.tStep + tR + opt.tHold];
vv = [0; -step/opt.tStep; -v; 0; 0];
vt = [t, vv];
p = pBase;  p.MgATP = opt.atp;
fs = fieldnames(opt.set);  for i = 1:numel(fs), p.(fs{i}) = opt.set.(fs{i}); end
p.SL0 = 2*L0;  p.Slim_l = 1.4;  p.Slim_r = 2.3;
p = getParams(p, [], true);
p.Velocity = vv;
[~, o] = evaluateModel(str2func(p.modelFcn), t, p);
o.tc = o.t(:);  o.tRampEnd = opt.tStep + tR;  o.tStep = opt.tStep;
end
