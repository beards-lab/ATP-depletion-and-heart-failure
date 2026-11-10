function [F, tSim] = simBeard(pBeard, e, rec)
% simBeard  Model force (kPa, single) on one Beard2022 record's own time grid,
% with the parameter(s) e.name multiplied by e.f (e.name empty = the incumbent).
p = pBeard;
if ~isempty(e.name)
    p.mods = cellstr(e.name);
    p = getParams(p, e.f, false, true);              % applies the modifier, clears mods
end
vt = rec.vt;
p.MgATP = rec.atp;
p.SL0 = 2 * vt(1, 4); p.Slim_l = 1.4; p.Slim_r = 2.3;
p = getParams(p, [], true);
p.Velocity = vt(:, 2);
t0 = tic;
[~, o] = evaluateModel(str2func(p.modelFcn), vt(:, 1), p);
tq = rec.t / 1000;
[tu, iu] = unique(o.t(:));
assert(tu(end) >= tq(end) - 1e-6, 'simulation stopped at %.3f s', tu(end));
F = single(interp1(tu, o.Force(iu), tq, 'linear'));
assert(all(isfinite(F)), 'non-finite force');
tSim = toc(t0);
end
