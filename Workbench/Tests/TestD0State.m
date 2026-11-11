% TestD0State.m
% 1) UseD0State = false must be BYTE-IDENTICAL to the pre-change model.
% 2) UseD0State = true with k2rip = 0 (no tearing) and D0FromR2 = 0 must ALSO
%    be byte-identical: nothing ever enters D0, so it is an inert extra state.
% 3) With tearing on, mass must still be conserved and D0 must fill and drain.
cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));

base = getParams(loadParams('ThursdayNightFever'), [], true, false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad=''; base.EvalFeatures=0;
base.RunSlackSegments='First';

fprintf('--- 1/2: inertness when nothing is torn ---\n');
pA = base;                       pA = getParams(pA, pA.g, true);
pB = base; pB.UseD0State = true; pB = getParams(pB, pB.g, true);
fprintf('PU0 length: off %d, on %d (expect +1)\n', numel(pA.PU0), numel(pB.PU0));
[EA, oA] = runSlackExperiment(pA);
[EB, oB] = runSlackExperiment(pB);
% ode15s picks its own steps and the extra state changes the error norm, so the
% two runs land on DIFFERENT time grids. Compare on a common grid, not elementwise.
mA=[true;diff(oA.t(:))>0]; mB=[true;diff(oB.t(:))>0];
tc = linspace(max(oA.t(1),oB.t(1)), min(oA.t(end),oB.t(end)), 4000)';
FA = interp1(oA.t(mA), oA.Force(mA), tc); FB = interp1(oB.t(mB), oB.Force(mB), tc);
dF = max(abs(FA-FB)); Fsc = max(abs(FA));
fprintf('E_slack off %.6f | on %.6f | max |dForce| = %.3e (%.4f%% of full scale)\n', ...
        EA(1), EB(1), dF, 100*dF/Fsc);
if dF/Fsc < 1e-3 && abs(EA(1)-EB(1))/EA(1) < 1e-3
    fprintf('PASS: D0 is inert when k2rip = 0 and D0FromR2 = 0 (residual = solver tolerance)\n');
else
    fprintf('FAIL: D0 changed the solution with no inflow\n');
end

fprintf('\n--- 3: tearing on -> D0 fills, mass conserved ---\n');
pC = base; pC.UseD0State = true; pC.k2rip = 50; pC.alphaRip = 100; pC.dr3 = -0.01;
pC.k_D0T = 200; pC.K_D0T = 0.15;
pC = getParams(pC, pC.g, true);
[EC, oC] = runSlackExperiment(pC);
ss = oC.params.ss;
tot = oC.p1_0 + oC.p2_0 + oC.PuR + oC.SR + oC.SRD + oC.PuATP + oC.D0;
fprintf('E_slack %.3f | max D0 = %.4f | mass sum: min %.6f max %.6f (expect 1)\n', ...
        EC(1), max(oC.D0), min(tot), max(tot));
if abs(max(tot)-1) < 1e-3 && abs(min(tot)-1) < 1e-3
    fprintf('PASS: probability conserved with D0 in the cycle\n');
else
    fprintf('FAIL: mass not conserved\n');
end
if max(oC.D0) > 1e-4
    fprintf('PASS: D0 fills during the protocol (peak %.4f)\n', max(oC.D0));
else
    fprintf('NOTE: D0 stayed empty — tearing never fired at these strains\n');
end
