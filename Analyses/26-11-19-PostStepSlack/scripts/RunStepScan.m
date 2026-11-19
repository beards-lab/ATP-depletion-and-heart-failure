%% RunStepScan  Post-step slack on the Beard2022 isovelocity cycle (#77): where it starts and what carries it.
%
% Synthetic cycle (simStepRamp): hold 1.10 Lo, step in 0.6 ms, ramp to 0.985 Lo. 8 mM.
% Incumbent: robust1_seed with the fit-A k2 gate (as RunFVBeard2022).
%   A  step scan: step 0.01-0.066 Lo x ramp 1, 2, 4, 6 ML/s -> slack time, lowest force
%   B  levers at the two measured cycles (1 ML/s / 0.021 Lo, 4 ML/s / 0.054 Lo), one at a time
% Output: results/step_scan.mat, results/step_scan.log

here = fileparts(mfilename('fullpath'));  ana = fileparts(here);  root = fileparts(fileparts(ana));
addpath(genpath(fullfile(root, 'Model')));  addpath(genpath(fullfile(root, 'Auxiliary')));  addpath(fullfile(root, 'params'));  addpath(here);
if ~isfolder(fullfile(ana, 'results')), mkdir(fullfile(ana, 'results')); end
logf = fullfile(ana, 'results', 'step_scan.log');  if isfile(logf), delete(logf); end
diary(logf);
pP = basePars();
fprintf('incumbent: kSE %.4g, ekSE %.4g, A0 %.3g, tau_reg %.3g s, v_ref_reg %.3g\n', pP.kSE, pP.ekSE, pP.A0, pP.tau_reg, pP.v_ref_reg);

%% A  step scan
steps = [0.010 0.021 0.030 0.042 0.054 0.066];  vv = [1 2 4 6];
A = struct();
fprintf('\n== A. step scan, 8 mM: slack time in the ramp (ms) / lowest force from 1 ms after the step to the ramp end (Fss) / force at the ramp end (Fss)\n');
for iv = 1:numel(vv)
    for is = 1:numel(steps)
        o = simStepRamp(pP, steps(is), vv(iv));  m = meas(o);
        A.slack(iv, is) = m.slack;  A.Fmin(iv, is) = m.Fmin;  A.Fend(iv, is) = m.Fend;  A.LSE0 = m.LSE0;  A.Fss = m.Fss;
        fprintf('v %g ML/s, step %.3f Lo: slack %5.2f ms, min %.3f, end %.3f\n', vv(iv), steps(is), m.slack, m.Fmin, m.Fend);
    end
end
fprintf('series-element extension at Fss: %.1f nm of sarcomere length = %.4f Lo; Fss %.1f kPa\n', 1e3*A.LSE0, A.LSE0/2, A.Fss);

%% B  levers, one at a time
n0 = pP.ekSE;  L0 = A.LSE0;  Fss = A.Fss;
lev = {'incumbent', struct()};
for n = [1.5 2 3 4]                                   % series element with a toe: same force and same stiffness at Fss
    Ln = L0*n/n0;  lev(end+1, :) = {sprintf('SE exponent %.3g (same stiffness at Fss)', n), struct('ekSE', n, 'kSE', Fss^(1/n)/Ln)}; %#ok<SAGROW>
end
for n = [2 3]                                         % same force and same extension at Fss (softer at Fss)
    lev(end+1, :) = {sprintf('SE exponent %.3g (same extension at Fss)', n), struct('ekSE', n, 'kSE', Fss^(1/n)/L0)}; %#ok<SAGROW>
end
lev(end+1, :) = {'gate fully open (A0 = 1)', struct('A0', 1)};
lev(end+1, :) = {'gate slower (tau_reg x 100)', struct('tau_reg', 100*pP.tau_reg)};
rn = {'kah', 'kamh', 'ka', 'kd', 'k1', 'k_1', 'k2'};  s2 = struct();  for i = 1:numel(rn), if isfield(pP, rn{i}) && isnumeric(pP.(rn{i})), s2.(rn{i}) = 2*pP.(rn{i}); end, end
lev(end+1, :) = {'cycle rates x 2', s2};
lev(end+1, :) = {'attachment ka x 2', struct('ka', 2*pP.ka)};
lev(end+1, :) = {'detachment k2 x 2', struct('k2', 2*pP.k2)};
cyc = [0.021 1; 0.054 4; 0.066 6];
B = struct('name', {}, 'set', {}, 'm', {});
fprintf('\n== B. levers. Columns per cycle (step Lo / ramp ML/s): slack ms, lowest force (Fss), force at ramp end (Fss). Data minima: 0.45 (1 ML/s), 0.09 (4 ML/s), 0.03 (6 ML/s)\n');
for i = 1:size(lev, 1)
    B(i).name = lev{i, 1};  B(i).set = lev{i, 2};  s = '';
    for c = 1:size(cyc, 1)
        try
            o = simStepRamp(pP, cyc(c, 1), cyc(c, 2), struct('set', lev{i, 2}));  m = meas(o);
        catch ME
            m = struct('slack', NaN, 'Fmin', NaN, 'Fend', NaN, 'Fss', NaN, 'LSE0', NaN);  fprintf('  (%s failed: %s)\n', lev{i, 1}, ME.message);
        end
        B(i).m(c) = m;  s = [s sprintf(' | %5.2f %6.3f %6.3f', m.slack, m.Fmin, m.Fend)]; %#ok<AGROW>
    end
    fprintf('%-42s Fss %5.1f%s\n', B(i).name, B(i).m(1).Fss, s);
end
diary off;
save(fullfile(ana, 'results', 'step_scan.mat'), 'A', 'B', 'steps', 'vv', 'cyc');

function m = meas(o)
t = o.t(:);  F = o.Force(:);  dt = gradient(t);
m.Fss = interp1(t, F, -0.001);  m.LSE0 = interp1(t, o.LSE(:), -0.001);
k = t > 0 & t < o.tRampEnd;
m.slack = 1e3*sum(dt(k & F < 0.01*m.Fss));
k2 = t > o.tStep + 1e-3 & t < o.tRampEnd;  if ~any(k2), k2 = k; end
m.Fmin = min(F(k2))/m.Fss;  m.Fend = interp1(t, F, o.tRampEnd - 1e-4)/m.Fss;
end
