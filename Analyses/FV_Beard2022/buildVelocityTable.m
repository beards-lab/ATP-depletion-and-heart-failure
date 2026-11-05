function vt = buildVelocityTable(t, L, tol, tWarm)
% buildVelocityTable  Piecewise-constant velocity table that reproduces a
% measured motor length record to within tol.
%
%   vt = buildVelocityTable(t, L, tol, tWarm)
%     t      (ms) record time, L (Lo) measured length
%     tol    (Lo) max deviation of the piecewise-linear fit (default 5e-4)
%     tWarm  (s)  isometric warm-up prepended at L(1) (default 2)
%
%   vt = [t(s), v(ML/s), v_um(bookkeeping), L(ML)], the evaluateModel /
%   Example/BuildExampleProtocol.m convention: row i's v holds over
%   [t_i, t_i+1]; the record starts at t = 0 and the warm-up at -tWarm.
%
% Breakpoints come from a Ramer-Douglas-Peucker simplification of the record
% (iterative, so it does not recurse 86k deep). Within each segment the speed is
% the chord slope, so the integrated length passes exactly through the measured
% length at every breakpoint.

if nargin < 3 || isempty(tol), tol = 5e-4; end
if nargin < 4 || isempty(tWarm), tWarm = 2; end
t = t(:); L = L(:);
% light de-noise so quantisation steps (1.4e-4 Lo) do not create breakpoints
Ls = movmean(L, 5);
keep = false(size(t)); keep([1 end]) = true;
stack = [1 numel(t)];
while ~isempty(stack)
    a = stack(end, 1); b = stack(end, 2); stack(end, :) = [];
    if b - a < 2, continue; end
    idx = (a+1:b-1)';
    chord = Ls(a) + (Ls(b) - Ls(a)) * (t(idx) - t(a)) / (t(b) - t(a));
    [dmax, j] = max(abs(Ls(idx) - chord));
    if dmax > tol
        m = idx(j); keep(m) = true;
        stack = [stack; a m; m b]; %#ok<AGROW>
    end
end
k = find(keep);
tb = t(k) / 1000; Lb = Ls(k);
v = [diff(Lb) ./ diff(tb); 0];
vt = [[-tWarm; tb], [0; v], [0; v] * 2, [Lb(1); Lb]];
end
