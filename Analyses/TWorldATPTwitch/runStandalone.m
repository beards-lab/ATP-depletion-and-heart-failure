function runStandalone(setName, idx)
%RUNSTANDALONE  Pace standalone T-World runs of a set (standaloneConditions) and cache each one
%   as results/sa_<label>.mat: the last two beats (V, Ca_i, Ca_SR, Na_i, Land_Ta, Jup, IKATP),
%   APD90, and the final state X (the 1 Hz steady state handed to the coupled runs).
%   IDX selects conditions (default all); separate processes can take disjoint IDX.
here = fileparts(mfilename('fullpath'));
C = standaloneConditions(setName);
if nargin < 2 || isempty(idx), idx = 1:numel(C); end
twPaths(C(1).flavor);
for i = idx
    c = C(i);  f = fullfile(here, 'results', ['sa_' c.label '.mat']);
    if exist(f, 'file'), fprintf('%s cached\n', c.label); continue; end
    param = struct('bcl', 1000);
    switch c.model
        case 'shock_ATP', param.model = @model_TWorld_shock_ATP;  param.extraParams = [0 310 c.pH c.ATP c.ADP c.gK];
        case 'plain',     param.model = @model_TWorld;
        case 'lowAtp',    param.model = @model_TWorld;  param.ATP_i = c.ATP;  param.ADP_i = c.ADP;  param.IKATP_Multiplier = c.gK;
    end
    X0 = getStartingState('TW_endo');
    tic;  [time, X] = modelRunner(X0, [], param, c.beats, c.beats - 2);  wall = toc;
    r = struct('c', c, 'wall', wall, 'Xend', X{end}(end, :));
    for k = 1:2
        s = getCurrentsStructure(time(k), X(k), param, 0);
        b = struct('t', s.time, 'V', s.V, 'Ca', 1e3*s.Ca_i, 'CaSR', s.Ca_SR, 'Na', s.Na_i, 'Ta', s.Land_Ta, ...
            'Jup', s.Jup, 'IKATP', s.IKATP, 'INaCa', s.INaCa, 'ICaL', s.ICaLtot, ...
            'APD90', apd90(s.time, s.V));
        for g = ["f_SERCA_ATP" "f_NaK_ATP" "fATP_SERCA" "fATP_NaK" "fATP_PCa" "fATP_KATP"]
            if isfield(s, g), b.(g) = s.(g)(1); end
        end
        r.beat(k) = b;
    end
    save(f, 'r');
    fprintf('%-24s %5.0f s | Ca %.4f..%.4f uM | APD90 %.1f ms\n', c.label, wall, min(r.beat(2).Ca), max(r.beat(2).Ca), r.beat(2).APD90);
end
end

function a = apd90(t, V)
% APD90 (ms): from the maximal upstroke slope to the first return below 90 % repolarisation
% (DataReporter.getAPD needs the Image Processing Toolbox, which is not installed).
[~, iu] = max(diff(V)./diff(t));  [Vp, ip] = max(V);
th = Vp - 0.9*(Vp - V(1));
j = find(V(ip:end) < th, 1) + ip - 1;
a = NaN;  if ~isempty(j), a = t(j) - t(iu); end
end
