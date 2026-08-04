function [o, k63] = seedWindow(p, PU0, LXBpivot, t0, dur, B, P)
%SEEDWINDOW Replay one recovery window from a captured state vector.
%
%   [O, K63] = seedWindow(P, PU0, LXBPIVOT, T0, DUR, B, P_PLATEAU)
%
%   Both recovery windows are holds, so a window is fully specified by its
%   initial state. This runs a pure hold from a supplied PU0 and measures k63
%   with the SAME crossing rule as recoveryWindows.m:118-125.
%
%   SEEDING CONTRACT (mirrors runSlackExperiment.m:225-233):
%     1. getParams(p, p.g, true) FIRST — it rebuilds the strain grid and PU0.
%     2. THEN assign p.PU0 and p.LXBpivot; evaluateModel calls
%        getParams(..., updateInit=false), so they survive.
%   getParams:571 silently rebuilds PU0 from defaults if its length does not
%   match Ns*ss + 7 + UseRegistrationAvailability + UseD0State, which would
%   discard the seed with no warning — so the length is asserted here.
%
%   B, P_PLATEAU are the baseline/plateau used for the k63 crossing. Pass the
%   captured window's own values so replays are comparable to the original.

    p = getParams(p, p.g, true);
    n = p.NumberOfStates*p.ss + 7 + p.UseRegistrationAvailability + p.UseD0State;
    assert(numel(PU0) == n, 'seedWindow:len', ...
           'PU0 has %d elements, this paramset needs %d — getParams would silently rebuild it.', ...
           numel(PU0), n);

    p.Velocity  = 0;                       % single hold segment
    p.PU0       = PU0(:)';
    p.LXBpivot  = LXBpivot;
    p.BreakOnODEUnstable = 0;
    p.PlotEachSeparately = 0;

    [~, o] = evaluateModel(str2func(p.modelFcn), [t0 t0+dur], p);

    if nargin < 6 || isempty(B); B = o.Force(1); end
    if nargin < 7 || isempty(P)
        tt = o.t(:); P = median(o.Force(tt >= tt(end)-0.15*(tt(end)-tt(1))));
    end
    k63 = crossRate(o.t(:), o.Force(:), B, P, 0.002);
end

function k = crossRate(t, F, B, P, SM)
% identical rule to recoveryWindows.m:118-125
    keep = [true; diff(t) > 0]; t = t(keep); F = F(keep);
    n  = max(3, round(SM/median(diff(t))));
    y  = (movmedian(F, n) - B)/max(P - B, eps);
    below = find(y < 0.05);
    if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(y >= 1-exp(-1) & (1:numel(y))' > i_on, 1);
    if isempty(i63) || t(i63) <= t(i_on); k = NaN; else; k = 1/(t(i63)-t(i_on)); end
end
