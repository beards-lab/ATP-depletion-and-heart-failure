function feats = extractKtrAttributes(t, F_norm, velocitytable)
%EXTRACTKTRATTRIBUTES Redevelopment rate of the ktr (Brenner) protocol.
%
%   FEATS = extractKtrAttributes(T, F_NORM, VELOCITYTABLE)
%
%   The ktr manoeuvre ends on a release back to the reference length, so force
%   redevelops isometrically from ~zero. The window is the final hold:
%   VELOCITYTABLE(end-1,1) .. VELOCITYTABLE(end,1).
%
%   Returns:
%     feats.ktrP_k63  - 1/(t63 - t_on), the model-free crossing rate (s^-1)
%     feats.ktrP_F20  - normalised force 20 ms after the release ends
%
%   WHY THE CROSSING RULE AND NOT A FREE-AMPLITUDE FIT. The window is only
%   ~58 ms and the trace has reached ~0.93 F_iso by its end, i.e. it is
%   truncated before the plateau. Letting A float on a span this short returns
%   ~86 % of the true rate (Analyses/RestretchVsKtrRecovery/conclusions.md,
%   "Estimator caveat for the ktr protocol"). The muscle is back at the
%   pre-release length, so the plateau is KNOWN -- it is F_iso, i.e. 1.0 in
%   normalised units. Fixing P = 1 makes the estimate span-insensitive because
%   t63 ~ 23 ms sits well inside the window.
%
%   The crossing rule is byte-identical to the one the analysis folder uses
%   (Analyses/RestretchVsKtrRecovery/recoveryWindows.m, crossRate), so model
%   and data go through the same estimator.
%
%   F_NORM must already be normalised by the pre-shortening isometric force.

    SMOOTH = 0.002;                          % median span for the crossing (s)
    feats = struct('ktrP_k63', NaN, 'ktrP_F20', NaN, 'ktrP_k80', NaN);

    t = t(:); F_norm = F_norm(:);
    keep = [true; diff(t) > 0];              % model output can repeat timestamps
    t = t(keep); F_norm = F_norm(keep);

    t_on  = velocitytable(end-1, 1);
    t_end = velocitytable(end,   1);
    m = t >= t_on & t <= t_end;
    if nnz(m) < 20; return; end
    tw = t(m) - t_on;  yw = F_norm(m);

    n = max(3, round(SMOOTH / median(diff(tw))));
    y = movmedian(yw, n);                    % P = 1.0 (F_iso), B = 0

    below = find(y < 0.05);
    if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(y >= 1 - exp(-1) & (1:numel(y))' > i_on, 1);
    if ~isempty(i63) && tw(i63) > tw(i_on)
        feats.ktrP_k63 = 1 / (tw(i63) - tw(i_on));
    end
    if tw(end) >= 0.020
        feats.ktrP_F20 = interp1(tw, yw, 0.020, 'linear');
    end

    % ktrP_k80 -- the SCORED rate. k63 measures from an onset (the last
    % sample below 0.05), and the onset of a noisy trace near zero is exactly
    % where a crossing is least well defined; it also carries the dead time.
    % ln(4)/(t80 - t20) needs no onset at all -- for P(1-exp(-k(t-t0))) the
    % t0 cancels -- and reads only the band where the rise is steepest. P is
    % still pinned at F_iso = 1 (the window is truncated; never float it).
    % 5 ms moving median: on a monotone rise it returns the centre sample,
    % so it removes noise without shifting either crossing.
    nx = max(3, round(0.005 / median(diff(tw))));
    zx = movmedian(yw, nx);
    i0 = find(zx < 0.05, 1, 'last'); if isempty(i0); i0 = 1; end
    t20 = crossUp(tw(i0:end), zx(i0:end), 0.20);
    t80 = crossUp(tw(i0:end), zx(i0:end), 0.80);
    if isfinite(t20) && isfinite(t80) && t80 > t20
        feats.ktrP_k80 = log(4) / (t80 - t20);
    end
end

function tx = crossUp(tt, z, lev)
% first upward crossing of LEV, linearly interpolated; NaN if none
    tx = NaN;
    i = find(z >= lev, 1);
    if isempty(i); return; end
    if i == 1; tx = tt(1); return; end
    z1 = z(i-1); z2 = z(i);
    if z2 == z1; tx = tt(i); else; tx = tt(i-1) + (lev - z1)/(z2 - z1)*(tt(i) - tt(i-1)); end
end
