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
    feats = struct('ktrP_k63', NaN, 'ktrP_F20', NaN);

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
end
