function f = protocolFeatures(P, t, F)
% PROTOCOLFEATURES  Per-protocol features for the PDE-vs-moments comparison.
% Rates use the ktr80 estimator ln4/(t80-t20); see CompareMomentsVsPDE.
% A few features per protocol, chosen to mirror what is scored or what a
% ventricle would feel. Rates use the ktr80 estimator ln4/(t80-t20).
tp = P.T(2);                                    % end of the 1 s activation
f.F_act50ms = interp1(t, F, 0.05);
f.F_pre = interp1(t, F, tp - 1e-4);
switch true
    case strcmp(P.name, 'iso')
        f.t50 = t(find(F >= 0.5*F(end), 1)) ;
    case any(strcmp(P.name, {'ktr', 'slack'}))
        tr = P.T(end-1);                        % start of the final hold
        f.F_min_rel = min(F(t > tp & t < tr));
        f.F_peak_post = max(F(t >= tr & t < tr + 0.05));
        f.k80 = rate80(t, F, tr);
        f.F_end = F(end);
    case strcmp(P.name, 'beat')
        f.F_max = max(F(t > tp)); f.F_min = min(F(t > tp));
        % length (ML) by integrating the commanded velocity
        L = cumtrapz(t, interp1(P.T(1:end-1), P.V, t, 'previous', 'extrap'));
        f.work = -trapz(L(t > tp), F(t > tp));             % kPa*ML, net shortening work
        f.F_mid_eject = interp1(t, F, tp + 0.1);
    case strcmp(P.name, 'stairs')
        for k = 1:4
            ts = P.T(2*k); te = P.T(2*k+1);
            f.(sprintf('peak%d', k)) = max(F(t >= ts & t <= ts + 0.02));
            f.(sprintf('hold%d', k)) = interp1(t, F, te - 1e-4);
        end
    otherwise   % fv
        f.F_ramp_end = interp1(t, F, P.T(3) - 1e-4);
        f.F_ramp_mid = interp1(t, F, (P.T(2) + P.T(3))/2);
end
end

function k = rate80(t, F, t0)
w = t >= t0; tt = t(w) - t0; FF = F(w);
[Fb, ib] = min(FF(1:max(2, find(tt < 0.03, 1, 'last')))); Fe = FF(end);
A = Fe - Fb; if A <= 0, k = NaN; return; end
i20 = find(FF(ib:end) >= Fb + 0.2*A, 1) + ib - 1; i80 = find(FF(ib:end) >= Fb + 0.8*A, 1) + ib - 1;
if isempty(i20) || isempty(i80) || i80 == i20 || i20 < 2, k = NaN; return; end
% interpolate the crossings: on a 0.5 ms output grid one sample is ~2 % of t80-t20
cross = @(i, lev) tt(i-1) + (lev - FF(i-1)) * (tt(i) - tt(i-1)) / (FF(i) - FF(i-1));
k = log(4)/(cross(i80, Fb + 0.8*A) - cross(i20, Fb + 0.2*A));
end
