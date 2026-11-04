% srxProtocol  run one candidate SRX-discriminating protocol (A activation step,
% B release-hold-restretch series, C paired release-restretch, D isometric SL) on params q.
% See RunSrxProtocols.m.
function r = srxProtocol(q, kind, x)
    mf = str2func(q.modelFcn); q.SL0 = 2.2; dL = 0.15;       % release amplitude (ML)
    switch kind
        case 'A'   % relaxed steady state, then attachment switched on
            qr = q; qr.ka = 0; qr.Velocity = [0; 0]; qr = getParams(qr, qr.g, true, false);
            [~, o0] = evaluateModel(mf, [-5 0], qr);
            qa = q; qa.Velocity = [0; 0]; qa = getParams(qa, qa.g, true, false);
            qa.PU0 = o0.PU(end, :);
            [~, o] = evaluateModel(mf, [0 0.4], qa);
            F0 = o.Force(1); Fi = o.Force(end); dF = Fi - F0;
            r.t50 = iCross(o.t, o.Force, F0 + 0.5*dF); r.t90 = iCross(o.t, o.Force, F0 + 0.9*dF);
            r.SRXrest = o0.SR(end) + o0.SRD(end); r.Frest = F0; r.Fiso = Fi;
        case 'B'   % release, hold x, restretch in 2 ms, hold 0.4 s
            vt = [-3, 0; 0, -dL/1e-3; 1e-3, 0; 1e-3 + x, dL/2e-3; 3e-3 + x, 0; 0.403 + x, 0];
            q.Velocity = vt(:, 2); q = getParams(q, q.g, true, false);
            [~, o] = evaluateModel(mf, vt(:, 1), q);
            tr = 3e-3 + x; r = iKtr(o, tr, vt(end, 1)); r.Fpre = iAt(o, -1e-3);
            m = o.t >= 0 & o.t <= tr; r.SRXmax = max(o.SR(m) + o.SRD(m));
            r.SRXatRestretch = iAt(o, tr, o.SR + o.SRD);
        case 'C'   % two ktr manoeuvres (T = 20 ms) separated by x after the first restretch
            T = 20e-3; t1 = 0; r1 = t1 + 1e-3 + T + 2e-3; t2 = r1 + x;
            vt = [-3, 0; t1, -dL/1e-3; t1 + 1e-3, 0; t1 + 1e-3 + T, dL/2e-3; r1, 0; ...
                  t2, -dL/1e-3; t2 + 1e-3, 0; t2 + 1e-3 + T, dL/2e-3; t2 + 3e-3 + T, 0; t2 + 0.403 + T, 0];
            q.Velocity = vt(:, 2); q = getParams(q, q.g, true, false);
            [~, o] = evaluateModel(mf, vt(:, 1), q);
            a = iKtr(o, r1, min(t2, r1 + 0.4)); b = iKtr(o, t2 + 3e-3 + T, vt(end, 1));
            r.k1 = a.k63; r.k2 = b.k63; r.ratio = b.k63 / a.k63;
            r.F2pre = iAt(o, t2 - 1e-4) / iAt(o, -1e-3);
            r.SRXpre2 = iAt(o, t2 - 1e-4, o.SR + o.SRD);
        case 'D'   % isometric at SL x
            q.SL0 = x; q.Velocity = [0; 0]; q = getParams(q, q.g, true, false);
            [~, o] = evaluateModel(mf, [-3 0], q);
            r.F = o.Force(end); r.SRX = o.SR(end) + o.SRD(end); r.att = o.p1_0(end) + o.p2_0(end);
    end
    r.kind = kind; r.x = x;
    if ismember(kind, {'A', 'B', 'C'})       % thinned trace for figures
        tq = linspace(o.t(find(o.t >= -0.01, 1)), o.t(end), 1500)'; k = makeMonotonous(o.t(:));
        r.tr = [tq, interp1(o.t(k), o.Force(k), tq), interp1(o.t(k), o.SR(k) + o.SRD(k), tq), ...
                interp1(o.t(k), o.SL(k), tq)];
    end
end

function r = iKtr(o, tr, tend)
    % k63 of the redevelopment after a restretch that ends at tr, timed from the
    % minimum that follows the restretch force spike (searched over 15 ms)
    m = find(o.t >= tr & o.t <= tr + 15e-3); [Fmin, j] = min(o.Force(m)); tmin = o.t(m(j));
    Fend = iAt(o, tend - 1e-3); lvl = Fmin + 0.63 * (Fend - Fmin);
    m2 = o.t >= tmin; t63 = iCross(o.t(m2), o.Force(m2), lvl);
    r.k63 = 1 / (t63 - tmin); r.Fmin = Fmin; r.Fend = Fend; r.tmin = tmin - tr;
    r.F20 = (iAt(o, tr + 20e-3) - Fmin) / (Fend - Fmin);
end

function t = iCross(tt, F, lvl)
    i = find(F >= lvl, 1);
    if isempty(i) || i == 1; t = NaN; return; end
    t = tt(i-1) + (lvl - F(i-1)) * (tt(i) - tt(i-1)) / (F(i) - F(i-1));
end

function v = iAt(o, t, y)
    if nargin < 3; y = o.Force; end
    k = makeMonotonous(o.t(:)); v = interp1(o.t(k), y(k), t, 'linear', 'extrap');
end
