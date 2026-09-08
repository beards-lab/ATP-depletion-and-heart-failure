function E = rsEstimators(t, y, smoothIt, opts)
%RSESTIMATORS  Seven ways to put a RATE on the post-re-stretch force recovery.
%
%   E = rsEstimators(T, Y, SMOOTHIT)
%
%   T, Y     the post-re-stretch window: T relative to anything (only
%            differences are used), Y in kPa. This is the same window
%            extractSlackAttributes uses, [velocity_segment(4) velocity_segment(5)].
%   SMOOTHIT true on the data path (a 5 ms movmedian is used to LOCATE the
%            anchor minimum and the level crossings; every fit is on raw Y).
%
%   All seven estimators share ONE preprocessing convention, which is the
%   "rescale it completely" option: locate the post-re-stretch minimum B
%   (the recovery starts there, not at the ramp end), take the plateau P as
%   the median of the last 15 % of the window, and work on the shape
%
%        z(t) = (y - B) / (P - B)        rising 0 -> 1, dimensionless
%
%   so no estimator can absorb a force-scale error. They differ only in WHICH
%   PART OF z THEY LOOK AT and in whether they are a fit or a functional.
%
%     .k_fit    CURRENT (extractSlackAttributes/fitRestretchRecovery): least
%               squares of z against 1-exp(-k*(t-t0)) with k and t0 both free,
%               over the whole window.
%     .k_t0z    the same with t0 == 0: the clock starts at the valley, so the
%               toe is charged to the rate instead of to a free delay.
%     .k_band   log-linear regression of log(1-z) on t over the BAND
%               z in [BAND(1) BAND(2)] only -- the "top of the rise": the toe
%               and the plateau are excluded by construction, no t0 needed.
%     .k_t50    model-free: ln(2) / (time from the valley to z = 0.5).
%     .k_2080   model-free, toe- and plateau-free: ln(4) / (t80 - t20).
%     .k_mrt    1 / MRT, MRT = INTEGRAL of (1 - z) dt over the window, with the
%               finite-window correction (1-exp(-k T))/k = MRT solved for k.
%               MRT is a LINEAR functional of the trace, so it is the only
%               one of the seven that commutes with trace averaging.
%     .k_area   MRT without the finite-T correction (1/integral) -- reported
%               to show the correction matters.
%
%   Also returned: .B .A .t0 .r2 (of k_fit), .tOn (z first leaves 5 %),
%   .nBand (points in the band), and .MRT.
%
%   Everything is NaN-safe: an estimator that cannot be formed returns NaN
%   rather than a wrong number.

if nargin < 3 || isempty(smoothIt); smoothIt = true; end
if nargin < 4; opts = struct(); end
if ~isfield(opts,'vallWin'); opts.vallWin = 0.080; end   % as extractSlackAttributes
if ~isfield(opts,'band');    opts.band    = [0.15 0.90]; end
if ~isfield(opts,'tailFrac');opts.tailFrac= 0.15;    end
if ~isfield(opts,'Aext');    opts.Aext    = [];      end   % external amplitude (kPa): bypasses the plateau estimate
if ~isfield(opts,'crossWin');opts.crossWin= 0.005;   end   % smoother width for the level crossings (s)

E = struct('k_fit',NaN,'k_t0z',NaN,'k_band',NaN,'k_t50',NaN,'k_2080',NaN, ...
           'k_mrt',NaN,'k_area',NaN,'B',NaN,'A',NaN,'t0',NaN,'r2',NaN, ...
           'tOn',NaN,'nBand',NaN,'MRT',NaN,'T',NaN);

t = t(:); y = y(:);
ok = isfinite(t) & isfinite(y); t = t(ok); y = y(ok);
if numel(t) < 20 || ~(t(end) > t(1)); return; end

dt = median(diff(t));
if smoothIt
    yl = smoothdata(y, 'movmedian', max(3, round(opts.crossWin/dt)));
else
    yl = y;
end

% ---- anchor: the minimum inside the first vallWin of the window ----------
nWin = nnz(t <= t(1) + opts.vallWin);
if nWin < 3; nWin = numel(t); end
[B, iB] = min(yl(1:nWin));

tt = t(iB:end) - t(iB);
yy = y(iB:end) - B;
yls = yl(iB:end) - B;
if numel(tt) < 20 || ~(tt(end) > 0); return; end

% ---- amplitude: plateau over the last tailFrac of what is left -----------
if isempty(opts.Aext)
    P = median(yy(tt >= tt(end) - opts.tailFrac*(tt(end)-tt(1))));
else
    P = opts.Aext;                        % amplitude supplied from outside
end
if ~isfinite(P) || P <= 0; return; end

A = P; T = tt(end);
z  = yy / A;              % raw shape
zs = yls / A;             % smoothed shape, for crossings only
E.B = B; E.A = A; E.T = T;

% ---- 1. the current estimator: k and t0 both free -----------------------
obj = @(q) sum((z - (1 - exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
q   = fminsearch(obj, [45 0.004], optimset('Display','off','MaxFunEvals',4e3,'MaxIter',4e3));
E.k_fit = abs(q(1)); E.t0 = max(q(2),0);
pred = 1 - exp(-E.k_fit.*max(tt-E.t0,0));
E.r2 = 1 - sum((z-pred).^2)/max(sum((z-mean(z)).^2), eps);

% ---- 2. the same with the clock started at the valley (t0 == 0) ---------
obj0 = @(k) sum((z - (1 - exp(-abs(k).*tt))).^2);
E.k_t0z = abs(fminbnd(obj0, 1, 600, optimset('Display','off')));

% ---- 3. log-linear on the band: only the top of the rise ---------------
% Take the FIRST contiguous stretch inside the band, so a noisy plateau
% excursion back through 0.9 cannot re-enter the regression.
inb = zs >= opts.band(1) & zs <= opts.band(2);
d   = diff([0; inb; 0]); rs = find(d==1); re = find(d==-1)-1;
if ~isempty(rs)
    [~, j] = max(re-rs);                 % the longest contiguous stretch
    idx = rs(j):re(j);
    if numel(idx) >= 5
        w = 1 - z(idx);                  % on the RAW shape
        gd = idx(w > 1e-3);
        if numel(gd) >= 5
            pf = polyfit(tt(gd), log(1-z(gd)), 1);
            E.k_band = -pf(1);
            E.nBand  = numel(gd);
        end
    end
end

% ---- 4/5. model-free level crossings -----------------------------------
E.tOn   = xLevel(tt, zs, 0.05);
t20     = xLevel(tt, zs, 0.20);
t50     = xLevel(tt, zs, 0.50);
t80     = xLevel(tt, zs, 0.80);
if isfinite(t50) && t50 > 0;              E.k_t50  = log(2)/t50;        end
if isfinite(t20) && isfinite(t80) && t80 > t20
    E.k_2080 = log(4)/(t80-t20);
end

% ---- 6/7. the mean recovery time -- a LINEAR functional ----------------
MRT = trapz(tt, 1 - z);                   % = t0 + 1/k for a first-order rise
E.MRT = MRT;
if MRT > 0
    E.k_area = 1/MRT;                     % uncorrected
    % finite window: INTEGRAL_0^T (1-z) = (1-exp(-kT))/k
    f = @(k) (1-exp(-k*T))/k - MRT;
    try
        E.k_mrt = fzero(f, [1/max(MRT,eps)*0.2, 5000]);
    catch
        E.k_mrt = E.k_area;
    end
end
end

% ------------------------------------------------------------------------
function tx = xLevel(tt, zs, lev)
%XLEVEL first upward crossing of LEV, linearly interpolated. NaN if none.
    tx = NaN;
    i = find(zs >= lev, 1);
    if isempty(i) || i == 1; if ~isempty(i); tx = tt(1); end; return; end
    z1 = zs(i-1); z2 = zs(i);
    if z2 == z1; tx = tt(i); else; tx = tt(i-1) + (lev-z1)/(z2-z1)*(tt(i)-tt(i-1)); end
end
