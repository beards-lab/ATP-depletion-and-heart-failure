function win = recoveryWindows(t, F, vt, kind, Fiso, opts)
%RECOVERYWINDOWS Extract and rate-fit the force-recovery windows of a protocol.
%
%   WIN = recoveryWindows(T, F, VT, KIND, FISO)
%
%   Shared by CompareRecoveryRates.m (data) and ModelVsDataRecovery.m (model)
%   so that model and experiment go through byte-identical window definitions
%   and estimators. See conclusions.md for why each choice is made.
%
%   Inputs:
%     T, F  - time (s) and force. F is divided by FISO internally, so pass
%             FISO = 1 for traces that are already normalised.
%     VT    - velocitytable of the protocol.
%     KIND  - 'slack' (5 slack/restretch cycles) or 'ktr' (single event).
%     FISO  - isometric reference force, in the units of F.
%     OPTS  - optional struct: .span (equal-span control, default 0.0585 s)
%                              .smooth (median span for crossings, 0.002 s)
%                              .vallWin (valley search span, 0.080 s)
%
%   Output WIN is a struct array with one element per window:
%     .type   'postSlack' | 'postRestretch' | 'ktr'
%     .cyc    slack cycle number (0 for ktr)
%     .t,.F   the window samples (F in units of F_iso)
%     .B      baseline: 0 for postSlack/ktr, the vall2 minimum for postRestretch
%     .P      MEASURED plateau: end-of-window median, or 1.0 for ktr (the
%             58.5 ms window ends before the plateau, and the muscle is back
%             at the pre-release length, so F_iso is the known endpoint)
%     .A      amplitude P-B, in units of F_iso
%     .k63    model-free 1/(t63-t_on)                      <- primary estimator
%     .kfixA  fit A*(1-exp(-k*(t-t0))) with A fixed at P-B
%     .kfit   same with A free
%     .kfit58 A free, restricted to 58.5 ms after onset (the ktr-equivalent)
%     .t0,.r2A,.r2

    if nargin < 6 || isempty(opts); opts = struct(); end
    if ~isfield(opts,'span');    opts.span    = 0.0585; end
    if ~isfield(opts,'smooth');  opts.smooth  = 0.002;  end
    if ~isfield(opts,'vallWin'); opts.vallWin = 0.080;  end

    t = t(:); F = F(:)/Fiso;
    keep = [true; diff(t) > 0];            % model output can repeat timestamps
    t = t(keep); F = F(keep);

    w = {};
    switch lower(kind)
        case 'slack'
            iS = find(vt(:,2) < -1);  iR = find(vt(:,2) > 1);
            nC = min(numel(iS), numel(iR));
            for c = 1:nC
                % (S) time origin = release ramp START, window = ramp END -> restretch
                t_org = vt(iS(c),1);  ta = vt(iS(c)+1,1);  tb = vt(iR(c),1);
                m = t >= ta & t <= tb;
                if nnz(m) > 20
                    w{end+1} = mkwin(sprintf('S%d',c),'postSlack',c, ...
                        t(m)-t_org, F(m), 0, vt(iS(c)+1,4), opts); %#ok
                end
                % (R) window starts at vall2 = min of the 80 ms after the ramp
                t_r1 = vt(iR(c)+1,1);
                if c < nC; td = vt(iS(c+1),1); else; td = max(t); end
                m = t >= t_r1 & t <= td;  tt = t(m); FF = F(m);
                if nnz(m) > 20
                    mv = tt <= tt(1)+opts.vallWin;  [Bv,iv] = min(FF(mv));
                    w{end+1} = mkwin(sprintf('R%d',c),'postRestretch',c, ...
                        tt(iv:end)-tt(iv), FF(iv:end), Bv, vt(iR(c)+1,4), opts); %#ok
                end
            end
        case 'ktr'
            % time origin = start of the final 1.05 -> 1.00 release ramp
            t_org = vt(end-2,1);  ta = vt(end-1,1);
            m = t >= ta & t <= vt(end,1);
            w{end+1} = mkwin('K','ktr',0, t(m)-t_org, F(m), 0, vt(end-1,4), opts);
        otherwise
            error('recoveryWindows:kind','KIND must be ''slack'' or ''ktr''.');
    end
    win = [w{:}];
end

% ------------------------------------------------------------------------
function w = mkwin(id,type,cyc,t,F,B,ML,opts)
    t=t(:); F=F(:); ok=isfinite(t)&isfinite(F); t=t(ok); F=F(ok);
    w.id=id; w.type=type; w.cyc=cyc; w.t=t; w.F=F; w.B=B; w.ML=ML;
    if strcmp(type,'ktr')
        w.P = 1.0;                                     % pre-release isometric force
    else
        w.P = median(F(t >= t(end)-0.15*(t(end)-t(1))));
    end
    Afix = w.P - B;
    [w.Afree, w.kfit,  w.t0,  w.r2 ] = fitKtr(t,F,B,[],  []);
    [~,       w.kfixA, w.t0A, w.r2A] = fitKtr(t,F,B,Afix,[]);
    m58 = t <= (t(1)+w.t0) + opts.span;
    if nnz(m58) < 20; m58 = t <= t(1)+opts.span; end
    [~, w.kfit58, ~, w.r258] = fitKtr(t(m58),F(m58),B,[],[w.Afree w.kfit w.t0]);
    w.A   = Afix;
    w.k63 = crossRate(t,F,B,w.P,opts.smooth);
end

function [A,k,t0,r2] = fitKtr(t,F,B,Afix,q0)
% F = B + A*(1-exp(-k*max(t-t0,0))). B always FIXED at the baseline.
% Afix = [] -> A free (3 params); Afix = value -> A fixed ([k t0] free).
    tt = t - t(1);  y = F - B;
    A0 = max(median(y(max(1,end-round(0.15*numel(y))):end)), 1e-6);
    opt = optimset('Display','off','MaxFunEvals',8e3,'MaxIter',8e3);
    if isempty(Afix)
        if isempty(q0); q0 = [A0, 45, 0.004]; end
        obj = @(q) sum((y - abs(q(1)).*(1-exp(-abs(q(2)).*max(tt-max(q(3),0),0)))).^2);
        q = fminsearch(obj, q0, opt);
        A = abs(q(1)); k = abs(q(2)); t0 = max(q(3),0);
    else
        A = Afix;
        obj = @(q) sum((y - A.*(1-exp(-abs(q(1)).*max(tt-max(q(2),0),0)))).^2);
        q = fminsearch(obj, [45, 0.004], opt);
        k = abs(q(1)); t0 = max(q(2),0);
    end
    pred = A.*(1-exp(-k.*max(tt-t0,0)));
    r2 = 1 - sum((y-pred).^2)/max(sum((y-mean(y)).^2),eps);
end

function k = crossRate(t,F,B,P,SM)
    n  = max(3, round(SM/median(diff(t))));
    y  = (movmedian(F,n) - B)/max(P-B,eps);
    below = find(y < 0.05);
    if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(y >= 1-exp(-1) & (1:numel(y))' > i_on, 1);
    if isempty(i63) || t(i63) <= t(i_on); k = NaN; else; k = 1/(t(i63)-t(i_on)); end
end
