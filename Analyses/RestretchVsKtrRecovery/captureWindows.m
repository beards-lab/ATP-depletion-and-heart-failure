% captureWindows.m
% =========================================================================
% Baseline substrate for the post-restretch excess diagnostic.
%
% Runs ONE slack + ONE ktr simulation on the current ThursdayNightFever and
% resolves the ABSOLUTE sample indices of each recovery window. recoveryWindows
% returns relative time only (the post-restretch window is shifted so t=0 is the
% vall2 minimum), so the absolute bounds are recomputed here by `windowBounds`,
% which mirrors recoveryWindows.m:44-95 exactly and is ASSERTED against it.
%
% Output: results/window_capture.mat with
%   os, ok      full model output structs (all out.R* fluxes + out.PU)
%   svt, kvt    velocity tables
%   p           the params actually used
%   W           struct array, one per analysed window:
%                 .name .kind .cyc .i0 .i1 .i63 .Fiso .k63 .B .P
%
% Analysis set: post-slack cycles 1-5 and post-restretch cycles 1-5 (cycle 5 is
% captured but flagged .excluded, it returns to 1.00 ML not 1.10) plus ktr.
%
% Run: cd(root); addpath(genpath('.')); Analyses/RestretchVsKtrRecovery/captureWindows
% =========================================================================

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','RestretchVsKtrRecovery','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

SLACK_F = 'protocol_03_27_2026_8mM_slack.mat';
KTR_F   = 'protocol_03_27_2026_velocitytable_ktr.mat';
S = load(fullfile(root,'data',SLACK_F)); svt = S.velocitytable;
K = load(fullfile(root,'data',KTR_F));   kvt = K.velocitytable;

%% ---------------- run the baseline ----------------
base = getParams(loadParams('ThursdayNightFever'), [], true, false);
base.PlotEachSeparately = 0; base.ghostSave = ''; base.ghostLoad = '';
base.EvalFeatures = 0;
% 'All', NOT 'AllPar' -- sequential, so cycle-to-cycle carryover is real and the
% captured state at each window start is the one the full protocol produces.
base.RunSlackSegments = 'All';
p = getParams(base, base.g, true);

fprintf('params: %s\n', which('ThursdayNightFever'));
fprintf('UseMaxwellTensionOnly=%d tau_reg=%g  numel(PU0)=%d\n', ...
        p.UseMaxwellTensionOnly, p.tau_reg, numel(p.PU0));
assert(numel(p.PU0) == p.NumberOfStates*p.ss + 7 + p.UseRegistrationAvailability + p.UseD0State, ...
       'captureWindows:PU0len','unexpected PU0 length');

tic; [~, os] = runSlackExperiment(p); fprintf('slack %.0f s\n', toc);

pk = p; pk.SL0 = 2*kvt(1,4); pk = getParams(pk, pk.g, true); pk.Velocity = kvt(:,2);
tic; [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk); fprintf('ktr   %.0f s\n', toc);

%% ---------------- reference estimator (the shared one) ----------------
iS = find(svt(:,2) < -1);
FisoS = interp1(os.t, os.Force, svt(iS(1),1)-0.02, 'linear', 'extrap');
FisoK = interp1(ok.t, ok.Force, kvt(2,1), 'linear', 'extrap');
wsRef = recoveryWindows(os.t, os.Force, svt, 'slack', FisoS);
wkRef = recoveryWindows(ok.t, ok.Force, kvt, 'ktr',   FisoK);

%% ---------------- absolute bounds + assertions ----------------
W = struct('name',{},'kind',{},'cyc',{},'i0',{},'i1',{},'i63',{}, ...
           'Fiso',{},'k63',{},'B',{},'P',{},'excluded',{},'src',{});

Ws = windowBounds(os.t, os.Force, svt, 'slack', FisoS);
Wk = windowBounds(ok.t, ok.Force, kvt, 'ktr',   FisoK);

fprintf('\n%-18s %6s %6s %6s %8s %8s %8s   assert\n','window','i0','i63','i1','B','P','k63');
for q = 1:numel(Ws)
    ref = wsRef(strcmp({wsRef.type}, Ws(q).kind) & [wsRef.cyc] == Ws(q).cyc);
    Ws(q) = checkAgainst(Ws(q), ref, os.t, os.Force, FisoS);
    Ws(q).src = 'os';
    W(end+1) = Ws(q); %#ok<SAGROW>
end
Wk.src = 'ok';
Wk = checkAgainst(Wk, wkRef, ok.t, ok.Force, FisoK);
W(end+1) = Wk;

save(fullfile(resdir,'window_capture.mat'), 'os','ok','svt','kvt','p','pk','W','FisoS','FisoK','-v7.3');
fprintf('\nsaved results/window_capture.mat  (%d windows)\n', numel(W));

%% ================= helpers =================
function W = windowBounds(t, F, vt, kind, Fiso)
% Mirrors recoveryWindows.m:44-95 but returns ABSOLUTE sample indices.
    t = t(:); Fn = F(:)/Fiso;
    keep = [true; diff(t) > 0];          % same duplicate-timestamp guard
    idx  = find(keep);                    % map filtered -> original index
    t = t(keep); Fn = Fn(keep);
    VALL = 0.080; SM = 0.002;
    w = {};
    switch lower(kind)
        case 'slack'
            iS = find(vt(:,2) < -1); iR = find(vt(:,2) > 1);
            nC = min(numel(iS), numel(iR));
            for c = 1:nC
                % (S) window = release-ramp END -> restretch ramp start
                ta = vt(iS(c)+1,1); tb = vt(iR(c),1);
                m = t >= ta & t <= tb;
                if nnz(m) > 20
                    j = find(m);
                    w{end+1} = mk(sprintf('postSlack%d',c),'postSlack',c, ...
                                  j(1), j(end), t, Fn, 0, SM, idx); %#ok<AGROW>
                end
                % (R) window starts at vall2 = min of the 80 ms after the ramp
                t_r1 = vt(iR(c)+1,1);
                if c < nC; td = vt(iS(c+1),1); else; td = max(t); end
                m = t >= t_r1 & t <= td;
                if nnz(m) > 20
                    j = find(m); tt = t(j); FF = Fn(j);
                    mv = tt <= tt(1)+VALL; [~,iv] = min(FF(mv));
                    w{end+1} = mk(sprintf('postRestretch%d',c),'postRestretch',c, ...
                                  j(iv), j(end), t, Fn, Fn(j(iv)), SM, idx); %#ok<AGROW>
                end
            end
        case 'ktr'
            ta = vt(end-1,1);
            m = t >= ta & t <= vt(end,1); j = find(m);
            w{end+1} = mk('ktr','ktr',0, j(1), j(end), t, Fn, 0, SM, idx);
    end
    W = [w{:}];
end

function w = mk(name, kind, cyc, j0, j1, t, Fn, B, SM, idx)
    w.name = name; w.kind = kind; w.cyc = cyc;
    w.B = B;
    if strcmp(kind,'ktr')
        w.P = 1.0;                                   % pre-release F_iso (known endpoint)
    else
        seg = j0:j1; ts = t(seg);
        w.P = median(Fn(seg(ts >= ts(end)-0.15*(ts(end)-ts(1)))));
    end
    % t63 with the SAME smoothing/crossing rule as recoveryWindows.m:118-125
    seg = j0:j1; y = (movmedian(Fn(seg), max(3, round(SM/median(diff(t(seg)))))) - B)/max(w.P-B,eps);
    below = find(y < 0.05); if isempty(below); i_on = 1; else; i_on = below(end); end
    i63 = find(y >= 1-exp(-1) & (1:numel(y))' > i_on, 1);
    if isempty(i63); w.k63 = NaN; w.i63 = j1; else
        w.k63 = 1/(t(seg(i63)) - t(seg(i_on)));  w.i63 = idx(seg(i63));
    end
    w.i0 = idx(j0); w.i1 = idx(j1);                  % ABSOLUTE indices into the raw out.*
    w.excluded = (cyc == 5);
    w.Fiso = NaN; w.src = '';
end

function w = checkAgainst(w, ref, t, F, Fiso)
% Assert the absolute bounds reproduce the shared estimator byte-for-byte.
    w.Fiso = Fiso;
    seg = w.i0:w.i1;
    % NB recoveryWindows measures postSlack/ktr time from the RELEASE-RAMP START
    % (t_org), not from the window start, so ref.t(1) > 0 there. The invariant is
    % the window DURATION, which is origin-independent.
    dur = t(w.i1) - t(w.i0);
    okDur = abs(dur - (ref.t(end) - ref.t(1))) < 1e-9;
    Fseg = F(seg)/Fiso;
    okF  = numel(Fseg) == numel(ref.F) && max(abs(Fseg(:) - ref.F(:))) < 1e-9;
    okK  = isnan(w.k63) == isnan(ref.k63) || abs(w.k63 - ref.k63) < 1e-6;
    tag = 'OK';
    if ~(okDur && okF && okK); tag = sprintf('*** MISMATCH dur=%d F=%d k=%d', okDur, okF, okK); end
    fprintf('%-18s %6d %6d %6d %8.3f %8.3f %8.1f   %s\n', ...
            w.name, w.i0, w.i63, w.i1, w.B, w.P, w.k63, tag);
    assert(okDur && okF && okK, 'captureWindows:mismatch', ...
           'window %s does not reproduce recoveryWindows', w.name);
end
