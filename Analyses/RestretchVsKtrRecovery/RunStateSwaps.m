% RunStateSwaps.m
% Causal test of what the flux ledger implicated. Both recovery windows are
% holds, so a window IS its initial state — transplant one component at a time
% from the post-slack start into the post-restretch start and see which one
% carries the 2x rate difference.
%
% H0a/H0b  unmodified replays              -> validation gate for everything else
% H1       p1,p2 <- donor (re-registered)  -> the surviving BOUND pool
% H2       PD    <- donor                  -> the primed free RESERVOIR
% H3       SR,SRD<- donor                  -> the SRX POOL
% H3p      LSE = 0 so Force(t0) ~ 0        -> the force GATE itself (decisive)
% H4       A_reg <- donor                  -> availability (expect no effect)
%
% Run: cd(root); addpath(genpath('.')); Analyses/RestretchVsKtrRecovery/RunStateSwaps

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','RestretchVsKtrRecovery','results');
C = load(fullfile(resdir,'window_capture.mat'));

p  = C.p;  os = C.os;
ss = p.ss; Ns = p.NumberOfStates;
iP1 = 1:ss;  iP2 = ss+1:2*ss;
iSR = Ns*ss+1; iSL = Ns*ss+3; iLSE = Ns*ss+4; iPD = Ns*ss+5; iSRD = Ns*ss+6;
iXV = Ns*ss+7; iAR = Ns*ss+8;

sel = @(nm) C.W(strcmp({C.W.name}, nm));
wD = sel('postSlack2');       % donor
wR = sel('postRestretch2');   % recipient
DUR = 0.25;

grab = @(w) struct('PU', os.PU(w.i0,:), 'piv', os.LXBPivot(w.i0), ...
                   't0', os.t(w.i0), 'B', w.B*w.Fiso, 'P', w.P*w.Fiso, 'k63', w.k63);
D = grab(wD);  R = grab(wR);

% strain origins, for re-registering the transplanted distributions
shift = @(S) (-(S.PU(iSL) - S.PU(iLSE)) + S.piv)/2;
sD = p.s(:) - shift(D);   sR = p.s(:) - shift(R);

fprintf('donor    postSlack2     : F=%.2f PD=%.3f SRX=%.3f p2=%.4f k63=%.1f\n', ...
        D.PU(iLSE)*p.kSE, D.PU(iPD), D.PU(iSR)+D.PU(iSRD), sum(D.PU(iP2))*p.dS, D.k63);
fprintf('recipient postRestretch2: F=%.2f PD=%.3f SRX=%.3f p2=%.4f k63=%.1f\n\n', ...
        R.PU(iLSE)*p.kSE, R.PU(iPD), R.PU(iSR)+R.PU(iSRD), sum(R.PU(iP2))*p.dS, R.k63);

%% ---- build the hybrids ----
H = {};
H{end+1} = mk('H0a donor replay',      D.PU, D, D);
H{end+1} = mk('H0b recip replay',      R.PU, R, R);

% H1: bound pool, re-registered onto the recipient's strain axis
v = R.PU;
v(iP1) = max(0, interp1(sD, D.PU(iP1), sR, 'linear', 0));
v(iP2) = max(0, interp1(sD, D.PU(iP2), sR, 'linear', 0));
H{end+1} = mk('H1 p1,p2 <- donor',     v, R, R);

% H1raw: same but raw grid-index copy (sensitivity to the re-registration)
v = R.PU; v(iP1) = D.PU(iP1); v(iP2) = D.PU(iP2);
H{end+1} = mk('H1raw grid copy',       v, R, R);

% H2: primed reservoir
v = R.PU; v(iPD) = D.PU(iPD);
H{end+1} = mk('H2 PD <- donor',        v, R, R);

% H3: SRX pool
v = R.PU; v(iSR) = D.PU(iSR); v(iSRD) = D.PU(iSRD);
H{end+1} = mk('H3 SRX <- donor',       v, R, R);

% H3p: the force GATE. LSE=0 -> Force=0 -> the SRX gate sees F=0 like post-slack.
% NB this also moves the strain origin (s depends on SL-LSE) — read with H1.
v = R.PU; v(iLSE) = 0; v(iXV) = 0;
H{end+1} = mk('H3p Force(t0)->0',      v, R, R);

% H4: availability
v = R.PU; v(iAR) = D.PU(iAR);
H{end+1} = mk('H4 A_reg <- donor',     v, R, R);

%% ---- run ----
% Each replay is self-normalised (B, P taken from the replay itself), exactly as
% recoveryWindows does for the data. Passing the ORIGINAL window's B/P instead
% breaks any hybrid whose force at t0 differs from the original valley: the
% normalised trace then starts above the 5% onset threshold, i_on collapses to
% sample 1 and the 63% crossing is found immediately, giving absurd rates.
fprintf('%-22s %8s %8s | %7s %7s %7s %7s | %s\n', ...
        'hybrid','k63','vs H0b','PT','PD','SRX','p1+p2','mass');
kH0b = NaN;
for q = 1:numel(H)
    h = H{q};
    p1_0 = sum(h.PU(iP1))*p.dS; p2_0 = sum(h.PU(iP2))*p.dS;
    SRX  = h.PU(iSR) + h.PU(iSRD); PD = h.PU(iPD);
    PT   = 1 - (p1_0 + p2_0 + PD + SRX);
    if PT < 0
        fprintf('%-22s %8s   REJECTED: implied PT = %.4f < 0\n', h.name, '--', PT);
        continue;
    end
    [o, k] = seedWindow(p, h.PU, h.ref.piv, h.ref.t0, DUR, [], []);
    if strncmp(h.name,'H0b',3); kH0b = k; end
    % VALIDITY: the post-restretch span is only ~0.17 F_iso, so a hybrid whose
    % force at t0 lands near its own plateau leaves no rise to measure and the
    % 63% crossing degenerates. Reject rather than report a spurious rate.
    span = max(o.Force) - o.Force(1);
    bad  = ~isfinite(k) || k > 5*h.base.k63 || span < 0.05*max(o.Force);
    if bad
        fprintf('%-22s %8s %8s | %7.4f %7.4f %7.4f %7.4f | rise=%.3f Fiso -- NOT MEASURABLE\n', ...
                h.name, 'n/a', '--', PT, PD, SRX, p1_0+p2_0, span/h.ref.P);
        H{q}.k = NaN; continue;
    end
    fprintf('%-22s %8.1f %8s | %7.4f %7.4f %7.4f %7.4f | %+8.1e\n', ...
            h.name, k, ternstr(isnan(kH0b), '--', sprintf('%+.1f', k-kH0b)), ...
            PT, PD, SRX, p1_0+p2_0, PT+p1_0+p2_0+PD+SRX-1);
    H{q}.k = k;
end

save(fullfile(resdir,'state_swaps.mat'),'H','D','R');
fprintf('\nbaselines: donor k63 %.1f | recipient k63 %.1f\n', D.k63, R.k63);
fprintf('saved results/state_swaps.mat\n');

function h = mk(name, PU, ref, base)
    h.name = name; h.PU = PU; h.ref = ref; h.base = base;
end

function s = ternstr(c, a, b); if c; s = a; else; s = b; end; end
