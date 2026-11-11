function y = sigLookup(sig, t)
% SIGLOOKUP  O(1) cubic Hermite interpolation on a uniform grid.
%
%   Y = SIGLOOKUP(SIG, T) evaluates the precomputed signal SIG at time T.
%   Clamps T to [t0, t0+(N-1)/Fs] — safe to call outside the signal range.
%
%   SIG struct fields (produced by preprocessDrivingSignal.m):
%     t0  — start time (s)
%     Fs  — sample rate (Hz)
%     N   — number of samples
%     y   — smoothed signal values (N×1)
%     m   — Catmull-Rom slopes (N×1), units/s
%
%   Uses cubic Hermite basis with Catmull-Rom slopes for C1 continuity.
%   Uniform grid means index k = floor((t-t0)*Fs)+1 — no binary search.
k   = max(1, min(sig.N - 1, 1 + round((t - sig.t0) * sig.Fs)));
y = sig.y(k);
return;


    k   = max(1, min(sig.N - 1, 1 + floor((t - sig.t0) * sig.Fs)));
    tau = (t - (sig.t0 + (k-1)/sig.Fs)) * sig.Fs;   % in [0, 1)

    t2  = tau*tau;  t3 = t2*tau;
    h00 =  2*t3 - 3*t2 + 1;
    h10 =    t3 - 2*t2 + tau;
    h01 = -2*t3 + 3*t2;
    h11 =    t3 -   t2;

    inv_Fs = 1 / sig.Fs;
    y = h00*sig.y(k) + h10*sig.m(k)*inv_Fs ...
      + h01*sig.y(k+1) + h11*sig.m(k+1)*inv_Fs;
end
