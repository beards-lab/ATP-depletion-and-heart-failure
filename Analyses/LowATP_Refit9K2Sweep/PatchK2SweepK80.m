% PatchK2SweepK80.m -- re-score the k2 sweep after extractKtrAttributes gained
% the exponential extrapolation of ktrP_k80 (the 2 mM ktr hold is 58.5 ms; at
% g <= 0.44 the 80 % crossing was never reached and k80 was NaN).
% Reads results/k2sweep_patched.mat (k80 re-extracted for the NaN runs by
% re-running only the ktr protocol), recomputes every score, overwrites
% results/k2sweep.mat and re-runs SummarizeK2Sweep.
rootP = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootP); addpath(genpath(rootP));
resP = fullfile(rootP,'Analyses','LowATP_Refit9K2Sweep','results');
load(fullfile(resP,'k2sweep_patched.mat')); %#ok<LOAD>
zFV = abs(FV2K2.r-1)/sqrt(FV2K2.s8r^2 + FV2K2.sr^2);
for kP = 1:numel(RK2)
    c   = evalFeatureCost(RK2(kP).fd, RK2(kP).fm, fn2K2);
    cFV = ((RK2(kP).fv2 - FV2K2.m2)/FV2K2.s2)^2;
    RK2(kP).cost = c;
    RK2(kP).score_plain = sum(c) + cFV;
    RK2(kP).score_rel   = sum(c(:)'.*aK2) + min(zFV,3)/3*cFV;
    fprintf('g = %.2f: plain %.4f  rel %.4f\n', RK2(kP).g, RK2(kP).score_plain, RK2(kP).score_rel);
end
save(fullfile(resP,'k2sweep.mat'), 'RK2','fm8K2','fd8K2','tr8K2','c8K2','fv8K2','fn8K2','fn2K2', ...
     'aK2','zK2','FV2K2','CORE','GLIST');
SummarizeK2Sweep;
