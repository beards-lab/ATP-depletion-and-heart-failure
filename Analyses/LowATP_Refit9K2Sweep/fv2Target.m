function T = fv2Target()
%FV2TARGET Pooled in-preparation force-velocity point (FV2) at 8 and 2 mM.
%
%   The 04-series slack recordings open with a 34 ms isovelocity ramp at
%   -2 ML/s (Analyses/DataPooling/RunFVConsistency.m, conclusions Sec. 6).
%   Two preparations have it (04/03, 04/10); 03/27 predates the manoeuvre.
%   P/P0 = mean force over [ramp start + 5 ms, ramp end] / pre-ramp isometric
%   (median over 2-40 ms), RAW (passive included -- the model carries passive
%   force too, so raw is the like-for-like number).
%
%   T.m8, T.s8   pooled P/P0 at 8 mM and its absolute SD (n = 2)
%   T.m2, T.s2   same at 2 mM
%   T.r,  T.sr   pooled 2/8 ratio of P/P0 and the SD of the per-prep ratios
%   T.s8r        relative 8 mM spread (for the ATP-effect z-score)
DAYS = {'04_03','04_10'};
V = nan(numel(DAYS), 2);
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
for d = 1:numel(DAYS)
    for c = 1:2
        suf = ''; if c == 2; suf = '_2mM'; end
        S = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_FV2%s.mat', DAYS{d}, suf)));
        V(d,c) = fv2Read(S.datatable(:,1), S.datatable(:,3), S.velocitytable);
    end
end
T.perPrep = V; T.days = DAYS; T.n = numel(DAYS);
T.m8 = mean(V(:,1)); T.s8 = std(V(:,1));
T.m2 = mean(V(:,2)); T.s2 = std(V(:,2));
rr = V(:,2)./V(:,1);
T.r = mean(rr); T.sr = std(rr); T.s8r = T.s8/T.m8;
end
