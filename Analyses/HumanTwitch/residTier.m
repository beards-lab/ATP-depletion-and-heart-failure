function [r, sg, fm] = residTier(x, cfg, tb, Fd)
%RESIDTIER  Normalised residual (sigma*F_model - F_data)/F_peak of tier CFG at parameter vector X.
[b, th] = buildModel(cfg, x);
fm = simHumanTwitch(b, th, tb);
sg = max((fm'*Fd)/(fm'*fm), 0);
r = (sg*fm - Fd)/max(Fd);
end
