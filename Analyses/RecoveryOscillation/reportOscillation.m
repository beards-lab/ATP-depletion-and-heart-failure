function reportOscillation(W, resDir)
%REPORTOSCILLATION  The five signatures, as a table.

fid = fopen(fullfile(resDir,'oscillation_id.txt'),'w');
p = @(varargin) tee(fid, varargin{:});

day  = {W.day};  cond = {W.cond};  type = {W.type};
f0   = [W.f0];   rel  = [W.relRMS]; bf = [W.bandFrac];
r2   = [W.r2];   A    = [W.A];      pk = [W.pkAmp];
f0L  = [W.f0L];  pkL  = [W.pkAmpL];
act  = ~strcmp(cond,'PNB');

p('RECOVERY RIPPLE -- rig resonance or biology?   (%d windows)\n', numel(W));
p('%s\n', repmat('=',1,86));

p('\nS1  ONE FREQUENCY?   f0 of the residual, Hz\n');
p('%-14s %8s %8s %8s %8s %8s\n','group','n','mean','sd','cv%','median');
grp = strcat(cond, {' '}, type);
u = unique(grp,'stable');
for i = 1:numel(u)
    s = strcmp(grp,u{i}) & isfinite(f0);
    p('%-14s %8d %8.1f %8.1f %8.1f %8.1f\n', u{i}, nnz(s), mean(f0(s)), ...
      std(f0(s)), 100*std(f0(s))/mean(f0(s)), median(f0(s)));
end
p('%-14s %8d %8.1f %8.1f %8.1f %8.1f\n','ALL', nnz(isfinite(f0)), ...
  mean(f0,'omitnan'), std(f0,'omitnan'), ...
  100*std(f0,'omitnan')/mean(f0,'omitnan'), median(f0,'omitnan'));

p('\nS2  ATP-INVARIANT?   f0 and ripple size by nucleotide\n');
p('%-8s %8s %10s %10s %10s\n','cond','n','f0 (Hz)','relRMS %','bandFrac');
for c = {'8mM','2mM','PNB'}
    s = strcmp(cond,c{1}) & isfinite(f0);
    if ~any(s); continue; end
    p('%-8s %8d %10.1f %10.1f %10.2f\n', c{1}, nnz(s), mean(f0(s)), ...
      100*mean(rel(s),'omitnan'), mean(bf(s),'omitnan'));
end
r8 = mean(f0(strcmp(cond,'8mM')),'omitnan'); r2m = mean(f0(strcmp(cond,'2mM')),'omitnan');
p('  f0 ratio 2mM/8mM = %.3f   (every measured RATE moves by 0.35-0.54)\n', r2m/r8);

p('\nS3  PRESENT WITHOUT CROSS-BRIDGES?   PNB+Mava, same rig, same manoeuvre\n');
s = strcmp(cond,'PNB') & isfinite(f0);
p('  PNB windows           n = %d\n', nnz(s));
p('  f0                    %.1f +/- %.1f Hz   (active: %.1f +/- %.1f)\n', ...
  mean(f0(s)), std(f0(s)), mean(f0(act&isfinite(f0))), std(f0(act&isfinite(f0))));
p('  ripple amplitude      %.2f kPa            (active: %.2f kPa)\n', ...
  mean(pk(s),'omitnan'), mean(pk(act),'omitnan'));

p('\nS4  AMPLITUDE TRACKS THE EXCITATION?   ripple vs rise amplitude / ramp speed\n');
sa = act & isfinite(A) & isfinite(pk);
if nnz(sa) > 3
    cA = corr(A(sa)', pk(sa)');
    cV = corr(abs([W(sa).vRamp])', pk(sa)');
    p('  corr(ripple amp, rise amplitude A) = %+.2f\n', cA);
    p('  corr(ripple amp, ramp velocity)    = %+.2f\n', cV);
    p('  ripple / rise amplitude            = %.1f %% (sd %.1f)\n', ...
      100*mean(pk(sa)./A(sa)), 100*std(pk(sa)./A(sa)));
end

p('\nS5  IN THE LENGTH CHANNEL?\n');
sL = isfinite(f0L) & isfinite(f0);
p('  length f0             %.1f +/- %.1f Hz\n', mean(f0L(sL)), std(f0L(sL)));
p('  matches force f0      %.0f %% of windows within 25 %%\n', ...
  100*mean(abs(f0L(sL)-f0(sL))./f0(sL) < 0.25));
p('  length ripple amp     %.2e ML\n', mean(pkL(sL)));

p('\nWHAT IT COSTS THE FIT\n');
for tp = {'postSlack','postRestretch'}
    s = act & strcmp(type,tp{1}) & isfinite(r2);
    p('  %-14s R2 of the pure exponential = %.3f +/- %.3f  (n=%d)\n', ...
      tp{1}, mean(r2(s)), std(r2(s)), nnz(s));
    p('  %-14s residual RMS = %.1f %% of the rise\n', '', 100*mean(rel(s),'omitnan'));
end
p('\n  A model with NO ripple can therefore never exceed R2 ~ %.2f on this data;\n', ...
  mean(r2(act & isfinite(r2))));
p('  it will be scored as if it were wrong by %.1f %% of the rise when it is not.\n', ...
  100*mean(rel(act),'omitnan'));

fclose(fid);
end

function tee(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end
