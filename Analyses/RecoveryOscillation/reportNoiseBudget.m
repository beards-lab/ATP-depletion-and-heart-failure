function reportNoiseBudget(B, resDir)
fid = fopen(fullfile(resDir,'noise_budget.txt'),'w');
p = @(varargin) tee(fid, varargin{:});

p('WHAT IS THE RESIDUAL OF A SINGLE EXPONENTIAL MADE OF?\n');
p('%s\n', repmat('=',1,104));
p('sigW  white transducer noise, measured on the plateau of the same window, above 60 Hz\n');
p('sigC  reproducible across the five cycles of a run (phase-locked to the ramp)\n');
p('sigS  everything else -- the only part a model could legitimately be asked to fit\n\n');

p('%-7s %-5s %-14s %3s | %6s %6s %6s %6s | %5s %5s %5s | %6s %6s %6s %6s\n', ...
  'day','cond','type','n','sigTot','sigW','sigC','sigS','%W','%C','%S','A(kPa)','R2','R2ceil','perC');
for i = 1:numel(B)
    b = B(i);
    p('%-7s %-5s %-14s %3d | %6.2f %6.2f %6.2f %6.2f | %5.0f %5.0f %5.0f | %6.1f %6.3f %6.3f %6.0f\n', ...
      b.day,b.cond,b.type,b.nc, b.sigTot,b.sigW,b.sigC,b.sigS, ...
      100*b.fracW,100*b.fracC,100*b.fracS, b.A, b.r2, b.r2ceil, b.perC);
end

act = ~strcmp({B.cond},'PNB');
for ty = {'postRestretch','postSlack'}
    s = act & strcmp({B.type},ty{1});
    p('\n%s, active runs (n = %d runs x 5 cycles)\n', ty{1}, nnz(s));
    p('  residual RMS            %.2f kPa   = %.1f %% of the rise (%.1f kPa)\n', ...
      mean([B(s).sigTot]), 100*mean([B(s).sigTot]./[B(s).A]), mean([B(s).A]));
    p('  of which WHITE NOISE    %.0f %% of the variance   (sigma = %.2f kPa)\n', ...
      100*mean([B(s).fracW]), mean([B(s).sigW]));
    p('  of which REPRODUCIBLE   %.0f %%                   (sigma = %.2f kPa, period ~ %.0f ms)\n', ...
      100*mean([B(s).fracC]), mean([B(s).sigC]), median([B(s).perC]));
    p('  UNEXPLAINED             %.0f %%                   (sigma = %.2f kPa)\n', ...
      100*mean([B(s).fracS]), mean([B(s).sigS]));
    p('  R2 of the exponential   %.3f;  a NOISE-FREE model can reach at most %.3f\n', ...
      mean([B(s).r2]), mean([B(s).r2ceil]));
end

sP = strcmp({B.cond},'PNB');
p('\nCONTROL -- PNB+Mava (cross-bridges off, same rig, same ramps)\n');
p('  sigW  %.2f kPa   (active %.2f)\n', mean([B(sP).sigW]), mean([B(act).sigW]));
p('  sigC  %.2f kPa   (active %.2f)  -- a reproducible ripple with NO cross-bridges\n', ...
  mean([B(sP).sigC]), mean([B(act).sigC]));

p('\nATP -- does the reproducible part know about nucleotide?\n');
for c = {'8mM','2mM'}
    s = strcmp({B.cond},c{1});
    p('  %-4s  sigC = %.2f kPa, period = %.0f ms\n', c{1}, mean([B(s).sigC]), median([B(s).perC]));
end
r8 = mean([B(strcmp({B.cond},'8mM')).perC]); r2 = mean([B(strcmp({B.cond},'2mM')).perC]);
p('  period ratio 2mM/8mM = %.2f   (ktr x1.9, rsK80 x2.6 slower at 2 mM)\n', r2/r8);

p('\nCONSEQUENCE FOR THE OBJECTIVE\n');
s = act & strcmp({B.type},'postRestretch');
p('  Any trace-wise least squares on the post-re-stretch window has a FLOOR of\n');
p('  %.2f kPa RMS (%.1f %% of the rise). A model sitting exactly on the truth is\n', ...
  sqrt(mean([B(s).sigW].^2 + [B(s).sigC].^2)), ...
  100*sqrt(mean([B(s).sigW].^2+[B(s).sigC].^2))/mean([B(s).A]));
p('  scored as that far wrong. Only %.0f %% of the residual variance is even in\n', ...
  100*mean([B(s).fracS]));
p('  principle reducible.\n');
fclose(fid);
end

function tee(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end
