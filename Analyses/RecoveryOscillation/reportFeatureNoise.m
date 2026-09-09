function reportFeatureNoise(N, FEATS, resDir)
%REPORTFEATURENOISE  s_noise against the between-preparation spread that the
%   objective's weights w = lambda/s^2 are built on, and what dropping the
%   fifth slack cycle does to it.

evalc('T = buildPooledTarget();');   % s8 = the moderated spread the weights use

fid = fopen(fullfile(resDir,'feature_noise.txt'),'w');
p = @(varargin) tee(fid, varargin{:});

p('HOW MUCH OF EACH FEATURE SPREAD IS THE TRANSDUCER?\n');
p('%s\n', repmat('=',1,104));
p('s_noise : SD across %d re-extractions of the SAME trace with its own measured\n', N(1).nrep);
p('          noise added -- biology held exactly fixed. Relative, per cent.\n');
p('s_prep  : the between-preparation spread the weight w = lambda/s^2 is built on.\n');
p('noise   : s_noise^2 / s_prep^2, the share of the weight that is rig, not biology.\n');
p('s_1-4   : the same, with the fifth slack cycle dropped. Cycle 5 re-stretches to a\n');
p('          DIFFERENT final length by the SMALLEST amplitude, so it is neither the\n');
p('          same manoeuvre nor as well conditioned (Analyses/RestretchMetricAudit).\n\n');
p('sigma per run (kPa): %s\n\n', strjoin(arrayfun(@(n) ...
    sprintf('%s %s %.2f', n.day, n.cond, n.sigma), N, 'UniformOutput', false), ' | '));

act  = strcmp({N.cond},'8mM');
name = {}; sn = []; s14 = []; sp = []; fr = [];
for k = 1:numel(FEATS)
    fk = FEATS{k};
    a = 100*mean([N(act).(fk)], 'omitnan');
    c = 100*mean([N(act).([fk '_s14'])], 'omitnan');
    j = find(strcmp(T.name, fk), 1);
    if isempty(j); b = NaN; else; b = 100*T.s8(j); end
    name{end+1} = fk; sn(end+1) = a; s14(end+1) = c; sp(end+1) = b;   %#ok<AGROW>
    fr(end+1) = 100*a^2/max(b,eps)^2;                                 %#ok<AGROW>
end
[~, o] = sort(fr, 'descend', 'MissingPlacement','last');

p('%-22s %9s %9s %9s %9s   %s\n','feature','s_noise %','s_prep %','noise %','s_1-4 %','verdict');
p('%s\n', repmat('-',1,104));
for i = o
    if ~isfinite(sp(i))
        p('%-22s %9.1f %9s %9s %9.1f   not in the pooled target\n', name{i}, sn(i), '-', '-', s14(i));
        continue
    end
    if     fr(i) >= 100; v = 'RIG -- s_prep is BELOW its own noise';
    elseif fr(i) >= 50;  v = 'mostly rig';
    elseif fr(i) >= 20;  v = 'part rig';
    else;                v = 'biology';
    end
    p('%-22s %9.1f %9.1f %9.0f %9.1f   %s\n', name{i}, sn(i), sp(i), fr(i), s14(i), v);
end

p('\nWHAT DROPPING CYCLE 5 BUYS  (features whose noise falls by more than a third)\n');
any5 = false;
for i = o
    if isfinite(s14(i)) && isfinite(sn(i)) && s14(i) < 0.67*sn(i)
        p('   %-20s  s_noise %6.1f %% -> %5.1f %%   (x %.2f)\n', name{i}, sn(i), s14(i), s14(i)/sn(i));
        any5 = true;
    end
end
if ~any5; p('   nothing -- cycle 5 is not the problem.\n'); end

p('\nPER-RUN, PER-CYCLE s_noise for the re-stretch rate (the term at issue)\n');
p('%-8s %-5s %7s | %s\n','day','cond','sigma','cycle 1..5, %');
for i = 1:numel(N)
    if ~isfield(N(i),'rsK80_seg') || all(isnan(N(i).rsK80_seg)); continue; end
    p('%-8s %-5s %7.2f | %s\n', N(i).day, N(i).cond, N(i).sigma, ...
      strjoin(arrayfun(@(x) sprintf('%6.1f', 100*x), N(i).rsK80_seg, 'UniformOutput', false), ' '));
end

p('\nFLOOR.  s must never be below s_noise. On the ALL-CYCLE feature:\n');
anyF = false;
for i = o
    if isfinite(sp(i)) && sn(i) > sp(i)
        p('   %-20s  s %.1f %% -> %.1f %%   (weight x %.2f)   [with cycle 5 dropped: %.1f %%]\n', ...
          name{i}, sp(i), sn(i), (sp(i)/sn(i))^2, s14(i)); anyF = true;
    end
end
if ~anyF; p('   none.\n'); end

p('\nATP CONTROL -- noise is a rig property and must not move with nucleotide:\n');
p('%-22s %10s %10s %8s\n','feature','8 mM','2 mM','ratio');
for k = 1:numel(FEATS)
    fk = FEATS{k};
    a = 100*mean([N(strcmp({N.cond},'8mM')).(fk)],'omitnan');
    b = 100*mean([N(strcmp({N.cond},'2mM')).(fk)],'omitnan');
    if isfinite(a) && a > 2
        p('%-22s %10.1f %10.1f %8.2f\n', fk, a, b, b/a);
    end
end
fclose(fid);
end

function tee(fid, varargin); fprintf(1, varargin{:}); fprintf(fid, varargin{:}); end
