function T = buildPooledTarget(opts)
% BUILDPOOLEDTARGET  Turn the 2026 protocol days into ONE fit target with an
% honest error bar on every number.
%
%   T = buildPooledTarget()                 % defaults, see below
%   T = buildPooledTarget(struct('shrink',false))
%
% This is the operational form of the DataPooling conclusions. Four steps, in
% this order:
%
%   1. SCALE.  Divide every force-dimension feature of preparation p by that
%      preparation's own pre-slack isometric force at SL 2.2 (`Fss22`), taken
%      from its 8 mM run. ONE divisor per preparation, shared by both ATP
%      conditions -- so the 2 mM / 8 mM ratio survives untouched and stays a
%      measurement rather than a definition.
%   2. RUNDOWN. Multiply the run recorded SECOND in each session by
%      1 + phi*|slope_first|*T_first / F_second, phi = 0.452. Without this the
%      ATP ratio depends on the running order, which differs across the days.
%   3. POOL.  Mean of the per-preparation FEATURES (never features of an
%      averaged trace -- see RunPoolingComparison: averaging traces biases the
%      re-stretch recovery rate by 24 %, more than a preparation width).
%   4. SPREAD. s = between-preparation relative SD, optionally shrunk toward
%      the median of its own class. With n = 3 a raw SD has 2 dof and a 95 %
%      CI spanning 12x; shrinkage is the standard low-replicate remedy.
%
% Returns
%   T.name{k}         feature name
%   T.cls{k}          'F' force-dimension (normalised) | '0' scale-free
%   T.m8(k,:)         pooled 8 mM value, one column per slack segment
%   T.s8(k)           relative between-preparation SD at 8 mM (the chi-2 scale)
%   T.s8raw(k)        the same before shrinkage
%   T.r(k)            pooled 2 mM / 8 mM ratio  (alpha-free)
%   T.sr(k)           its between-preparation relative SD
%   T.z(k)            |r-1| / sqrt(s8^2 + sr^2) -- ATP discriminability
%   T.scale           absolute anchor: mean Fss22 in kPa, its CV, per-prep values
%   T.opts            what was actually applied
%
% The one absolute number the normalisation removes is kept in T.scale: fit the
% shape in normalised units and the LEVEL through that single term, rather than
% dropping it and leaving the objective scale-blind.

if nargin < 1; opts = struct(); end
if ~isfield(opts,'anchor');  opts.anchor  = 'Fss22'; end
if ~isfield(opts,'rundown'); opts.rundown = true;    end
if ~isfield(opts,'shrink');  opts.shrink  = true;    end
if ~isfield(opts,'phi');     opts.phi     = 0.452;   end
% Which slack cycles each feature's SPREAD is measured over. The pooled VALUES
% always carry all five; this sets only the segments behind s, and s must be
% measured on exactly the segments the objective scores. Cycle 5 re-stretches
% to a DIFFERENT final length by the SMALLEST amplitude, so for the re-stretch
% rate it is both a different manoeuvre (Analyses/RestretchMetricAudit sec. 3)
% and 3-5x noisier (Analyses/RecoveryOscillation) -- there, s(1:5) is entirely
% accounted for by the transducer and the feature carries no resolvable
% between-preparation information at all.
%   opts.segments = 1:4                      applies to every feature
%   opts.segments = struct('rsK80', 1:4)     per-feature override, rest 1:5
if ~isfield(opts,'segments'); opts.segments = 1:5;   end
if ~isfield(opts,'segDefault'); opts.segDefault = 1:5; end

R    = poolLoadRuns();
days = unique({R.day},'stable'); nP = numel(days);
COND = {'8mM','2mM'}; nEl = 5;

FEAT = { 'A','F'; 'Am','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'steady','F'; 'vall2_dy','F'; 'rsA','F'; 'restretchSlopeStart','F'
         'ktr','0'; 'rsK','0'; 'rsT0','0'; 'rsR2','0'; 't0','0'
         'peak1_t','0'; 'vall_t','0'; 'vall2_t','0'; 'peak1_dSL','0'
         'rsK80','0'; 'rsT20','0'; 'ktr80','0'
         'rsKfit','0'; 'rsTfit','0'; 'ktrfit','0' };
% *fit = the same rates and delay from a weighted log-linear regression over
% the whole 20-80 % band instead of two crossings (Analyses/KtrSegmentSlope,
% bandRateFit): on the model's numerically rough traces a two-point rate
% inherits both points' jitter.
% ktr80 is the audited HOLD-redevelopment rate (Analyses/KtrSegmentSlope):
% rsK80's construction on the slack hold -- plateau amplitude, 20-80 % band,
% no force threshold. It replaces the free-amplitude `ktr`, whose threshold
% cut biased it 9-11 % on model traces and 1.6 % on data.
% rsK80 / rsT20 are the audited re-stretch rate and delay
% (Analyses/RestretchMetricAudit): the amplitude comes from `steady` instead of
% from the tail of the fit window, and the rate is read off the 20-80 % band of
% the rise. They supersede rsK / rsT0, which are kept for continuity -- compare
% their s8 and z columns in the report below before using either.
nF = size(FEAT,1);

%% ---- gather ------------------------------------------------------------
Y = nan(nP,2,nF,nEl); anc = nan(nP,2); ord = nan(nP,2); slp = ord; tac = ord;
for p = 1:nP
    for c = 1:2
        i = find(strcmp({R.day},days{p}) & strcmp({R.cond},COND{c}));
        for k = 1:nF
            v = double(R(i).feat.(FEAT{k,1})(:))'; m = min(nEl,numel(v));
            Y(p,c,k,1:m) = v(1:m);
        end
        anc(p,c) = R(i).A.(opts.anchor);
        ord(p,c) = R(i).order; slp(p,c) = R(i).slope; tac(p,c) = R(i).Tact;
    end
end

%% ---- 2. rundown --------------------------------------------------------
lift = ones(nP,2);
if opts.rundown
    for p = 1:nP
        e = ord(p,:)==1; l = ord(p,:)==2;
        lift(p,l) = 1 + opts.phi*abs(slp(p,e))*tac(p,e) / anc(p,l);
    end
end
isF = strcmp(FEAT(:,2),'F'); kF = find(isF);
for c = 1:2
    Y(:,c,kF,:) = bsxfun(@times, Y(:,c,kF,:), lift(:,c));
end
ancC = anc .* lift;

%% ---- 1. scale (one per preparation, from its 8 mM run) -----------------
alpha = ancC(:,1);
for c = 1:2
    Y(:,c,kF,:) = bsxfun(@rdivide, Y(:,c,kF,:), alpha);
end

%% ---- 3-4. pool and spread ---------------------------------------------
T = struct();
T.name = FEAT(:,1)'; T.cls = FEAT(:,2)';
[T.m8, T.m2] = deal(nan(nF,nEl));
[T.s8, T.s8raw, T.r, T.sr, T.z] = deal(nan(1,nF));
T.seg = repmat({1:nEl}, 1, nF);   % segments behind each s, see segmentsFor
Rp = nan(nP,nF);                                % per-prep 2 mM / 8 mM ratio
for k = 1:nF
    V8 = squeeze(Y(:,1,k,:)); V2 = squeeze(Y(:,2,k,:));
    T.m8(k,:)  = mean(V8,1,'omitnan');
    T.m2(k,:)  = mean(V2,1,'omitnan');
    sg         = segmentsFor(opts, T.name{k}, nEl);
    T.seg{k}   = sg;
    T.s8raw(k) = relSD(V8(:,sg));
    for p = 1:nP
        qq = V2(p,:)./V8(p,:); qq = qq(isfinite(qq) & qq > 0);
        if ~isempty(qq); Rp(p,k) = exp(mean(log(qq))); end
    end
    T.r(k)  = mean(Rp(:,k),'omitnan');
    T.sr(k) = std(Rp(:,k),'omitnan')/abs(T.r(k));
end
T.s8 = T.s8raw;

% LEVEL vs SHAPE. Under a per-preparation scale every force feature carries the
% same ATP effect on the LEVEL, so nine force ratios are nine copies of one
% number. Dividing each by the level (taken from `steady`) leaves the SHAPE
% ratio, which FeatureWeighting S3 found flat. The objective therefore wants
% ONE level term, not nine.
kS = find(strcmp(T.name,'steady'),1);
T.level  = mean(Rp(:,kS));  T.sLevel = std(Rp(:,kS))/T.level;
Sp = Rp ./ Rp(:,kS);
T.rShape = mean(Sp,1,'omitnan');
T.srShape = std(Sp,0,1,'omitnan')./abs(T.rShape);

%% ---- 4b. empirical-Bayes shrinkage of the spreads ---------------------
%   s^2_moderated = (d0*s0^2 + d*s^2) / (d0 + d),   d = n-1 = 2
% The prior (d0, s0^2) is ESTIMATED from the spread of the observed s^2 rather
% than fixed. This matters: a fixed prior pulls every feature toward the
% middle, which would inflate `ktr` (11 %) toward the median of a class that
% also holds `vall_t` (84 %) -- features that plainly are NOT exchangeable.
% Matching the observed dispersion of log(s^2) against its chi-square sampling
% value drives d0 -> 0 exactly when the features genuinely differ, so the
% shrinkage disables itself instead of manufacturing agreement. (This is
% limma's moderated variance / fitFDist, the standard low-replicate remedy.)
% THE OTHER END OF THE SAME ESTIMATOR. fitFDist can also return d0 = Inf --
% "the observed variances are no MORE dispersed than sampling alone explains",
% i.e. shrink everything to the common prior. That happened the moment the
% peak-detector fix (Analyses/RestretchPeakAudit) removed the two wild
% features: with peak1_dSL and vall_t down from 52 %/57 % to 18 %/19 %, the
% remaining 20 spreads span only 4-38 % and read as exchangeable. The formula
% has no limit at d0 = Inf and silently produced s = NaN for every feature,
% which would have poisoned every floor and weight downstream.
%
% Two things are needed. (1) The correct limit: as d0 -> Inf the moderated
% variance is exactly the prior s0. (2) A CAP, because collapsing all twenty
% features onto one common variance -- making every weight equal -- is far too
% strong a claim from n = 3, and it would undo the whole point of weighting by
% measured precision. D0_CAP = 5 (2.5x the data's own d = 2) keeps at least
% ~29 % of each feature's own variance while still compressing the range. This
% is a judgement call; it is reported below so it can be argued with.
D0_CAP = 5;
T.d0 = 0; T.s0 = NaN; T.d0raw = 0;
if opts.shrink
    dd = nP-1;
    m  = isfinite(T.s8raw) & T.s8raw > 0;
    [T.d0, T.s0] = fitFDist(T.s8raw(m).^2, dd);
    T.d0raw = T.d0;
    if ~isfinite(T.d0) || T.d0 > D0_CAP
        T.d0 = D0_CAP;
        warning('buildPooledTarget:d0capped', ...
                ['fitFDist returned d0 = %g (variances look exchangeable); ' ...
                 'capped at %g so the weights do not collapse to uniform.'], T.d0raw, D0_CAP);
    end
    T.s8(m) = sqrt((T.d0*T.s0 + dd*T.s8raw(m).^2)/(T.d0 + dd));
end
T.z = abs(T.r-1) ./ sqrt(T.s8.^2 + T.sr.^2);

%% ---- how far is each preparation from the pool it belongs to? ----------
% The objective currently fits ONE preparation (03/27). If that preparation
% sits at chi-2 ~ 1 against the pooled mean it is a fair stand-in; the further
% it is, the more of the fit is that fibre's idiosyncrasy.
T.chi2Prep = nan(nP,1); T.chi2PrepK = nan(nP,nF);
for p = 1:nP
    for k = 1:nF
        v = squeeze(Y(p,1,k,:))'; m = T.m8(k,:);
        g = isfinite(v) & isfinite(m) & m ~= 0;
        if ~any(g); continue; end
        T.chi2PrepK(p,k) = mean(((v(g)-m(g))./m(g)).^2) / T.s8(k)^2;
    end
    T.chi2Prep(p) = mean(T.chi2PrepK(p,:),'omitnan');
end

%% ---- the absolute level, kept ------------------------------------------
T.scale = struct('anchor',opts.anchor, 'perPrep',ancC(:,1)', 'raw',anc(:,1)', ...
                 'mean',mean(ancC(:,1)), 'cv',std(ancC(:,1))/mean(ancC(:,1)), ...
                 'cvRaw',std(anc(:,1))/mean(anc(:,1)), 'lift',lift, 'days',{days});
T.opts = opts; T.nP = nP;

%% ---- report ------------------------------------------------------------
if nargout == 0 || true
    fprintf('=== pooled target: %d preparations, anchor %s, rundown %s, shrink %s ===\n', ...
        nP, opts.anchor, onoff(opts.rundown), onoff(opts.shrink));
    fprintf('\nabsolute level (the one number normalisation removes, kept here):\n');
    fprintf('  per prep %s : %s kPa\n', opts.anchor, strtrim(sprintf('%.1f ', anc(:,1))));
    fprintf('  after rundown correction     : %s kPa\n', strtrim(sprintf('%.1f ', ancC(:,1))));
    fprintf('  CV  %.1f%% raw  ->  %.1f%% corrected   (mean %.1f kPa)\n', ...
        100*T.scale.cvRaw, 100*T.scale.cv, T.scale.mean);
    fprintf('  -> rundown explains %.0f%% of the apparent between-prep force variance\n', ...
        100*(1 - T.scale.cv^2/T.scale.cvRaw^2));

    if opts.shrink
        fprintf('\nvariance shrinkage: prior d0 = %.2f against d = %d from the data\n', T.d0, nP-1);
        if T.d0 < 0.5
            fprintf('  -> near zero: the features are NOT exchangeable, so the estimator\n');
            fprintf('     correctly declines to shrink. s (used) ~ s (raw).\n');
        end
    else
        fprintf('\nvariance shrinkage: OFF -- s (used) = s (raw), 2 dof each.\n');
    end
    fprintf(['CAVEAT on `steady`: it is the post-re-stretch isometric at ML 1.10 and\n' ...
             '  the anchor is the PRE-slack isometric at the same length, so their\n' ...
             '  ratio is close to a self-normalisation. Its raw 4 %% is a within-run\n' ...
             '  repeatability, not a between-preparation spread -- do not hand it the\n' ...
             '  largest weight in the objective. Shrinkage corrects this on its own.\n']);

    fprintf('\nATP effect, split into the one LEVEL term and the SHAPE residuals:\n');
    fprintf('  level (from `steady`) : x%.3f +- %.1f%%  (per prep: %s)\n', ...
        T.level, 100*T.sLevel, strtrim(sprintf('%.3f ', Rp(:,kS))));

    fprintf('\n%-22s %4s %10s %8s %8s %9s %9s %7s %9s\n', ...
        'feature','cls','pooled 8mM','s (raw)','s (used)','2/8 ratio','shape r','z','w=1/s^2');
    for k = 1:nF
        if strcmp(T.cls{k},'F'); sh = sprintf('%9.3f', T.rShape(k)); else; sh = sprintf('%9s','-'); end
        fprintf('%-22s %4s %10.4g %7.1f%% %7.1f%% %9.3f %s %7.2f %9.1f\n', ...
            T.name{k}, T.cls{k}, mean(T.m8(k,:),'omitnan'), 100*T.s8raw(k), ...
            100*T.s8(k), T.r(k), sh, T.z(k), 1/T.s8(k)^2);
    end
    fprintf(['\nz >= 2 : the ATP effect on this feature is resolved by the data and\n' ...
             '         belongs in the objective as a RATIO term.\n' ...
             'z <  1 : flat to within noise -- a ratio term there only adds a constant.\n' ...
             'shape r near 1.000 : the force feature carries no ATP information beyond\n' ...
             '         the level, so its ratio term would be redundant with `level`.\n']);

    fprintf('\nweight dynamic range 1/s^2 : %.0fx raw  ->  %.0fx moderated\n', ...
        max(1./T.s8raw.^2)/min(1./T.s8raw.^2), max(1./T.s8.^2)/min(1./T.s8.^2));

    fprintf('\n=== is fitting ONE preparation the same as fitting the pool? ===\n');
    fprintf('chi-2 of each preparation against the pooled mean (per feature, 8 mM):\n');
    for p = 1:nP
        [~, wk] = max(T.chi2PrepK(p,:));
        fprintf('  %-7s  chi2 = %5.2f    worst feature: %s (%.1f)\n', ...
            days{p}, T.chi2Prep(p), T.name{wk}, T.chi2PrepK(p,wk));
    end
    fprintf(['  chi2 ~ 0.67 is what three preparations scattering around their own\n' ...
             '  mean must give ((n-1)/n). Above that, the preparation is not a fair\n' ...
             '  stand-in for the pool and fitting it alone imports its idiosyncrasy.\n']);
end
end

%% ======================= helpers ========================================
function sg = segmentsFor(opts, name, nEl)
%SEGMENTSFOR  Which slack cycles this feature's spread is measured over.
    s = opts.segments;
    if isstruct(s)
        if isfield(s, name); sg = s.(name); else; sg = opts.segDefault; end
    else
        sg = s;
    end
    sg = sg(sg >= 1 & sg <= nEl);
    if isempty(sg); sg = 1:nEl; end
end

function s = relSD(V)
    mu = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
    g = isfinite(mu) & abs(mu) > 0 & isfinite(sd);
    if ~any(g); s = NaN; return; end
    s = sqrt(mean((sd(g)./abs(mu(g))).^2));
end

function [d0, s0] = fitFDist(s2, df)
% Method-of-moments prior for a scaled-inverse-chi-square on s2 with df dof
% (Smyth 2004, limma::fitFDist). Returns d0 = prior dof, s0 = prior variance.
% d0 -> 0 when the observed s2 scatter exceeds what sampling alone explains,
% i.e. when the features are not exchangeable and shrinkage is unwarranted.
    z    = log(s2(:));
    e    = z - psi(df/2) + log(df/2);          % unbiased for log(sigma^2)
    ez   = mean(e);
    vz   = var(e)*numel(e)/(numel(e)-1) - psi(1, df/2);   % excess dispersion
    if vz <= 0
        d0 = Inf; s0 = exp(ez); return          % no excess: shrink fully
    end
    d0 = 2*trigammaInv(vz);
    s0 = exp(ez + psi(d0/2) - log(d0/2));
end

function y = trigammaInv(x)
% Solve psi'(y) = x by Newton iteration on 1/y (Smyth's initialisation).
    y = 0.5 + 1/x;
    for it = 1:50
        t  = psi(1,y) - x;
        dy = t / psi(2,y);
        y  = y - dy;
        if abs(dy) < 1e-10*max(y,1); break; end
    end
    y = max(y, 1e-8);
end

function o = onoff(b); if b; o = 'on'; else; o = 'off'; end; end
