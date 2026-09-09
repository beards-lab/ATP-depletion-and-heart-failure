% RunFeatureNoise.m -- how much of each feature's "between-preparation
% spread" is actually the transducer?
%
% The objective's weights are w = lambda / s^2, where s is the between-
% preparation relative SD of the feature (Analyses/FeatureWeighting). That
% rule is right only if s measures how much PREPARATIONS differ. It does not:
%
%       s^2  =  s_bio^2  +  s_noise^2
%
% and s_noise is the part contributed by the 1.2-1.6 kPa of white transducer
% noise that RunNoiseBudget measured in situ (and found unchanged with the
% cross-bridges chemically switched off). A feature located as the position of
% an extremum on a noisy trace -- peak1_t, vall_t, peak1_dSL -- can have
% s_noise comparable to s, in which case its weight is not describing biology
% and the same preparation would not reproduce the number twice.
%
% MEASUREMENT. Take each run's real trace, add NREP independent realisations
% of ITS OWN measured noise (sigma read off that run's own plateaus, so no
% noise model is assumed), re-extract every feature through the production
% extractor, and take the SD across realisations. That is s_noise, per
% feature, per run, with the biology held exactly fixed.
%
% WHAT IT IS FOR. Two things, both conservative:
%   (1) a FLOOR on s. No feature can be known better than its own noise, so a
%       weight built on s < s_noise is too large. Floor s at s_noise.
%   (2) a RETIREMENT test. A feature whose s_noise is most of its s is not
%       measuring a preparation, it is measuring the rig.
%
% Output -> results/feature_noise.mat / .txt
%
% Usage:  NREP = 200; run('Analyses/RecoveryOscillation/RunFeatureNoise.m')

clearvars -except NREP; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('NREP','var') || isempty(NREP); NREP = 200; end
rng(7);

CFG = { '03/27','8mM','protocol_03_27_2026_8mM_slack.mat'
        '04/03','8mM','protocol_04_03_2026_8mM_slack.mat'
        '04/10','8mM','protocol_04_10_2026_8mM_slack.mat'
        '03/27','2mM','protocol_03_27_2026_2mM_slack.mat'
        '04/03','2mM','protocol_04_03_2026_2mM_slack.mat'
        '04/10','2mM','protocol_04_10_2026_2mM_slack.mat'};

% every feature the objective reads off the slack protocol
FEATS = {'ktr','A','t0','Am','restretchSlopeStart','peak1_y','peak1_t','peak1_dSL', ...
         'vall_y','vall_t','peak2','steady','vall2_dy','vall2_t','rsK','rsA','rsT0', ...
         'rsR2','rsK80','rsT20','ktr_rmse','coolDownLS','t0_crossing','SLdiff'};

N = struct([]);
for i = 1:size(CFG,1)
    S  = load(fullfile('data', CFG{i,3}));
    dt = S.datatable;  vt = S.velocitytable;
    t  = dt(:,1);  L = dt(:,2);  F = dt(:,3);

    sig = traceSigma(t, F, vt);
    f0  = extractSlackAttributes(t, F, L, vt, [], [], false, true);

    Fmat = repmat(F, 1, NREP) + sig*randn(numel(F), NREP);
    rep  = cell(1, NREP);
    parfor r = 1:NREP
        rep{r} = extractSlackAttributes(t, Fmat(:,r), L, vt, [], [], false, true);
    end

    n = struct('day',CFG{i,1}, 'cond',CFG{i,2}, 'sigma',sig, 'nrep',NREP);
    for k = 1:numel(FEATS)
        fk = FEATS{k};
        if ~isfield(f0, fk); n.(fk) = NaN; n.([fk '_base']) = NaN; ...
                             n.([fk '_s14']) = NaN; n.([fk '_seg']) = NaN; continue; end
        b0 = double(f0.(fk));
        nS = numel(b0);
        % per-realisation, per-segment
        M = nan(NREP, nS);
        for r = 1:NREP
            if ~isfield(rep{r}, fk); continue; end
            v = double(rep{r}.(fk));
            M(r, 1:min(nS,numel(v))) = v(1:min(nS,numel(v)));
        end
        n.([fk '_base']) = mean(b0, 'omitnan');
        n.(fk)           = relSD(mean(M, 2, 'omitnan'));              % all segments
        if nS >= 5
            % segment 5 re-stretches to a DIFFERENT final length by the
            % SMALLEST amplitude (Analyses/RestretchMetricAudit sec. 3), so it
            % is both a different manoeuvre and the worst-conditioned one.
            n.([fk '_s14']) = relSD(mean(M(:,1:4), 2, 'omitnan'));
        else
            n.([fk '_s14']) = n.(fk);
        end
        n.([fk '_seg']) = arrayfun(@(c) relSD(M(:,c)), 1:nS);         % per segment
    end
    if isempty(N); N = n; else; N(end+1) = n; end %#ok<AGROW>
    fprintf('  %s %s : sigma = %.2f kPa, %d reps done\n', CFG{i,1}, CFG{i,2}, sig, NREP);
end

save(fullfile(resDir,'feature_noise.mat'), 'N', 'FEATS', 'NREP');
reportFeatureNoise(N, FEATS, resDir);

% ------------------------------------------------------------------------
function s = relSD(v)
%RELSD  Relative SD of a sample. The SD, not a robust surrogate: the cost is
%   quadratic, so a feature whose estimator blows up on 1 % of realisations
%   really is that unreliable to a least-squares objective, and hiding the
%   tail behind a MAD would understate exactly the failure mode being hunted.
    v = v(isfinite(v));
    if numel(v) < 5; s = NaN; return; end
    s = std(v) / max(abs(mean(v)), eps);
end

function s = traceSigma(t, F, vt)
%TRACESIGMA  White noise of a run, read off the quasi-isometric plateaus (the
%   last 60 % of each hold-long segment) by the second-difference estimator,
%   which is immune to any drift or ramp inside the window.
    seg = find(vt(:,2) < -1)';
    v = [];
    for i = 1:numel(seg)
        s4 = min(seg(i)+3, size(vt,1)); s5 = min(seg(i)+4, size(vt,1));
        ta = vt(s4,1); tb = vt(s5,1);
        w = t >= ta + 0.4*(tb-ta) & t <= tb;
        if nnz(w) > 50; v(end+1) = std(diff(F(w),2))/sqrt(6); end %#ok<AGROW>
    end
    s = median(v);
end
