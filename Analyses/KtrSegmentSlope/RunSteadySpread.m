% RunSteadySpread.m -- what is the 17.1 % deadband on `steady` actually made of?
%
% The pooled target gives every force-class feature a between-preparation
% spread of 17-18 %, and the floor rule fl = (1/3)*nEl*s^2 turns that into a
% deadband the model currently sits well inside. But the report also says
%   Fref = 69.35 kPa, CV 15.5 %, per prep 79.8 / 70.0 / 58.3
% so the suspicion is that the deadband is the between-HEART difference in
% absolute force, not any uncertainty in a plateau within one heart.
%
% This script separates the two, with no simulation:
%   raw s      : spread of the feature across preps, as pooled
%   level s    : spread of each prep's own mean (the absolute force level)
%   shape s    : spread after dividing each prep by its own mean
% and reports the same for the other force-class features, because if the
% diagnosis holds it is not a `steady` problem, it is a whole-block problem.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunSteadySpread.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end

DAYS = {'03_27','04_03','04_10'};
PT   = pooledTarget8mM();

% what the objective actually scores, in the order it appears in fn
FEATS = {'steady','A','peak1_y','vall_y','peak2','restretchSlopeStart', ...
         'vall2_dy','ktr','rsK80','t0','vall_t','peak1_dSL'};
FCLASS = {'steady','A','peak1_y','vall_y','peak2','restretchSlopeStart','vall2_dy'};

% per-prep normalisation exactly as buildPooledTarget does it: each prep by its
% own Fss22, then the block rescaled onto Fref
F = struct();
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    F(d).f = S.features_data;
end
fprintf('=== per-preparation `steady` (as curated, kPa) ===\n');
fprintf('%-10s %8s %8s %8s %8s %8s   %8s\n','prep','seg1','seg2','seg3','seg4','seg5','mean');
St = nan(3,5);
for d = 1:3
    St(d,:) = F(d).f.steady(:)';
    fprintf('%-10s %8.2f %8.2f %8.2f %8.2f %8.2f   %8.2f\n', strrep(DAYS{d},'_','/'), St(d,:), mean(St(d,:)));
end
lev = mean(St,2);
fprintf('level CV across preps: %.1f %%\n', 100*std(lev)/mean(lev));
Sh = St ./ lev;
fprintf('\n=== the same, each prep divided by its OWN mean (pure shape) ===\n');
for d = 1:3
    fprintf('%-10s %8.4f %8.4f %8.4f %8.4f %8.4f\n', strrep(DAYS{d},'_','/'), Sh(d,:));
end
fprintf('%-10s %8.4f %8.4f %8.4f %8.4f %8.4f   <- CV per segment\n','shape CV', std(Sh,0,1)./mean(Sh,1));

%% ---- the decomposition, for every scored feature -----------------------
fprintf('\n=== spread decomposition (relative, n = 3) ===\n');
fprintf('%-22s %9s %9s %9s %9s   %s\n','feature','s_used','s_raw','s_level','s_shape','what the deadband is');
rows = {};
for i = 1:numel(FEATS)
    nm = FEATS{i};
    X = nan(3, numel(F(1).f.(nm)));
    ok = true;
    for d = 1:3
        v = F(d).f.(nm);
        if numel(v) ~= size(X,2); ok = false; break; end
        X(d,:) = v(:)';
    end
    if ~ok; continue; end
    mu  = mean(X, 1, 'omitnan');
    s_raw = mean(std(X,0,1,'omitnan') ./ abs(mu), 'omitnan');       % per-element CV, averaged
    lv  = mean(X, 2, 'omitnan');
    s_level = std(lv) / abs(mean(lv));
    Xs  = X ./ lv;
    mus = mean(Xs,1,'omitnan');
    s_shape = mean(std(Xs,0,1,'omitnan') ./ abs(mus), 'omitnan');
    j = find(strcmp(PT.name, nm), 1);
    s_used = NaN; if ~isempty(j); s_used = PT.s(j); end
    verdict = 'shape';
    if s_level > 2*s_shape; verdict = 'LEVEL (between-heart force)'; end
    if s_used > 1.8*s_raw;  verdict = [verdict ' + moderation inflated']; end
    fprintf('%-22s %9.4f %9.4f %9.4f %9.4f   %s\n', nm, s_used, s_raw, s_level, s_shape, verdict);
    rows(end+1,:) = {nm, s_used, s_raw, s_level, s_shape}; %#ok<SAGROW>
end

%% ---- what that means for the floors -------------------------------------
fprintf('\n=== floors: current vs shape-based, for the force block ===\n');
fprintf('  rule fl = (1/3) * nEl * s^2 ;  w = lambda / s^2 ;  lambda = 0.0106294\n');
LAMBDA = 0.0106294;
fprintf('%-22s %5s %11s %11s %9s %9s\n','feature','nEl','fl_now','fl_shape','w_now','w_shape');
for i = 1:size(rows,1)
    nm = rows{i,1};
    if ~ismember(nm, FCLASS); continue; end
    nEl = 5; if strcmp(nm,'rsK80'); nEl = 4; end
    su = rows{i,2}; ss = max(rows{i,5}, 0.01);
    fprintf('%-22s %5d %11.5f %11.5f %9.4f %9.4f\n', nm, nEl, (1/3)*nEl*su^2, (1/3)*nEl*ss^2, ...
            LAMBDA/su^2, LAMBDA/ss^2);
end

save(fullfile(resdir,'steady_spread.mat'), 'rows','St','Sh','FEATS','FCLASS');
fprintf('\nDONE\n');
