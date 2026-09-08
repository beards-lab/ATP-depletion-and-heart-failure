function R = poolLoadRuns(forceReload)
% POOLLOADRUNS  One code path that loads every 2026 protocol run and returns
% the raw material every DataPooling script needs: the traces, the candidate
% normalisation anchors, and a fresh feature extraction.
%
%   R = poolLoadRuns()       % uses results/pool_raw.mat if present
%   R = poolLoadRuns(true)   % force re-read from the raw files
%
% R is a struct array, one element per RUN (= one day x one condition), with
%
%   .day .cond .order      day tag, '8mM'|'2mM'|'PNB', 1 or 2 = order in session
%   .feat                  extractSlackAttributes on the curated slack .mat
%   .A                     struct of candidate ANCHORS (kPa), see below
%   .slope .Tact           within-run rundown slope (kPa/s) and activation time
%   .t .L .F               the merged trace (decimated to 1 kHz)
%   .fv2                   struct of FV2 ramp features, or [] if the day has none
%
% ANCHORS (.A) -- every one of these is a force level that could stand in for
% "how strong is this preparation", i.e. a candidate estimator of the
% cross-sectional-area nuisance:
%
%   Fss22    isometric plateau at ML 1.10 (SL 2.2) measured in the high-rate
%            burst file immediately BEFORE the slack train. "the 8 mM steady
%            state before the slack".
%   Fiso20a  quiescent plateau at ML 1.00, after the ktr burst, before the
%            stairs  (t in [71.50 72.35] of the merged file)
%   Fiso20b  quiescent plateau at ML 1.00, after the slack train
%            (t in [77.70 78.30])
%   Fpeak    peak of the activation transient (t in [64 70])
%   Fsteady  mean of the five slack-segment quasi-isometric steady levels
%            (this is the anchor RunFeatureSpread/SlackDataAnalysis used)
%
% The passive (PNB+Mava) run of each day is loaded as its own element so the
% passive level can be subtracted or scored on the same footing.

if nargin < 1 || isempty(forceReload); forceReload = false; end

here   = fileparts(mfilename('fullpath'));
root   = fullfile(here, '..', '..');
resDir = fullfile(here, 'results');
cache  = fullfile(resDir, 'pool_raw.mat');
if ~isfolder(resDir); mkdir(resDir); end

if ~forceReload && isfile(cache)
    S = load(cache, 'R'); R = S.R; return
end

oldDir = pwd; cd(root); c = onCleanup(@() cd(oldDir));

% day, cond, order, merged trace, curated slack .mat, raw burst .txt
CFG = { ...
 '03/27','8mM',1,'03 27 2026 M/02_Merged_8mM_Active.txt',        'protocol_03_27_2026_8mM_slack.mat',        '03 27 2026 M/8mM_slack.txt'
 '03/27','2mM',2,'03 27 2026 M/03_Merged_2mM_Active.txt',        'protocol_03_27_2026_2mM_slack.mat',        '03 27 2026 M/2mM_slack.txt'
 '03/27','PNB',3,'03 27 2026 M/06_Merged_8mM_Active_PNB_Mava.txt','protocol_03_27_2026_ActivePNBMava_slack.mat','03 27 2026 M/8mM_slack_PNB_Mava.txt'
 '04/03','8mM',1,'04 03 2026 F/02_Merged_8mM_Active.txt',        'protocol_04_03_2026_8mM_slack.mat',        '04 03 2026 F/8mM_slack.txt'
 '04/03','2mM',2,'04 03 2026 F/03_Merged_2mM_Active.txt',        'protocol_04_03_2026_2mM_slack.mat',        '04 03 2026 F/2mM_slack.txt'
 '04/03','PNB',3,'04 03 2026 F/04_Merged_8mM_Active_PNB_Mava.txt','protocol_04_03_2026_ActivePNBMava_slack.mat','04 03 2026 F/8mM_slack_PNB_Mava.txt'
 '04/10','2mM',1,'04 10 2026 Male 2-8/02_Merged_2mM_Active.txt',  'protocol_04_10_2026_2mM_slack.mat',        '04 10 2026 Male 2-8/2mM_slack.txt'
 '04/10','8mM',2,'04 10 2026 Male 2-8/03_Merged_8mM_Active.txt',  'protocol_04_10_2026_8mM_slack.mat',        '04 10 2026 Male 2-8/8mM_slack.txt'
 '04/10','PNB',3,'04 10 2026 Male 2-8/04_Merged_8mM_Active_PNB_Mava.txt','protocol_04_10_2026_ActivePNBMava_slack.mat','04 10 2026 Male 2-8/8mM_slack_PNB_Mava.txt'};

WIN = struct('zA',[71.50 72.35], 'zB',[77.70 78.30], 'pk',[64 70]);

R = struct([]);
for i = 1:size(CFG,1)
    r = struct('day',CFG{i,1}, 'cond',CFG{i,2}, 'order',CFG{i,3});

    % ---- merged trace: anchors on the ML 1.00 plateaus + rundown slope ----
    M = readmatrix(fullfile('data', CFG{i,4}));
    mono = [true; diff(M(:,1)) > 0]; M = M(mono,:);
    t = M(:,1); L = M(:,2); F = M(:,3);
    keep = t > 60 & t < 85;                      % the protocol window
    r.t = t(keep); r.L = L(keep); r.F = F(keep);

    wA = t >= WIN.zA(1) & t <= WIN.zA(2);
    wB = t >= WIN.zB(1) & t <= WIN.zB(2);
    A.Fiso20a = median(F(wA));
    A.Fiso20b = median(F(wB));
    A.Fpeak   = max(movmean(F(t > WIN.pk(1) & t < WIN.pk(2)), 201));

    % within-run rundown: one line through both quiescent ML 1.00 plateaus
    pf = polyfit(t(wA|wB), F(wA|wB), 1);
    r.slope = pf(1);                             % kPa/s (negative)

    Fs = movmean(F, 201); thr = 0.2*prctile(Fs, 99);
    on = find(Fs > thr, 1); off = find(Fs > thr, 1, 'last');
    r.Tact = t(off) - t(on);                     % activation duration (s)

    % ---- burst file: the pre-slack isometric at ML 1.10 (SL 2.2) ---------
    B = readmatrix(fullfile('data', CFG{i,6}));
    i0 = find(B(:,1) >= -0.1, 1);
    tb = B(i0:end,1)/1000; Lb = B(i0:end,2); Fb = B(i0:end,3);
    Lbase = median(Lb(tb >= 0 & tb < 0.040));
    iEv   = find(Lb < Lbase - 0.005, 1);         % first length change of any kind
    if isempty(iEv); iEv = numel(tb); end
    wPre  = tb >= 0.002 & tb < tb(iEv) - 0.002;
    A.Fss22 = median(Fb(wPre));
    r.tPre  = tb(iEv);
    r.burst = [tb Lb Fb];

    % ---- curated slack file: fresh feature extraction, one code path -----
    S = load(fullfile('data', CFG{i,5}));
    r.feat = extractSlackAttributes(S.datatable(:,1), S.datatable(:,3), ...
                                    S.datatable(:,2), S.velocitytable, [], [], false, true);
    A.Fsteady = mean(double(r.feat.steady), 'omitnan');
    r.A = A;

    % ---- FV2 isovelocity ramp, when the day has one ----------------------
    r.fv2 = extractFV2(tb, Lb, Fb);

    if isempty(R); R = r; else; R(end+1) = r; end %#ok<AGROW>
end

save(cache, 'R');
fprintf('poolLoadRuns: %d runs -> %s\n', numel(R), cache);
end

% ------------------------------------------------------------------------
function fv = extractFV2(t, L, F)
% The 04-series slack bursts open with a load-clamp manoeuvre: fast step
% release, then a 34 ms ISOVELOCITY RAMP at -2 ML/s, then the slack. The ramp
% plateau is a genuine sub-Vmax force-velocity point measured in the SAME
% activation as the slack train. The 03/27 files have no such ramp.
    fv = [];
    dt = median(diff(t));
    v  = gradient(movmean(L, max(1, round(5e-4/dt))), t);
    clip   = t >= 0 & t <= 0.45;
    isRamp = clip & v < -0.5 & v > -5;           % between the holds and the step
    d  = diff([0; isRamp; 0]); rs = find(d==1); re = find(d==-1)-1;
    if isempty(rs); return; end
    [dur, iL] = max(t(re) - t(rs));
    if dur < 0.010; return; end                  % no sustained ramp -> no FV2
    t1 = t(rs(iL)); t2 = t(re(iL));
    Fss = median(F(t >= 0.002 & t < 0.040));
    w   = t >= t1 + 0.005 & t <= t2;             % skip the step transient
    fv = struct('v', (L(re(iL))-L(rs(iL)))/(t2-t1), 'dur', dur, ...
                't1', t1, 't2', t2, 'Fss', Fss, ...
                'Fplateau', mean(F(w)), 'Fnorm', mean(F(w))/Fss, ...
                'L1', L(rs(iL)), 'L2', L(re(iL)));
end
