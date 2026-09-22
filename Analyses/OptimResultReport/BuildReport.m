% BuildReport.m — assemble the single-file HTML report.
%
% Consumes what RunDataEnvelope.m and RunFitReport.m already wrote, so it is
% cheap to re-run: no model is evaluated here. Every figure is inlined as a
% base64 data URI, so results/report_<tag>.html is ONE self-contained file that
% can be mailed or dropped in a shared folder with nothing else attached.
%
%   run('Analyses/OptimResultReport/RunDataEnvelope.m')   % ~40 s
%   run('Analyses/OptimResultReport/RunFitReport.m')      % ~4 min
%   run('Analyses/OptimResultReport/BuildReport.m')       % ~2 s

clearvars -except TAG_NEW TAG_BASE; clc;          % keep caller-set knobs
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');

if ~exist('TAG_NEW','var');  TAG_NEW  = 'spreadw_opt';  end
if ~exist('TAG_BASE','var'); TAG_BASE = 'realign4_opt'; end

L  = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
EV = load(fullfile(resDir, 'data_envelope.mat')); ENV = EV.ENV;
R  = L.R; nm = L.nm; c1 = L.c1; c2 = L.c2;

%% ---- table 1: per-feature cost ----------------------------------------
gain = sum(c1) - sum(c2);
[~, ord] = sort(c2 - c1, 'ascend');            % biggest improvement first
rows = strings(0,1);
for i = ord
    d = c2(i) - c1(i);
    if abs(d) < 5e-4 && max(c1(i),c2(i)) < 5e-3; continue; end   % inert terms
    if d < 0; cls = ' class="good"'; else; cls = ' class="bad"'; end
    rows(end+1) = sprintf( ...
        '<tr><td><code>%s</code></td><td>%.4f</td><td>%.4f</td><td%s>%+.4f</td><td>%.1f %%</td></tr>', ...
        nm{i}, c1(i), c2(i), cls, d, -100*d/gain);                      %#ok<SAGROW>
end
rows(end+1) = sprintf( ...
    '<tr><td><b>total</b></td><td><b>%.4f</b></td><td><b>%.4f</b></td><td><b>%+.4f</b></td><td></td></tr>', ...
    sum(c1), sum(c2), sum(c2)-sum(c1));
TABLE_COST = strjoin(rows, newline);

iFV = find(strcmp(nm, 'FV_normpowerAvg'), 1);
FV_SHARE_OF_GAIN = sprintf('%.1f %%', 100*(c1(iFV)-c2(iFV))/gain);

%% ---- table 2: distance from the envelope in SDs ------------------------
NMZ = {'rsT0','vall_t','peak1_y','peak1_dSL','ktr','Am','rsK','steady','A','rsA', ...
       'vall2_dy','vall_y','restretchSlopeStart','peak2'};
rows = strings(0,1); Z = nan(1,numel(NMZ));
iR = ENV.iRef;
for j = 1:numel(NMZ)
    f = NMZ{j};
    if ~isfield(ENV.feat, f); continue; end
    E = ENV.feat.(f); v = nan(1,2); z = NaN; zr = NaN;
    for k = 1:2
        if ~isfield(R(k).fmodel, f); continue; end
        x = double(R(k).fmodel.(f)(:))'; x = x(1:numel(E.mean));
        if E.cls == 'F'; x = x / mean(double(R(k).fmodel.steady),'omitnan'); end
        v(k) = mean(x);
        if k == 2
            z  = mean((x - E.mean)      ./E.sd, 'omitnan');   % vs the 3-prep mean
            zr = mean((x - E.raw(iR,:)) ./E.sd, 'omitnan');   % vs 03/27, the fitted prep
        end
    end
    Z(j) = z;
    if abs(z)  >= 2; cls  = ' class="bad"'; else; cls  = ''; end
    if abs(zr) >= 2; clsr = ' class="bad"'; else; clsr = ''; end
    rows(end+1) = sprintf( ...
        ['<tr><td><code>%s</code></td><td>%.4g</td><td class="ref">%.4g</td><td>%.4g</td>' ...
         '<td>%.4g</td><td%s>%+.1f</td><td%s>%+.1f</td></tr>'], ...
        f, mean(E.mean), mean(E.raw(iR,:)), v(1), v(2), cls, z, clsr, zr);   %#ok<SAGROW>
end
[~, oz] = sort(abs(Z), 'descend', 'MissingPlacement','last');
TABLE_Z = strjoin(rows(oz), newline);

%% ---- predictions: ktr protocol and staircase ---------------------------
% Recomputed here rather than written into the template: these change with every
% refit, and a stale hardcoded number in a report is worse than no number.
K = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
dtk = K.datatable; vtk = K.velocitytable; t0k = vtk(2,1);
wk  = dtk(:,1) > 1.5965 & dtk(:,1) < 1.6005;  ovD = max(dtk(wk,3));
tr  = 1.6035;
% Guarded: `end-100` underflows when the recovery window has <101 samples, and
% find() returns empty when F never reaches 90 % of its own tail mean. Both
% assign empty into a scalar and abort the whole report. NaN instead -- a
% missing t90 is a blank cell, not a lost report.
t90 = @(t,F) min([reshape(t(find(F >= 0.9*mean(F(max(1,end-100):end)), 1)), 1, []), NaN]) - tr;
mdk = dtk(dtk(:,1) > tr, :);  t90D = t90(mdk(:,1), movmean(mdk(:,3),21));

S  = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_stairs.mat'));
ds = S.datatable; vs = S.velocitytable;
tst = vs(vs(:,2)~=0, 1); tst = tst(diff([0;tst]) > 0.02);
t0s = vs(find(vs(:,2)~=0,1), 1);

ovM = nan(1,2); t90M = nan(1,2); stM = nan(1,2);
for k = 1:2
    o = R(k).out_ktr; F0 = interp1(o.t,o.Force,t0k,'linear','extrap'); Fn = o.Force/F0;
    ovM(k)  = max(Fn(o.t > 1.5965 & o.t < 1.6005));
    t90M(k) = t90(o.t(o.t>tr), Fn(o.t>tr));
    o = R(k).out_stairs; F0 = interp1(o.t,o.Force,t0s,'linear','extrap'); Fn = o.Force/F0;
    e = nan(1,numel(tst));
    for j = 1:numel(tst)
        w2 = ds(:,1) > tst(j) & ds(:,1) < tst(j)+0.06;  if ~any(w2); continue; end
        mm = o.t > tst(j) & o.t < tst(j)+0.06;
        e(j) = max(Fn(mm))/max(movmean(ds(w2,3),9)) - 1;
    end
    stM(k) = 100*mean(e,'omitnan');
end
bad = @(tf) repmat(' class="bad"', 1, tf);
TABLE_PRED = strjoin({ ...
  sprintf(['<tr><td><code>ktr</code> protocol, re-stretch overshoot</td><td>%.3f</td>' ...
           '<td>%.3f (%+.0f %%)</td><td>%.3f (%+.0f %%)</td></tr>'], ovD, ...
           ovM(1), 100*(ovM(1)/ovD-1), ovM(2), 100*(ovM(2)/ovD-1)), ...
  sprintf(['<tr><td><code>ktr</code> protocol, t to 90 %% of final</td><td>%.1f ms</td>' ...
           '<td%s>%.1f ms (%.1f&times; fast)</td><td%s>%.1f ms (%.1f&times; fast)</td></tr>'], ...
           1000*t90D, bad(true), 1000*t90M(1), t90D/t90M(1), bad(true), 1000*t90M(2), t90D/t90M(2)), ...
  sprintf(['<tr><td>staircase, mean step-peak overshoot</td><td>&mdash;</td>' ...
           '<td>%+.1f %%</td><td%s>%+.1f %%</td></tr>'], stM(1), bad(stM(2) > stM(1)), stM(2)) ...
  }, newline);

%% ---- envelope band width ----------------------------------------------
rel = ENV.tc.Fsd(:,1) ./ abs(ENV.tc.Fmean(:,1));
BAND_MED = sprintf('%.1f %%', 100*median(rel,'omitnan'));
BAND_P90 = sprintf('%.1f %%', 100*prctile(rel,90));

%% ---- fill and write ----------------------------------------------------
s = struct();
s.TITLE      = 'Cross-bridge model — spread-weighted refit';
s.TAG_NEW    = TAG_NEW;
s.TAG_BASE   = TAG_BASE;
s.DATE       = datestr(now, 'yyyy-mm-dd');                              %#ok<TNOW1,DATST>
s.COST_BASE  = sprintf('%.4f', sum(c1));
s.COST_NEW   = sprintf('%.4f', sum(c2));
s.COST_PCT   = sprintf('%+.1f %%', 100*(sum(c2)-sum(c1))/sum(c1));
s.FV_SHARE_OF_GAIN = FV_SHARE_OF_GAIN;
s.SCALES     = strjoin(arrayfun(@(x) sprintf('%.1f',x), ENV.scale, 'uni',0), ' / ');
s.BAND_MED   = BAND_MED;
s.BAND_P90   = BAND_P90;
s.TABLE_COST = TABLE_COST;
s.TABLE_Z    = TABLE_Z;
s.TABLE_PRED = TABLE_PRED;
% headline feature numbers, so no claim in the prose is hand-typed
zget = @(f) iZ(ENV, R(2), f);
[zM,zR,pc] = zget('rsK');     s.RSK_Z = sprintf('%+.1f',zM); s.RSK_ZREF = sprintf('%+.1f',zR);
                              s.RSK_PCT = sprintf('%+.0f %%',pc);
[zM,~,pc]  = zget('peak1_y'); s.P1_Z  = sprintf('%+.1f',zM); s.P1_PCT  = sprintf('%+.0f %%',pc);
[zM,~,pc]  = zget('ktr');     s.KTR_Z = sprintf('%+.1f',zM); s.KTR_PCT = sprintf('%+.0f %%',pc);
[zM,~,~]   = zget('rsT0');    s.RST0_Z = sprintf('%+.1f',zM);
s.KTR_FAST   = sprintf('%.1f', t90D/t90M(2));
s.KTR_OVER   = sprintf('%+.0f %%', 100*(ovM(2)/ovD-1));
s.STAIRS_OVER= sprintf('%+.1f %%', stM(2));
s.STAIRS_BASE= sprintf('%+.1f %%', stM(1));
s.FIG_COSTSHIFT    = fullfile(resDir,'fit_costshift.png');
s.FIG_FV           = fullfile(resDir,'fit_fv.png');
s.FIG_ENVELOPE     = fullfile(resDir,'fit_envelope_overlay.png');
s.FIG_FEATURES     = fullfile(resDir,'fit_features.png');
s.FIG_SLACK_FULL   = fullfile(resDir,'fit_slack_full.png');
s.FIG_SLACK_DETAIL = fullfile(resDir,'fit_slack_detail.png');
s.FIG_KTR          = fullfile(resDir,'fit_ktr_pred.png');
s.FIG_STAIRS       = fullfile(resDir,'fit_stairs_pred.png');
% the slide block (MakeSlideFigures.m) -- 03/27 only, oversized type
s.SLIDE_FV     = fullfile(resDir,'slide_fv.png');
s.SLIDE_STAIRS = fullfile(resDir,'slide_stairs.png');
s.SLIDE_RSALL  = fullfile(resDir,'slide_restretch_all.png');
s.SLIDE_RSDET  = fullfile(resDir,'slide_restretch_det.png');
s.PPTX_NAME    = sprintf('slide_%s.pptx', TAG_NEW);

missing = {}; needSlides = false;
f = fieldnames(s);
for i = 1:numel(f)
    if ~(startsWith(f{i},'FIG_') || startsWith(f{i},'SLIDE_')); continue; end
    if ~isfile(s.(f{i}))
        missing{end+1} = f{i};                                      %#ok<SAGROW>
        needSlides = needSlides || startsWith(f{i},'SLIDE_');
    end
end
if ~isempty(missing)
    hint = 'run RunFitReport.m first';
    if needSlides; hint = [hint ' (and MakeSlideFigures.m for the slide block)']; end
    error('BuildReport:missingFigures', '%s — missing: %s', hint, strjoin(missing, ', '));
end

out = buildHtmlReport(fullfile(here,'report_template.html'), ...
                      fullfile(resDir, sprintf('report_%s.html', TAG_NEW)), s);
fprintf('open: %s\n', out);

function [zMean, zRef, pctOverMean] = iZ(ENV, Rk, f)
% Model's mean distance from the 3-prep envelope, in envelope SDs, and from the
% fitted preparation; plus the plain percentage over the envelope mean.
    E = ENV.feat.(f);
    x = double(Rk.fmodel.(f)(:))'; x = x(1:numel(E.mean));
    if E.cls == 'F'; x = x / mean(double(Rk.fmodel.steady),'omitnan'); end
    zMean = mean((x - E.mean)              ./ E.sd, 'omitnan');
    zRef  = mean((x - E.raw(ENV.iRef,:))   ./ E.sd, 'omitnan');
    pctOverMean = 100*(mean(x)/mean(E.mean) - 1);
end
