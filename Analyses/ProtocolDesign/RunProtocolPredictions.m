% RunProtocolPredictions.m -- what the CURRENT model predicts for each proposed
% new measurement: commanded length (SL) trace + expected force read-out.
% Feeds the data-needs deck (2026-10). Core = robust1_seed with the k2 power
% gate R2 x (ATP/8)^p, p = log(0.54)/log(1/4) (set A of the 2026-09 status deck).
%
%   matlab -batch "run('Analyses/ProtocolDesign/RunProtocolPredictions.m')"
%
% Outputs -> results/pred_*.png, results/pred_summary.txt
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
res = fullfile(root, 'Analyses', 'ProtocolDesign', 'results'); if ~isfolder(res), mkdir(res); end
fid = fopen(fullfile(res, 'pred_summary.txt'), 'w'); say = @(varargin) fprintf(fid, varargin{:});
say('RunProtocolPredictions  %s\n', datestr(now));

CORE = 'robust1_seed'; PK2 = log(0.54) / log(1/4);
p0 = getParams(loadParams(CORE), [], true, false);
p0.UseAtpPowerGates = true; p0.MgATP_ref = 8; p0.UseAtpDetachR2 = 0; p0.UseAtpKmsrd = 0; p0.pAtp_k2 = PK2;
p0.ghostLoad = ''; p0.ghostSave = ''; p0.SaveBest = 0; p0.MaxRunTime = 900; p0.BreakOnODEUnstable = false;
p0.mods = {}; p0.g = ones(1, 0);
mf = str2func(p0.modelFcn);

C8 = [0 39 76]/255; C4 = [0.85 0.45 0.10]; C2 = [0.75 0.20 0.12]; CR = [0.45 0.45 0.45]; CT = [0 0.5 0.5];
set(groot, 'defaultAxesFontSize', 15, 'defaultLineLineWidth', 2);

%% ---- 1. The session battery: slack train at 8 / 4 / 2 mM + rundown repeat --
S = load('data/protocol_03_27_2026_8mM_slack.mat'); vtS = S.velocitytable; vtS(1,1) = vtS(2,1) - 2;
t0S = vtS(2,1);
cases = {'8 mM', 8, C8, '-', 0; '4 mM (prediction)', 4, C4, '-', 0; '2 mM', 2, C2, '-', 0; ...
         '8 mM repeat, run-down (kSE x0.65, SL -0.098 um)', 8, CR, '--', 1};
OS = cell(1, 4);
for i = 1:4
    p = p0; p.MgATP = cases{i, 2};
    OS{i} = runVT(p, vtS, mf, cases{i, 5});
    [Fss, k63, pk] = slackReads(OS{i}, vtS);
    say('SLACK %-48s Fss22 %6.2f kPa  slack k63(c2) %5.1f /s  restretch peak1 %6.2f kPa\n', cases{i,1}, Fss, k63, pk);
    OS{i}.Fss = Fss; OS{i}.k63 = k63;
end
say('  ratios vs 8 mM: 4 mM force %.3f ktr %.3f | 2 mM force %.3f ktr %.3f | rundown force %.3f ktr %.3f\n', ...
    OS{2}.Fss/OS{1}.Fss, OS{2}.k63/OS{1}.k63, OS{3}.Fss/OS{1}.Fss, OS{3}.k63/OS{1}.k63, OS{4}.Fss/OS{1}.Fss, OS{4}.k63/OS{1}.k63);
f = newFig(); tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(1); plotL(ax1, vtS, t0S, 1); ylabel('SL (\mum)'); title('Slack train (the standard battery): commanded length');
ax2 = nexttile([2 1]); hold on;
for i = 1:4, plot(OS{i}.t - t0S, OS{i}.Force, cases{i, 4}, 'Color', cases{i, 3}); end
xlabel('t from first release (s)'); ylabel('force (kPa)'); legend(cases(:, 1), 'Location', 'southeast'); grid on;
title('Expected force (current model; 4 mM is a prediction the fit never saw)');
linkaxes([ax1 ax2], 'x'); xlim([-0.2 3.6]);
print(f, fullfile(res, 'pred_slack_8_4_2.png'), '-dpng', '-r110'); close(f);

%% ---- 2. ktr protocol: 8 / 4 / 2 mM + rundown; protocol-ktr vs slack-ktr under rundown --
K = load('data/protocol_03_27_2026_velocitytable_ktr.mat'); vtK = K.velocitytable; vtK(1,1) = vtK(2,1) - 2;
t0K = vtK(end-1, 1);
OK = cell(1, 4);
for i = 1:4
    p = p0; p.MgATP = cases{i, 2};
    OK{i} = runVT(p, vtK, mf, cases{i, 5});
    Fss = interp1(OK{i}.t, OK{i}.Force, vtK(2,1) - 1e-3);
    fe = extractKtrAttributes(OK{i}.t, OK{i}.Force / Fss, vtK);
    OK{i}.Fss = Fss; OK{i}.k80 = fe.ktrP_k80; OK{i}.k63 = fe.ktrP_k63;
    say('KTR   %-48s Fss20 %6.2f kPa  k80 %5.1f /s  k63 %5.1f /s\n', cases{i,1}, Fss, fe.ktrP_k80, fe.ktrP_k63);
end
say('  BRACKET TEST: rundown/fresh  protocol k80 %.3f  k63 %.3f  |  slack k63 %.3f\n', ...
    OK{4}.k80/OK{1}.k80, OK{4}.k63/OK{1}.k63, OS{4}.k63/OS{1}.k63);
say('  ATP 2/8:   protocol k80 %.3f  k63 %.3f  |  slack k63 %.3f\n', OK{3}.k80/OK{1}.k80, OK{3}.k63/OK{1}.k63, OS{3}.k63/OS{1}.k63);
f = newFig(); tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(1); plotL(ax1, vtK, t0K, 1e3); ylabel('SL (\mum)'); title('ktr protocol (hi-res burst): commanded length');
ax2 = nexttile([2 1]); hold on;
for i = 1:4, plot((OK{i}.t - t0K)*1e3, OK{i}.Force / OK{i}.Fss, cases{i, 4}, 'Color', cases{i, 3}); end
xlabel('t from final release (ms)'); ylabel('F / F_{iso}'); grid on; legend(cases(:, 1), 'Location', 'southeast');
title('Expected redevelopment (current model)');
linkaxes([ax1 ax2], 'x'); xlim([-25 60]);
print(f, fullfile(res, 'pred_ktr_8_4_2.png'), '-dpng', '-r110'); close(f);

%% ---- 3. FV2 ramp at 8 / 4 / 2 mM -------------------------------------------
Fv = load('data/protocol_04_10_2026_velocitytable_FV2.mat'); vtF = Fv.velocitytable; vtF(1,1) = max(vtF(1,1), -1);
OF = cell(1, 3); PP = nan(1, 3);
for i = 1:3
    p = p0; p.MgATP = cases{i, 2};
    OF{i} = runVT(p, vtF, mf, 0);
    tq = (0:1e-4:0.12)'; Fq = interp1(OF{i}.t, OF{i}.Force, tq, 'linear');
    PP(i) = fv2Read(tq, Fq, vtF); OF{i}.Fss = median(Fq(tq >= 0.002 & tq < 0.040));
    say('FV2   %-20s P/P0 %.3f  (Fiso %.2f kPa)\n', cases{i,1}, PP(i), OF{i}.Fss);
end
say('  FV2 ratios: 4/8 %.3f  2/8 %.3f  (data 2/8 = 0.71 in-prep; Baker borrowed 0.49)\n', PP(2)/PP(1), PP(3)/PP(1));
f = newFig(); tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(1); plotL(ax1, vtF, 0, 1e3); ylabel('SL (\mum)'); title('FV2 ramp (-2 ML/s, 35 ms): commanded length');
ax2 = nexttile([2 1]); hold on;
for i = 1:3, plot(OF{i}.t*1e3, OF{i}.Force / OF{i}.Fss, 'Color', cases{i, 3}); end
xlabel('t (ms)'); ylabel('F / F_{iso}'); grid on;
legend(arrayfun(@(i) sprintf('%s  P/P_0 = %.2f', cases{i,1}, PP(i)), 1:3, 'uni', 0), 'Location', 'northeast');
title('Expected force on the ramp (current model)');
linkaxes([ax1 ax2], 'x'); xlim([-20 160]);
print(f, fullfile(res, 'pred_fv2_8_4_2.png'), '-dpng', '-r110'); close(f);

%% ---- 4. SL offset sensitivity (why L0 criterion / optical SL matters) -------
SLs = 2.10:0.05:2.30; Fiso = nan(size(SLs));
for i = 1:numel(SLs)
    p = p0; p.SL0 = SLs(i); p.Velocity = [0; 0]; p = getParams(p, p.g, true, false);
    [~, o] = evaluateModel(mf, [-3 0], p); Fiso(i) = o.Force(end);
end
say('SLOFF  SL %s\n       F/F(2.2) %s\n', mat2str(SLs), mat2str(Fiso / Fiso(SLs == 2.2), 3));
f = newFig(); tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot([-1 0 1], [2.2 2.2 2.2], 'Color', C8); hold on; plot([-1 0 1], [2.15 2.15 2.15], '--', 'Color', CR); plot([-1 0 1], [2.25 2.25 2.25], '--', 'Color', CR);
ylim([2.05 2.35]); xlabel('t (s)'); ylabel('SL (\mum)'); title('Isometric at nominal SL 2.2 \pm 0.05 (the measured per-prep offset)'); grid on;
nexttile; plot(SLs, Fiso / Fiso(SLs == 2.2), '-o', 'Color', C8); hold on; xline(2.2, ':'); yline(1, ':');
xlabel('true SL (\mum)'); ylabel('F_{iso} / F_{iso}(2.2)'); grid on; title('Expected isometric force vs true SL (current model)');
print(f, fullfile(res, 'pred_sl_offset.png'), '-dpng', '-r110'); close(f);

%% ---- 5. Dose–response: a third [ATP]; power-law vs Hill gate -------------
ATP = [0.5 1 2 4 8 12 20]; Km = 2/(0.54*(8+2)/8 - 1) * (1 - 0.54*(8+2)/8) ; % placeholder, solved below
% Hill gate normalised to 8 mM: g = (A/(A+Km)) / (8/(8+Km)), with g(2) = 0.54  -> Km
Km = fzero(@(k) (2/(2+k))/(8/(8+k)) - 0.54, 3);
gP = (ATP/8).^PK2; gH = (ATP./(ATP+Km)) ./ (8/(8+Km));
FP = nan(size(ATP)); kP = FP; FH = FP; kH = FP;
for i = 1:numel(ATP)
    for form = 1:2
        p = p0; p.MgATP = ATP(i);
        if form == 2 && ATP(i) ~= 8, p.pAtp_k2 = log(gH(i)) / log(ATP(i)/8); end
        o = runVT(p, vtK, mf, 0);
        Fss = interp1(o.t, o.Force, vtK(2,1) - 1e-3); fe = extractKtrAttributes(o.t, o.Force / Fss, vtK);
        if form == 1, FP(i) = Fss; kP(i) = fe.ktrP_k80; else, FH(i) = Fss; kH(i) = fe.ktrP_k80; end
    end
    say('DOSE  %5.1f mM  gate power %.3f Hill %.3f | Fiso %6.2f / %6.2f kPa | k80 %5.1f / %5.1f /s\n', ATP(i), gP(i), gH(i), FP(i), FH(i), kP(i), kH(i));
end
say('  Hill Km = %.2f mM\n', Km);
f = newFig(); tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; semilogx(ATP, FP / FP(ATP == 8), '-o', 'Color', C8); hold on; semilogx(ATP, FH / FH(ATP == 8), '--s', 'Color', CT);
xline(2, ':'); xline(8, ':'); xline(20, ':', 'Color', CT); xline(4, ':', 'Color', C4);
xlabel('[MgATP] (mM)'); ylabel('F_{iso} / F_{iso}(8 mM)'); grid on; legend({'power-law gate (current)', 'Hill gate, same 2 mM value'}, 'Location', 'northeast');
title('Isometric force vs [ATP]: the two gate forms');
nexttile; semilogx(ATP, kP / kP(ATP == 8), '-o', 'Color', C8); hold on; semilogx(ATP, kH / kH(ATP == 8), '--s', 'Color', CT);
xline(2, ':'); xline(8, ':'); xline(20, ':', 'Color', CT); xline(4, ':', 'Color', C4);
xlabel('[MgATP] (mM)'); ylabel('k_{tr} / k_{tr}(8 mM)'); grid on; title('ktr (protocol k80) vs [ATP]');
print(f, fullfile(res, 'pred_dose_forms.png'), '-dpng', '-r110'); close(f);

%% ---- 6. Release–hold–restretch, holds 5–320 ms (SRX refill window) ----------
dL = 0.15; T = [5 10 20 40 80 160 320] * 1e-3; k63T = nan(size(T)); OB = cell(1, 3); showT = [5 20 80] * 1e-3;
for i = 1:numel(T)
    x = T(i);
    vt = [-2, 0, 0, 1.1; 0, -dL/1e-3, 0, 1.1; 1e-3, 0, 0, 1.1-dL; 1e-3 + x, dL/2e-3, 0, 1.1-dL; 3e-3 + x, 0, 0, 1.1; 0.403 + x, 0, 0, 1.1];
    o = runVT(p0, vt, mf, 0);
    tr = 3e-3 + x; m = find(o.t >= tr & o.t <= tr + 15e-3); [Fmin, j] = min(o.Force(m)); tmin = o.t(m(j));
    Fend = interp1(o.t, o.Force, vt(end,1) - 1e-3); lvl = Fmin + 0.63*(Fend - Fmin);
    m2 = o.t >= tmin; tt = o.t(m2); FF = o.Force(m2); ii = find(FF >= lvl, 1);
    t63 = tt(ii-1) + (lvl - FF(ii-1)) * (tt(ii) - tt(ii-1)) / (FF(ii) - FF(ii-1)); k63T(i) = 1/(t63 - tmin);
    say('RHR   hold %4.0f ms  k63 %5.1f /s\n', x*1e3, k63T(i));
    j = find(abs(showT - x) < 1e-6); if ~isempty(j), OB{j} = o; OB{j}.vt = vt; end
end
L10 = [43.0 47.4 59.2 68.5 58.7 51.7 53.1]; S3 = [37.0 38.2 40.1 41.0 43.8 46.7 46.9];   % SRXDuringContraction/results/protocols.txt
f = newFig(); tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = {C8, C4, C2};
ax1 = nexttile(1); hold on; for j = 1:3, plotL(ax1, OB{j}.vt, 0, 1e3, cols{j}); end
ylabel('SL (\mum)'); title('Release 0.15 ML, hold 5 / 20 / 80 ms, restretch in 2 ms'); xlim([-5 120]);
ax2 = nexttile(3); hold on; for j = 1:3, plot(OB{j}.t*1e3, OB{j}.Force, 'Color', cols{j}); end
xlabel('t from release (ms)'); ylabel('force (kPa)'); grid on; xlim([-5 120]); legend({'5 ms', '20 ms', '80 ms'}, 'Location', 'southeast');
title('Expected force (current model)');
nexttile([2 1]); semilogx(T*1e3, k63T, '-o', 'Color', C8); hold on;
semilogx(T*1e3, L10, '--s', 'Color', CT); semilogx(T*1e3, S3, '--^', 'Color', C2);
xlabel('hold duration (ms)'); ylabel('k_{63} after restretch (s^{-1})'); grid on;
legend({'incumbent (SRX_{iso} 0.03)', 'SRX pool x10 (0.14)', 'slow-exchange pool x3 (0.22)'}, 'Location', 'southeast');
title('Read-out: redevelopment rate vs hold');
print(f, fullfile(res, 'pred_release_hold.png'), '-dpng', '-r110'); close(f);

%% ---- 7. Symmetric sawtooth (+4 % L0 at v, back at -v) ----------------------
vs = [0.1 0.4]; OW = cell(2, 2);
for i = 1:2
    v = vs(i); dt = 0.04 / v;
    vt = [-1, 0, 0, 1; 0, v, 0, 1; dt, -v, 0, 1.04; 2*dt, 0, 0, 1; 2*dt + 0.4, 0, 0, 1];
    for onoff = 1:2
        p = p0; p.UseRegistrationAvailability = 2 - onoff; o = runVT(p, vt, mf, 0); o.vt = vt;
        o.F0 = interp1(o.t, o.Force, -1e-3); OW{i, onoff} = o;
        [mx, im] = max(o.Force(o.t > 2*dt)); tt = o.t(o.t > 2*dt);
        say('SAW   v %.2f ML/s  regAvail %d  post-release overshoot %.3f F0 at %.0f ms after return\n', v, 2-onoff, mx/o.F0, (tt(im)-2*dt)*1e3);
    end
end
f = newFig(); tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:2
    ax = nexttile(i); plotL(ax, OW{i,1}.vt, 2*0.04/vs(i), 1e3); ylabel('SL (\mum)'); title(sprintf('Sawtooth +4 %% L_0 at %.1f ML/s', vs(i)));
    nexttile(i+2); hold on; plot((OW{i,1}.t - 2*0.04/vs(i))*1e3, OW{i,1}.Force/OW{i,1}.F0, 'Color', C8);
    plot((OW{i,2}.t - 2*0.04/vs(i))*1e3, OW{i,2}.Force/OW{i,2}.F0, '--', 'Color', CR); yline(1, ':');
    xlabel('t from return to L_0 (ms)'); ylabel('F / F_0'); grid on; legend({'current model', 'registration availability OFF'}, 'Location', 'southwest');
    xlim([-2*0.04/vs(i)*1e3 - 20, 300]);
end
print(f, fullfile(res, 'pred_sawtooth.png'), '-dpng', '-r110'); close(f);

%% ---- 8. Force–pCa (probe gate: pCa50 5.33, nH 3.33; ATP shift not computed) --
pCa = 7:-0.02:4.4; Fp = 1 ./ (1 + 10.^(3.33*(pCa - 5.33)));
f = newFig(); tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot([0 1], [2.2 2.2], 'Color', C8); ylim([2.05 2.35]); xlabel('t (s)'); ylabel('SL (\mum)'); grid on;
title('Isometric at SL 2.2; steady force at each pCa, 8 and 2 mM');
nexttile; plot(pCa, Fp, 'Color', C8); hold on; set(gca, 'XDir', 'reverse'); xlabel('pCa'); ylabel('F / F_{max}'); grid on;
xline(5.33, ':'); xline(4.51, ':', 'Color', CR);
annotation(f, 'textarrow', [0.80 0.74], [0.5 0.5], 'String', '2 mM: pCa_{50} predicted to shift LEFT (bound heads hold sites on); size not computed', 'FontSize', 13, 'Color', C2);
title('Ca probe on this core (literature thin-filament constants): pCa_{50} 5.33, n_H 3.33');
print(f, fullfile(res, 'pred_force_pca.png'), '-dpng', '-r110'); close(f);

%% ---- 9. Mavacamten + mant-ATP (numbers from SRXDuringContraction §10) --------
f = newFig(); tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; b = bar([0.27 0.46 0.51]); b.FaceColor = 'flat'; b.CData = [C8; CT; C2]; set(gca, 'XTickLabel', {'control', 'recruitment x3', 'retention /3'});
ylabel('resting SRX fraction (mant-ATP)'); title('mant chase: both drug models raise rest SRX'); grid on;
nexttile; b = bar([44 36 22]); b.FaceColor = 'flat'; b.CData = [C8; CT; C2]; set(gca, 'XTickLabel', {'control', 'recruitment x3', 'retention /3'});
ylabel('slack k_{tr} (s^{-1})'); title('ktr under the same dose separates them'); grid on;
nexttile; b = bar([0 -2 -2]); b.FaceColor = 'flat'; b.CData = [C8; CT; C2]; set(gca, 'XTickLabel', {'control', 'recruitment x3', 'retention /3'});
ylabel('max isometric force change (%)'); ylim([-10 2]); title('Prediction: force barely moves (-1 to -3 %)'); grid on;
print(f, fullfile(res, 'pred_mava_mant.png'), '-dpng', '-r110'); close(f);

fclose(fid); fprintf('DONE -> %s\n', res); type(fullfile(res, 'pred_summary.txt'));

%% ======================================================================
function o = runVT(p, vt, mf, rundown)
    p.SL0 = 2 * vt(1, 4);
    if rundown, p.kSE = p.kSE * 0.65; p.SL0 = p.SL0 - 0.098; end
    p.Velocity = vt(:, 2); p = getParams(p, p.g, true, false);
    [~, out] = evaluateModel(mf, vt(:, 1), p);
    k = makeMonotonous(out.t(:)); t = out.t(:); F = out.Force(:); o.t = t(k); o.Force = F(k);
end
function [Fss, k63, pk] = slackReads(o, vt)
    rel = find(vt(:, 2) < -20);                 % release rows
    Fss = interp1(o.t, o.Force, vt(rel(1), 1) - 1e-3);
    r = rel(2); tOn = vt(r+1, 1); tEnd = vt(r+2, 1);   % cycle 2: hold = [end of release, start of restretch]
    m = o.t >= tOn & o.t <= tEnd; tt = o.t(m); FF = o.Force(m);
    [Fmin, j] = min(FF(tt < tOn + 0.02)); Fend = FF(end); lvl = Fmin + 0.63*(Fend - Fmin);
    ii = find(FF(:) >= lvl & (1:numel(FF))' > j, 1);
    if isempty(ii), k63 = NaN; else, k63 = 1 / (tt(ii) - tt(j)); end
    m2 = o.t >= tEnd & o.t <= tEnd + 0.05; pk = max(o.Force(m2));
end
function plotL(ax, vt, t0, scale, col)
    if nargin < 5, col = [0.2 0.2 0.2]; end
    L = vt(1, 4) + [0; cumsum(vt(1:end-1, 2) .* diff(vt(:, 1)))];
    plot(ax, (vt(:, 1) - t0) * scale, 2 * L, 'Color', col); grid(ax, 'on');
end
function f = newFig()
    f = figure('Color', 'w', 'Position', [20 20 1600 760], 'Visible', 'off');
end
