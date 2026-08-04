% RunMechanismSimulation.m
% Which mechanical lesion reproduces the OBSERVED rundown signature?
%
% The data (RunRundownMechanisms.m) say that 694 s of rundown does three
% measurable things to the 03/27 8 mM fibre:
%     (i)   force            x0.829
%     (ii)  ktr              x0.884
%     (iii) the length-tension curve BENDS -- a horizontal shift fits it better
%           than a vertical scale at equal degrees of freedom
%
% Here each candidate lesion is imposed ON THE MODEL, driven through the same
% slack protocol, and read out through the same three observables. Because the
% mechanisms are not equally "strong" per unit of their own parameter, each is
% swept and then interpolated to the level that reproduces the OBSERVED FORCE
% LOSS exactly. With force matched by construction, ktr and the length-tension
% shape become clean, independent discriminators.
%
%   M1  fewer attached heads (myosin loss / density loss)   kstiff1, kstiff2 x s
%   M2  creep in the SERIAL elastic (tearing of the ends)   kSE x c
%   M3  lost sarcomere length at fixed ML (slack)           SL0 - 2d
%   M4  uniform kinetic slowdown (motor damage)             xrate x x
%   M5  reduced attachment / "less power" per sarcomere     ka x a
%
% M6, a dead region in PARALLEL (torn myofibrils still mechanically attached),
% is not simulated: on the active side it is indistinguishable from M1, and it
% differs only by adding parallel passive force, which is ~5 kPa here and cannot
% produce the observed bend (see conclusions.md).
%
% Outputs -> results/ : mechanism_simulation.png / .mat

clear; close all; clc;
cd(fileparts(which('RunMechanismSimulation')));
addpath(genpath('../..'));
resDir = 'results';  if ~isfolder(resDir), mkdir(resDir); end

OBS_F   = 0.829;    % observed force ratio  (RunRundownMechanisms.m)
OBS_KTR = 0.884;    % observed ktr   ratio

%% ── Baseline ───────────────────────────────────────────────────────────────
p0 = getParams(loadParams('ThursdayNightFever'), [], true, false);
p0.PlotEachSeparately = 0;  p0.PlotFullscreen = 0;  p0.PlotFeatureFitting = 0;
p0.ghostSave = '';  p0.ghostLoad = '';  p0.EvalFeatures = 1;
p0.RunSlackSegments = 'All';
p0.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';

S  = load('../../data/protocol_03_27_2026_8mM_slack.mat');
vt = S.velocitytable;  is = find(vt(:,2) < -1);
holdEnd = vt(is+2,1);  MLpts = [vt(is+1,4); vt(1,4)]';   % 5 slack holds + ML 1.10
tPts    = [holdEnd; 74.20]';

fprintf('=== baseline ===\n');
tic; [F0, k0] = iRun(p0, tPts); fprintf('  %.1f s | force %s | ktr %.2f\n', toc, mat2str(round(F0,1)), k0);

%% ── Mechanism sweep ────────────────────────────────────────────────────────
MECH = struct( ...
 'id',    {'M1','M2','M3','M4','M5'}, ...
 'name',  {'fewer heads (kstiff)','serial creep (kSE)','lost SL (SL0)', ...
           'uniform slowdown (xrate)','less attachment (ka)'}, ...
 'lev',   {[1 0.90 0.80],  [1 0.40 0.20],  [0 0.05 0.10], ...
           [1 0.85 0.70],  [1 0.80 0.60]}, ...
 'apply', {@(p,v) setfield(setfield(p,'kstiff1',p.kstiff1*v),'kstiff2',p.kstiff2*v), ...
           @(p,v) setfield(p,'kSE',p.kSE*v), ...
           @(p,v) setfield(p,'SL0',p.SL0-2*v), ...
           @(p,v) setfield(p,'xrate',v), ...
           @(p,v) setfield(p,'ka',p.ka*v)});

for m = 1:numel(MECH)
    L = MECH(m).lev;  nF = nan(numel(L),numel(tPts));  nK = nan(1,numel(L));
    fprintf('\n=== %s : %s ===\n', MECH(m).id, MECH(m).name);
    for j = 1:numel(L)
        p = MECH(m).apply(p0, L(j));
        try
            tic; [Fv, kk] = iRun(p, tPts);
            nF(j,:) = Fv;  nK(j) = kk;
            fprintf('  level %6.3f | %5.1f s | force x%.3f | ktr x%.3f\n', ...
                    L(j), toc, mean(Fv(1:5))/mean(F0(1:5)), kk/k0);
        catch ME
            fprintf('  level %6.3f | FAILED: %s\n', L(j), ME.message);
        end
    end
    MECH(m).F = nF;  MECH(m).K = nK;
    MECH(m).fr = mean(nF(:,1:5),2)'/mean(F0(1:5));   % force ratio per level
    MECH(m).kr = nK/k0;                              % ktr ratio per level
end

%% ── Interpolate each mechanism to the OBSERVED force loss ──────────────────
fprintf('\n\n============ MATCHED-FORCE COMPARISON (all set to force x%.3f) ============\n', OBS_F);
fprintf('%-4s %-28s %10s %10s %10s %12s\n','id','mechanism','level','force x','ktr x','ktr error');
for m = 1:numel(MECH)
    fr = MECH(m).fr; kr = MECH(m).kr; L = MECH(m).lev;
    good = isfinite(fr);
    if sum(good) < 2 || min(fr(good)) > OBS_F
        MECH(m).lvlStar = NaN; MECH(m).krStar = NaN; MECH(m).FStar = nan(1,numel(tPts));
        fprintf('%-4s %-28s %10s %10s %10s %12s\n', MECH(m).id, MECH(m).name, ...
                '-','out of range','-','-');
        continue;
    end
    MECH(m).lvlStar = interp1(fr(good), L(good), OBS_F, 'linear', 'extrap');
    MECH(m).krStar  = interp1(fr(good), kr(good), OBS_F, 'linear', 'extrap');
    MECH(m).FStar   = interp1(fr(good), MECH(m).F(good,:), OBS_F, 'linear', 'extrap');
    fprintf('%-4s %-28s %10.3f %10.3f %10.3f %+11.1f%%\n', MECH(m).id, MECH(m).name, ...
            MECH(m).lvlStar, OBS_F, MECH(m).krStar, 100*(MECH(m).krStar-OBS_KTR)/OBS_KTR);
end
fprintf('%-4s %-28s %10s %10.3f %10.3f %12s\n','--','OBSERVED (03/27 rundown)','',OBS_F,OBS_KTR,'--');

%% ── Length-tension shape: does the mechanism bend the curve? ───────────────
% Same test the data got: fit the perturbed curve as a pure vertical SCALE of
% the baseline, and as a pure horizontal SHIFT, and see which wins.
pq = polyfit(MLpts, F0, 2);  Fmod = @(x) polyval(pq, x);
dd = linspace(0, 0.20, 2001);
fprintf('\n\n============ LENGTH-TENSION SHAPE (at matched force) ============\n');
fprintf('%-4s %-28s %11s %11s %14s\n','id','mechanism','RMSE scale','RMSE shift','shape verdict');
for m = 1:numel(MECH)
    Fs = MECH(m).FStar;
    if any(~isfinite(Fs)); MECH(m).rScale=NaN; MECH(m).rShift=NaN; continue; end
    s  = sum(F0.*Fs)/sum(F0.^2);        rs = sqrt(mean((Fs - s*F0).^2));
    [~,i2] = min(arrayfun(@(d) sum((Fs - Fmod(MLpts-d)).^2), dd));
    rh = sqrt(mean((Fs - Fmod(MLpts-dd(i2))).^2));
    MECH(m).rScale = rs;  MECH(m).rShift = rh;
    if rh < rs*0.8; v = 'SHIFT (bends)'; elseif rs < rh*0.8; v = 'scale'; else; v = 'ambiguous'; end
    fprintf('%-4s %-28s %11.3f %11.3f %14s\n', MECH(m).id, MECH(m).name, rs, rh, v);
end
fprintf('%-4s %-28s %11.3f %11.3f %14s\n','--','OBSERVED (03/27 rundown)',2.110,1.238,'SHIFT (bends)');

%% ── Combinations: no single lesion hits all three observables ──────────────
% The sweep splits the observables cleanly:
%   force + the L-T BEND  <- only M3 (lost SL)
%   ktr                   <- only M2 (serial creep) or M4 (uniform slowdown)
% M2 and M3 are the SAME physical lesion: a series elastic element that creeps
% longer and softer simultaneously lets the sarcomeres sit shorter at a given ML
% (M3) and adds compliance (M2). So test that pairing against the alternative
% pairing (fewer heads + slower motors), which reaches the same (force, ktr)
% point but by a route that should NOT bend the length-tension curve.
COMB = struct( ...
 'id',   {'C1','C2'}, ...
 'name', {'series-elastic damage: kSE x0.65, SL -0.098um','fewer heads + slower motors (kstiff + xrate)'}, ...
 'apply',{@(p) setfield(setfield(p,'kSE',p.kSE*0.65),'SL0',p.SL0-2*0.049), ...
          @(p) setfield(setfield(setfield(p,'kstiff1',p.kstiff1*0.84),'kstiff2',p.kstiff2*0.84),'xrate',0.85)});

fprintf('\n\n============ COMBINATIONS ============\n');
fprintf('%-4s %-46s %9s %9s %11s %11s %s\n','id','combination','force x','ktr x','RMSE scale','RMSE shift','shape');
for c = 1:numel(COMB)
    [Fv, kk] = iRun(COMB(c).apply(p0), tPts);
    COMB(c).F = Fv;  COMB(c).fr = mean(Fv(1:5))/mean(F0(1:5));  COMB(c).kr = kk/k0;
    s  = sum(F0.*Fv)/sum(F0.^2);   rs = sqrt(mean((Fv - s*F0).^2));
    [~,i2] = min(arrayfun(@(d) sum((Fv - Fmod(MLpts-d)).^2), dd));
    rh = sqrt(mean((Fv - Fmod(MLpts-dd(i2))).^2));
    COMB(c).rScale = rs; COMB(c).rShift = rh;
    if rh < rs*0.8; v='SHIFT (bends)'; elseif rs < rh*0.8; v='scale'; else; v='ambiguous'; end
    fprintf('%-4s %-46s %9.3f %9.3f %11.3f %11.3f %s\n', COMB(c).id, COMB(c).name, ...
            COMB(c).fr, COMB(c).kr, rs, rh, v);
end
fprintf('%-4s %-46s %9.3f %9.3f %11.3f %11.3f %s\n','--','OBSERVED (03/27 rundown)', ...
        OBS_F, OBS_KTR, 2.110, 1.238, 'SHIFT (bends)');
fprintf('\n  Both combinations can reach the observed (force, ktr) point -- that pair of\n');
fprintf('  numbers alone does NOT identify the lesion. The length-tension SHAPE does.\n');

%% ── Figure ─────────────────────────────────────────────────────────────────
fig = figure(940); clf; set(fig,'Position',[20 20 1500 760]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cm = lines(numel(MECH));

nexttile([1 2]); hold on; box on;
for m = 1:numel(MECH)
    plot(MECH(m).fr, MECH(m).kr, 'o-', 'Color', cm(m,:), 'LineWidth', 1.8, ...
         'MarkerFaceColor', cm(m,:), 'DisplayName', [MECH(m).id ' ' MECH(m).name]);
end
plot(OBS_F, OBS_KTR, 'kp', 'MarkerSize', 20, 'MarkerFaceColor','y', 'DisplayName','OBSERVED rundown');
plot(1,1,'k*','MarkerSize',12,'HandleVisibility','off');
for c = 1:numel(COMB)
    plot(COMB(c).fr, COMB(c).kr, 'kd', 'MarkerSize', 11, 'LineWidth', 1.5, ...
         'MarkerFaceColor', [0.6 0.6 0.6]*c/numel(COMB), 'DisplayName', [COMB(c).id ' ' COMB(c).name]);
end
xline(OBS_F,':','HandleVisibility','off'); yline(OBS_KTR,':','HandleVisibility','off');
xlabel('force  (x baseline)'); ylabel('ktr  (x baseline)');
title('Each mechanism traces a different path in (force, ktr)');
legend('Location','southwest','FontSize',7); grid on;

nexttile; hold on; box on;
kk = [MECH.krStar];
b = bar([kk OBS_KTR]); b.FaceColor='flat';
b.CData = [cm; 1 0.85 0];
set(gca,'XTick',1:numel(kk)+1,'XTickLabel',[{MECH.id} {'obs'}]);
yline(OBS_KTR,'k--','LineWidth',1.5); ylabel('ktr x baseline (at matched force)');
title('ktr at the observed force loss');

nexttile; hold on; box on;
plot(MLpts, F0, 'k^-','MarkerFaceColor','k','LineWidth',1.5,'DisplayName','baseline');
for m = 1:numel(MECH)
    if any(~isfinite(MECH(m).FStar)); continue; end
    plot(MLpts, MECH(m).FStar, 'o--','Color',cm(m,:),'DisplayName',MECH(m).id);
end
xlabel('ML'); ylabel('force (kPa)'); title('Length-tension at matched force loss');
legend('Location','northwest','FontSize',7);

nexttile; hold on; box on;
rr = [[MECH.rScale] [COMB.rScale] 2.110; [MECH.rShift] [COMB.rShift] 1.238]';
b = bar(rr,'grouped'); b(1).FaceColor=[.4 .4 .8]; b(2).FaceColor=[.3 .7 .3];
set(gca,'XTick',1:size(rr,1),'XTickLabel',[{MECH.id} {COMB.id} {'OBS'}]);
ylabel('RMSE (kPa)'); title('Does the mechanism BEND the L-T curve?');
legend({'pure scale','pure shift'},'Location','northwest','FontSize',7);
text(numel(MECH)+1, max(rr(:))*0.9, 'shift < scale = bends', 'FontSize',7);

nexttile; hold on; box on;
score = abs([MECH.krStar]-OBS_KTR)/OBS_KTR*100;
shapeOK = double([MECH.rShift] < [MECH.rScale]*0.8);
b = bar([score(:) 20*shapeOK(:)],'grouped');
b(1).FaceColor=[.85 .33 .1]; b(2).FaceColor=[.2 .6 .3];
set(gca,'XTick',1:numel(MECH),'XTickLabel',{MECH.id});
ylabel('ktr error (%)   /   shape match (20 = yes)');
title('Verdict per mechanism'); legend({'|ktr error|','L-T shape matches'},'FontSize',7);

sgtitle('Which lesion reproduces rundown? Model perturbations at matched force loss');
exportgraphics(fig, fullfile(resDir,'mechanism_simulation.png'), 'Resolution', 150);
save(fullfile(resDir,'mechanism_simulation.mat'), 'MECH','F0','k0','MLpts','OBS_F','OBS_KTR');
fprintf('\nSaved %s\n', fullfile(resDir,'mechanism_simulation.png'));

% ---------------------------------------------------------------------------
function [Fv, ktr] = iRun(p, tPts)
%IRUN Drive the slack protocol and read force at the hold ends plus mean ktr.
    [~, out, fm] = runSlackExperiment(p);
    mono = makeMonotonous(out.t);
    Fv   = interp1(out.t(mono), out.Force(mono), tPts, 'linear', 'extrap');
    ktr  = mean(fm.ktr, 'omitnan');
end
