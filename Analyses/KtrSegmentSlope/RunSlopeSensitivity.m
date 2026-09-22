% RunSlopeSensitivity.m -- why does the model's slack ktr not FALL across the
% five slack segments, the way the data does, and is anything in the model
% sensitive to that slope at all?
%
% WHAT THE SEGMENTS ACTUALLY ARE.  The slack protocol is not five repeats: it
% is a LADDER of increasing slack depth, all five with the same 2.2 ms release
% duration, so both the shortening velocity and the held length change:
%
%     seg      SL_hold (um)      v_release (um/s)
%      1          2.041              -36
%      2          2.001              -45
%      3          1.960              -53
%      4          1.920              -62
%      5          1.879              -70
%
% and `ktr` is fitted to force redevelopment during the 297 ms isometric hold
% AT THAT SHORTENED LENGTH.  So the feature is ktr(SL), sampled at five lengths
% between 2.04 and 1.88 um.  The data falls ~14 % over that span; the model is
% flat or rises.
%
% WHAT THIS SCRIPT DOES
%   1. the data slope, pooled and per preparation, so we know it is not noise
%   2. the baseline model slope, plus the two geometric SL factors that
%      multiply the attachment flux (N_overlap and f_lattice), evaluated
%      analytically at each segment's held length
%   3. structural variants -- one switch or one parameter at a time -- scored
%      on BOTH the ktr slope and the full objective
%   4. a one-at-a-time sweep of the whole 77-parameter lean pool at +15 %,
%      ranked by |d(slope)|, as a net for anything not thought of above
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunSlopeSensitivity.m')"
clc;
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrSegmentSlope','results');
if ~isfolder(resdir); mkdir(resdir); end

cl = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_kss_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(cl, 6); end

BASE = 'refit4_opt';
p0 = getParams(loadParams(BASE), [], true, false);
p0.RunSlackSegments = 'All';        % serial inside -- the outer parfor owns the pool
p0.EvalFeatures = 1; p0.PlotEachSeparately = 0; p0.PlotFeatureFitting = 0;
p0.ghostSave = ''; p0.ghostLoad = ''; p0.SaveBest = 0;
p0.BreakOnODEUnstable = 0; p0.MaxRunTime = 600;
p0.mods = {}; p0.g = ones(1,0);

%% ---- 1. the data ---------------------------------------------------------
DAYS = {'03_27','04_03','04_10'};
PT = pooledTarget8mM();
jk = find(strcmp(PT.name,'ktr'),1);
ktr_pool = PT.val(jk,:);
fprintf('=== DATA ===\n');
fprintf('%-10s %7s %7s %7s %7s %7s   %10s %8s\n','source','seg1','seg2','seg3','seg4','seg5','slope/seg','%%span');
K = nan(3,5);
for d = 1:3
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    if isfield(S,'features_data') && isfield(S.features_data,'ktr')
        K(d,:) = S.features_data.ktr(:)';
    end
    b = polyfit(1:5, K(d,:), 1);
    fprintf('%-10s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f %7.1f%%\n', strrep(DAYS{d},'_','/'), ...
            K(d,:), b(1), 100*(K(d,5)-K(d,1))/K(d,1));
end
b_pool = polyfit(1:5, ktr_pool, 1);
SLOPE_DATA = b_pool(1);
fprintf('%-10s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f %7.1f%%\n','POOLED', ktr_pool, ...
        SLOPE_DATA, 100*(ktr_pool(5)-ktr_pool(1))/ktr_pool(1));
fprintf('  pooled between-prep s for ktr = %.4f (relative)\n', PT.s(jk));

%% ---- 2. the geometry, analytically ---------------------------------------
SL_hold = [2.0406 2.0006 1.9602 1.9202 1.8792];
No = zeros(1,5); fl = zeros(1,5);
for i = 1:5
    SL = SL_hold(i);
    L1 = min(p0.L_thick*0.5, SL*0.5);
    L2 = max(SL*0.5 - (SL - p0.L_thin), p0.L_hbare*0.5);
    No(i) = (L1-L2)*2/(p0.L_thick - p0.L_hbare);
    d10 = p0.d10_ref * sqrt(p0.SL_ref_lattice / SL);
    ex  = max(0, d10 - p0.R_thick - p0.R_thin - p0.d_optimal);
    fl(i) = exp(-(ex^2)/(2*p0.sigma_lattice^2));
end
fprintf('\n=== the SL factors that multiply the attachment flux RD1 ===\n');
fprintf('%-14s %8s %8s %8s %8s %8s   %8s\n','factor','seg1','seg2','seg3','seg4','seg5','seg5/seg1');
fprintf('%-14s %8.4f %8.4f %8.4f %8.4f %8.4f   %8.3f\n','N_overlap',No,No(5)/No(1));
fprintf('%-14s %8.4f %8.4f %8.4f %8.4f %8.4f   %8.3f\n','f_lattice',fl,fl(5)/fl(1));
fprintf('%-14s %8.4f %8.4f %8.4f %8.4f %8.4f   %8.3f\n','product',No.*fl,(No(5)*fl(5))/(No(1)*fl(1)));
fprintf(['  If ktr tracked the attachment prefactor, it would fall %.0f %% across the ladder.\n' ...
         '  The data falls %.0f %%. UseOverlapFactor = %d.\n'], ...
        100*(1-(No(5)*fl(5))/(No(1)*fl(1))), 100*(1-ktr_pool(5)/ktr_pool(1)), p0.UseOverlapFactor);

%% ---- 3. variants ---------------------------------------------------------
% {label, field, value}  -- value [] means "toggle the logical"
V = {
  'BASELINE',                  '',                  []
  'UseOverlapFactor=0',        'UseOverlapFactor',  0
  'UsePassiveForSR=1',         'UsePassiveForSR',   1
  'useHalfActiveForSR=1',      'useHalfActiveForSR',1
  'sigma_lattice x0.7',        'sigma_lattice',     p0.sigma_lattice*0.7
  'sigma_lattice x1.5',        'sigma_lattice',     p0.sigma_lattice*1.5
  'd_optimal x0.97',           'd_optimal',         p0.d_optimal*0.97
  'd10_ref x1.03',             'd10_ref',           p0.d10_ref*1.03
  'P_bound_max x0.25',         'P_bound_max',       p0.P_bound_max*0.25
  'L_thin x1.03',              'L_thin',            p0.L_thin*1.03
  'L_thick x0.97',             'L_thick',           p0.L_thick*0.97
  'tau_reg = 20 ms',           'tau_reg',           0.020
  'RegAvailShorteningOnly=0',  'RegAvailShorteningOnly', 0
  'sigma1 x0.5',               'sigma1',            p0.sigma1*0.5
  'sigma1 x2',                 'sigma1',            p0.sigma1*2
  'ksr2srd x0.5',              'ksr2srd',           p0.ksr2srd*0.5
  'ksrd2sr x0.5',              'ksrd2sr',           p0.ksrd2sr*0.5
  'kmsrd x0.5',                'kmsrd',             p0.kmsrd*0.5
  'ksr0 x0.5',                 'ksr0',              p0.ksr0*0.5
  'Lsc0 x1.03',                'Lsc0',              p0.Lsc0*1.03
};
nV = size(V,1);
Vk = nan(nV,5); Vs = nan(nV,1); Verr = cell(nV,1);
parfor i = 1:nV
    p = p0;
    if ~isempty(V{i,2}); p.(V{i,2}) = V{i,3}; p = getParams(p, p.g, true, false); end
    try
        [~,~,fm] = runSlackExperiment(p);
        k = fm.ktr(:)';
        Vk(i,:) = k; bb = polyfit(1:5,k,1); Vs(i) = bb(1);
    catch ME
        Verr{i} = ME.message;
    end
end
fprintf('\n=== VARIANTS: per-segment ktr and the slope ===\n');
fprintf('%-26s %7s %7s %7s %7s %7s   %10s %10s\n','variant','seg1','seg2','seg3','seg4','seg5','slope','d(slope)');
s_base = Vs(1);
for i = 1:nV
    if ~isempty(Verr{i}); fprintf('%-26s  FAILED: %s\n', V{i,1}, Verr{i}); continue; end
    fprintf('%-26s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f %+10.3f\n', ...
            V{i,1}, Vk(i,:), Vs(i), Vs(i)-s_base);
end
fprintf('%-26s %7.2f %7.2f %7.2f %7.2f %7.2f   %+10.3f   <- TARGET\n','POOLED DATA',ktr_pool,SLOPE_DATA);

%% ---- 4. whole lean pool, one at a time ----------------------------------
pool = p0.leanPool;
np = numel(pool);
Ps = nan(np,1); Pk = nan(np,5);
parfor i = 1:np
    p = p0; nmp = pool{i}; us = strfind(nmp,'__');
    try
        if isempty(us)
            if ~isfield(p,nmp) || ~isscalar(p.(nmp)) || p.(nmp)==0; continue; end
            p.(nmp) = p.(nmp)*1.15;
        else
            b = nmp(1:us(1)-1); ix = str2double(nmp(us(1)+2:end));
            if ~isfield(p,b) || numel(p.(b))<ix || p.(b)(ix)==0; continue; end
            v = p.(b); v(ix) = v(ix)*1.15; p.(b) = v;
        end
        p = getParams(p, p.g, true, false);
        [~,~,fm] = runSlackExperiment(p);
        k = fm.ktr(:)'; Pk(i,:) = k; bb = polyfit(1:5,k,1); Ps(i) = bb(1);
    catch
    end
end
d = Ps - s_base;
[~,ord] = sort(abs(d),'descend','MissingPlacement','last');
fprintf('\n=== LEAN POOL, +15 %% one at a time: effect on the ktr slope ===\n');
fprintf('  (baseline slope %+0.3f /segment, data %+0.3f -- a fix needs d(slope) <= %+0.3f)\n', ...
        s_base, SLOPE_DATA, SLOPE_DATA-s_base);
fprintf('%-34s %10s %10s %9s %9s\n','param','slope','d(slope)','ktr(1)','ktr(5)');
for i = ord(1:min(25,numel(ord)))'
    if ~isfinite(d(i)); continue; end
    fprintf('%-34s %+10.3f %+10.3f %9.2f %9.2f\n', pool{i}, Ps(i), d(i), Pk(i,1), Pk(i,5));
end
nneg = sum(d < (SLOPE_DATA - s_base) & isfinite(d));
fprintf('\n  parameters that on their own reach the data slope: %d of %d\n', nneg, sum(isfinite(d)));

save(fullfile(resdir,'slope_sensitivity.mat'), 'V','Vk','Vs','pool','Ps','Pk','d', ...
     'ktr_pool','K','SLOPE_DATA','s_base','No','fl','SL_hold','BASE');

%% ---- plot ---------------------------------------------------------------
f = figure('Position',[30 30 1500 780],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on
for dd = 1:3; plot(1:5,K(dd,:),'-','Color',[0.75 0.75 0.75],'LineWidth',0.8); end
plot(1:5,ktr_pool,'ko-','LineWidth',2,'MarkerFaceColor','k');
plot(1:5,Vk(1,:),'o-','Color',[0 0.45 0.74],'LineWidth',1.8);
xlabel('slack segment'); ylabel('ktr (s^{-1})'); xticks(1:5);
title(sprintf('A  data %+0.2f/seg  vs  model %+0.2f/seg', SLOPE_DATA, s_base));
legend({'','','','pooled data','refit4\_opt'},'Location','best','Box','off');

nexttile; hold on; box on; grid on
yyaxis left;  plot(SL_hold,No,'o-','LineWidth',1.6); ylabel('N\_overlap');
yyaxis right; plot(SL_hold,fl,'s-','LineWidth',1.6); ylabel('f\_lattice');
set(gca,'XDir','reverse'); xlabel('SL during the hold (\mum)');
title('B  the two SL factors on the attachment flux');

nexttile; hold on; box on; grid on
plot(1:5,(No.*fl)/(No(1)*fl(1)),'k^-','LineWidth',1.8);
plot(1:5,ktr_pool/ktr_pool(1),'ko-','LineWidth',1.8,'MarkerFaceColor','k');
plot(1:5,Vk(1,:)/Vk(1,1),'o-','Color',[0 0.45 0.74],'LineWidth',1.8);
xlabel('slack segment'); ylabel('normalised to segment 1'); xticks(1:5);
title('C  the prefactor falls 40 %, the model does not');
legend({'N\_overlap \times f\_lattice','data ktr','model ktr'},'Location','best','Box','off');

nexttile([1 2]); hold on; box on; grid on
dv = Vs(2:end) - s_base; lab = V(2:end,1);
[~,o2] = sort(dv,'ascend','MissingPlacement','last');
for i = 1:numel(o2)
    if ~isfinite(dv(o2(i))); continue; end
    c = [0.8 0.3 0.2]; if dv(o2(i)) < 0; c = [0.2 0.6 0.3]; end
    barh(i, dv(o2(i)), 0.8, 'FaceColor', c, 'EdgeColor','none');
end
xline(SLOPE_DATA - s_base,'k--','LineWidth',1.5,'Label','data slope');
xline(0,'k-');
set(gca,'YTick',1:numel(o2),'YTickLabel',lab(o2),'YDir','reverse','FontSize',9, ...
        'TickLabelInterpreter','none'); grid on
xlabel('\Delta (ktr slope) per segment'); title('D  structural variants');

nexttile; hold on; box on; grid on
[~,o3] = sort(d,'ascend','MissingPlacement','last'); o3 = o3(1:min(12,numel(o3)));
for i = 1:numel(o3)
    if ~isfinite(d(o3(i))); continue; end
    c = [0.8 0.3 0.2]; if d(o3(i))<0; c=[0.2 0.6 0.3]; end
    barh(i, d(o3(i)), 0.8,'FaceColor',c,'EdgeColor','none');
end
xline(SLOPE_DATA - s_base,'k--','LineWidth',1.5);
set(gca,'YTick',1:numel(o3),'YTickLabel',pool(o3),'YDir','reverse','FontSize',8, ...
        'TickLabelInterpreter','none');
xlabel('\Delta slope at +15 %'); title('E  lean pool, most negative 12');

sgtitle(sprintf('slack ktr across the slack ladder (SL 2.04 -> 1.88 \\mum): %s', BASE), ...
        'FontSize',12,'Interpreter','tex');
exportgraphics(f, fullfile(resdir,'ktr_slope.png'), 'Resolution',110);
fprintf('\nwrote %s\n', fullfile(resdir,'ktr_slope.png'));
fprintf('\nDONE\n');
