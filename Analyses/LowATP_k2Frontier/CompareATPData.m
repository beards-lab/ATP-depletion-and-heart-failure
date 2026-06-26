%% CompareATPData.m
% =====================================================================
% Compares the 8 mM and 2 mM slack feature data (03/27/2026) segment by
% segment, verifies the 5 slack points are paired, and EMITS the low-ATP
% optimization target object used by the lever fits.
%
% Findings that set the target shape:
%   * The 5 slack segments are the SAME protocol in both files (SLslack,
%     SLdiff, t_seg identical; v_restretch within ~1-2%) -> paired, so
%     per-segment ratios are valid.
%   * Force-amplitude ratios are FLAT across segments (CV 1-3%) -> a single
%     scalar per feature is sufficient.
%   * ktr and t0 ratios SLOPE with slack (CV 12-15%) -> per-segment VECTOR.
%   * restretchSlopeStart / ktr2_overshoot are noisy (CV 30-52%) -> de-weight.
%
% Scoring is RELATIVE: model(2mM)/model(8mM) vs these data ratios.
%
% Out: results/atp_target.mat (struct `target`) + results/atp_per_segment.png
% =====================================================================
clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
addpath(genpath(root));
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

f2 = load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;
f8 = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;

% --- consistency check (paired protocol) ---
chk = {'t_seg','SLslack','SLdiff','v_restretch'};
fprintf('Protocol consistency (8mM vs 2mM):\n');
for i=1:numel(chk)
    d = max(abs(double(f8.(chk{i}))-double(f2.(chk{i}))));
    fprintf('  %-12s max|Δ| = %.4g\n', chk{i}, d);
end
SLslack = double(f8.SLslack);

% --- per-segment ratios + constancy ---
forceFeats = {'steady','A','Am','peak2','peak1_y','vall_y'};   % flat -> scalar
kinFeats   = {'ktr','t0'};                                     % sloped -> vector
noisyFeats = {'restretchSlopeStart','ktr2_overshoot'};         % de-weight
allF = [forceFeats kinFeats noisyFeats];

target = struct('scalarFeats',{forceFeats},'vectorFeats',{kinFeats}, ...
                'deweightFeats',{noisyFeats},'SLslack',SLslack);
fprintf('\n%-20s  per-segment ratio (2mM/8mM)        meanR   CV%%   target\n','feature');
for i=1:numel(allF)
    r = double(f2.(allF{i}))./double(f8.(allF{i}));
    cv = 100*std(r)/mean(r);
    if any(strcmp(allF{i},forceFeats)); typ='scalar'; target.scalar.(allF{i})=mean(r);
    elseif any(strcmp(allF{i},kinFeats)); typ='VECTOR'; target.vector.(allF{i})=r;
    else; typ='(noisy)'; target.deweight.(allF{i})=mean(r); end
    fprintf('%-20s  %s  %5.3f  %4.1f  %s\n', allF{i}, num2str(r,'%6.3f '), mean(r), cv, typ);
end

% recommended weights: force scalars high, kinetic vectors high, noisy low
target.weights = struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10, ...
                        'ktr',2,'t0',1,'restretchSlopeStart',0.1,'ktr2_overshoot',0.1);

%% --- figure: per-segment ratio profiles (flat force vs sloped kinetics) ---
fig=figure('Position',[100 100 820 520]); hold on; box on;
xs = 1:5;
for i=1:numel(forceFeats)
    r=double(f2.(forceFeats{i}))./double(f8.(forceFeats{i}));
    h1=plot(xs,r,'-o','Color',[0.10 0.35 0.80],'LineWidth',1.4,'MarkerFaceColor',[0.10 0.35 0.80]);
end
for i=1:numel(kinFeats)
    r=double(f2.(kinFeats{i}))./double(f8.(kinFeats{i}));
    h2=plot(xs,r,'-s','Color',[0.85 0.30 0.10],'LineWidth',1.8,'MarkerFaceColor',[0.85 0.30 0.10]);
    text(5.05,r(end),kinFeats{i},'Color',[0.85 0.30 0.10],'FontSize',9);
end
yline(1,'k:'); yline(1.18,'--','Color',[0.10 0.35 0.80]);
text(1.05,1.20,'force gain \approx 1.18 (flat \rightarrow scalar)','Color',[0.10 0.35 0.80],'FontSize',9);
xlim([0.8 5.8]); xticks(1:5);
xticklabels(arrayfun(@(s) sprintf('%d\\newlineSL%.2f',s,SLslack(s)),1:5,'UniformOutput',false));
xlabel('slack segment'); ylabel('ratio  (2 mM / 8 mM)');
title('Low-ATP effect per slack segment: force-amplitude flat, ktr/t0 slope with slack');
legend([h1 h2],{'force-amplitude features','kinetic features (ktr, t0)'},'Location','east');
exportgraphics(fig, fullfile(resdir,'atp_per_segment.png'),'Resolution',130);

save(fullfile(resdir,'atp_target.mat'),'target');
fprintf('\nSaved results/atp_target.mat (target object) + results/atp_per_segment.png\n');
