% RunAmplitudeVsData.m — the falsification figure.
%
% RunAmplitudeTest showed the MODEL's post-restretch rate rises steeply with
% the size of the re-stretch (68 s^-1 in the small-signal limit -> 137 s^-1 at
% the protocol's own depth). If the real muscle did that too, the model would
% merely be mis-parameterised. So: does it?
%
% The protocol itself sweeps re-stretch amplitude — cycles 1..5 re-stretch by
% 0.080, 0.100, 0.121, 0.141 and 0.061 ML — so the answer is already in the
% data, on every protocol day, for free.
%
% Writes results/amplitude_vs_data.png.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

FILES = { 'protocol_03_27_2026_8mM_slack.mat', '03/27 8 mM'; ...
          'protocol_04_03_2026_8mM_slack.mat', '04/03 8 mM'; ...
          'protocol_04_10_2026_8mM_slack.mat', '04/10 8 mM'; ...
          'protocol_03_27_2026_2mM_slack.mat', '03/27 2 mM'; ...
          'protocol_04_03_2026_2mM_slack.mat', '04/03 2 mM'; ...
          'protocol_04_10_2026_2mM_slack.mat', '04/10 2 mM'};

D = struct('label',{},'amp',{},'rsK',{},'atp',{});
for i = 1:size(FILES,1)
    F = fullfile('data', FILES{i,1});
    if ~isfile(F); fprintf('missing %s\n', F); continue; end
    ds = load(F);
    vt = ds.velocitytable;
    iR = find(vt(:,2) > 1);
    amp = vt(iR+1, 4) - vt(iR, 4);          % ML gained during the re-stretch
    n = min(numel(amp), numel(ds.features_data.rsK));
    D(end+1) = struct('label', FILES{i,2}, 'amp', amp(1:n)', ...
                      'rsK', ds.features_data.rsK(1:n), ...
                      'atp', contains(FILES{i,2},'8 mM')*8 + contains(FILES{i,2},'2 mM')*2); %#ok<SAGROW>
    fprintf('%-12s amp %s\n              rsK %s\n', FILES{i,2}, ...
        mat2str(round(amp(1:n)',3)), mat2str(round(ds.features_data.rsK(1:n),1)));
end

% model, real protocol (from the baseline run) + the synthetic amplitude sweep
B = load(fullfile(here,'results','baseline.mat'));
A = load(fullfile(here,'results','amplitude_test.mat'));
ds1 = load(fullfile('data','protocol_03_27_2026_8mM_slack.mat'));
iR1 = find(ds1.velocitytable(:,2) > 1);
amp_m = (ds1.velocitytable(iR1+1,4) - ds1.velocitytable(iR1,4))';
rsK_m = B.features_model.rsK;

%% ---- slope test: is the trend real? ------------------------------------
fprintf('\n%-14s %8s %10s %10s\n','series','n','slope','r');
allA = []; allK = [];
for i = 1:numel(D)
    if D(i).atp ~= 8; continue; end
    a = D(i).amp(:); k = D(i).rsK(:);
    p = polyfit(a, k, 1); r = corr(a, k);
    fprintf('%-14s %8d %10.0f %10.2f\n', D(i).label, numel(a), p(1), r);
    allA = [allA; a]; allK = [allK; k]; %#ok<AGROW>
end
pD = polyfit(allA, allK, 1); rD = corr(allA, allK);
fprintf('%-14s %8d %10.0f %10.2f\n', 'DATA pooled 8mM', numel(allA), pD(1), rD);

pM = polyfit(amp_m(:), rsK_m(:), 1); rM = corr(amp_m(:), rsK_m(:));
fprintf('%-14s %8d %10.0f %10.2f\n', 'MODEL protocol', numel(amp_m), pM(1), rM);

useS = [A.res.okRestretch];
aS = [A.res(useS).depth]'; kS = [A.res(useS).kPostRestretch]';
pS = polyfit(aS, kS, 1); rS = corr(aS, kS);
fprintf('%-14s %8d %10.0f %10.2f\n', 'MODEL synthetic', numel(aS), pS(1), rS);

%% ---- figure ------------------------------------------------------------
fig = figure(703); clf; set(fig,'Position',[70 70 1080 440],'Color','w');

subplot(1,2,1); hold on; box on;
co = lines(6);
for i = 1:numel(D)
    if D(i).atp ~= 8; continue; end
    plot(D(i).amp, D(i).rsK, 'o', 'Color', co(i,:), 'MarkerFaceColor', co(i,:), ...
         'MarkerSize', 7, 'DisplayName', ['data ' D(i).label]);
end
plot(amp_m, rsK_m, 'rs-', 'LineWidth', 1.8, 'MarkerFaceColor','r', ...
     'MarkerSize', 8, 'DisplayName', 'MODEL (same protocol)');
plot(aS, kS, 'r^--', 'LineWidth', 1.4, 'DisplayName', 'MODEL (synthetic sweep)');
xl = xlim; plot(xl, polyval(pD, xl), 'k-', 'LineWidth', 1.2, 'DisplayName', 'data trend');
xlabel('re-stretch amplitude (ML)'); ylabel('post-restretch rate  rsK (s^{-1})');
title('8 mM ATP: rate vs re-stretch size');
legend('Location','northwest','FontSize',7.5); ylim([0 150]);

subplot(1,2,2); hold on; box on;
lbl = {}; vals = [];
for i = 1:numel(D)
    vals(end+1) = mean(D(i).rsK, 'omitnan'); lbl{end+1} = D(i).label; %#ok<SAGROW>
end
vals(end+1) = mean(rsK_m, 'omitnan'); lbl{end+1} = 'MODEL 8 mM';
b = bar(vals, 'FaceColor', 'flat');
for i = 1:numel(D); b.CData(i,:) = [0.2 0.4 0.8]*(1 - 0.4*(D(i).atp==2)) + 0.4*(D(i).atp==2); end
b.CData(end,:) = [0.85 0.2 0.2];
set(gca,'XTick',1:numel(lbl),'XTickLabel',lbl,'XTickLabelRotation',35,'FontSize',8);
ylabel('mean rsK (s^{-1})'); title('all preps + model');
yline(mean(allK), 'k--', '8 mM data mean');

exportgraphics(fig, fullfile(here,'results','amplitude_vs_data.png'), 'Resolution',140);
fprintf('\nwrote results/amplitude_vs_data.png\n');
