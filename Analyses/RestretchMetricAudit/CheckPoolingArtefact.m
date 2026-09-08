% CheckPoolingArtefact.m -- was the -24 % "trace averaging destroys rsK"
% result real, or a splice artefact?
%
% RunMetricAudit's T5 averages the post-re-stretch window directly and finds a
% gap of only -2.0 % for the CURRENT estimator -- not -24 %. The two runs
% differ in ONE thing: DataPooling/RunPoolingComparison builds the grand
% average on a relative grid TG = -0.020 : 0.570 s and then SPLICES it into
% 03/27's own datatable, leaving the rest of 03/27's (merely rescaled) trace in
% place. The slack segment is ~0.58 s long after the release and the
% post-re-stretch window alone is 0.280 s, so the LAST ~10 ms of every recovery
% window is not the average at all -- it is one preparation, rescaled.
%
% That matters because the recovery amplitude A is the median of the LAST 15 %
% of the window (42 ms). If the splice leaves a step there, A is wrong and the
% rate inherits it.
%
% This script measures the step, reproduces the -24 %, and re-runs the same
% comparison with the grid extended past the end of the segment.
%
% Output -> results/pooling_artefact.mat / .png

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root)); addpath(here);
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

days = {'03/27','04/03','04/10'}; nP = 3; nSeg = 5;
R = poolLoadRuns();
Fref = nan(nP,1);
for p = 1:nP
    Fref(p) = R(find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'),1)).A.Fss22;
end

Sm  = load(fullfile('data', slackFile(days{1},'8mM')));
t0v = Sm.datatable(:,1); L0 = Sm.datatable(:,2); vt0 = Sm.velocitytable;
rel0 = find(vt0(:,2) < -20);
seg0 = find(vt0(:,2) < -1)';

% how long IS a segment, measured from the release?
segLen = nan(nSeg,1);
for g = 1:nSeg
    vs = vt0(seg0(g):seg0(g)+4,1);
    segLen(g) = vs(5) - vt0(rel0(g),1);
end
fprintf('=== segment length after the release ===\n');
fprintf('  %.3f s (min) .. %.3f s (max);  DataPooling grid ends at 0.570 s\n', ...
        min(segLen), max(segLen));
fprintf('  -> the last %.0f ms of the longest segment is NOT averaged\n', ...
        1000*(max(segLen)-0.570));

%% ---- build the grand average on two grids ----------------------------
GRIDS = {0.570, 'DataPooling grid (-0.020 .. 0.570)'
         max(segLen)+0.02, 'extended past the segment end'};
out = struct('lbl',{},'rsK',{},'rsA',{},'rsT0',{},'jump',{},'Fgrand',{},'tt',{});
for q = 1:size(GRIDS,1)
    TG = (-0.020 : 2e-4 : GRIDS{q,1})';
    Tr = nan(nP,nSeg,numel(TG));
    for p = 1:nP
        S = load(fullfile('data', slackFile(days{p},'8mM')));
        t = S.datatable(:,1); F = S.datatable(:,3); vt = S.velocitytable;
        rel = find(vt(:,2) < -20);
        for g = 1:nSeg
            Tr(p,g,:) = interp1(t - vt(rel(g),1), F, TG, 'linear', NaN);
        end
    end
    Fbar0 = mean(Fref);
    Fbar  = squeeze(mean(Tr ./ Fref, 1, 'omitnan'));
    Fgrand = Sm.datatable(:,3) * (Fbar0/Fref(1));
    jump = nan(nSeg,1);
    for g = 1:nSeg
        tt = t0v - vt0(rel0(g),1);
        m  = tt >= TG(1) & tt <= TG(end);
        newv = Fbar0 * interp1(TG, Fbar(g,:), tt(m), 'linear', 'extrap');
        % the step left at the far edge of the splice
        iL = find(m); iL = iL(end);
        if iL < numel(Fgrand)
            jump(g) = Fgrand(iL+1) - newv(end);
        end
        Fgrand(m) = newv;
    end
    fP = extractSlackAttributes(t0v, Fgrand, L0, vt0, [], [], false, true);
    out(q) = struct('lbl',GRIDS{q,2}, 'rsK',mean(double(fP.rsK),'omitnan'), ...
        'rsA',mean(double(fP.rsA),'omitnan'), 'rsT0',mean(double(fP.rsT0),'omitnan'), ...
        'jump',mean(abs(jump),'omitnan'), 'Fgrand',Fgrand, 'tt',t0v);
end

%% ---- the mean-of-features reference ----------------------------------
per = nan(nP,nSeg); perA = per;
for p = 1:nP
    i = find(strcmp({R.day},days{p}) & strcmp({R.cond},'8mM'),1);
    v = double(R(i).feat.rsK(:))';  per(p,1:numel(v))  = v(1:min(nSeg,numel(v)));
    a = double(R(i).feat.rsA(:))';  perA(p,1:numel(a)) = a(1:min(nSeg,numel(a)));
end
mof = mean(per(:),'omitnan');

fprintf('\n=== rsK: mean-of-features vs features-of-mean-trace ===\n');
fprintf('%-42s %10s %9s %10s %9s\n','grid','rsK','gap','rsA (kPa)','splice step');
fprintf('%-42s %10.1f %8s %10.2f %9s\n','mean-of-features (per prep, then average)', ...
        mof, '--', mean(perA(:),'omitnan'), '--');
for q = 1:numel(out)
    fprintf('%-42s %10.1f %+7.1f%% %10.2f %8.2f kPa\n', out(q).lbl, out(q).rsK, ...
            100*(out(q).rsK/mof-1), out(q).rsA, out(q).jump);
end

fprintf('\nInterpretation: if extending the grid moves rsK back onto the\n');
fprintf('mean-of-features value, the -24 %% was the splice, not the averaging.\n');

%% ---- figure ----------------------------------------------------------
f = figure('Position',[80 80 1100 400],'Color','w');
tiledlayout(f,1,2,'Padding','compact');
nexttile; hold on; box on
g = 1; vs = vt0(seg0(g):seg0(g)+4,1); rt = vt0(rel0(g),1);
for q = 1:numel(out)
    m = out(q).tt >= vs(4)-0.01 & out(q).tt <= vs(5)+0.01;
    plot(1000*(out(q).tt(m)-vs(4)), out(q).Fgrand(m), 'LineWidth',1.3, ...
         'DisplayName', out(q).lbl);
end
xline(1000*(rt+0.570-vs(4)),'k--','DisplayName','splice edge (0.570 s)');
xlabel('time from ramp end (ms)'); ylabel('grand-average force (kPa)');
legend('Location','southeast','FontSize',7);
title('(a) the two grand-average traces in the recovery window');

nexttile; hold on; box on
bar([mof out(1).rsK out(2).rsK]);
set(gca,'XTick',1:3,'XTickLabel',{'mean-of-feat','DataPooling grid','extended grid'}, ...
    'XTickLabelRotation',20);
ylabel('rsK (1/s)'); title('(b) the -24 % and where it went');
exportgraphics(f, fullfile(resDir,'pooling_artefact.png'), 'Resolution',150);
save(fullfile(resDir,'pooling_artefact.mat'), 'out','mof','per','segLen');

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
