% RunFinalMetricChoice.m -- the decision table.
%
% Established by the three scripts before this one:
%   * the "-24 % trace averaging" result was ONE broken segment on the pooled
%     trace (segment 5, whose window is not covered by the splice grid), not a
%     systematic smearing: segments 1-4 pool to within 2 % (CheckPoolingArtefact);
%   * segment 5 is a DIFFERENT MANOEUVRE -- re-stretch to ML 1.00 not 1.10, the
%     SMALLEST re-stretch (0.061 ML) at the slowest ramp, window 3x longer -- and
%     it supplies 76 % of the apparent dk/dSLslack (RunSegment5Check);
%   * every estimator's span- and dilation-sensitivity came from the shared
%     amplitude convention A = median(last 15 % of the window), not from the
%     rate formula. Taking A from `steady` instead collapses the 60 ms span
%     error from +77 % to -6.7 % and the identity error from 24 % to 4 %
%     (RunAmplitudeIsTheWeakLink).
%
% So the real choice is a 5 x 2 grid: rate formula x amplitude source. This
% scores all ten on what a fit target actually needs.
%
% Output -> results/final_choice.mat, results/final_choice.png
%           results/rsk_candidates.csv  (the numbers, for the write-up)

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root)); addpath(here);
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

EST = {'k_fit','k_t0z','k_band','k_t50','k_2080'};
LBL = {'current','t0=0','band15-90','t50','t20-t80'};
AMP = {'tail','steady'};
nE = numel(EST); nA = 2;
days = {'03/27','04/03','04/10'}; COND = {'8mM','2mM'};
nP = 3; nC = 2; nS = 5;
CROSSWIN = 0.005;      % fixed a priori, not tuned: 5 ms, as extractSlackAttributes

%% ---- windows -----------------------------------------------------------
W = cell(nP,nC,nS); Aext = nan(nP,nC,nS); SLs = nan(nP,nC,nS); dRS = nan(nS,1);
for p = 1:nP, for c = 1:nC
    S = load(fullfile('data', slackFile(days{p},COND{c})));
    t = S.datatable(:,1); F = S.datatable(:,3); L = S.datatable(:,2);
    v = S.velocitytable; sg = find(v(:,2) < -1)';
    for g = 1:nS
        r = sg(g):sg(g)+4; vs = v(r,1);
        w = t >= vs(4) & t <= vs(5);
        W{p,c,g} = struct('t',t(w)-vs(4),'y',F(w));
        E0 = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true, struct('crossWin',CROSSWIN));
        Aext(p,c,g) = median(F(t >= vs(5)-0.02 & t <= vs(5))) - E0.B;
        SLs(p,c,g)  = median(L(t >= vs(2) & t <= vs(3)));
        if p==1 && c==1; dRS(g) = v(r(4),4) - v(r(2),4); end
    end
end, end

%% ---- every estimator x amplitude source -------------------------------
K = nan(nP,nC,nS,nE,nA);
for p = 1:nP, for c = 1:nC, for g = 1:nS
    for a = 1:nA
        o = struct('crossWin',CROSSWIN);
        if a == 2; o.Aext = Aext(p,c,g); end
        E0 = rsEstimators(W{p,c,g}.t, W{p,c,g}.y, true, o);
        for e = 1:nE; K(p,c,g,e,a) = E0.(EST{e}); end
    end
end, end, end

%% ---- identity + span, per combination ---------------------------------
CF = [0.6 0.8 1.25 1.6667]; SPAN = 0.060;
ID = nan(nE,nA); SP = nan(nE,nA);
for a = 1:nA
    ri = nan(nP,nC,nS,nE,numel(CF)); rs = nan(nP,nC,nS,nE);
    for p = 1:nP, for c = 1:nC, for g = 1:nS
        o = struct('crossWin',CROSSWIN); if a==2; o.Aext = Aext(p,c,g); end
        for j = 1:numel(CF)
            [td,yd] = dilate(W{p,c,g}.t, W{p,c,g}.y, CF(j));
            Ed = rsEstimators(td, yd, true, o);
            for e = 1:nE; ri(p,c,g,e,j) = Ed.(EST{e})/K(p,c,g,e,a)/CF(j) - 1; end
        end
        m = W{p,c,g}.t <= SPAN;
        if nnz(m) >= 25
            Ed = rsEstimators(W{p,c,g}.t(m), W{p,c,g}.y(m), true, o);
            for e = 1:nE; rs(p,c,g,e) = Ed.(EST{e})/K(p,c,g,e,a) - 1; end
        end
    end, end, end
    for e = 1:nE
        v = ri(:,:,:,e,:); ID(e,a) = max(median(abs(reshape(v,[],numel(CF))),1,'omitnan'));
        v = rs(:,:,:,e);   SP(e,a) = median(abs(v(:)),'omitnan');
    end
end

%% ---- target quality, per combination ----------------------------------
% Reported on segments 1-4 (the true amplitude ladder: 0.080 -> 0.141 ML
% re-stretch at a common final length of ML 1.10) AND on all five.
fprintf('=== re-stretch amplitude per segment (ML) : %s\n', sprintf('%.3f ',dRS));
fprintf('    segment 5 is the SMALLEST re-stretch, to a different final length\n\n');

RES = struct('lbl',{},'amp',{},'segs',{},'val',{},'sB',{},'sW',{},'bio',{}, ...
             'r',{},'sr',{},'z',{},'slope',{},'id',{},'sp',{});
for a = 1:nA, for e = 1:nE, for sset = 1:2
    if sset == 1; gg = 1:4; tag = '1-4'; else; gg = 1:5; tag = '1-5'; end
    pm = nan(nP,1); pw = nan(nP,1); rt = nan(nP,1); sl = nan(nP,1);
    for p = 1:nP
        v8 = squeeze(K(p,1,gg,e,a)); v2 = squeeze(K(p,2,gg,e,a));
        x  = squeeze(SLs(p,1,gg));
        pm(p) = mean(v8,'omitnan'); rt(p) = mean(v2,'omitnan')/pm(p);
        gd = isfinite(v8) & isfinite(x);
        if nnz(gd) >= 3
            pf = polyfit(x(gd), v8(gd), 1);
            pw(p) = std(v8(gd)-polyval(pf,x(gd)))/abs(pm(p));
            sl(p) = pf(1);
        end
    end
    sB = std(pm,'omitnan')/mean(pm,'omitnan');
    sW = mean(pw,'omitnan'); sMeas = sW/sqrt(numel(gg));
    RES(end+1) = struct('lbl',LBL{e},'amp',AMP{a},'segs',tag, ...
        'val',mean(pm,'omitnan'),'sB',sB,'sW',sW, ...
        'bio',max(0,1-(sMeas/sB)^2),'r',mean(rt,'omitnan'), ...
        'sr',std(rt,'omitnan')/abs(mean(rt,'omitnan')), ...
        'z',abs(mean(rt,'omitnan')-1)/sqrt(sB^2+(std(rt,'omitnan')/mean(rt,'omitnan'))^2), ...
        'slope',mean(sl,'omitnan'),'id',ID(e,a),'sp',SP(e,a)); %#ok<SAGROW>
end, end, end

fprintf('=== the decision table (segments 1-4; 8 mM value, spreads, ATP, robustness) ===\n');
fprintf('%-11s %-7s %8s %8s %8s %6s %8s %7s %7s %9s %8s %8s\n', ...
        'rate','A from','k(1/s)','s_betw','s_with','bio%','r 2/8','s_r','z','dk/dSL','ident','span60');
for i = 1:numel(RES)
    if ~strcmp(RES(i).segs,'1-4'); continue; end
    fprintf('%-11s %-7s %8.1f %7.1f%% %7.1f%% %5.0f%% %8.3f %6.1f%% %7.2f %9.0f %7.1f%% %7.1f%%\n', ...
        RES(i).lbl, RES(i).amp, RES(i).val, 100*RES(i).sB, 100*RES(i).sW, ...
        100*RES(i).bio, RES(i).r, 100*RES(i).sr, RES(i).z, RES(i).slope, ...
        100*RES(i).id, 100*RES(i).sp);
end
fprintf('\n=== the same on all five segments (for comparison) ===\n');
fprintf('%-11s %-7s %8s %8s %7s %9s\n','rate','A from','k(1/s)','s_betw','z','dk/dSL');
for i = 1:numel(RES)
    if ~strcmp(RES(i).segs,'1-5'); continue; end
    fprintf('%-11s %-7s %8.1f %7.1f%% %7.2f %9.0f\n', RES(i).lbl, RES(i).amp, ...
        RES(i).val, 100*RES(i).sB, RES(i).z, RES(i).slope);
end

%% ---- a composite that says what it is optimising ---------------------
% A target feature wants: small between-prep spread (so a residual means
% something), high ATP z (so it carries the effect we are after), and small
% identity/span distortion (so it measures a RATE and not the window).
fprintf('\n=== ranking (segments 1-4) ===\n');
fprintf('  score = s_betw + ident + span60 - z/20   (lower better)\n');
sc = nan(numel(RES),1);
for i = 1:numel(RES)
    if ~strcmp(RES(i).segs,'1-4'); continue; end
    sc(i) = RES(i).sB + RES(i).id + RES(i).sp - RES(i).z/20;
end
[~, ord] = sort(sc);
for i = ord(:)'
    if isnan(sc(i)); continue; end
    fprintf('  %6.3f  %-11s A=%-7s  s=%.1f%% z=%.2f id=%.1f%% sp=%.1f%%\n', ...
        sc(i), RES(i).lbl, RES(i).amp, 100*RES(i).sB, RES(i).z, ...
        100*RES(i).id, 100*RES(i).sp);
end

%% ---- csv ------------------------------------------------------------
T = struct2table(RES);
writetable(T, fullfile(resDir,'rsk_candidates.csv'));
save(fullfile(resDir,'final_choice.mat'),'RES','K','ID','SP','Aext','SLs','EST','LBL','AMP','dRS');

%% ---- figure ---------------------------------------------------------
f = figure('Position',[60 60 1250 460],'Color','w');
tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
i14 = strcmp({RES.segs},'1-4');
sB = 100*[RES(i14).sB]; zz = [RES(i14).z]; idv = 100*[RES(i14).id];
nm = strcat({RES(i14).lbl},'/',{RES(i14).amp});

nexttile; hold on; box on
scatter(sB, zz, 40+3*idv, idv, 'filled');
for i = 1:numel(nm); text(sB(i)+0.3, zz(i), nm{i}, 'FontSize',7); end
xlabel('between-preparation SD (%)'); ylabel('ATP discriminability z');
cb = colorbar; cb.Label.String = 'identity error (%)'; grid on
title('(a) spread vs ATP signal (size/colour = distortion)');

nexttile; hold on; box on
bar([[RES(i14).slope]' [RES(strcmp({RES.segs},'1-5')).slope]']);
set(gca,'XTick',1:nnz(i14),'XTickLabel',nm,'XTickLabelRotation',40,'FontSize',7);
ylabel('dk/dSLslack (1/s per ML)');
legend({'segments 1-4','segments 1-5'},'Location','southwest','FontSize',8);
title('(b) the slope the model must match, and segment 5''s share');

nexttile; hold on; box on
cols = lines(nP);
[~, bi] = min(sc);
be = find(strcmp(LBL, RES(bi).lbl)); ba = find(strcmp(AMP, RES(bi).amp));
for p = 1:nP
    plot(dRS, squeeze(K(p,1,:,1,1)), 'o--','Color',cols(p,:),'LineWidth',1.1, ...
         'DisplayName',[days{p} ' current/tail']);
    plot(dRS, squeeze(K(p,1,:,be,ba)), 's-','Color',cols(p,:),'LineWidth',1.6, ...
         'DisplayName',[days{p} ' ' RES(bi).lbl '/' RES(bi).amp]);
end
xlabel('re-stretch amplitude (ML)'); ylabel('rate (1/s)');
legend('Location','northwest','FontSize',6);
title('(c) against the RIGHT x-axis: re-stretch amplitude');
exportgraphics(f, fullfile(resDir,'final_choice.png'), 'Resolution',150);
fprintf('\nwrote %s\n', fullfile(resDir,'final_choice.png'));

function [td, yd] = dilate(t, y, c)
    t0 = t(1); tsrc = t0 + (t - t0)/c; td = t;
    yd = interp1(tsrc, y, td, 'linear');
    keep = isfinite(yd); td = td(keep); yd = yd(keep);
end

function f = slackFile(day, cond)
    tag = struct('x03_27','03_27_2026','x04_03','04_03_2026','x04_10','04_10_2026');
    f = sprintf('protocol_%s_%s_slack.mat', tag.(['x' strrep(day,'/','_')]), cond);
end
