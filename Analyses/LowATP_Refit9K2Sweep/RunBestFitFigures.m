% RunBestFitFigures.m -- full-slide figure of the best high + low ATP fit(s):
% slack ladder (pooled band), ktr protocol (pooled band) and force-velocity,
% 8 mM and 2 mM, one figure per parameter set.
%
%   matlab -batch "run('Analyses/LowATP_Refit9K2Sweep/RunBestFitFigures.m')"
%
% Parameter sets (8 mM = the incumbent in both; gates act at 2 mM only):
%   A  k2only   R2 x0.54                 best k2-only (RunRefit9K2Sweep)
%   B  k2ksr0   R2 x0.54 + ksr0 x0.4     best of the k2 x SRX-entry grid (RunKsr0Grid)
if ~exist('CORE','var'); CORE = 'robust1_seed'; end
if ~exist('CFG','var')
    CFG = struct('tag', {'k2only', 'k2ksr0'}, ...
        'label', {'A: R_2 \times0.54 at 2 mM', 'B: R_2 \times0.54 + SRX entry k_{sr0} \times0.4 at 2 mM'}, ...
        'gates', {struct('pAtp_k2', log(0.54)/log(1/4)), ...
                  struct('pAtp_k2', log(0.54)/log(1/4), 'pAtp_ksr0', log(0.4)/log(1/4))}, ...
        'plaus', {'ATP binding inside A2 detachment: the canonical ATP step (plausible, tier A)', ...
                  '+ fewer ATP-loaded heads parked in SRX (plausible, unproven; tier A)'}, ...
        'score', {'2 mM \chi^2 (rel / plain) 4.52 / 7.86', '2 mM \chi^2 (rel / plain) 4.13 / 6.73'});
end
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resB = fullfile(root,'Analyses','LowATP_Refit9K2Sweep','results');
diary(fullfile(resB,'bestfit.log')); diary on;
clB = parcluster('local'); jslB = fullfile(tempdir, sprintf('mljob_bf_%d', feature('getpid')));
if ~isfolder(jslB); mkdir(jslB); end; clB.JobStorageLocation = jslB;
if ~exist('PLOTONLY','var') || ~PLOTONLY; if isempty(gcp('nocreate')); parpool(clB, 5); end; end

pB = getParams(loadParams(CORE), [], true, false);
pB.RunForceVelocity=false; pB.RunSlack=true; pB.RunStairs=false; pB.RunKtr=false;
pB.RunForceVelocityTime=false; pB.RunSlackPassive=false; pB.EvalFeatures=true;
pB.BreakOnODEUnstable=false; pB.PlotEachSeparately=0; pB.PlotFeatureFitting=0;
pB.RunSlackSegments='AllPar'; pB.MaxRunTime=400; pB.OptimizeOn='Feats';
pB.ghostSave=''; pB.ghostLoad=''; pB.SaveBest=0; pB.mods={}; pB.g=ones(1,0);
pB.recalculateDataFeats=false; pB.UsePooledTarget=true;
pB.UseAtpPowerGates=true; pB.MgATP_ref=8; pB.UseAtpDetachR2=0; pB.UseAtpKmsrd=0;
p8B = pB; p8B.MgATP = 8;
p8B.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
p8B.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr.mat';

LoadData;                                   % Data_ATP (Baker F-V, 8 and 2 mM)
DataFV = Data_ATP;
[BAND, TG, REL, Fref] = iSlackBands();
[KB, KT] = iKtrBands();
FV2 = fv2Target();

% 8 mM legs (shared). PLOTONLY = true re-draws from the saved bestfit_<tag>.mat.
if ~exist('PLOTONLY','var'); PLOTONLY = false; end
if ~PLOTONLY
    fprintf('8 mM ...\n');
    os8 = iSlack(p8B); ok8 = iKtr(p8B, 'protocol_03_27_2026_velocitytable_ktr.mat');
    [~,~,fv8] = runFVExperiment(p8B, [8 4 2], DataFV); fvp8 = fv2Model(p8B);
end

for iC = 1:numel(CFG)
    fprintf('%s ...\n', CFG(iC).tag);
    p2B = p8B; p2B.MgATP = 2;
    p2B.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
    gn = fieldnames(CFG(iC).gates);
    for j = 1:numel(gn); p2B.(gn{j}) = CFG(iC).gates.(gn{j}); end
    if PLOTONLY
        load(fullfile(resB, ['bestfit_' CFG(iC).tag '.mat'])); %#ok<LOAD>
    else
        os2 = iSlack(p2B); ok2 = iKtr(p2B, 'protocol_03_27_2026_velocitytable_ktr_2mM.mat');
        [~,~,fv2] = runFVExperiment(p2B, [2 4 8], DataFV); fvp2 = fv2Model(p2B);
        save(fullfile(resB, ['bestfit_' CFG(iC).tag '.mat']), 'os8','os2','ok8','ok2','fv8','fv2','fvp8','fvp2');
    end
    iFigure(CFG(iC), BAND, TG, REL, Fref, KB, KT, {os8, os2}, {ok8, ok2}, fv8, fv2, fvp8, fvp2, FV2, DataFV, ...
            fullfile(resB, ['bestfit_' CFG(iC).tag '.png']));
end
fprintf('DONE\n'); diary off;

%% ======================================================================
function o = iSlack(p)
    params0 = p;                                               %#ok<NASGU>
    LoadData; RunBakersExp;                                    %#ok<*ASGLU>
    o = struct('t', out_slack.t(:), 'F', out_slack.Force(:));
end
function o = iKtr(p, file)
    p.ktr_velocitytableonfile = file;
    [~, out] = runKtrExperiment(p);
    S = load(fullfile('data', file)); vt = S.velocitytable;
    k = makeMonotonous(out.t(:)); t = out.t(:); t = t(k); F = out.Force(:); F = F(k);
    Fss = interp1(t, F, vt(2,1), 'linear', 'extrap');
    o = struct('t', t - vt(end-1,1), 'F', F/Fss);              % t = 0 at the final release
end
function [BAND, TG, REL, Fref] = iSlackBands()
    evalc('T = buildPooledTarget(struct());');
    Fref = T.scale.mean; lift = T.scale.lift; alpha = T.scale.perPrep(:);
    DAYS = {'03_27_2026','04_03_2026','04_10_2026'}; CONDS = {'8mM','2mM'};
    TG = (-0.02 : 2e-4 : 0.57)'; BAND = nan(5, numel(TG), 3, 2); REL = cell(1,2);
    for c = 1:2
        for g = 1:5
            V = nan(3, numel(TG));
            for p = 1:3
                S = load(fullfile('data', sprintf('protocol_%s_%s_slack.mat', DAYS{p}, CONDS{c})));
                rl = S.velocitytable(S.velocitytable(:,2) < -20, 1);
                if p == 1; REL{c} = rl; end
                V(p,:) = interp1(S.datatable(:,1) - rl(g), S.datatable(:,3)*lift(p,c)/alpha(p)*Fref, TG, 'linear', NaN);
            end
            BAND(g,:,1,c) = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
            BAND(g,:,2,c) = BAND(g,:,1,c) - sd; BAND(g,:,3,c) = BAND(g,:,1,c) + sd;
        end
    end
end
function [KB, KT] = iKtrBands()
    DAYS = {'03_27','04_03','04_10'}; SUF = {'', '_2mM'};
    KT = (-0.03 : 2e-4 : 0.058)'; KB = nan(numel(KT), 3, 2);
    for c = 1:2
        V = nan(3, numel(KT));
        for p = 1:3
            S = load(fullfile('data', sprintf('protocol_%s_2026_velocitytable_ktr%s.mat', DAYS{p}, SUF{c})));
            V(p,:) = interp1(S.datatable(:,1) - S.velocitytable(end-1,1), S.datatable(:,3), KT, 'linear', NaN);
        end
        KB(:,1,c) = mean(V,1,'omitnan'); sd = std(V,0,1,'omitnan');
        KB(:,2,c) = KB(:,1,c) - sd(:); KB(:,3,c) = KB(:,1,c) + sd(:);
    end
end
function iFigure(C, BAND, TG, REL, Fref, KB, KT, OS, OK, fv8, fv2, fvp8, fvp2, FV2, D, outf)
    navy = [0 39 76]/255; red = [0.75 0.2 0.12]; CC = {navy, red};
    f = figure('Color','w','Position',[20 20 2000 1000],'Visible','off');
    tl = tiledlayout(f, 2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
    for g = 1:5
        nexttile(g); hold on; box on; grid on;
        for c = 1:2
            fill(1000*[TG; flipud(TG)], [BAND(g,:,2,c)'; flipud(BAND(g,:,3,c)')], CC{c}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(1000*TG, BAND(g,:,1,c), '-', 'Color', CC{c}, 'LineWidth', 1.5);
            o = OS{c}; r = REL{c}(g); m = o.t - r >= 0 & o.t - r <= TG(end);
            plot(1000*(o.t(m) - r), o.F(m), '--', 'Color', CC{c}, 'LineWidth', 2);
        end
        ylim([-5 130]); xlim(1000*[TG(1) TG(end)]); set(gca, 'FontSize', 11);
        title(sprintf('slack cycle %d', g), 'FontSize', 12); xlabel('ms from release');
        if g == 1; ylabel('force (kPa, pooled level)'); end
    end
    % ktr protocol
    nexttile(6, [1 2]); hold on; box on; grid on;
    for c = 1:2
        okb = all(isfinite(KB(:,2:3,c)), 2); kt = KT(okb);   % fill() draws nothing if any vertex is NaN
        fill(1000*[kt; flipud(kt)], [KB(okb,2,c); flipud(KB(okb,3,c))], CC{c}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(1000*KT, KB(:,1,c), '-', 'Color', CC{c}, 'LineWidth', 1.5);
        o = OK{c}; m = o.t >= KT(1) & o.t <= KT(end);
        plot(1000*o.t(m), o.F(m), '--', 'Color', CC{c}, 'LineWidth', 2);
    end
    xlim(1000*[KT(1) KT(end)]); ylim([-0.1 1.4]); set(gca, 'FontSize', 11);
    xlabel('ms from final release'); ylabel('F / own pre-shortening isometric');
    title('k_{tr} protocol (pooled, 3 preps)', 'FontSize', 12);
    % F-V normalised
    nexttile(8); hold on; box on; grid on;
    plot(D(:,1), D(:,2)/D(1,2), 'o', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 6);
    plot(D(:,1), D(:,4)/D(1,4), 'o', 'Color', red, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.3);
    plot(fv8.FV_v, fv8.FV_fnorm, '--', 'Color', navy, 'LineWidth', 2);
    plot(fv2.FV_v, fv2.FV_fnorm, '--', 'Color', red, 'LineWidth', 2);
    errorbar(2, FV2.m8, FV2.s8, 'd', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 9, 'LineWidth', 1.5);
    errorbar(2.12, FV2.m2, FV2.s2, 'd', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 9, 'LineWidth', 1.5);
    plot(2, fvp8, 'x', 'Color', navy, 'MarkerSize', 13, 'LineWidth', 2.5); plot(2.12, fvp2, 'x', 'Color', red, 'MarkerSize', 13, 'LineWidth', 2.5);
    xlim([0 6]); ylim([0 1.05]); set(gca, 'FontSize', 11); xlabel('shortening velocity (ML/s)'); ylabel('F / F_0');
    title('force–velocity', 'FontSize', 12);
    % power
    nexttile(9); hold on; box on; grid on;
    plot(D(:,1), D(:,1).*D(:,2)/D(1,2), 'o', 'Color', navy, 'MarkerFaceColor', navy, 'MarkerSize', 6);
    plot(D(:,1), D(:,1).*D(:,4)/D(1,4), 'o', 'Color', red, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.3);
    plot(fv8.FV_v, fv8.FV_v(:).*fv8.FV_fnorm(:), '--', 'Color', navy, 'LineWidth', 2);
    plot(fv2.FV_v, fv2.FV_v(:).*fv2.FV_fnorm(:), '--', 'Color', red, 'LineWidth', 2);
    xlim([0 6]); set(gca, 'FontSize', 11); xlabel('shortening velocity (ML/s)'); ylabel('normalised power (F/F_0 \cdot v)');
    title('power', 'FontSize', 12);
    % legend / label tile
    nexttile(10); axis off; hold on;
    h = gobjects(0);
    h(end+1) = plot(nan, nan, '-', 'Color', navy, 'LineWidth', 2); h(end+1) = plot(nan, nan, '-', 'Color', red, 'LineWidth', 2);
    h(end+1) = plot(nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    h(end+1) = plot(nan, nan, 'o', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3]);
    h(end+1) = plot(nan, nan, 'o', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', 'w');
    h(end+1) = plot(nan, nan, 'd', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3]);
    h(end+1) = plot(nan, nan, 'x', 'Color', [0.3 0.3 0.3], 'MarkerSize', 10, 'LineWidth', 2);
    legend(h, {'8 mM data (band = \pm SD)', '2 mM data', 'model', 'Baker F–V 8 mM (scored)', 'Baker F–V 2 mM (not scored)', ...
               'FV2 in-prep, 2 ML/s (n=2)', 'model FV2'}, 'Location', 'north', 'FontSize', 10);
    title(tl, ['Best high + low ATP fit — ' C.label], 'FontSize', 15);
    exportgraphics(f, outf, 'Resolution', 170); close(f);
end
