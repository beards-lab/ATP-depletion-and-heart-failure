% RunMechanismTest.m — what in the model makes a symmetric input asymmetric?
%
% The sawtooth commands a length trajectory that is exactly symmetric in time
% about the apex, yet the force response is not: at slow ramp speeds the model
% returns MORE work than it absorbs (negative hysteresis), the attached fraction
% at the end of the release exceeds its isometric value, and force OVERSHOOTS F0
% after the muscle is back at L0.
%
% The suspect is directional registration availability. params/estiff125_opt.m
% has UseRegistrationAvailability = 1 with RegAvailShorteningOnly = 1, i.e.
% sliding refreshes myosin/actin registration and so raises attachment, but ONLY
% while SHORTENING, half-recovering at v_ref_reg = 0.63. That makes the two
% halves of a symmetric sawtooth run at different availability by construction.
% It was introduced to fix the steady-state force-velocity shoulder
% (Analyses/FV_Shoulder), so a transient test of it is close to independent
% evidence rather than a restatement of the fit.
%
% Test: rerun the sawtooth with the mechanism ON and OFF and see whether the
% asymmetry survives. Nothing else is changed.
%
% Writes results/mechanism.mat

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..'); cd(root); addpath(genpath(root));
resdir = fullfile(here, 'results');
C = load(fullfile(resdir, 'characterize.mat'));

base = getParams(loadParams('params/estiff125_opt.m'), [], true, false);
base.ghostLoad = ''; base.ghostSave = ''; base.SaveBest = false;
base.PlotEachSeparately = false; base.PlotProbsOnFig = false;
base.ValuesInTime = true; base.MaxRunTime = 320;

L0 = 1.00; FRAC = 0.04; T_WARM = 2.00; T_POST = 1.00; SRS_EPS = 0.0045;
V_VERB = FRAC / 0.100;
vlist  = [0.1, 0.2, V_VERB, 0.35*C.Vmax];
onoff  = [true false];

M = struct();
fprintf('%-6s %-5s | %8s %8s %8s %8s %8s %8s\n', 'v', 'reg', 'F0kPa', ...
        'hyst', 'elRet%', 'att(t2)', 'postMax', 'SRSeps');
for i = 1:numel(vlist)
    [vt, meta] = buildSawtoothProtocol(L0, FRAC, vlist(i), T_WARM, T_POST, base.ML);
    for k = 1:2
        p = base;
        p.UseRegistrationAvailability = onoff(k);
        out = iRun(p, vt);
        a   = extractSawtoothAttributes(out, meta, SRS_EPS);
        M(i,k).v = vlist(i); M(i,k).reg = onoff(k); M(i,k).a = a; M(i,k).meta = meta;
        M(i,k).t = out.t; M(i,k).Fn = out.Force / a.F0_kPa;
        M(i,k).att = out.p1_0 + out.p2_0;
        fprintf('%-6.3f %-5d | %8.2f %8.3f %8.0f %8.3f %8.3f %8.2f\n', ...
                vlist(i), onoff(k), a.F0_kPa, a.W_hyst, 100*a.elasticReturn, ...
                a.att_t2, a.Fn_postMax, a.SRS_eps);
    end
end

fprintf('\n--- effect of switching the mechanism OFF ---\n');
for i = 1:numel(vlist)
    on = M(i,1).a; off = M(i,2).a;
    fprintf('v=%5.3f ML/s: hysteresis %+.3f -> %+.3f | att(t2) %.3f -> %.3f (iso %.3f -> %.3f) | postMax %.3f -> %.3f\n', ...
            vlist(i), on.W_hyst, off.W_hyst, on.att_t2, off.att_t2, ...
            on.att_t0, off.att_t0, on.Fn_postMax, off.Fn_postMax);
end

%% Figure -----------------------------------------------------------------
f = figure('Visible', 'off', 'Position', [50 50 1180 700]);
tiledlayout(2, numel(vlist), 'Padding', 'compact');
for row = 1:2
    for i = 1:numel(vlist)
        nexttile; hold on; box on;
        for k = 1:2
            m = M(i,k).meta;
            if row == 1
                plot(1000*(M(i,k).t - m.t2), M(i,k).Fn, '-', 'LineWidth', 1.6);
            else
                plot(1000*(M(i,k).t - m.t2), M(i,k).att, '-', 'LineWidth', 1.6);
            end
        end
        if row == 1
            yline(1, ':k'); ylabel('F / F_0');
        else
            ylabel('attached fraction');
        end
        xlim([-1000*2*M(i,1).meta.dt, 300]);
        xlabel('t - t(back at L_0) (ms)');
        title(sprintf('v = %.2f ML/s', vlist(i)), 'FontSize', 9);
        if i == 1 && row == 1
            legend({'reg. availability ON', 'OFF'}, 'Location', 'best', 'FontSize', 7);
        end
    end
end
sgtitle('Directional registration availability is what makes the symmetric sawtooth asymmetric', 'FontSize', 10);
try
    exportgraphics(f, fullfile(resdir, 'mechanism.png'), 'Resolution', 150);
catch
    print(f, fullfile(resdir, 'mechanism.png'), '-dpng', '-r150');
end
close(f);

save(fullfile(resdir, 'mechanism.mat'), 'M', 'vlist', 'onoff');
fprintf('\nSaved results/mechanism.mat and results/mechanism.png\n');


% =========================================================================
function out = iRun(p, vt)
    p.SL0 = 2 * vt(1, 4);
    p.Slim_l = 1.4; p.Slim_r = 2.3;
    if isfield(p, 'PU0'); p = rmfield(p, 'PU0'); end
    p = getParams(p, p.g, true);
    p.Velocity = vt(:, 2);
    lastwarn('');
    [~, out] = evaluateModel(str2func(p.modelFcn), vt(:, 1), p);
    if ~isempty(lastwarn); fprintf('    [warn] %s\n', lastwarn); end
end
