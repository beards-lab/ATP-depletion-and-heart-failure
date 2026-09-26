%% AnalyzeModelNoise.m -- per feature definition: numerical noise, bias, phase sensitivity.
% Reads results/model_noise_<PARAMSET>.mat (RunModelNoise). For every solver
% configuration and every read of the identical traces (see RunModelNoise):
%   num   : relative SD over the 1 + NREP points 1e-6 apart     (the noise)
%   bias  : unperturbed value vs the reference configuration E  (the error)
%   phase : relative SD over 10 sub-millisecond shifts of the 1 kHz read grid
% each as an RMS over the cycles the feature is scored on.
if ~exist('PARAMSET', 'var'), PARAMSET = 'refit9_opt'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(here, '..', '..')));
res = fullfile(here, 'results');
L = load(fullfile(res, ['model_noise_' PARAMSET '.mat'])); M = L.M; vt = L.vt;
cf = fieldnames(M); cf = cf(:)';
PIPE = {'rsK80','rsKfit','rsT20','rsTfit','peak1_t','peak1_dSL','vall_t','ktr80','ktrfit','t0_crossing','t0','restretchSlopeStart','steady'};
ROBF = {'rsK80','rsKfit','rsT20','rsT20e','rsTfit','rsTfitE','peak1_t','peak1_tq','peak1_dSL','peak1_dSLq', ...
        'vall_t','vall_tq','tFall50','ktr80','ktrfit','tOn20','t0x','restretchSlopeStart','rsSlope1ms','steady','rsKls','rsT0ls','peak1_y','peak1_yq','vall_y','vall_yq'};
SEG14 = {'rsK80','rsKfit','rsT20','rsT20e','rsTfit','rsTfitE','rsKls','rsT0ls'};
rd = @(tr, dt, ph) readAt(tr, vt, dt, ph);

X = struct();
for c = cf
    O = M.(c{1}); n = numel(O.TR);
    for r = 1:n
        X.(c{1}).pipe(r) = O.FM{r};
        X.(c{1}).a1k(r)  = rd(O.TR{r}, 1e-3, 0);
        X.(c{1}).a4k(r)  = rd(O.TR{r}, 2.5e-4, 0);
    end
    for ph = 0:9
        X.(c{1}).phase(ph+1)   = rd(O.TR{1}, 1e-3, ph*1e-4);
        X.(c{1}).phase4k(ph+1) = rd(O.TR{1}, 2.5e-4, ph*2.5e-5);
    end
end
ref = X.(cf{end});                   % the last configuration is the reference (E)

rows = {};
for c = cf
    for mode = {'pipe', 'a1k', 'a4k'}
        FL = PIPE; if ~strcmp(mode{1}, 'pipe'), FL = ROBF; end
        for f = FL
            sg = 1:5; if any(strcmp(f{1}, SEG14)), sg = 1:4; end
            V = cell2mat(arrayfun(@(s) take(s, f{1}), X.(c{1}).(mode{1})(:), 'uni', 0));
            mu = mean(V, 1, 'omitnan');
            num = rmsv(std(V(:, sg), 0, 1, 'omitnan') ./ abs(mu(sg)));
            rf = take(ref.a4k(1), f{1}); if all(isnan(rf)), rf = take(ref.pipe(1), f{1}); end
            bias = rmsv((V(1, sg) - rf(sg)) ./ abs(rf(sg)));
            ph = NaN;
            if ~strcmp(mode{1}, 'pipe')
                pf = 'phase'; if strcmp(mode{1}, 'a4k'), pf = 'phase4k'; end
                Vp = cell2mat(arrayfun(@(s) take(s, f{1}), X.(c{1}).(pf)(:), 'uni', 0));
                ph = rmsv(std(Vp(:, sg), 0, 1, 'omitnan') ./ abs(mean(Vp(:, sg), 1, 'omitnan')));
            end
            rows(end+1, :) = {M.(c{1}).name(1), mode{1}, f{1}, num, bias, ph, mat2str(round(V(1, :), 4))}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'cfg','read','feature','num','bias','phase','value'});
T.cfg = string(T.cfg); T.read = string(T.read); T.feature = string(T.feature);
save(fullfile(res, ['model_noise_summary_' PARAMSET '.mat']), 'T', 'X');
writetable(T, fullfile(res, ['model_noise_summary_' PARAMSET '.csv']));

% ---- compact print: one row per feature, columns = configs ----------------
for mode = {'pipe', 'a1k', 'a4k'}
    fprintf('\n==== read: %s   (num = relSD over 1e-6 points | bias vs E | phase = 1 kHz grid shift) ====\n', mode{1});
    fprintf('%-20s', 'feature'); for c = cf, fprintf('  %-22s', M.(c{1}).name(1:min(22, end))); end; fprintf('\n');
    FL = unique(T.feature(T.read == mode{1}), 'stable');
    for f = FL'
        fprintf('%-20s', f);
        for c = cf
            k = T.read == mode{1} & T.feature == f & T.cfg == string(M.(c{1}).name(1));
            if ~strcmp(mode{1}, 'pipe')
                fprintf('  %5.2f %6.2f %6.2f %%  ', 100*T.num(k), 100*T.bias(k), 100*T.phase(k));
            else
                fprintf('  %5.2f %6.2f %%        ', 100*T.num(k), 100*T.bias(k));
            end
        end
        fprintf('\n');
    end
end
fprintf('\nslack cost (current definitions, pipeline):\n');
for c = cf, O = M.(c{1}); fprintf('  %-44s base %.4f  SD %.4f (%.1f %%)\n', O.name, O.cost(1), std(O.cost), 100*std(O.cost)/mean(O.cost)); end

function F = readAt(tr, vt, dt, ph)
    [t, iu] = unique(tr.t); y = tr.F(iu); s = tr.SL(iu);
    g = ((ceil(t(1)/dt) : floor(t(end)/dt) - 1)' + ph/dt) * dt;
    g = g(g >= t(1) & g <= t(end));
    F = robustSlackFeatures(g, interp1(t, y, g), interp1(t, s, g), vt, false, t, tr.LSE(iu));
end

function v = take(s, f)
    if isfield(s, f), v = double(s.(f)(:))'; else, v = nan(1, 5); end
    if numel(v) < 5, v(end+1:5) = NaN; end; v = v(1:5);
end

function r = rmsv(x), r = sqrt(mean(x(isfinite(x)).^2)); end
