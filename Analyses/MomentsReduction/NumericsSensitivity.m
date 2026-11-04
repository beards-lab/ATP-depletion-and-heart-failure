% NumericsSensitivity.m
% How much of the ktr / slack "truth" is numerics? Feature table across PDE
% variants (grid spacing, window width, kernel extrapolation clamp) and the
% Gaussian moments model with and without the clamp. Errors relative to the
% production PDE (dS 2 nm, 120 nm window, unclamped).

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
P = momentsProtocols();
variants = {'ref_', 'production (dS 2, W 120 nm)'; 'refdS1_', 'dS 1 nm'; 'refdS0.5_', 'dS 0.5 nm'; ...
            'ref_W2_', 'window 240 nm'; 'ref_clamp_', 'clamped kernels'; 'ref_W2_clamp_', 'window 240 + clamp'};
rows = {};
for name = {'ktr', 'slack'}
    Pp = P(strcmp({P.name}, name{1}));
    for iv = 1:size(variants, 1)
        fnm = fullfile(here, 'results', [variants{iv,1} name{1} '.mat']);
        if ~exist(fnm, 'file'), continue; end
        S = load(fnm); [t, i] = unique(S.ref.t);
        f = protocolFeatures(Pp, t, S.ref.Force(i));
        rows(end+1,:) = {name{1}, variants{iv,2}, f.F_peak_post, f.k80, S.ref.wall}; %#ok<SAGROW>
    end
    for clampK = [false true]
        pm = params0; pm.ClampStrainDepExtrapolation = clampK;
        o = evaluateMoments(Pp.T, Pp.V, pm, 'gauss', 3);
        [t, i] = unique(o.t); f = protocolFeatures(Pp, t, o.Force(i));
        rows(end+1,:) = {name{1}, sprintf('moments gauss3%s', repmat(' + clamp', 1, clampK)), f.F_peak_post, f.k80, o.wall}; %#ok<SAGROW>
    end
end
Tn = cell2table(rows, 'VariableNames', {'protocol','variant','F_peak_post','k80','wall_s'});
for name = {'ktr','slack'}
    w = strcmp(Tn.protocol, name{1}); b = find(w, 1);
    Tn.k80_vs_prod_pct(w) = 100*(Tn.k80(w) - Tn.k80(b))/Tn.k80(b);
    Tn.peak_vs_prod_pct(w) = 100*(Tn.F_peak_post(w) - Tn.F_peak_post(b))/Tn.F_peak_post(b);
end
disp(Tn); writetable(Tn, fullfile(here, 'results', 'numerics_sensitivity.csv'));
