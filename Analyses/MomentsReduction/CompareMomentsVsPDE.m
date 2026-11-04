% CompareMomentsVsPDE.m
% A-posteriori test: run the closed moments model (dMdT_Moments) on every
% protocol and compare with the PDE reference (results/ref_*.mat) at the
% level of traces AND the features the objective/physiology care about.
% All errors are in % of the PDE isometric force (one common scale).

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
P = momentsProtocols();
variants = {'gauss', 3; 'qmom', 4; 'qmom', 6};
S = load(fullfile(here, 'results', 'ref_iso.mat')); Fiso = S.ref.Force(end);

R = struct(); rows = {};
for ip = 1:numel(P)
    S = load(fullfile(here, 'results', ['ref_' P(ip).name '.mat'])); ref = S.ref;
    [tu, iu] = unique(ref.t); Fr = ref.Force(iu);
    fr = protocolFeatures(P(ip), tu, Fr);
    for iv = 1:size(variants, 1)
        tag = sprintf('%s%d', variants{iv,1}, variants{iv,2});
        try
            o = evaluateMoments(P(ip).T, P(ip).V, params0, variants{iv,1}, variants{iv,2});
            Fm = interp1(o.t, o.Force, tu, 'linear', 'extrap');
            fm = protocolFeatures(P(ip), tu, Fm);
            win = tu >= 1.0;                                   % protocol proper (after activation)
            rmsE = 100*sqrt(mean((Fm(win) - Fr(win)).^2))/Fiso;
            maxE = 100*max(abs(Fm(win) - Fr(win)))/Fiso;
            R.(matlab.lang.makeValidName(P(ip).name)).(tag) = struct('o', o, 'fm', fm, 'fr', fr, 'rms', rmsE, 'max', maxE);
            fn = fieldnames(fr);
            for k = 1:numel(fn)
                rows(end+1, :) = {P(ip).name, tag, fn{k}, fr.(fn{k}), fm.(fn{k}), o.wall, ref.wall, rmsE, maxE}; %#ok<SAGROW>
            end
            fprintf('%-7s %-7s  rms %5.2f%%  max %6.2f%%  wall %5.2f s (PDE %5.1f s)  nf %d\n', ...
                P(ip).name, tag, rmsE, maxE, o.wall, ref.wall, o.nfevals);
        catch e
            fprintf('%-7s %-7s  FAILED: %s\n', P(ip).name, tag, e.message);
        end
    end
end
Tf = cell2table(rows, 'VariableNames', {'protocol','closure','feature','pde','mom','wall_mom','wall_pde','rms_pct','max_pct'});
Tf.rel_err_pct = 100*(Tf.mom - Tf.pde)./max(abs(Tf.pde), 1e-9);
save(fullfile(here, 'results', 'aposteriori.mat'), 'R', 'Tf', 'Fiso');
writetable(Tf, fullfile(here, 'results', 'aposteriori_features.csv'));
disp(Tf(strcmp(Tf.closure, 'gauss3'), {'protocol','feature','pde','mom','rel_err_pct'}));

