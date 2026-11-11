% PDEGridConvergence.m
% Is the moments error bigger than the PDE's OWN discretisation error?
% Compares features of the dS = 2 nm reference (the production grid), the
% dS = 1 nm rerun, and the Gaussian moments model, protocol by protocol.

here = fileparts(mfilename('fullpath'));
params0 = momentsSetup('fastgate_best');
P = momentsProtocols();
A = load(fullfile(here, 'results', 'aposteriori.mat'));
names = {'iso','ktr','slack','fv2','fv4'};
rows = {};
for ip = find(ismember({P.name}, names))
    S2 = load(fullfile(here, 'results', ['ref_' P(ip).name '.mat']));
    S1 = load(fullfile(here, 'results', ['refdS1_' P(ip).name '.mat']));
    [t2, i2] = unique(S2.ref.t); [t1, i1] = unique(S1.ref.t);
    f2 = protocolFeatures(P(ip), t2, S2.ref.Force(i2));
    f1 = protocolFeatures(P(ip), t1, S1.ref.Force(i1));
    fm = A.R.(matlab.lang.makeValidName(P(ip).name)).gauss3.fm;
    fn = fieldnames(f2);
    for k = 1:numel(fn)
        rows(end+1,:) = {P(ip).name, fn{k}, f1.(fn{k}), f2.(fn{k}), fm.(fn{k}), ...
            100*(f2.(fn{k}) - f1.(fn{k}))/max(abs(f1.(fn{k})),1e-9), ...
            100*(fm.(fn{k}) - f1.(fn{k}))/max(abs(f1.(fn{k})),1e-9)}; %#ok<SAGROW>
    end
end
Tg = cell2table(rows, 'VariableNames', {'protocol','feature','pde_dS1','pde_dS2','gauss3','dS2_vs_dS1_pct','gauss3_vs_dS1_pct'});
disp(Tg); writetable(Tg, fullfile(here, 'results', 'grid_convergence.csv'));
