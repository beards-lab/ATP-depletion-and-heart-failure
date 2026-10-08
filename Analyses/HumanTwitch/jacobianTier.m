function Jr = jacobianTier(name)
%JACOBIANTIER  Forward-difference Jacobian, local 95 % intervals and correlations at a finished tier's optimum.
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');
F = load(fullfile(res, [name '_fit.mat']));  R = F.R;  cfg = F.cfg;
L = load(fullfile(here, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);
x = R.x;  r0 = residTier(x, cfg, tb, Fd);  J = zeros(numel(r0), numel(x));
for i = 1:numel(x)
    xi = x;  h = 2e-3*max(1, abs(x(i)));  xi(i) = xi(i) + h;  J(:, i) = (residTier(xi, cfg, tb, Fd) - r0)/h;
end
s2 = sum(r0.^2)/(numel(r0) - numel(x));  Cv = s2*pinv(J'*J);  sd = sqrt(diag(Cv))';  Rho = Cv./(sd'*sd);  sv = svd(J);
sem = sqrt(mean(D.sdFn(sel).^2))/sqrt(nnz(D.isCtrl));  ac1 = corr(r0(1:end-1), r0(2:end));
names = {cfg.par.name};
fprintf('%s: rms %.3f %% | SEM %.2f %% (ratio %.1f) | lag-1 autocorr %.2f | cond(J) %.0f\n', name, 100*sqrt(mean(r0.^2)), 100*sem, sqrt(mean(r0.^2))/sem, ac1, sv(1)/sv(end));
ci = cell(numel(x), 1);
for i = 1:numel(x)
    if cfg.par(i).log, ci{i} = sprintf('x/ %.3f', 10^(1.96*sd(i)));  fprintf('%-10s %9.4g  %s (95 %%)\n', names{i}, 10^x(i), ci{i});
    else, ci{i} = sprintf('+- %.3g', 1.96*sd(i));  fprintf('%-10s %9.4g  %s (95 %%)\n', names{i}, x(i), ci{i}); end
end
disp(array2table(round(Rho, 2), 'VariableNames', names, 'RowNames', names));
Jr = struct('J', J, 'Cv', Cv, 'Rho', Rho, 'sd', sd, 'r0', r0, 'sem', sem, 'ac1', ac1, 'names', {names}, 'ci', {ci}, 'cond', sv(1)/sv(end));
save(fullfile(res, [name '_jac.mat']), '-struct', 'Jr');
end
