% BuildCombined.m -- realign4 baseline + the low-ATP ATP-gate block, at full precision.
%
% "Update the baseline from the best fit in Further-optim, keep the low-ATP specifics."
% The low-ATP specifics are a clean separable block: UseAtpPowerGates + 21 pAtp_*
% exponents + the 2 mM condition itself. Everything else in the incumbent snapshot IS
% the old 8 mM baseline, so the transfer is: take realign4_opt.m verbatim, append the
% block. Exponents come from the optimiser state at FULL precision (%.17g) rather than
% from the committed snapshot's %.6g -- that rounding alone is worth ~0.60 of the 2 mM
% objective (Analyses/AtpGateHillForm/conclusions.md 3a).

here = fileparts(mfilename('fullpath')); root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

S = load(fullfile(root,'Analyses','LowATP_AbsoluteFit','results','powergate2_state.mat'));

base = fileread(fullfile(root,'params','realign4_opt.m'));
out  = fullfile(root,'params','lowATP_realign4_gates.m');

L = {};
L{end+1} = sprintf('\n%%%% ==== Low-ATP block, transferred verbatim from lowATP_powergate2_opt (fit 2026-08-10) ====');
L{end+1} = sprintf('%%%% Baseline above is params/realign4_opt.m (Further-optim 8 mM best, 2026-08-18 merge).');
L{end+1} = sprintf('%%%% Exponents are the optimiser''s own best.x at FULL precision (its recorded 2 mM cost');
L{end+1} = sprintf('%%%% on the OLD baseline/feature list was %.6f); NOT refitted on this baseline.', S.best.c);
L{end+1} = sprintf('%%%% Gates are exactly 1 at MgATP_ref, so the 8 mM fit above is untouched.');
L{end+1} = 'params0.UseAtpPowerGates = 1;';
L{end+1} = 'params0.MgATP_ref = 8;';
L{end+1} = 'params0.UseAtpDetachR2 = 0;';
L{end+1} = 'params0.UseAtpKmsrd = 0;';
for i = 1:numel(S.P)
    L{end+1} = sprintf('params0.pAtp_%s = %.17g;', S.P{i}, S.best.x(i)); %#ok<SAGROW>
end
L{end+1} = 'params0.MgATP = 2;';
L{end+1} = 'params0.MgADP = 0;';
L{end+1} = 'params0.Pi = 0.8;';
L{end+1} = 'params0.velocitytableonfile = ''protocol_03_27_2026_2mM_slack.mat'';';
L{end+1} = '';

fid = fopen(out,'w');
fwrite(fid, base);
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('written: %s\n', out);

% ---- the incumbent's feature list, for the record ------------------------------
fprintf('\nincumbent fit was scored on %d features:\n', numel(S.FN));
fprintf('  %s\n', strjoin(S.FN, ' | '));
fprintf('\nrecorded 8 mM floor on that list: %.6f\n', S.FLOOR);

% ---- sanity: both snapshots load and differ only where expected ----------------
a = loadParams('lowATP_powergate2_opt');
b = loadParams('lowATP_realign4_gates');
fa = fieldnames(a); fb = fieldnames(b);
fprintf('\nfields: incumbent %d, combined %d, only-in-combined: %s\n', ...
    numel(fa), numel(fb), strjoin(setdiff(fb,fa)', ', '));
nd = 0;
for i = 1:numel(fa)
    f = fa{i};
    if ~isfield(b,f), continue; end
    if isnumeric(a.(f)) && isnumeric(b.(f)) && isequal(size(a.(f)),size(b.(f))) ...
            && ~isequaln(a.(f), b.(f))
        nd = nd + 1;
    end
end
fprintf('numeric fields differing between incumbent and combined: %d\n', nd);
fprintf('\nBUILD DONE\n');
