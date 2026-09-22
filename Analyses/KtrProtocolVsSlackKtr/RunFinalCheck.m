% RunFinalCheck.m -- regression + validation for params/fastgate_best.m.
%
% fastgate_best = estiff125_opt with tau_reg 26.6487 ms -> 0.2 ms, nothing else.
% The ktr protocol stays OUT of the objective (RunKtr = 0, no ktrP_* in fn) so it
% remains a held-out validation; validateKtrProtocol is how it is measured.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));

q = getParams(loadParams('fastgate_best'), [], true, false);
b = getParams(loadParams('estiff125_opt'), [], true, false);
for s = {'PlotEachSeparately','BreakOnODEUnstable','SaveBest'}
    q.(s{1}) = 0; b.(s{1}) = 0;
end
q.ghostSave=''; q.ghostLoad=''; b.ghostSave=''; b.ghostLoad='';

fprintf('\n===== 1. the objective is the SAME 24 features, ktr excluded =====\n');
assert(q.RunKtr == 0 && ~any(contains(q.fn,'ktrP_')), 'ktr must stay out of the objective');
assert(isequal(q.fn, b.fn), 'fn must be identical to the incumbent''s');
fprintf('RunKtr = %d, fn entries = %d, ktrP_* in fn = %d, fn identical to incumbent = %d\n', ...
        q.RunKtr, numel(q.fn), any(contains(q.fn,'ktrP_')), isequal(q.fn,b.fn));
[Etq,Eq] = evaluateBakersExp([], q, false, false);
[Etb,Eb] = evaluateBakersExp([], b, false, false);
fprintf('  estiff125_opt  Et = %.4f   raw E = %s\n', Etb, mat2str(round(Eb,4)));
fprintf('  fastgate_best  Et = %.4f   raw E = %s\n', Etq, mat2str(round(Eq,4)));
fprintf('  -> %+.4f (%.1f%% better)\n', Etq-Etb, 100*(1-Etq/Etb));

fprintf('\n===== 2. only tau_reg differs =====\n');
diff_fields = {};
for f = intersect(fieldnames(q), fieldnames(b))'
    if ~isequaln(q.(f{1}), b.(f{1})); diff_fields{end+1} = f{1}; end %#ok
end
fprintf('fields that differ: %s   (%.6g ms vs %.6g ms)\n', ...
        strjoin(diff_fields,', '), q.tau_reg*1e3, b.tau_reg*1e3);
assert(isequal(diff_fields, {'tau_reg'}), 'expected tau_reg to be the only change');

fprintf('\n===== 3. staircase, also held out (RunStairs = 0) =====\n');
NM = {'fastgate_best','estiff125_opt'};
for v = 1:2
    p = q; if v==2; p = b; end
    ps = p; ps.RunStairs = 1; ps = getParams(ps, ps.g, true);
    fprintf('  %-14s stairs E = %.4f\n', NM{v}, runStairsExperiment(ps));
end

fprintf('\n===== 4. the held-out ktr protocol =====');
validateKtrProtocol('estiff125_opt');
validateKtrProtocol('fastgate_best');
