% BuildFastGateBest.m -- promote the fast gate into the 8 mM best paramset.
%
% params/fastgate_best.m = params/estiff125_opt.m with ONE change:
%     tau_reg  26.6487 ms  ->  0.2 ms
%
% The ktr protocol is deliberately NOT scored: RunKtr stays 0 and ktrP_k63 stays
% out of fn, so it remains a clean held-out validation. The feature extractor
% (Model/extractKtrAttributes.m) is still there -- it is how the validation is
% MEASURED, and it costs nothing while the feature is absent from fn.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));

base = getParams(loadParams('estiff125_opt'), [], true, false);
assert(base.RunKtr == 0, 'expected the incumbent to have RunKtr = 0');
assert(~any(contains(base.fn, 'ktrP_')), 'expected no ktr-protocol feature in fn');

p = base;
p.tau_reg = 2e-4;
p.mods    = {};                    % refit drivers set their own tunable list
writeParamsToMFile(fullfile(root,'params','fastgate_best.m'), p, [], ...
    ['8 mM best. estiff125_opt with tau_reg 26.6487 ms -> 0.2 ms (fast registration-' ...
     'availability gate). ktr protocol deliberately NOT scored (RunKtr = 0, no ktrP_* in fn) ' ...
     '- it is the held-out validation. See Analyses/KtrProtocolVsSlackKtr/conclusions.md']);

%% ---- verify the round trip against the incumbent, same objective ----
q = getParams(loadParams('fastgate_best'), [], true, false);
q.PlotEachSeparately=0; q.ghostSave=''; q.ghostLoad=''; q.BreakOnODEUnstable=0;
b = base; b.PlotEachSeparately=0; b.ghostSave=''; b.ghostLoad=''; b.BreakOnODEUnstable=0;

fprintf('\n=========== fastgate_best round trip ===========\n');
fprintf('tau_reg      %.6g ms  (incumbent %.6g ms)\n', q.tau_reg*1e3, b.tau_reg*1e3);
fprintf('RunKtr       %d        fn entries %d   (ktrP_* present: %d)\n', ...
        q.RunKtr, numel(q.fn), any(contains(q.fn,'ktrP_')));
d = setdiff(fieldnames(q), fieldnames(b));
fprintf('new fields vs incumbent: %s\n', strjoin([{'(none)'} d'],' '));
diff_fields = {};
for f = intersect(fieldnames(q), fieldnames(b))'
    a1=q.(f{1}); a2=b.(f{1});
    if isnumeric(a1)&&isnumeric(a2)&&isequal(size(a1),size(a2))
        if ~isequaln(a1,a2); diff_fields{end+1}=f{1}; end %#ok
    elseif ~isequaln(a1,a2); diff_fields{end+1}=f{1}; end %#ok
end
fprintf('fields that differ: %s\n', strjoin(diff_fields, ', '));

[Etq,~,fcq,fn] = evaluateBakersExp([], q, false, false);
[Etb,~,fcb,~ ] = evaluateBakersExp([], b, false, false);
fprintf('\nOBJECTIVE (unchanged 24-feature set, ktr protocol excluded)\n');
fprintf('  estiff125_opt  Et = %.4f\n', Etb);
fprintf('  fastgate_best  Et = %.4f   (%+.4f, %.1f%% better)\n', Etq, Etq-Etb, 100*(1-Etq/Etb));
dd = fcq-fcb; [~,o]=sort(abs(dd),'descend');
fprintf('  movers: ');
for i=o(1:5); fprintf('%s %+.4f  ', regexprep(fn{i},'\|.*$',''), dd(i)); end
fprintf('\n');
