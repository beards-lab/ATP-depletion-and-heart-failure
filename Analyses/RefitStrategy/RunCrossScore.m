% RunCrossScore.m -- two lineages, two objectives. Score every combination so
% the refit decision is made on numbers, not on which file was edited last.
%
%   A  params/fastgate_best.m  = estiff125_opt + fast gate. Objective: 24 terms,
%      NO floors, coolDownLS, rsK, FV weight 125, single-prep 03/27 target.
%   B  params/lean_opt.m       = the MechanismAblation output. Objective: floored
%      at 1 prep-width, rsK80, no coolDownLS, FV share 0.25, POOLED target,
%      and three mechanisms deleted (A2Reattaching, MaxwellTensionOnly,
%      DynamicPassive).
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

TAGS = {'fastgate_best','estiff125_opt','lean_opt'};
P = struct();
for i=1:numel(TAGS)
    q = getParams(loadParams(TAGS{i}), [], true, false);
    q.PlotEachSeparately=0; q.ghostSave=''; q.ghostLoad=''; q.BreakOnODEUnstable=0; q.SaveBest=0;
    P.(TAGS{i}) = q;
end
% the two objectives = fn + which target the slack block is scored against
OBJ = { 'A unfloored / 03-27',  P.fastgate_best.fn, P.fastgate_best.UsePooledTarget
        'B floored / pooled',   P.lean_opt.fn,      P.lean_opt.UsePooledTarget };

fprintf('\n================ CROSS-SCORING ================\n');
fprintf('%-16s %24s %24s\n','paramset',OBJ{1,1},OBJ{2,1});
M = nan(numel(TAGS), 2);
for i=1:numel(TAGS)
    for o=1:2
        p = P.(TAGS{i}); p.fn = OBJ{o,2}; p.UsePooledTarget = OBJ{o,3};
        p = getParams(p, p.g, true);
        try
            M(i,o) = evaluateBakersExp([], p, false, false);
        catch ME
            fprintf('  %s under %s FAILED: %s\n', TAGS{i}, OBJ{o,1}, ME.message);
        end
    end
    fprintf('%-16s %24.4f %24.4f\n', TAGS{i}, M(i,1), M(i,2));
end
fprintf('\n(columns are different objectives -- compare DOWN a column, never across)\n');

%% ---- the mechanism/flag differences that separate the two lineages ----
fprintf('\n================ WHAT DIFFERS ================\n');
a = P.fastgate_best; b = P.lean_opt;
F = {'tau_reg','UseMaxwellTensionOnly','UseA2Reattaching','UseDynamicPassive', ...
     'UsePooledTarget','estiff','A0','v_ref_reg','kstiff2','ka','kd','k1','k2','kah'};
fprintf('%-24s %14s %14s\n','field','fastgate_best','lean_opt');
for f = F
    va = a.(f{1}); vb = b.(f{1});
    fprintf('%-24s %14.6g %14.6g%s\n', f{1}, va, vb, repmat('   <-- differs', 1, ~isequaln(va,vb)));
end

%% ---- held-out ktr for all three ----
fprintf('\n================ HELD-OUT ktr PROTOCOL ================\n');
K = nan(numel(TAGS),3);
for i=1:numel(TAGS)
    V = validateKtrProtocol(TAGS{i}, false);
    K(i,:) = [mean([V.k63]), mean([V.Ektr]), mean([V.rmse]./[V.noise])];
    fprintf('%-16s k63 %6.1f (data 43.2)   traceMSE %6.4f   %5.2fx the data noise floor\n', ...
            TAGS{i}, K(i,1), K(i,2), K(i,3));
end
save(fullfile(resdir,'crossscore.mat'),'M','K','TAGS','OBJ');
