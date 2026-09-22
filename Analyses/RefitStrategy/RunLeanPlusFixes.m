% RunLeanPlusFixes.m -- can the lean lineage carry the two ktr fixes?
%   fix 1  UseMaxwellTensionOnly = 1   (the titin-dashpot force floor that only
%          exists in the ktr protocol; the ablation deleted it as "+0.037 free",
%          but the ktr protocol was not in the objective that measured that)
%   fix 2  tau_reg -> 0.2 ms           (the fast registration-availability gate)
% Priced on the LEAN objective (floored, pooled), plus the held-out ktr.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');

b = getParams(loadParams('lean_opt'), [], true, false);
b.PlotEachSeparately=0; b.ghostSave=''; b.ghostLoad=''; b.BreakOnODEUnstable=0; b.SaveBest=0;

V = { 'L0  lean_opt as is',        @(p) p
      'L1  + MaxwellTensionOnly',  @(p) setfield(p,'UseMaxwellTensionOnly',1)
      'L2  + fast gate',           @(p) setfield(p,'tau_reg',2e-4)
      'L3  + both',                @(p) setfield(setfield(p,'UseMaxwellTensionOnly',1),'tau_reg',2e-4) };

fprintf('\n%-28s %10s %10s %12s %12s\n','variant','Et(lean)','dEt','ktr k63','rmse/noise');
R = struct('name',{},'Et',{},'fc',{},'k63',{},'rn',{});
for v=1:size(V,1)
    p = V{v,2}(b); p = getParams(p, p.g, true);
    [Et,~,fc,fn] = evaluateBakersExp([], p, false, false);
    % held-out ktr, inline (validateKtrProtocol loads from disk, this is in-memory)
    kvt = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_ktr.mat'));
    pk = p; pk.SL0 = 2*kvt.velocitytable(1,4); pk = getParams(pk,pk.g,true);
    pk.Velocity = kvt.velocitytable(:,2);
    [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt.velocitytable(:,1), pk);
    FisoK = interp1(ok.t, ok.Force, kvt.velocitytable(2,1),'linear','extrap');
    fm = extractKtrAttributes(ok.t, ok.Force/FisoK, kvt.velocitytable);
    Fi = interp1(ok.t, ok.Force/FisoK, kvt.datatable(:,1),'linear','extrap');
    m  = kvt.datatable(:,1) >= kvt.velocitytable(end-1,1);
    rmse = sqrt(mean((Fi(m)-kvt.datatable(m,3)).^2));
    R(v) = struct('name',V{v,1},'Et',Et,'fc',fc,'k63',fm.ktrP_k63,'rn',rmse/0.0222);
    fprintf('%-28s %10.4f %10s %12.1f %11.2fx\n', V{v,1}, Et, ...
        iif(v==1,'--',sprintf('%+.4f',Et-R(1).Et)), fm.ktrP_k63, R(v).rn);
end
fprintf('%-28s %10s %10s %12.1f\n','DATA 03/27','--','--',49.4);

fprintf('\n---- per-feature, L3 vs L0 ----\n');
d = R(4).fc - R(1).fc; [~,o] = sort(abs(d),'descend');
lab = cellfun(@(s) regexprep(s,'\|.*$',''), fn, 'UniformOutput',false);
for i = o(1:8); fprintf('%-24s %8.4f -> %8.4f  (%+.4f)\n', lab{i}, R(1).fc(i), R(4).fc(i), d(i)); end
save(fullfile(resdir,'leanplusfixes.mat'),'R','fn','lab');

function o = iif(c,a,b); if c; o=a; else; o=b; end; end
