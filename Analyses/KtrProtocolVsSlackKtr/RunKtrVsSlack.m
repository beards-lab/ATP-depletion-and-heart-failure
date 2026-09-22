% RunKtrVsSlack.m
% =========================================================================
% The slack protocol and the ktr protocol both end with force redeveloping
% isometrically at SL 2.0 um from ~zero.  In the DATA the two rates agree
% (Analyses/RestretchVsKtrRecovery/conclusions.md Part 1: ktr/post-slack = 1.04).
%
% But the CURRENT best 8 mM paramset (params/estiff125_opt.m) has
% params0.RunKtr = 0 -- the ktr protocol is not simulated and contributes
% nothing to the feature objective.  Every `ktr` feature in params0.fn comes
% from extractSlackAttributes, i.e. from the SLACK protocol.
%
% This script asks the held-out question: with the model fitted only to the
% slack redevelopment, what does it predict for the ktr protocol?
%
% Reproduce:
%   cd(root); addpath(genpath('.'));
%   Analyses/KtrProtocolVsSlackKtr/RunKtrVsSlack
% =========================================================================

root   = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

DAYS = {'03_27','04_03','04_10'};

%% ---------------- DATA: both windows, all three days ----------------
D = struct('day',{},'kS2',{},'kK',{},'kS',{},'AS2',{},'AK',{});
for d = 1:numel(DAYS)
    S = load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
    svt = S.velocitytable;  kvt = K.velocitytable;
    iS  = find(svt(:,2) < -1);
    % F_iso of the slack file: median of the 300 ms isometric run-up before slack 1
    FisoD = median(S.datatable(S.datatable(:,1) > svt(iS(1),1)-0.30 & ...
                               S.datatable(:,1) < svt(iS(1),1)-0.01, 3));
    ws = recoveryWindows(S.datatable(:,1), S.datatable(:,3), svt, 'slack', FisoD);
    wk = recoveryWindows(K.datatable(:,1), K.datatable(:,3), kvt, 'ktr',   1.0);
    ps = ws(strcmp({ws.type},'postSlack'));
    D(d).day = DAYS{d};
    D(d).kS  = [ps.k63];
    D(d).kS2 = ps(2).k63;          % cycle 2 == SL 2.0 um == ktr's length
    D(d).AS2 = ps(2).A;
    D(d).kK  = wk(1).k63;
    D(d).AK  = wk(1).A;
    D(d).wsS2 = ps(2); D(d).wkK = wk(1);
end

fprintf('\n===================== DATA, 8 mM =====================\n');
fprintf('%-8s %10s %10s %8s | %9s %9s\n','day','postSlack2','ktr','ratio','A_slack2','A_ktr');
for d=1:numel(D)
    fprintf('%-8s %10.1f %10.1f %8.2f | %9.3f %9.3f\n', ...
        D(d).day, D(d).kS2, D(d).kK, D(d).kK/D(d).kS2, D(d).AS2, D(d).AK);
end
fprintf('%-8s %10.1f %10.1f %8.2f\n','MEAN',mean([D.kS2]),mean([D.kK]),mean([D.kK])/mean([D.kS2]));

%% ---------------- MODEL: current best paramset ----------------
base = getParams(loadParams('estiff125_opt'), [], true, false);
base.PlotEachSeparately = 0; base.ghostSave = ''; base.ghostLoad = '';
base.EvalFeatures = 1;
p = getParams(base, base.g, true);

% --- slack protocol ---
[E_slack, os, fm_s, fd_s] = runSlackExperiment(p);
svt = load(fullfile(root,'data',p.velocitytableonfile)).velocitytable;
iS  = find(svt(:,2) < -1);
FisoM = interp1(os.t, os.Force, svt(iS(1),1)-0.02, 'linear','extrap');
wsM = recoveryWindows(os.t, os.Force, svt, 'slack', FisoM);
psM = wsM(strcmp({wsM.type},'postSlack'));

% --- ktr protocol ---
kvt = load(fullfile(root,'data',p.ktr_velocitytableonfile)).velocitytable;
pk  = p; pk.SL0 = 2*kvt(1,4); pk = getParams(pk, pk.g, true);
pk.Velocity = kvt(:,2);
[~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk);
FisoK = interp1(ok.t, ok.Force, kvt(2,1), 'linear','extrap');
wkM = recoveryWindows(ok.t, ok.Force, kvt, 'ktr', FisoK);

fprintf('\n============ MODEL estiff125_opt (ktr is HELD OUT) ============\n');
fprintf('F_iso  slack=%.2f kPa   ktr-protocol=%.2f kPa\n', FisoM, FisoK);
fprintf('%-16s %8s %8s %8s %8s %8s\n','window','k63','kfixA','kfree','A/Fiso','r2');
for c=1:numel(psM)
    fprintf('postSlack %-6d %8.1f %8.1f %8.1f %8.3f %8.3f\n', c, ...
        psM(c).k63, psM(c).kfixA, psM(c).kfit, psM(c).A, psM(c).r2);
end
fprintf('%-16s %8.1f %8.1f %8.1f %8.3f %8.3f\n','ktr', ...
    wkM(1).k63, wkM(1).kfixA, wkM(1).kfit, wkM(1).A, wkM(1).r2);
fprintf('\nMODEL ratio ktr/postSlack2 = %.2f     (DATA %.2f)\n', ...
    wkM(1).k63/psM(2).k63, mean([D.kK])/mean([D.kS2]));

fprintf('\n-- the objective''s own estimator (extractSlackAttributes) --\n');
fprintf('model ktr per slack : %s\n', mat2str(round(fm_s.ktr(:)',1)));
fprintf('data  ktr per slack : %s\n', mat2str(round(fd_s.ktr(:)',1)));

save(fullfile(resdir,'baseline.mat'),'D','psM','wkM','FisoM','FisoK','fm_s','fd_s','E_slack','-v7.3');
save(fullfile(resdir,'baseline_outs.mat'),'os','ok','svt','kvt','-v7.3');
