% RunSrxDrug.m -- does a drug experiment (mavacamten-like) identify the ACTIVATED SRX pool?
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxDrug.m')"
%
% Idea: the control data cannot tell fibres with different activated pools apart (sec. 8,
% twin). But a drug that acts on the SRX is measured WITHIN a preparation (control vs drug).
% If fibres that are identical at control respond DIFFERENTLY to the same drug, then
% control + drug data do pin the activated pool (given the drug's mechanism).
% Fibres (all 13-knob parameter sets, same knob list as RunSrxCompensation / RunSrxTwin):
%   twin truth (pool x1), twin battery-refits at x0.8 and x1.2   (identical twin control data)
%   real-data refits sec. 8 set B at x0.5 and x2, and the x1 control refit
% Drug models (on top of each fibre):
%   D1 recruitment   SRX inflow  x3   (ksr0 & ksrd)       "stabilises SRX from DRX"
%   D2 retention     SRX outflow /3   (kmsrd)             "holds heads in SRX, also under load"
% Read-outs: isometric force, activated SRX, battery features, protocol B k63(T)/F20(T).
% The drug EFFECT (drug - control of the same fibre) is compared across fibres in noise
% units (within-prep difference SD = sqrt(2) x one-prep SD).
% Output -> results/drug.mat (PlotSrxDrug.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsD = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
cdirC = fullfile(rsD, 'comp'); cdirT = fullfile(rsD, 'twin'); cdirD = fullfile(rsD, 'drug');
if ~isfolder(cdirD); mkdir(cdirD); end
diary(fullfile(rsD, 'drug.log')); diary on;
fprintf('SRX drug test, started %s\n', datestr(now));
clD = parcluster('local'); jslD = fullfile(tempdir, sprintf('mljob_srxd_%d', feature('getpid')));
if ~isfolder(jslD); mkdir(jslD); end
clD.JobStorageLocation = jslD;
if isempty(gcp('nocreate')); parpool(clD, 5); end

pD = getParams(loadParams('robust1_seed'), [], true, false);
pD.RunForceVelocity=true; pD.RunSlack=true; pD.RunStairs=false; pD.RunKtr=true;
pD.RunForceVelocityTime=false; pD.RunSlackPassive=true; pD.EvalFeatures=true;
pD.BreakOnODEUnstable=false; pD.PlotEachSeparately=0; pD.PlotFeatureFitting=0;
pD.RunSlackSegments='AllPar'; pD.MaxRunTime=400; pD.OptimizeOn='Feats';
pD.ghostSave=''; pD.ghostLoad=''; pD.SaveBest=0; pD.mods={}; pD.g=ones(1,0);

W = load(fullfile(rsD, 'twin.mat')); C = load(fullfile(rsD, 'compensation.mat'));
KN = W.KN; nK = size(KN, 1); iIn = find(strcmp(KN(:, 1), 'inflow')); iOut = find(strcmp(KN(:, 1), 'kmsrd'));
TB = W.TB;
pr = @(t, f) W.PR([W.PR.target] == t & strcmp({W.PR.fit}, f)).theta;
cr = @(t, s) C.R([C.R.target] == t & strcmp({C.R.set}, s)).theta;
FB = struct('name', {'twin truth x1', 'twin refit x0.8', 'twin refit x1.2', 'real refit x0.5', 'real refit x1 (control)', 'real refit x2'}, ...
            'theta', {zeros(nK, 1), pr(0.8, 'bat'), pr(1.2, 'bat'), cr(0.5, 'B'), cr(1, 'B'), cr(2, 'B')}, ...
            'group', {'twin', 'twin', 'twin', 'real', 'real', 'real'});
DR = struct('name', {'control', 'D1 inflow x3', 'D2 outflow /3'}, 'k', {[], iIn, iOut}, 'f', {1, 3, 1/3});

R = struct('fibre', {}, 'drug', {}, 'E', {});
for i = 1:numel(FB)
    for d = 1:numel(DR)
        th = FB(i).theta; if ~isempty(DR(d).k); th(DR(d).k) = th(DR(d).k) + log(DR(d).f); end
        tag = sprintf('drug_f%d_d%d', i, d);
        E = iEval(pD, KN, th, tag, cdirC, cdirT, cdirD, TB);
        R(end+1) = struct('fibre', i, 'drug', d, 'E', E); %#ok<SAGROW>
        fprintf('%-24s %-14s SRX %.4f  F %.2f  k63(5ms) %.1f  ktrfit %.1f\n', FB(i).name, DR(d).name, E.srx, E.F, ...
                E.P.k63(1), mean(E.fm.ktrfit));
    end
end
save(fullfile(rsD, 'drug.mat'), 'FB', 'DR', 'R', 'KN', 'TB');
fprintf('DONE %s\n', datestr(now)); diary off;

function E = iEval(p, KN, th, tag, cdirC, cdirT, cdirD, TB)
    % reuse any cached run with this exact theta (compensation, twin or drug caches)
    E = [];
    for dd = {cdirC, cdirT, cdirD}
        L = dir(fullfile(dd{1}, '*.mat'));
        for k = 1:numel(L)
            if endsWith(L(k).name, '_proto.mat'); continue; end
            S = load(fullfile(dd{1}, L(k).name));
            if isfield(S, 'E') && isfield(S.E, 'theta') && isequal(size(S.E.theta), size(th)) && max(abs(S.E.theta - th)) < 1e-12
                E = S.E; fP = fullfile(cdirT, [erase(L(k).name, '.mat') '_proto.mat']);
                if isfile(fP); P = load(fP); E.P = P.P; end
                break
            end
        end
        if ~isempty(E); break; end
    end
    q = srxApplyTheta(p, KN, th);
    if isempty(E)
        t0 = tic; s = srxIsoState(q); q2 = q; q2.SRXT_0 = s.SR; q2.SRXD_0 = s.SRD;
        [fm, fd] = iRun(q2);
        E = struct('tag', tag, 'theta', th, 'fm', fm, 'fd', fd, 'srx', s.SRX, 'F', s.F, 'sec', toc(t0));
        save(fullfile(cdirD, [tag '.mat']), 'E');
    end
    if ~isfield(E, 'P')
        fP = fullfile(cdirD, [tag '_proto.mat']);
        if isfile(fP); P = load(fP); E.P = P.P; else; E.P = srxProtoFeatures(q, TB); P = E.P; save(fP, 'P'); end
    end
    if ~isfield(E, 'srx') || ~isfield(E, 'F'); s = srxIsoState(q); E.srx = s.SRX; E.F = s.F; end
end

function [fm, fd] = iRun(p)
    params0 = p; %#ok<NASGU>
    LoadData; RunBakersExp; %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end
