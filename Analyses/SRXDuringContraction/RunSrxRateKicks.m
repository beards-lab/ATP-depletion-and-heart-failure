% RunSrxRateKicks.m -- what does each SRX rate do? Kick it +1000 % (x11) and watch the
% features. UNCOMPENSATED: force is allowed to move, so this is the raw effect of the rate,
% not an identifiability test (that is RunSrxFamily, which re-solves ka).
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxRateKicks.m')"
%
% Rates of the SRX block (dPUdT_CombinedTransitions):
%   in  : ksr0    PT  -> P_SR   (x e^(-F/sigma2))      DRX.ATP   -> SRX.ATP
%         ksrd    PD  -> P_SRD  (x e^(-F/sigma_srd2))  DRX.ADPPi -> SRX.ADPPi
%   out : kmsrd   P_SRD -> PD   (x e^(+F/sigma_srd1))  SRX.ADPPi -> DRX.ADPPi
%         (kmsr, P_SR -> PT, is 0 in the incumbent: nothing to scale)
%   internal: ksr2srd / ksrd2sr  hydrolysis inside SRX
% Each kick: the full 8 mM battery (slack-restretch x5, ktr protocol, F-V, passive),
% the isometric state at SL 2.2 (force, SRX, attached, net ATPase R2T) and one traced
% slack-restretch cycle (cycle 1 of the 03/27 protocol) and one traced ktr protocol.
% Output -> results/kick_run_<k>.mat, results/kicks.mat (SummarizeSrxRateKicks.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsK = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
diary(fullfile(rsK, 'kicks.log')); diary on;
fprintf('SRX rate kicks on robust1_seed, started %s\n', datestr(now));
clK = parcluster('local');
jslK = fullfile(tempdir, sprintf('mljob_srxk_%d', feature('getpid')));
if ~isfolder(jslK); mkdir(jslK); end
clK.JobStorageLocation = jslK;
if isempty(gcp('nocreate')); parpool(clK, 5); end

pK = getParams(loadParams('robust1_seed'), [], true, false);
pK.RunForceVelocity=true; pK.RunSlack=true; pK.RunStairs=false; pK.RunKtr=true;
pK.RunForceVelocityTime=false; pK.RunSlackPassive=true; pK.EvalFeatures=true;
pK.BreakOnODEUnstable=false; pK.PlotEachSeparately=0; pK.PlotFeatureFitting=0;
pK.RunSlackSegments='AllPar'; pK.MaxRunTime=400; pK.OptimizeOn='Feats';
pK.ghostSave=''; pK.ghostLoad=''; pK.SaveBest=0; pK.mods={}; pK.g=ones(1,0);

K = 11;   % +1000 %
KV = struct( ...
  'name',  {'ref', 'in x11', 'ksr0 x11', 'ksrd x11', 'out x11', 'in+out x11', 'ksr2srd x11', 'ksrd2sr x11', 'in /11', 'out /11', ...
            'in x3', 'in /3', 'out x3', 'out /3'}, ...
  'rates', {{},    {'ksr0','ksrd'}, {'ksr0'}, {'ksrd'}, {'kmsrd'}, {'ksr0','ksrd','kmsrd'}, {'ksr2srd'}, {'ksrd2sr'}, {'ksr0','ksrd'}, {'kmsrd'}, ...
            {'ksr0','ksrd'}, {'ksr0','ksrd'}, {'kmsrd'}, {'kmsrd'}}, ...
  'x',     {1,     K, K, K, K, K, K, K, 1/K, 1/K, 3, 1/3, 3, 1/3});
% The analysis reads SRX as ONE pool: SRX+ = inflow (ksr0 + ksrd together: PT -> SRX.ATP and
% PD -> SRX.ADP.Pi), SRX- = outflow (kmsrd: SRX.ADP.Pi -> PD), each as a dose series
% /11, /3, x3, x11. The single-rate and SR<->SRD kicks (3,4,6-8) are kept as side runs.
nV = numel(KV);

%% ---- isometric state + traced slack cycle 1 and ktr protocol per kick (parallel) --
dP = load(fullfile(root, 'data', 'protocol_03_27_2026_8mM_slack.mat'));
vt0 = dP.velocitytable; iRel = find(vt0(:, 2) < -1, 2);        % first two releases
tR = vt0(iRel(1), 1);
vtC = [-2, 0; vt0(iRel(1):iRel(2)-1, 1:2) - [tR 0]; vt0(iRel(2), 1) - tR - 0.05, 0];
vtC(end-1:end, :) = [vtC(end-1, 1), 0; vtC(end, 1), 0];         % end on an isometric hold
mD = dP.datatable(:, 1) >= tR - 0.1 & dP.datatable(:, 1) <= vt0(iRel(2), 1) - 0.05;
dataC = [dP.datatable(mD, 1) - tR, dP.datatable(mD, 2:3)];
dK = load(fullfile(root, 'data', 'protocol_03_27_2026_velocitytable_ktr.mat'));
vtK = dK.velocitytable(:, 1:2); tK = vtK(2, 1); vtK(:, 1) = vtK(:, 1) - tK;   % release at t = 0
dataK = [dK.datatable(:, 1) - tK, dK.datatable(:, 2:3)];                       % force normalised
ISO = zeros(nV, 6); TRC = cell(1, nV); TRK = cell(1, nV);
parfor iV = 1:nV
    q = srxScaleRates(pK, KV(iV).rates, KV(iV).x);
    [s, ~] = srxIsoState(q);
    ISO(iV, :) = [s.F, s.SRX, s.att, s.R2T, s.RTD, s.SR];
    qc = q; qc.SRXT_0 = s.SR; qc.SRXD_0 = s.SRD; qc.SL0 = 2 * vt0(1, 4);
    qc.Velocity = vtC(:, 2); qc = getParams(qc, qc.g, true, false);
    [~, o] = evaluateModel(str2func(qc.modelFcn), vtC(:, 1), qc);
    k = makeMonotonous(o.t(:)); tq = (-0.05:5e-4:vtC(end, 1))';
    TRC{iV} = [tq, interp1(o.t(k), o.Force(k), tq), interp1(o.t(k), o.SR(k) + o.SRD(k), tq), ...
               interp1(o.t(k), o.SL(k), tq), interp1(o.t(k), o.p1_0(k) + o.p2_0(k), tq)];
    qk = q; qk.SRXT_0 = s.SR; qk.SRXD_0 = s.SRD; qk.SL0 = 2 * dK.velocitytable(1, 4);
    qk.Velocity = vtK(:, 2); qk = getParams(qk, qk.g, true, false);
    [~, o] = evaluateModel(str2func(qk.modelFcn), vtK(:, 1), qk);
    k = makeMonotonous(o.t(:)); tq = (-0.01:2.5e-4:vtK(end, 1))';
    TRK{iV} = [tq, interp1(o.t(k), o.Force(k), tq), interp1(o.t(k), o.SR(k) + o.SRD(k), tq), ...
               interp1(o.t(k), o.SL(k), tq), interp1(o.t(k), o.p1_0(k) + o.p2_0(k), tq)];
end
for iV = 1:nV
    fprintf('%-12s F_iso %6.2f  SRX %.4f  att %.4f  R2T %.2f\n', KV(iV).name, ISO(iV, 1:4));
end

%% ---- full battery per kick --------------------------------------------------------
RK = struct('name',{},'rates',{},'x',{},'iso',{},'trace',{},'traceK',{},'fm',{},'fd',{},'cost',{},'sec',{});
for iV = 1:nV
    fR = fullfile(rsK, sprintf('kick_run_%d.mat', iV));
    if isfile(fR)                                   % resume: keep the battery, refresh the traces
        L = load(fR);
        if strcmp(L.R1.name, KV(iV).name)
            R1 = L.R1; R1.iso = ISO(iV, :); R1.trace = TRC{iV}; R1.traceK = TRK{iV};
            RK(iV) = orderfields(R1, RK); save(fR, 'R1'); continue
        end
    end
    if isempty(gcp('nocreate')); parpool(clK, 5); end
    tV = tic;
    q = srxScaleRates(pK, KV(iV).rates, KV(iV).x);
    q.SRXT_0 = ISO(iV, 6); q.SRXD_0 = ISO(iV, 2) - ISO(iV, 6);
    try
        [fm, fd] = iRun(q); c = evalFeatureCost(fd, fm, pK.fn);
    catch eV
        warning('kick %d failed: %s', iV, eV.message); fm = struct(); fd = struct(); c = NaN;
    end
    RK(iV).name = KV(iV).name; RK(iV).rates = KV(iV).rates; RK(iV).x = KV(iV).x;
    RK(iV).iso = ISO(iV, :); RK(iV).trace = TRC{iV}; RK(iV).traceK = TRK{iV}; RK(iV).fm = fm; RK(iV).fd = fd;
    RK(iV).cost = c; RK(iV).sec = toc(tV);
    fprintf('%-12s cost %.4f  (%.0f s)\n', KV(iV).name, sum(c), RK(iV).sec);
    R1 = RK(iV); save(fR, 'R1'); %#ok<NASGU>
end
save(fullfile(rsK, 'kicks.mat'), 'RK', 'KV', 'dataC', 'vtC', 'dataK', 'vtK');
fprintf('DONE %s\n', datestr(now)); diary off;

function [fm, fd] = iRun(p)
    % RunBakersExp is a SCRIPT; a function workspace keeps it from clobbering ours.
    params0 = p; %#ok<NASGU>
    LoadData; RunBakersExp; %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end
