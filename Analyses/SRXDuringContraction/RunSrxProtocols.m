% RunSrxProtocols.m -- which experiment would tell a low-SRX from a high-SRX fibre?
%
%   matlab -batch "run('Analyses/SRXDuringContraction/RunSrxProtocols.m')"
%
% Variants (from RunSrxFamily.m, each ka-compensated to the incumbent's isometric force):
%   ref   robust1_seed                       SRX_iso ~0.03, rest ~0.27
%   L10   recruitment x10                    SRX_iso ~0.15, rest ~0.70
%   S3    load sensitivity / 3               SRX_iso ~0.2,  rest ~0.33
%   M1e-3 entry/escape x1e-3                 SRX_iso ~0.008 (SR empties), frozen
%   MA1e-3 every SRX rate x1e-3              SRX_iso ~0.03, frozen (mant-like speed)
% Candidate protocols, all at 8 mM, run directly through evaluateModel:
%   A  activation step from rest (caged-Ca / rapid solution-switch analogue):
%      relaxed steady state (ka = 0), then ka on at t = 0 -> force rise t50, t90
%   B  release-hold series: release 0.15 ML, hold T, restretch in 2 ms, ktr after
%      restretch vs T in {5 10 20 40 80 160 320} ms
%   C  paired release-restretch: two identical ktr manoeuvres (T = 20 ms) separated
%      by D in {30 60 120 250 500} ms; second/first k63
%   D  isometric SL series 1.9-2.3 um: force and SRX (length-dependent activation)
% Output -> results/protocols.mat, results/protocols.txt (SummarizeSrxProtocols.m)
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsP = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
diary(fullfile(rsP, 'protocols.log')); diary on;
fprintf('SRX protocols, started %s\n', datestr(now));
clP = parcluster('local');
jslP = fullfile(tempdir, sprintf('mljob_srxp_%d', feature('getpid')));
if ~isfolder(jslP); mkdir(jslP); end
clP.JobStorageLocation = jslP;
if isempty(gcp('nocreate')); parpool(clP, 2); end

p0 = getParams(loadParams('robust1_seed'), [], true, false);
VP = struct('name', {'ref', 'L10', 'S3', 'M1e-3', 'MA1e-3'}, 'fam', {'ref', 'L', 'S', 'M', 'MA'}, 'x', {1, 10, 3, 1e-3, 1e-3});
nV = numel(VP);

%% ---- ka compensation (same recipe as RunSrxFamily) ------------------------------
[Ft, ssr] = srxIsoForce(p0, p0.ka, [], 3);
kaV = zeros(1, nV); icV = zeros(nV, 2);
parfor iV = 1:nV
    q = srxApplyVariant(p0, VP(iV)); ic = [];
    if startsWith(VP(iV).fam, 'M'); ic = ssr; end
    ka = q.ka; [F1, s1] = srxIsoForce(q, ka, ic, 3); ic = s1;
    ka2 = ka * Ft / F1; [F2, s2] = srxIsoForce(q, ka2, ic, 3); ic = s2;
    for it = 1:4
        if abs(F2 - Ft) < 0.05; break; end
        ka3 = ka2 + (Ft - F2) * (ka2 - ka) / (F2 - F1);
        ka = ka2; F1 = F2; ka2 = max(ka3, 1); [F2, s2] = srxIsoForce(q, ka2, ic, 3); ic = s2;
    end
    kaV(iV) = ka2; icV(iV, :) = s2;
end
for iV = 1:nV; fprintf('%-6s ka x%.3f  SRX_iso %.4f\n', VP(iV).name, kaV(iV)/p0.ka, sum(icV(iV, :))); end

%% ---- job list -------------------------------------------------------------------
TB = [5 10 20 40 80 160 320] * 1e-3;  DC = [30 60 120 250 500] * 1e-3;  SLs = [1.9 2.0 2.1 2.2 2.3];
J = struct('v', {}, 'kind', {}, 'x', {});
for iV = 1:nV
    J(end+1) = struct('v', iV, 'kind', 'A', 'x', 0); %#ok<SAGROW>
    for T = TB;  J(end+1) = struct('v', iV, 'kind', 'B', 'x', T);  end %#ok<SAGROW>
    for D = DC;  J(end+1) = struct('v', iV, 'kind', 'C', 'x', D);  end %#ok<SAGROW>
    for L = SLs; J(end+1) = struct('v', iV, 'kind', 'D', 'x', L);  end %#ok<SAGROW>
end
nJ = numel(J); RES = cell(1, nJ);
parfor iJ = 1:nJ
    q = srxApplyVariant(p0, VP(J(iJ).v)); q.ka = kaV(J(iJ).v);
    q.SRXT_0 = icV(J(iJ).v, 1); q.SRXD_0 = icV(J(iJ).v, 2);
    try
        RES{iJ} = srxProtocol(q, J(iJ).kind, J(iJ).x);
    catch e
        RES{iJ} = struct('err', e.message);
    end
end
save(fullfile(rsP, 'protocols.mat'), 'VP', 'kaV', 'icV', 'J', 'RES', 'TB', 'DC', 'SLs', 'Ft');
fprintf('DONE %s\n', datestr(now)); diary off;
