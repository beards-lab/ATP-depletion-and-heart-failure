% RunSrxEnergetics.m -- net ATP turnover of the ka-compensated SRX variants.
% The objective's XTOR is out.RTD = kah*PT, the GROSS DRX hydrolysis flux. It misses
% reverse hydrolysis (RDT) and hydrolysis inside SRX (SR -> SRD), which grows when SRX
% recruitment grows. At steady state the net ATPase per head is the strong-detachment
% flux R2T (the only step that releases ADP here; D0 is off); check:
% RTD - RDT + (RSR2SRD - RSRD2SR) = R2T.
% Needs results/family.mat (RunSrxFamily). Output: results/energetics.txt
root = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
rsE = fullfile(root, 'Analyses', 'SRXDuringContraction', 'results');
Fm = load(fullfile(rsE, 'family.mat'));
p0 = getParams(loadParams('robust1_seed'), [], true, false);
clE = parcluster('local'); jslE = fullfile(tempdir, sprintf('mljob_srxe_%d', feature('getpid')));
if ~isfolder(jslE); mkdir(jslE); end
clE.JobStorageLocation = jslE;
if isempty(gcp('nocreate')); parpool(clE, 4); end
nV = numel(Fm.RF); E = zeros(nV, 8);
parfor iV = 1:nV
    r = Fm.RF(iV); q = srxApplyVariant(p0, struct('fam', r.fam, 'x', r.x));
    q.ka = r.ka; q.SRXT_0 = r.ic(1); q.SRXD_0 = r.ic(2); q.Velocity = [0; 0]; q.SL0 = 2.2;
    q = getParams(q, q.g, true, false);
    [~, o] = evaluateModel(str2func(q.modelFcn), [-3 0], q);
    E(iV, :) = [o.Force(end), o.SR(end) + o.SRD(end), o.p1_0(end) + o.p2_0(end), o.RTD(end), ...
                o.RDT(end), o.RSR2SRD(end) - o.RSRD2SR(end), o.R2T(end), o.PuATP(end)];
end
fid = fopen(fullfile(rsE, 'energetics.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Isometric SL 2.2, ka-compensated: net ATPase R2T vs the objective''s XTOR (gross RTD)\n\n');
out('%-10s %7s %7s %7s %8s %8s %8s %8s %8s %8s\n', 'variant', 'F', 'SRX', 'att', 'RTD', 'RDT', 'SRXhyd', 'R2T', 'check', 'R2T/ref');
for iV = 1:nV
    e = E(iV, :); chk = e(4) - e(5) + e(6);
    out('%-10s %7.2f %7.4f %7.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.3f\n', sprintf('%s x%g', Fm.RF(iV).fam, Fm.RF(iV).x), ...
        e(1:7), chk, e(7)/E(1, 7));
end
fclose(fid);
save(fullfile(rsE, 'energetics.mat'), 'E');
fprintf('DONE\n');
