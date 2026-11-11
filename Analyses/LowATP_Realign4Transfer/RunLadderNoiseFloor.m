% RunLadderNoiseFloor.m -- how big is the noise floor on the NEW 2 mM objective?
%
% The 0.25 accept threshold used throughout was measured on the slack-only 2 mM
% objective (sigma ~ 0.076, cost ~10). This objective now also carries the 2 mM
% force-velocity and power terms and sits at ~14-20, so that threshold is not
% transferable and the small marginal gains in the pathway ladder cannot be
% called real until this is measured.
%
% Method (as Analyses/AtpGateHillForm/RunCostNoiseAmplitude.m): nudge ONE gate
% exponent along a ladder of deltas far too small to mean anything physically,
% and look at the spread. If the cost moves non-monotonically with delta, it is
% scatter, not sensitivity.
%
% Run: matlab -batch "run('Analyses/LowATP_Realign4Transfer/RunLadderNoiseFloor.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
if isempty(gcp('nocreate')); parpool('local',5); end
resdir = fullfile(here,'results');

L  = load(fullfile(resdir,'atp_ladder.mat'));
L2 = load(fullfile(resdir,'atp_ladder2.mat'));
fn2 = L.fn2;

% start from the refitted H2 (the physiologically defensible rung)
iH2 = find(contains({L2.OUT.name},'H2'), 1);
p = L2.OUT(iH2).p;

D = [0 1e-9 -1e-9 1e-8 -1e-8 1e-7 -1e-7 1e-6 -1e-6 1e-5 -1e-5];
c = nan(size(D));
fprintf('\n%12s %14s\n','delta','2 mM cost');
for i = 1:numel(D)
    q = p; q.pAtp_kmsrd = p.pAtp_kmsrd + D(i);
    c(i) = iCost(q, fn2, ATP_c, Data_ATP);
    fprintf('%12.0e %14.6f\n', D(i), c(i));
end

fprintf('\nmean          %.4f\n', mean(c));
fprintf('sigma         %.4f\n', std(c));
fprintf('peak-to-peak  %.4f\n', max(c)-min(c));
fprintf('\nmonotone in delta? (a real sensitivity would be): %s\n', ...
        string(issorted(c(2:end)) || issorted(flip(c(2:end)))));

fprintf('\n--- the ladder''s marginal gains against this floor ---\n');
g = struct('n',{'H1->H2 (4 SRX DOF)','H2->H3 (3 mech DOF)','H3->H4 (11 other DOF)'}, ...
           'v',{L2.OUT(1).c-L2.OUT(2).c, L2.OUT(2).c-L2.OUT(3).c, L2.OUT(3).c-L2.OUT(4).c});
for i = 1:numel(g)
    fprintf('  %-24s %+7.3f  = %5.1f x peak-to-peak  %s\n', g(i).n, -g(i).v, ...
            abs(g(i).v)/max(max(c)-min(c),eps), ...
            string(abs(g(i).v) > (max(c)-min(c))));
end

save(fullfile(resdir,'ladder_noise.mat'), 'D','c');
fprintf('\nNOISE FLOOR DONE\n');

function c = iCost(q, fn2, ATP_c, Data_ATP)
[~,~,fm,fd] = runSlackExperiment(q);
qf = q; qf.RunSlack = false;
[~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
f = fieldnames(fmv); for i=1:numel(f); fm.(f{i}) = fmv.(f{i}); end
f = fieldnames(fdv); for i=1:numel(f); fd.(f{i}) = fdv.(f{i}); end
c = sum(evalFeatureCost(fd, fm, fn2, 2));
end
