% RunFvidProfile.m -- stage 2: do the "badly identified by the Beard records"
% verdicts survive a real, finite parameter change?
%
% AnalyzeFvid picks the three parameters the slack battery identifies and the
% Beard2022 records constrain worst, and for each writes two points: the parameter
% moved x2 and x1/2, every other parameter compensating as the Beard-record Fisher
% information says it can (results/profile_plan.mat). Each point is simulated here
% on the three Beard2022 records and on the slack battery. If the linear verdict
% is right, the Beard-record features barely move (a few chi2 units) while the
% slack-battery features move a lot.
%
%   matlab -batch "NW=6; run('Analyses/FV_Identifiability/RunFvidProfile.m')"
% Output -> results/profile/ (checkpoints), results/profile.mat, table below.

here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
cd(root); addpath(genpath(root));
wt = strsplit(genpath(fullfile(root, '.claude')), pathsep); wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
if ~exist('NW', 'var'), NW = 6; end
[pBat, pBeard, ~, recs] = fvidBase();
% the plan and the incumbent's vectors: from AnalyzeFvid.m, or from its Python twin
if isfile(fullfile(here, 'results', 'fvid.mat'))
    L = load(fullfile(here, 'results', 'profile_plan.mat')); A = load(fullfile(here, 'results', 'fvid.mat'));
else
    A = load(fullfile(here, 'results', 'fvid_py.mat')); L = A;
end
plan = L.plan; nP = numel(plan);
outP = fullfile(here, 'results', 'profile'); if ~isfolder(outP), mkdir(outP); end

% tasks: (point, record 1..3) on the Beard protocol, (point, 4) = the slack battery
[ip, ir] = ndgrid(1:nP, 1:4); T = [ip(:), ir(:)];
files = arrayfun(@(i) fullfile(outP, sprintf('p%02d_r%d.mat', T(i, 1), T(i, 2))), (1:size(T, 1))', 'uni', 0);
todo = find(~cellfun(@isfile, files));
fprintf('profile: %d points, %d tasks, %d to run\n', nP, size(T, 1), numel(todo));
if ~isempty(todo)
    pl = gcp('nocreate'); if isempty(pl) || pl.NumWorkers ~= NW, delete(pl); parpool('local', NW); end
    Tt = T(todo, :); Ft = files(todo);
    parfor ii = 1:numel(todo)
        e = plan(Tt(ii, 1)); r = struct('ok', false, 'msg', '', 'F', [], 'fm', [], 'fd', []); %#ok<PFBNS>
        try
            if Tt(ii, 2) <= 3
                r.F = simBeard(pBeard, e, recs(Tt(ii, 2))); %#ok<PFBNS>
            else
                p = pBat; k = ~strcmp(e.name, 'pAtp_k2');   % the battery is 8 mM only
                p.mods = e.name(k); p.g = e.f(k);
                [r.fm, r.fd] = batteryFeatures(p);
            end
            r.ok = true;
        catch err
            r.msg = err.message;
        end
        parsaveR(Ft{ii}, r);
        fprintf('%s point %d task %d: ok=%d %s\n', datestr(now, 'HH:MM:SS'), Tt(ii, 1), Tt(ii, 2), r.ok, r.msg); %#ok<TNOW1,DATST>
    end
end

% --- compare the real move with the linear prediction
seB = [A.semPts, A.semTrc]; datB = [A.dPts, A.dTrc]; v0B = [A.Vp(1, :), A.Vt(1, :)];
v0C = A.Vc(1, A.okB); seC = A.semBat(A.okB); datC = A.dBat(A.okB);
chi = @(v, d, se) sum(((v - d) ./ se).^2, 'omitnan');
PR = table();
for i = 1:nP
    Fe = cell(1, 3); okB = true;
    for a = 1:3, r = load(fullfile(outP, sprintf('p%02d_r%d.mat', i, a))); okB = okB && r.r.ok; Fe{a} = r.r.F; end
    mvB = NaN; dChiB = NaN; mvC = NaN; dChiC = NaN; lost = NaN;
    if okB
        [vp, vt] = beardVector(recs, Fe, [], A.S); v = [vp, vt];
        mvB = chi(v, v0B, seB); dChiB = chi(v, datB, seB) - chi(v0B, datB, seB);
    end
    r = load(fullfile(outP, sprintf('p%02d_r4.mat', i)));
    if r.r.ok
        v = batteryVector(r.r.fm, r.r.fd, pBat.fn); v = v(A.okB); fe = isfinite(v); lost = sum(~fe);
        mvC = chi(v(fe), v0C(fe), seC(fe)); dChiC = chi(v(fe), datC(fe), seC(fe)) - chi(v0C(fe), datC(fe), seC(fe));
    end
    PR = [PR; table(string(plan(i).target), exp(plan(i).delta), plan(i).predChi2B, mvB, dChiB, plan(i).predChi2C, mvC, dChiC, lost, ...
        'VariableNames', {'target', 'factor', 'predMoveB', 'realMoveB', 'dChi2DataB', 'predMoveC', 'realMoveC', 'dChi2DataC', 'batteryFeaturesLost'})]; %#ok<AGROW>
end
fprintf('\nProfile check. "Move" = squared distance of the model features from the incumbent''s, in standard errors (chi2 units):\n');
fprintf('pred = from the Fisher information, real = from the simulation. dChi2Data = change of the misfit to the data.\n');
fprintf('B = Beard2022 records 8/4/2 mM (%d elements), C = slack battery (%d elements).\n', numel(v0B), numel(v0C));
disp(PR);
save(fullfile(here, 'results', 'profile.mat'), 'PR', 'plan');

function parsaveR(f, r) %#ok<INUSD>
save(f, 'r');
end
