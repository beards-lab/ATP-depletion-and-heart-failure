% RunMaxwellTensionOnly.m — is `UseMaxwellTensionOnly` still the free win it was
% called, at the CURRENT parameter point and the CURRENT objective?
%
% Background. RestretchVsKtrRecovery Part 2 "Defect 1": the flag is unset, so the
% Maxwell (titin) element only ever CHARGES —
%       default   dX/dt = kSE_M*max(0,dSLc) - X/eta_M
%       with flag dX/dt = kSE_M*max(0,dSLc) - X/eta_M - X*max(0,-dSLc)/x_M_slack
% i.e. by default it gains stress on lengthening and can shed it only on the eta_M
% CLOCK; with the flag, shortening pays out slack and dumps the stored stress over a
% LENGTH constant x_M_slack. The ktr protocol ends on a 1.05 -> 1.00 ML release,
% which by default cannot discharge it, so the redevelopment window started at
% 0.18 F_iso with an apparent rate of 192 s^-1. Setting the flag gave ktr 192 -> 85.5,
% E_ktr 0.98 -> 0.225 AND E_slack 183 -> 173, and was called "strictly an
% improvement — take it before anything else". It was never applied:
% params/rskR2_w025_opt.m still has UseMaxwellTensionOnly = 0.
%
% Why the answer may have changed at this snapshot:
%   (a) the current objective sets RunKtr = 0, so the protocol the fix was FOR is not
%       scored at all;
%   (b) the release term only fires while SHORTENING, and in the slack protocol every
%       shortening step is preceded by a ~275 ms hold with eta_M = 27 ms — about ten
%       time constants, so there should be nothing left to dump.
% Prediction from (b): nearly inert here. This script measures it rather than
% assuming it, and reads out how much stress is actually standing at each release.
%
% SLACK PROTOCOL ONLY and SEQUENTIAL (no parpool), so it can be run beside a
% live optimiser without either job tripping its MaxRunTime guard. FV/passive are
% excluded from the scoring here for the same reason — see the note printed below.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

S  = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
vt = S.velocitytable; iR = find(vt(:,2) > 1); iS = find(vt(:,2) < -1);

CASES = { 'flag OFF (current)',  struct('UseMaxwellTensionOnly', false)
          'flag ON  x_M=0.01',   struct('UseMaxwellTensionOnly', true)
          'flag ON  x_M=0.003',  struct('UseMaxwellTensionOnly', true, 'x_M_slack', 0.003) };

R = struct('name',{},'ct2',{},'L2',{},'feat',{},'out',{});
for i = 1:size(CASES,1)
    p = getParams(loadParams('params/rskR2_w025_opt.m'), [], true, false);
    f = fieldnames(CASES{i,2}); for k=1:numel(f); p.(f{k}) = CASES{i,2}.(f{k}); end
    p = getParams(p, [], true, false);
    p.MaxRunTime = 600; p.RunSlackSegments = 'All';   % sequential, 1 core
    p.EvalFeatures = 1; p.PlotEachSeparately = 0; p.PlotFeatureFitting = 0;

    fprintf('[%d/%d] %-22s ', i, size(CASES,1), CASES{i,1}); tic;
    [~, o, fm, fd] = runSlackExperiment(p);
    el = toc;

    if i == 1
        prim = cellfun(@(s) regexprep(strtok(s,'|'), '\[.*$',''), p.fn, 'uni', 0);
        scorable = cellfun(@(q) isfield(fm,q) && isfield(fd,q), prim);
        FN = p.fn(scorable); FD = fd;
    end
    ct2 = evalFeatureCost(fd, fm, FN, 2);
    R(end+1) = struct('name', CASES{i,1}, 'ct2', ct2(:), 'L2', sum(ct2), ...
                      'feat', fm, 'out', o); %#ok<SAGROW>
    fprintf('slack L2 %8.4f   (%.0f s)\n', sum(ct2), el);
end

nmf = cellfun(@(x) regexprep(strtok(x,'|'), '\[.*$',''), FN, 'uni', 0);
fprintf('\nNOTE: slack-protocol features only (%d terms). FV and passive are not run\n', numel(FN));
fprintf('here so as not to compete with the optimiser for cores.\n');

for i = 2:numel(R)
    d = R(i).ct2 - R(1).ct2; [~, o] = sort(abs(d), 'descend');
    fprintf('\n%s  vs  %s   (L2 %.4f -> %.4f)\n', R(i).name, R(1).name, R(1).L2, R(i).L2);
    moved = false;
    for j = o(1:min(8,numel(o)))
        if abs(d(j)) < 1e-4; break; end
        fprintf('    %-20s %8.4f -> %8.4f  %+8.4f\n', nmf{j}, R(1).ct2(j), R(i).ct2(j), d(j));
        moved = true;
    end
    if ~moved
        fprintf('    (nothing moves by more than 1e-4 — the flag is INERT on this protocol)\n');
    end
end

fprintf('\n%-20s %10s %10s %10s\n', 'observable', 'DATA', R(1).name, R(2).name);
for k = {'rsK','ktr','steady','vall_y','peak1_y'}
    if ~isfield(FD, k{1}); continue; end
    fprintf('%-20s %10.2f %10.2f %10.2f\n', k{1}, mean(FD.(k{1}),'omitnan'), ...
        mean(R(1).feat.(k{1}),'omitnan'), mean(R(2).feat.(k{1}),'omitnan'));
end

%% ---- the physical question: is there any stress left to dump? ------------
% The flag can only act while shortening. In this protocol that is the release at
% the START of each cycle. If F_Maxwell has already decayed to ~0 by then, the flag
% has nothing to do — which is the (b) prediction above.
fprintf('\n---- standing Maxwell stress at each event (flag OFF) ----\n');
o = R(1).out;
fprintf('%-34s %10s %10s %8s\n', 'event', 'F_total', 'F_Maxwell', 'share');
ev = {};
for c = 1:4
    ev(end+1,:) = {sprintf('cyc%d restretch ramp END', c), vt(iR(c)+1,1)}; %#ok<SAGROW>
    ev(end+1,:) = {sprintf('cyc%d  +40 ms', c),           vt(iR(c)+1,1)+0.040}; %#ok<SAGROW>
    if c < 4
        ev(end+1,:) = {sprintf('cyc%d NEXT RELEASE (flag acts here)', c+1), vt(iS(c+1),1)}; %#ok<SAGROW>
    end
end
for e = 1:size(ev,1)
    [~, iq] = min(abs(o.t - ev{e,2}));
    fprintf('%-34s %10.2f %10.3f %7.1f%%\n', ev{e,1}, o.Force(iq), o.F_Maxwell(iq), ...
        100*o.F_Maxwell(iq)/max(o.Force(iq),eps));
end
fprintf('\neta_M = %.4f s  =>  time constant %.1f ms;  hold between restretch and\n', ...
        R(1).out.t(1)*0 + getParams(loadParams('params/rskR2_w025_opt.m'),[],true,false).eta_M, ...
        1000*getParams(loadParams('params/rskR2_w025_opt.m'),[],true,false).eta_M);
fprintf('the next release is ~%.0f ms = %.1f time constants.\n', ...
        1000*(vt(iS(2),1) - vt(iR(1)+1,1)), ...
        (vt(iS(2),1) - vt(iR(1)+1,1))/getParams(loadParams('params/rskR2_w025_opt.m'),[],true,false).eta_M);

save(fullfile(here,'results','maxwell_tension.mat'), 'R', 'CASES', 'FN', 'FD');
