% probeP1dSL.m
% =========================================================================
% Test the sensitivity-surfaced peak1_dSL levers kstiff1 (S=2.69) and k_1
% (S=2.32) -- both flagged as reducing peak1_dSL, and kstiff1 also reduces
% peak1_y (both currently too high). Question: do they fix peak1_dSL WITHOUT
% the ktr penalty that kSE carried (kSE x2 -> peak1_dSL 0.028 but ktr 52->61)?
% If kstiff1/k_1 down move peak1_dSL toward data (0.0256) while holding ktr(~49-52)
% and steady, that's a clean 2-state win on the dominant restretch residual.
% Base = params/params_2state_a2hop.m (2.95).
%
% Writes Analyses/FeatureSensitivity/probeP1dSL_report.txt
% Run:  cd(root); addpath(genpath('.')); probeP1dSL
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP = 'params/params_2state_a2hop.m';
REPORT = fullfile('Analyses','FeatureSensitivity','probeP1dSL_report.txt');
MRT = 90;

params0 = getParams(); run(SNAP); FN = params0.fn;
ks1 = params0.kstiff1; k_1 = params0.k_1; kSE0 = params0.kSE;

L = {};
L{end+1} = {'BASE',              struct()};
L{end+1} = {'kstiff1 x0.7',      struct('kstiff1',ks1*0.7)};
L{end+1} = {'kstiff1 x0.5',      struct('kstiff1',ks1*0.5)};
L{end+1} = {'k_1 x0.6',          struct('k_1',k_1*0.6)};
L{end+1} = {'k_1 x0.3',          struct('k_1',k_1*0.3)};
L{end+1} = {'ks1x0.6 k_1x0.5',   struct('kstiff1',ks1*0.6,'k_1',k_1*0.5)};
L{end+1} = {'ks1x0.5 kSEx1.4',   struct('kstiff1',ks1*0.5,'kSE',kSE0*1.4)};

fid = fopen(REPORT,'w');
fprintf(fid, '==== probeP1dSL on %s ====\n', SNAP);
fprintf(fid, 'Data: peak1_dSL~0.0256, peak1_y~89.3, vall_y~71, peak2~78, steady~77, ktr~49\n');
fprintf(fid, 'GOAL: peak1_dSL down toward 0.0256 while HOLDING ktr(~52) and steady.\n\n');
fprintf(fid, '%-16s %8s %7s %7s %7s %7s %6s | %-18s %7s\n', ...
    'lever','pk1dSL','peak1y','vall_y','peak2','steady','ktr','FV(1,2)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc,~,fm,~] = costOfSnap(SNAP, FN, MRT, extra);
        g = @(n) lm(fm,n); fv = lfv(fm);
        fprintf(fid, '%-16s %8.4f %7.1f %7.1f %7.1f %7.1f %6.1f | %-18s %7.3f\n', ...
            label, g('peak1_dSL'), g('peak1_y'), g('vall_y'), g('peak2'), g('steady'), g('ktr'), ...
            sprintf('%.2f %.2f', fv(2), fv(3)), tc);
    catch e
        fprintf(fid, '%-16s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end
fclose(fid); type(REPORT); disp('DONE probeP1dSL');

function v = lm(fm,n); if isfield(fm,n)&&~isempty(fm.(n)); v=mean(fm.(n),'omitnan'); else; v=NaN; end; end
function fv = lfv(fm); fv=[NaN NaN NaN NaN NaN]; if isfield(fm,'FV_fnorm')&&numel(fm.FV_fnorm)>=5; fv=fm.FV_fnorm(:)'; end; end
