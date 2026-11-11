% probeP1dSLcomp.m
% =========================================================================
% Deciding test: is peak1_dSL threadable in 2-state, or doubly-walled?
% kstiff1 down fixes peak1_dSL + peak1_y and HOLDS ktr, but collapses the valley
% (vall_y, peak2) by removing p1 force. The A2 hop adds valley force via p2
% INDEPENDENTLY -> kstiff1 down + kA2hop up (refill valley) [+ kstiff2 up to hold
% peak2/steady] may recover a point with peak1_dSL~0.026 AND valley/ktr held.
% If yes -> peak1_dSL is a (multi-lever) 2-state win, not a wall.
% Base = params/params_2state_a2hop.m (2.95).
%
% Writes Analyses/FeatureSensitivity/probeP1dSLcomp_report.txt
% Run:  cd(root); addpath(genpath('.')); probeP1dSLcomp
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP = 'params/params_2state_a2hop.m';
REPORT = fullfile('Analyses','FeatureSensitivity','probeP1dSLcomp_report.txt');
MRT = 110;

params0 = getParams(); run(SNAP); FN = params0.fn;
ks1 = params0.kstiff1; ks2 = params0.kstiff2; kh = params0.kA2hop; kSE0 = params0.kSE;

L = {};
L{end+1} = {'BASE',                    struct()};
L{end+1} = {'ks1x0.7',                 struct('kstiff1',ks1*0.7)};
L{end+1} = {'ks1x0.7 hopx1.6',         struct('kstiff1',ks1*0.7,'kA2hop',kh*1.6)};
L{end+1} = {'ks1x0.7 hopx1.6 ks2x1.08',struct('kstiff1',ks1*0.7,'kA2hop',kh*1.6,'kstiff2',ks2*1.08)};
L{end+1} = {'ks1x0.7 hopx2.2 ks2x1.12',struct('kstiff1',ks1*0.7,'kA2hop',kh*2.2,'kstiff2',ks2*1.12)};
L{end+1} = {'ks1x0.6 hopx2 ks2x1.1',   struct('kstiff1',ks1*0.6,'kA2hop',kh*2.0,'kstiff2',ks2*1.10)};

fid = fopen(REPORT,'w');
fprintf(fid, '==== probeP1dSLcomp on %s ====\n', SNAP);
fprintf(fid, 'Data: peak1_dSL~0.0256, peak1_y~89.3, vall_y~71, peak2~78, steady~77, ktr~49\n');
fprintf(fid, 'GOAL: peak1_dSL~0.026 AND vall_y~71 AND peak2~78 AND ktr~52 together.\n\n');
fprintf(fid, '%-24s %8s %7s %7s %7s %7s %6s | %-12s %7s\n', ...
    'lever','pk1dSL','peak1y','vall_y','peak2','steady','ktr','FV(1,2)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc,~,fm,~] = costOfSnap(SNAP, FN, MRT, extra);
        g = @(n) lm(fm,n); fv = lfv(fm);
        fprintf(fid, '%-24s %8.4f %7.1f %7.1f %7.1f %7.1f %6.1f | %-12s %7.3f\n', ...
            label, g('peak1_dSL'), g('peak1_y'), g('vall_y'), g('peak2'), g('steady'), g('ktr'), ...
            sprintf('%.2f %.2f', fv(2), fv(3)), tc);
    catch e
        fprintf(fid, '%-24s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end
fclose(fid); type(REPORT); disp('DONE probeP1dSLcomp');

function v = lm(fm,n); if isfield(fm,n)&&~isempty(fm.(n)); v=mean(fm.(n),'omitnan'); else; v=NaN; end; end
function fv = lfv(fm); fv=[NaN NaN NaN NaN NaN]; if isfield(fm,'FV_fnorm')&&numel(fm.FV_fnorm)>=5; fv=fm.FV_fnorm(:)'; end; end
