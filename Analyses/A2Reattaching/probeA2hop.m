% probeA2hop.m
% =========================================================================
% Extended A2 mechanism: does the restretch HOP (high-strain p2 heads detach
% faster and re-land at low strain, kA2hop) cap the restretch PEAK (peak1_y,
% peak1_dSL -- the major cost 0.74+0.39), alone and combined with the valley-
% fill (kA2re)? Base = params/opt2state_v2_opt.m (2-state, cost 3.224).
%
% peak1 forms over a few ms, so the hop rate must be O(100-500/s) at peak
% strains: R_hop = kA2hop*(s-sA2hop)*p2, so at s-sA2hop~0.015 um, kA2hop~2e4
% gives ~300/s. Sweep a wide range to find where it bites.
%
% Writes Analyses/A2Reattaching/probeA2hop_report.txt
% Run:  cd(root); addpath(genpath('.')); probeA2hop
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_v2_opt.m';
REPORT = fullfile('Analyses','A2Reattaching','probeA2hop_report.txt');
MRT    = 60;

params0 = getParams(); run(SNAP); FN = params0.fn;

on = @(varargin) struct('UseA2Reattaching',1, varargin{:});
L = {};
L{end+1} = {'BASE off',          struct()};
L{end+1} = {'hop2k',             on('kA2hop',2000,'sA2hop',0.005)};
L{end+1} = {'hop8k',             on('kA2hop',8000,'sA2hop',0.005)};
L{end+1} = {'hop30k',            on('kA2hop',30000,'sA2hop',0.005)};
L{end+1} = {'hop80k',            on('kA2hop',80000,'sA2hop',0.005)};
L{end+1} = {'re10+hop8k',        on('kA2re',10,'kA2hop',8000,'sA2hop',0.005)};
L{end+1} = {'re10+hop30k',       on('kA2re',10,'kA2hop',30000,'sA2hop',0.005)};
L{end+1} = {'re10+hop30k s.002', on('kA2re',10,'kA2hop',30000,'sA2hop',0.002)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probeA2hop on %s ====\n', SNAP);
fprintf(fid, 'Data targets: peak1_y~96, peak1_dSL~0.025, vall_y~71, peak2~81, steady~[80..64], FV~[.92 .66 .32 .11], ktr~49\n');
fprintf(fid, 'Goal: hop should DROP peak1_y and/or peak1_dSL (the major cost) -- ideally with valley-fill (kA2re) holding vall_y.\n\n');
fprintf(fid, '%-16s %7s %8s %7s %7s %7s %7s %6s | %-22s %7s\n', ...
    'lever','peak1y','pk1dSL','vall_y','peak2','ovrsht','steady','ktr','FV(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        g  = @(nm) local_mean(fm, nm);
        fv = local_fv(fm);
        fprintf(fid, '%-16s %7.1f %8.4f %7.1f %7.1f %7.2f %7.1f %6.1f | %-22s %7.3f\n', ...
            label, g('peak1_y'), g('peak1_dSL'), g('vall_y'), g('peak2'), ...
            g('ovrsht_dy'), g('steady'), g('ktr'), ...
            sprintf('%.2f %.2f %.2f %.2f', fv(1),fv(2),fv(3),fv(4)), tc);
    catch e
        fprintf(fid, '%-16s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end

fclose(fid);
type(REPORT);
disp('DONE probeA2hop');

function v = local_mean(fm, name)
    if isfield(fm, name) && ~isempty(fm.(name)); v = mean(fm.(name), 'omitnan'); else; v = NaN; end
end
function fv = local_fv(fm)
    fv = [NaN NaN NaN NaN];
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5; x = fm.FV_fnorm(:)'; fv = x(2:5); end
end
