% probeA2retune.m
% =========================================================================
% Best A2 config is hop-only (kA2hop=30000, sA2hop=0.005 -> cost 2.952), which
% fixes vall_y/peak2 but leaves the DOMINANT restretch residual peak1_dSL (0.036
% vs data 0.025) untouched -- because peak1_dSL is a series-compliance/geometry
% quantity, not a p2-population one. The lever for peak1_dSL is kSE (stiffer
% series element -> peak reached at smaller SL excursion). Test hop + kSE retune
% (+ kstiff1 for the p1 spike height) to see if the peak cluster finally moves
% and whether total cost drops below 2.95. Base = params/opt2state_v2_opt.m.
%
% Writes Analyses/A2Reattaching/probeA2retune_report.txt
% Run:  cd(root); addpath(genpath('.')); probeA2retune
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_v2_opt.m';
REPORT = fullfile('Analyses','A2Reattaching','probeA2retune_report.txt');
MRT    = 60;

params0 = getParams(); run(SNAP); FN = params0.fn;
kSE0 = params0.kSE; ks1 = params0.kstiff1;

% explicit structs (struct() rejects duplicate field names, so no override helper)
b = {'UseA2Reattaching',1,'sA2hop',0.005};
L = {};
L{end+1} = {'hop30k',              struct(b{:},'kA2hop',30000)};
L{end+1} = {'hop30k kSEx1.3',      struct(b{:},'kA2hop',30000,'kSE',kSE0*1.3)};
L{end+1} = {'hop30k kSEx1.6',      struct(b{:},'kA2hop',30000,'kSE',kSE0*1.6)};
L{end+1} = {'hop30k kSEx2.0',      struct(b{:},'kA2hop',30000,'kSE',kSE0*2.0)};
L{end+1} = {'hop30k kSEx1.6 ks1x1.1', struct(b{:},'kA2hop',30000,'kSE',kSE0*1.6,'kstiff1',ks1*1.1)};
L{end+1} = {'hop20k kSEx1.6',      struct(b{:},'kA2hop',20000,'kSE',kSE0*1.6)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probeA2retune on %s ====  (hop30k alone = 2.952)\n', SNAP);
fprintf(fid, 'Targets: peak1_y~96, peak1_dSL~0.025, vall_y~71, peak2~81, steady~[80..64], ktr~49\n\n');
fprintf(fid, '%-22s %7s %8s %7s %7s %7s %6s | %-22s %7s\n', ...
    'lever','peak1y','pk1dSL','vall_y','peak2','steady','ktr','FV(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        g  = @(nm) local_mean(fm, nm);
        fv = local_fv(fm);
        fprintf(fid, '%-22s %7.1f %8.4f %7.1f %7.1f %7.1f %6.1f | %-22s %7.3f\n', ...
            label, g('peak1_y'), g('peak1_dSL'), g('vall_y'), g('peak2'), ...
            g('steady'), g('ktr'), sprintf('%.2f %.2f %.2f %.2f', fv(1),fv(2),fv(3),fv(4)), tc);
    catch e
        fprintf(fid, '%-22s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end

fclose(fid);
type(REPORT);
disp('DONE probeA2retune');

function v = local_mean(fm, name)
    if isfield(fm, name) && ~isempty(fm.(name)); v = mean(fm.(name), 'omitnan'); else; v = NaN; end
end
function fv = local_fv(fm)
    fv = [NaN NaN NaN NaN];
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5; x = fm.FV_fnorm(:)'; fv = x(2:5); end
end
