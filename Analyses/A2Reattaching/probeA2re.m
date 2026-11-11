% probeA2re.m
% =========================================================================
% Analysis: does restretch-gated A2 re-attachment (UseA2Reattaching) improve the
% 2-state restretch fit? Base = params/opt2state_v2_opt.m (2-state, cost ~3.22).
%
% Mechanism (fixed + gated this campaign): during lengthening (vel>0), ATP-cocked
% heads (PT) bind DIRECTLY into the strong post-stroke state p2 at strain dA2re,
% then are carried to +strain by the ongoing stretch -> add force through the
% restretch/valley phase. Rate RT2 = kA2re*PT. Unlike catch bond (freezes existing
% strained bridges -> couples valley^ with peak^), this adds NEW low-strain heads,
% so the hypothesis is: raise vall_y toward data WITHOUT inflating peak1_y.
%
% Target residuals (base 3.22): vall_y 66.9 -> data ~71 (too deep); peak1_y ~100
% -> data ~96 (must NOT rise); peak1_dSL, ovrsht. steady/FV must stay put (the
% mechanism is vel>0-gated so isometric+shortening are untouched by construction).
%
% Writes Analyses/A2Reattaching/probeA2re_report.txt
% Run:  cd(root); addpath(genpath('.')); probeA2re
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_v2_opt.m';
REPORT = fullfile('Analyses','A2Reattaching','probeA2re_report.txt');
MRT    = 60;

params0 = getParams(); run(SNAP); FN = params0.fn;

% lever list: {label, extra-struct}. kA2re = re-attach rate; dA2re = landing strain.
L = {};
L{end+1} = {'BASE (off)',        struct()};
L{end+1} = {'kA2re10 d0',        struct('UseA2Reattaching',1,'kA2re',10,'dA2re',0)};
L{end+1} = {'kA2re30 d0',        struct('UseA2Reattaching',1,'kA2re',30,'dA2re',0)};
L{end+1} = {'kA2re80 d0',        struct('UseA2Reattaching',1,'kA2re',80,'dA2re',0)};
L{end+1} = {'kA2re30 d.005',     struct('UseA2Reattaching',1,'kA2re',30,'dA2re',0.005)};
L{end+1} = {'kA2re30 d.010',     struct('UseA2Reattaching',1,'kA2re',30,'dA2re',0.010)};
L{end+1} = {'kA2re80 d.010',     struct('UseA2Reattaching',1,'kA2re',80,'dA2re',0.010)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probeA2re on %s ====\n', SNAP);
fprintf(fid, 'Data targets: vall_y~71, peak1_y~96, peak1_dSL~0.025, steady~[80..64], FV~[.92 .66 .32 .11], ktr~49\n');
fprintf(fid, 'Hypothesis: raise vall_y toward 71 WITHOUT raising peak1_y (the catch-bond failure mode).\n\n');
fprintf(fid, '%-14s %7s %7s %8s %7s %7s %7s %6s | %-22s %7s\n', ...
    'lever','vall_y','peak1y','pk1dSL','peak2','ovrsht','steady','ktr','FV(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        g  = @(nm) local_mean(fm, nm);
        fv = local_fv(fm);
        fprintf(fid, '%-14s %7.1f %7.1f %8.4f %7.1f %7.2f %7.1f %6.1f | %-22s %7.3f\n', ...
            label, g('vall_y'), g('peak1_y'), g('peak1_dSL'), g('peak2'), ...
            g('ovrsht_dy'), g('steady'), g('ktr'), ...
            sprintf('%.2f %.2f %.2f %.2f', fv(1),fv(2),fv(3),fv(4)), tc);
    catch e
        fprintf(fid, '%-14s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end

fclose(fid);
type(REPORT);
disp('DONE probeA2re');

function v = local_mean(fm, name)
    if isfield(fm, name) && ~isempty(fm.(name)); v = mean(fm.(name), 'omitnan'); else; v = NaN; end
end
function fv = local_fv(fm)
    fv = [NaN NaN NaN NaN];
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5; x = fm.FV_fnorm(:)'; fv = x(2:5); end
end
