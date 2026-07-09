% probe3state.m — calibrate the 3-state p3 population.
% Seed params_3state_seedB.m over-accumulates p3 (steady ~400-600 vs ~80) because
% k3 is too slow relative to the p2->p3 feed. Sweep k3 (turnover), kstiff3 (p3
% force) and alpha3 (negative-strain detachment enhancement) to find the basin
% where steady~80 and FV holds (tail not collapsed). costOfSnap `extra` overrides.
% Writes Analyses/RestretchFeatureFit/probe3state_report.txt (incremental).
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));
SNAP   = 'params/params_3state_seedB.m';
REPORT = fullfile('Analyses','RestretchFeatureFit','probe3state_report.txt');
MRT    = 60;
params0 = getParams(); run(SNAP); FN = params0.fn;

% {label, extra}
L = {};
L{end+1} = {'BASE k3=44',        struct()};
L{end+1} = {'k3=200',            struct('k3',200)};
L{end+1} = {'k3=500',            struct('k3',500)};
L{end+1} = {'k3=1000',           struct('k3',1000)};
L{end+1} = {'k3=500 ks12k',      struct('k3',500,'kstiff3',12000)};
L{end+1} = {'k3=1000 ks12k',     struct('k3',1000,'kstiff3',12000)};
L{end+1} = {'k3=500 a250',       struct('k3',500,'alpha3',250)};
L{end+1} = {'k3=500 ks12k a250', struct('k3',500,'kstiff3',12000,'alpha3',250)};
L{end+1} = {'k3=2000 ks12k',     struct('k3',2000,'kstiff3',12000)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probe3state on %s ====\n', SNAP);
fprintf(fid, 'Goal: steady~[80..64], FV~[.92 .66 .32 .11], ktr~49, peak1y~[96..77]. seed base steady~400 (p3 over-accum).\n\n');
fprintf(fid, '%-20s %8s %7s %7s %7s %8s | %-24s %9s\n', ...
    'config','steady','peak1y','vall_y','ktr','att_ss','FV_fnorm(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        fprintf(fid, '%-20s %8.1f %7.1f %7.1f %7.1f %8.3f | %-24s %9.2f\n', ...
            label, mf(fm,'steady'), mf(fm,'peak1_y'), mf(fm,'vall_y'), mf(fm,'ktr'), mf(fm,'attached_ss'), ...
            fvstr(fm), tc);
    catch e
        fprintf(fid, '%-20s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end
fclose(fid);
type(REPORT);
disp('DONE probe3state');

function v = mf(fm, name)
    v = NaN;
    if isfield(fm, name) && ~isempty(fm.(name))
        v = mean(fm.(name), 'omitnan');
    end
end

function s = fvstr(fm)
    s = 'NaN';
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5
        x = fm.FV_fnorm(:)';
        s = sprintf('%.2f %.2f %.2f %.2f', x(2), x(3), x(4), x(5));
    end
end
