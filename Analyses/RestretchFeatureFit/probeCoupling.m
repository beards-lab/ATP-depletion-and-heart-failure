% probeCoupling.m
% =========================================================================
% Coupling probe for the current best 2-state snapshot (params/opt2state_opt.m).
% Question the user posed: are the ktr and restretch-peak residuals genuine
% structural "iron walls", or can a single lever lower them WITHOUT breaking
% force / FV / physiology?
%
% Strategy: perturb one candidate lever at a time (absolute overrides via
% costOfSnap's `extra` struct), and record the outcome vector
%   [ktr, steady, peak1_y, peak2, vall_y, XTOR_vmax, attached_ss, SRX_ss,
%    FV_fnorm(2:5), total_cost].
% A lever that moves ktr (or peak1_y) but LEAVES steady/FV/physiology put is a
% viable decoupling lever -> the wall is soft. A lever that drags everything
% together is the wall.
%
% Candidate levers (mechanistic rationale):
%   k1x/k2x/cyclex  - Brenner f+g: does ktr track the force-generating cycle?
%   kmsrd           - SRX recruitment speed (currently 60, 3x its soft bound):
%                     does slowing SRX return slow ktr at fixed force?
%   dr1             - Route B: p1 pre-stroke force-duty, the code's own
%                     documented ktr<->FV decoupler (currently ~0.0047).
%   kstiff2/kSE/Pbm - restretch peak height vs steady coupling.
%   catchbond       - built-in restretch detachment-suppression shaper (off).
%
% Writes Analyses/RestretchFeatureFit/probeCoupling_report.txt
% Run:  cd(root); addpath(genpath('.')); probeCoupling
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_opt.m';
REPORT = fullfile('Analyses','RestretchFeatureFit','probeCoupling_report.txt');
MRT    = 60;

% --- read fn + base rate values off the snapshot ---
params0 = getParams(); run(SNAP);
FN = params0.fn;
b.k1=params0.k1; b.k2=params0.k2; b.ka=params0.ka; b.kd=params0.kd;
b.k_1=params0.k_1; b.k_2=params0.k_2; b.kstiff2=params0.kstiff2; b.kSE=params0.kSE;
b.Pbm=params0.P_bound_max; b.kmsrd=params0.kmsrd; b.dr1=params0.dr1;

% --- build lever list: {label, extra-struct} ---
L = {};
L{end+1} = {'BASE',            struct()};
L{end+1} = {'k1 x0.75',        struct('k1', b.k1*0.75)};
L{end+1} = {'k2 x0.75',        struct('k2', b.k2*0.75)};
L{end+1} = {'cycle x0.80',     struct('k1',b.k1*0.8,'k2',b.k2*0.8,'ka',b.ka*0.8,'kd',b.kd*0.8,'k_1',b.k_1*0.8,'k_2',b.k_2*0.8)};
L{end+1} = {'kmsrd 30',        struct('kmsrd', 30)};
L{end+1} = {'kmsrd 15',        struct('kmsrd', 15)};
L{end+1} = {'dr1 0.010',       struct('dr1', 0.010)};
L{end+1} = {'dr1 0.015',       struct('dr1', 0.015)};
L{end+1} = {'kstiff2 x0.85',   struct('kstiff2', b.kstiff2*0.85)};
L{end+1} = {'kSE x0.70',       struct('kSE', b.kSE*0.70)};
L{end+1} = {'P_bound x0.85',   struct('P_bound_max', b.Pbm*0.85)};
L{end+1} = {'catchbond 0.3',   struct('UseCatchBond',1,'k_catch_bond',0.3)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probeCoupling on %s ====\n', SNAP);
fprintf(fid, 'Each row: mean-over-segments of features; FV = fnorm(v=.5,1,2,4); COST = total feature cost.\n');
fprintf(fid, 'Reference data targets: ktr~49, steady~[80 80 81 81 64], peak1_y~96, FV~[.92 .66 .32 .11]\n\n');
fprintf(fid, '%-14s %6s %7s %7s %7s %7s %8s %8s %7s | %-26s %7s\n', ...
    'lever','ktr','steady','peak1y','peak2','vall_y','XTORvmx','att_ss','SRX_ss','FV_fnorm(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        g  = @(nm) local_meanfeat(fm, nm);
        fv = local_fvtail(fm);
        fprintf(fid, '%-14s %6.1f %7.1f %7.1f %7.1f %7.1f %8.2f %8.3f %7.3f | %-26s %7.3f\n', ...
            label, g('ktr'), g('steady'), g('peak1_y'), g('peak2'), g('vall_y'), ...
            g('XTOR_vmax'), g('attached_ss'), g('SRX_ss'), ...
            sprintf('%.2f %.2f %.2f %.2f', fv(1),fv(2),fv(3),fv(4)), tc);
    catch e
        fprintf(fid, '%-14s  ERROR: %s\n', label, e.message);
    end
    fprintf('done lever %d/%d: %s\n', i, numel(L), label);
end

fclose(fid);
type(REPORT);
disp('DONE probeCoupling');

% ---- local functions ----
function v = local_meanfeat(fm, name)
    if isfield(fm, name) && ~isempty(fm.(name))
        v = mean(fm.(name), 'omitnan');
    else
        v = NaN;
    end
end

function fv = local_fvtail(fm)
    fv = [NaN NaN NaN NaN];
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5
        x = fm.FV_fnorm(:)';
        fv = x(2:5);   % drop v=0 (=1 by normalisation)
    end
end
