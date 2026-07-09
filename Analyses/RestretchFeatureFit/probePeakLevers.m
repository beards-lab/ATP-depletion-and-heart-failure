% probePeakLevers.m
% =========================================================================
% Probe two under-explored lever families against the new best 2-state
% snapshot (params/opt2state_v2_opt.m, cost ~3.22). Targets the now-dominant
% residuals: peak1_dSL (0.74), FV mid-tail (0.61), ovrsht (0.40), peak1_y (0.39).
%
% FAMILY A -- strain-limited CATCH BOND (restretch-only detachment suppression).
%   Coded but OFF; blew up in probeCoupling at k=0.3 because CatchBondStrainMax
%   defaults to Inf (ALL strains suppressed -> force runaway). Here it is made
%   strain-limited (finite CatchBondStrainMax) and swept at small k, plus the
%   milder R1D-only variant. Expected target: vall_y (raise the too-deep valley)
%   and the valley->peak2 recovery, WITHOUT touching isometric steady (vel>0 gate).
%
% FAMILY B -- NEGATIVE-STRAIN STIFFNESS (kstiff1_n, kstiff2_n).
%   The force from compressed heads (s<0) is a single linear coefficient -- a
%   simplification of nonlinear XB elasticity, never tuned/pooled. During
%   shortening, compressed p2 heads contribute RESISTIVE force
%   kstiff2_n*p2_1_neg (<0); lowering kstiff2_n reduces drag -> should LIFT the
%   FV mid-tail (the residual ktr<->FV-tail tension). Check it holds steady/peaks.
%
% Writes Analyses/RestretchFeatureFit/probePeakLevers_report.txt
% Run:  cd(root); addpath(genpath('.')); probePeakLevers
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/opt2state_v2_opt.m';
REPORT = fullfile('Analyses','RestretchFeatureFit','probePeakLevers_report.txt');
MRT    = 60;
SMAX   = 0.008;   % strain-limited catch-bond cutoff (um): s<=SMAX suppressed, high-strain slip detaches normally

params0 = getParams(); run(SNAP);
FN = params0.fn;
b.k1n = params0.kstiff1_n; b.k2n = params0.kstiff2_n;

L = {};
L{end+1} = {'BASE',              struct()};
% --- Family A: strain-limited catch bond ---
L{end+1} = {'cb k.005 sm.008',   struct('UseCatchBond',1,'k_catch_bond',0.005,'CatchBondStrainMax',SMAX,'UseCatchBondR1DOnly',0)};
L{end+1} = {'cb k.01 sm.008',    struct('UseCatchBond',1,'k_catch_bond',0.01, 'CatchBondStrainMax',SMAX,'UseCatchBondR1DOnly',0)};
L{end+1} = {'cb k.02 sm.008',    struct('UseCatchBond',1,'k_catch_bond',0.02, 'CatchBondStrainMax',SMAX,'UseCatchBondR1DOnly',0)};
L{end+1} = {'cb k.02 R1Donly',   struct('UseCatchBond',1,'k_catch_bond',0.02, 'UseCatchBondR1DOnly',1)};
L{end+1} = {'cb k.05 R1Donly',   struct('UseCatchBond',1,'k_catch_bond',0.05, 'UseCatchBondR1DOnly',1)};
% --- Family B: negative-strain stiffness ---
L{end+1} = {'k2n x0.3',          struct('kstiff2_n', b.k2n*0.3)};
L{end+1} = {'k2n x0.5',          struct('kstiff2_n', b.k2n*0.5)};
L{end+1} = {'k2n x0.7',          struct('kstiff2_n', b.k2n*0.7)};
L{end+1} = {'k2n x1.5',          struct('kstiff2_n', b.k2n*1.5)};
L{end+1} = {'k1n x0.5',          struct('kstiff1_n', b.k1n*0.5)};
L{end+1} = {'k1n x2.0',          struct('kstiff1_n', b.k1n*2.0)};
L{end+1} = {'k1n+k2n x0.5',      struct('kstiff1_n', b.k1n*0.5,'kstiff2_n', b.k2n*0.5)};

fid = fopen(REPORT, 'w');
fprintf(fid, '==== probePeakLevers on %s ====\n', SNAP);
fprintf(fid, 'Targets: peak1_dSL(0.74), FV mid-tail(0.61), ovrsht(0.40), peak1_y(0.39).\n');
fprintf(fid, 'Data: peak1_dSL~0.025, peak1_y~[96 94 90 89 77], vall_y~[77 75 72 70 61], FV~[.92 .66 .32 .11], steady~[80 80 81 81 64], ktr~49\n\n');
fprintf(fid, '%-16s %7s %7s %7s %7s %7s %6s | %-24s %7s\n', ...
    'lever','pk1dSL','peak1y','vall_y','peak2','steady','ktr','FV_fnorm(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc, ~, fm, ~] = costOfSnap(SNAP, FN, MRT, extra);
        g  = @(nm) local_meanfeat(fm, nm);
        fv = local_fvtail(fm);
        pk2 = local_meanfeat(fm,'peak2');
        flag = '';
        if pk2 > 200 || isnan(local_meanfeat(fm,'peak1_y')); flag = '  <-- UNSTABLE'; end
        fprintf(fid, '%-16s %7.4f %7.1f %7.1f %7.1f %7.1f %6.1f | %-24s %7.3f%s\n', ...
            label, g('peak1_dSL'), g('peak1_y'), g('vall_y'), pk2, g('steady'), g('ktr'), ...
            sprintf('%.2f %.2f %.2f %.2f', fv(1),fv(2),fv(3),fv(4)), tc, flag);
    catch e
        fprintf(fid, '%-16s  ERROR: %s\n', label, e.message);
    end
    fprintf('done lever %d/%d: %s\n', i, numel(L), label);
end
fclose(fid);
type(REPORT);
disp('DONE probePeakLevers');

function v = local_meanfeat(fm, name)
    if isfield(fm, name) && ~isempty(fm.(name)); v = mean(fm.(name), 'omitnan'); else; v = NaN; end
end
function fv = local_fvtail(fm)
    fv = [NaN NaN NaN NaN];
    if isfield(fm, 'FV_fnorm') && numel(fm.FV_fnorm) >= 5; x = fm.FV_fnorm(:)'; fv = x(2:5); end
end
