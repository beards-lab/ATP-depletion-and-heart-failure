% RunEstiffScan.m — does flipping the p2 force-strain exponent past 1 save the
% staircase prediction, and what does it cost the fit that is already paid for?
%
% estiff is the exponent e in the post-stroke force sum
%     F_2 = kstiff2 * sum_s sign(s+dr) |s+dr|^e p2(s) dS
% e = 1 is a Hookean cross-bridge; e < 1 softens at high strain, e > 1 stiffens.
%
% *** COMPENSATION ***  |s+dr| ~ 0.0095 um at the isometric operating point, so
% |x|^e is a number far below 1 and changing e alone rescales the whole active
% force by ~x_ref^(e0-e) (a factor of 12 between e=0.965 and e=1.5). A raw scan
% would therefore measure "wrong isometric force", not "wrong curvature". For
% every e we rescale kstiff2/kstiff2_n by the factor that holds the isometric
% active force at its baseline value, computed on the frozen isometric p2
% distribution. UNCOMP=true reruns the scan without that rescale as a control.
%
% Outputs -> results/estiff_scan[_uncomp].mat + .png

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end

if ~exist('SNAP','var');    SNAP    = 'spreadw_opt'; end
if ~exist('EGRID','var');   EGRID   = [0.5 0.7 0.85 0.965322 1.10 1.25 1.50 1.75 2.00]; end
if ~exist('UNCOMP','var');  UNCOMP  = false; end
if ~exist('DOFULL','var');  DOFULL  = true;  end   % also evaluate the paid-for objective

base = getParams(loadParams(SNAP), [], true, false);
d    = load('data/protocol_03_27_2026_velocitytable_stairs.mat');
dt   = d.datatable; vt = d.velocitytable;
fD   = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);

%% ---- reference: isometric distribution of the baseline ------------------
pr = base; pr.RunStairs = 1; pr.PlotEachSeparately = false;
pr.ghostSave=''; pr.ghostLoad='';
[E0_lsq, o0] = runStairsExperiment(pr);
ss   = numel(o0.params.s);
tMov = vt(find(vt(:,2)~=0,1),1);
[~,i0] = min(abs(o0.t - (tMov-0.005)));
p2iso = o0.PU(i0, ss+1:2*ss)';
p1iso = o0.PU(i0, 1:ss)';
cShift = (-(o0.SL(i0)-o0.LSE(i0)) + o0.LXBPivot(i0))/2;   % the ODE's OWN shift
xiso   = o0.params.s(:) - cShift + base.dr;
dS     = o0.params.s(2)-o0.params.s(1);

pw   = @(e) dS*sum(sign(xiso).*abs(xiso).^e.*p2iso);      % p2_1 at exponent e
e0   = base.estiff;
F2_0 = base.kstiff2*pw(e0);
F_iso_target = o0.Force(i0);
F_p1 = F_iso_target - F2_0;                                % p1 + passive share
fprintf('=== estiff scan on %s ===\n', SNAP);
fprintf('isometric: F=%.4g kPa, p2 share=%.1f%%, x_ref(p2-weighted)=%.5g um\n', ...
    F_iso_target, 100*F2_0/F_iso_target, dS*sum(xiso.*p2iso)/(dS*sum(p2iso)));
fprintf('baseline E_stairs(LSQ)=%.5f  featCost=%.4f\n', E0_lsq, stairsCost(o0.t, o0.Force/interp1(o0.t,o0.Force,tMov), vt, fD));

%% ---- scan ---------------------------------------------------------------
n = numel(EGRID);
R = struct('estiff',num2cell(EGRID), 'scale',[], 'Elsq',[], 'Cfeat',[], ...
           'parts',[], 'Fiso',[], 'Et',[], 'featCost',[], 'fM',[]);
featNames = {};
for i = 1:n
    e = EGRID(i);
    p = base;
    p.estiff = e;
    if UNCOMP
        sc = 1;
    else
        sc = (F_iso_target - F_p1) / (base.kstiff2 * pw(e));   % hold F2 constant
    end
    p.kstiff2   = base.kstiff2   * sc;
    p.kstiff2_n = base.kstiff2_n * sc;
    R(i).scale  = sc;

    % --- staircase prediction ---
    ps = p; ps.RunStairs = 1; ps.PlotEachSeparately=false; ps.ghostSave=''; ps.ghostLoad='';
    try
        [Elsq, os] = runStairsExperiment(ps);
        Fm_ss = interp1(os.t, os.Force, tMov, 'linear','extrap');
        [Cf, parts, fM] = stairsCost(os.t, os.Force/Fm_ss, vt, fD);
        R(i).Elsq = Elsq; R(i).Cfeat = Cf; R(i).parts = parts; R(i).fM = fM;
        R(i).Fiso = Fm_ss;
    catch ME
        fprintf('  stairs FAILED at e=%.3f: %s\n', e, ME.message);
        R(i).Elsq = NaN; R(i).Cfeat = NaN; R(i).Fiso = NaN;
    end

    % --- the objective this parameter set was actually bought with ---
    if DOFULL
        pf = p; pf.RunStairs = 0;
        try
            [Et, ~, fc, fnm] = evaluateBakersExp(pf.g, pf, true, false);
            R(i).Et = Et; R(i).featCost = fc; featNames = fnm;
        catch ME
            fprintf('  objective FAILED at e=%.3f: %s\n', e, ME.message);
            R(i).Et = NaN;
        end
    end
    fprintf('e=%5.3f  kstiff2x%7.3f  Fiso=%6.1f | stairs LSQ=%.5f feat=%.4f | Et=%.4f\n', ...
        e, sc, R(i).Fiso, R(i).Elsq, R(i).Cfeat, R(i).Et);
end

tag = 'estiff_scan'; if UNCOMP; tag = [tag '_uncomp']; end
save(fullfile(resdir,[tag '.mat']), 'R','EGRID','fD','featNames','SNAP','UNCOMP','E0_lsq');
fprintf('saved -> %s\n', fullfile(resdir,[tag '.mat']));
