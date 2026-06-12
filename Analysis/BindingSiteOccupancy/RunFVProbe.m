% RunFVProbe
% =====================================================================
% Probe FV-shaping mechanisms on the BoundedFit basis. The model FV is too
% steep (normalised force falls too fast during shortening). Candidates:
%   mu                          : viscous half-sarcomere drag (lower => flatter FV)
%   PieceWiseStrainDep2Params   : R2 detachment burst at NEGATIVE strain
%                                 (lower => heads stay attached during shortening)
%   global occupancy + kstiff2  : occupancy at MATCHED isometric force (fair FV test)
%
% Captures, per config: total feature cost, the FV_fnorm sub-cost (cost_vec(1),
% weighted), the normalised FV vector, and collateral slack/physiology metrics.
% Saves fv_probe_results.mat.
%
% Run: cd repo root; addpath(genpath('.')); RunFVProbe
% =====================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(root));

%% ---- build BoundedFit base (mirrors DriverBoundedFit) ----
params0 = getParams();
run(fullfile(root, 'params', 'ModelOptParams_TL3_iter_17.m'));
params0.mods = {}; params0.g = [];
params0 = getParams(params0, [], true, false);
params0.RunForceVelocity=true; params0.RunSlack=true; params0.RunKtr=false; params0.RunStairs=false;
params0.RunForceVelocityTime=false; params0.RunForceLengthEstim=false;
params0.EvalFeatures=true; params0.BreakOnODEUnstable=false; params0.MaxRunTime=120;
params0.RunSlackSegments='All';
try
    if isempty(gcp('nocreate')); parpool('Threads',5); end
    params0.RunSlackSegments='AllPar';
catch ME
    warning(ME.identifier,'no pool (%s)',ME.message);
end
params0.fn = boundedOutputFn(0.001);
params0.ksr2srd="=kah"; params0.ksrd2sr="=kah/10"; params0.kah=100; params0.kamh="=kah/10";
params0.ksr0=500; params0.kmsr=0; params0.kmsrd=40; params0.ksrd=0;
params0.ksr2srd="=kah"; params0.ksr0=190; params0.kmsrd=60;
params0.kstiff2=6.9794e+04*0.5;
params0.UseSuperRelaxed=true; params0.UseSuperRelaxedADP=true;
params0.UseVelGaussAttachment=false; params0.v_att_amplitude=0.5; params0.v_att_center=0; params0.v_att_sigma=0.15;
params0.FV_velocities = -[0,0.5,1,2,4];
run(fullfile(root,'params','params_reseeded.m'));
params0.justPlotStateTransitionsFlag=false; params0.PlotFeatureFitting=false;
params_base = params0;

negmask  = params_base.PieceWiseStrainDep2X < 0;          % negative-strain (shortening) knots
zeromask = abs(params_base.PieceWiseStrainDep2X) < 1e-9;  % knot at s2=0 (onset of shortening)
fprintf('R2 neg-strain knots (s2<0): vals = %s ; knot@0 = %s\n', ...
    num2str(params_base.PieceWiseStrainDep2Params(negmask),'%.3g '), ...
    num2str(params_base.PieceWiseStrainDep2Params(zeromask),'%.3g '));

%% ---- config list (mu fixed at base: low mu destabilises the ODE) ----
C = {};
C{end+1} = struct('name','base');
C{end+1} = struct('name','R2reshape(k0=5,neg0.4)','r2neg',0.4,'k0',5);
C{end+1} = struct('name','occ0.30 +kstiff1.6','occ',true,'cap',0.30,'kstiffx',1.6);
C{end+1} = struct('name','occ0.30 kst1.6 +R2reshape','occ',true,'cap',0.30,'kstiffx',1.6,'r2neg',0.4,'k0',5);

R = struct([]);
for k = 1:numel(C)
    c = C{k};
    params0 = params_base;
    if isfield(c,'r2neg'), params0.PieceWiseStrainDep2Params(negmask) = params_base.PieceWiseStrainDep2Params(negmask)*c.r2neg; end
    if isfield(c,'k0'),    params0.PieceWiseStrainDep2Params(zeromask) = c.k0; end
    if isfield(c,'occ') && c.occ
        params0.UseGlobalOccupancySaturation=true; params0.OccupancyForm='linear'; params0.P_bound_max=c.cap;
    end
    if isfield(c,'kstiffx'), params0.kstiff2 = params_base.kstiff2*c.kstiffx; end

    features_data=struct(); features_model=struct();
    fprintf('\n===== %s =====\n', c.name); t0=tic;
    figure(91000+k); clf;
    RunBakersExp;
    cv = plotFeatures(features_data, features_model, features_model, params0.fn, params0);
    R(k).name=c.name; R(k).total=sum(cv,'omitnan'); R(k).FVcost=cv(1);
    R(k).FV_fnorm = getv(features_model,'FV_fnorm');
    R(k).XTOR=getf(features_model,'XTOR'); R(k).attached_ss=getf(features_model,'attached_ss');
    R(k).A=getm(features_model,'A'); R(k).steady=getm(features_model,'steady'); R(k).peak1_y=getm(features_model,'peak1_y');
    R(k).ktr=getm(features_model,'ktr');
    fprintf('%s: total=%.2f  FVcost=%.2f  XTOR=%.3g  att=%.3f  A=%.1f steady=%.1f ktr=%.3g (%.0fs)\n', ...
        c.name, R(k).total, R(k).FVcost, R(k).XTOR, R(k).attached_ss, R(k).A, R(k).steady, R(k).ktr, toc(t0));
end

%% ---- summary ----
if isfield(features_data,'FV_fnorm'), dtarget = features_data.FV_fnorm(:)'; else, dtarget = []; end
fprintf('\n\n============ FV PROBE SUMMARY ============\n');
fprintf('%-24s %7s %7s %6s %6s %6s %6s %6s\n','config','total','FVcost','XTOR','att','A','steady','ktr');
for k=1:numel(R)
    fprintf('%-24s %7.2f %7.2f %6.3g %6.3f %6.1f %6.1f %6.3g\n', ...
        R(k).name,R(k).total,R(k).FVcost,R(k).XTOR,R(k).attached_ss,R(k).A,R(k).steady,R(k).ktr);
end
fprintf('\nNormalised FV (model) vs data target:\n');
if ~isempty(dtarget), fprintf('  %-22s %s\n','DATA', num2str(dtarget,'%6.3f ')); end
for k=1:numel(R), fprintf('  %-22s %s\n', R(k).name, num2str(R(k).FV_fnorm,'%6.3f ')); end
save(fullfile(fileparts(mfilename('fullpath')),'fv_probe_results.mat'),'R','C');
fprintf('\nSaved fv_probe_results.mat\n');

%% helpers
function v=getf(fm,n), if isfield(fm,n)&&~isempty(fm.(n)), x=fm.(n)(:); v=x(find(~isnan(x),1)); if isempty(v),v=NaN;end, else, v=NaN; end, end
function v=getm(fm,n), if isfield(fm,n)&&~isempty(fm.(n)), v=mean(fm.(n),'omitnan'); else, v=NaN; end, end
function v=getv(fm,n), if isfield(fm,n)&&~isempty(fm.(n)), v=fm.(n)(:)'; else, v=NaN; end, end
