% unusedMechScreen.m
% =========================================================================
% Screen currently-OFF mechanisms on the residual features. Base =
% params/params_2state_a2hop.m (2-state, cost ~2.95). Focus on restretch-transient
% and FV levers not yet tested this campaign:
%   UseViscoelasticSE (c_SE_visc)  - Kelvin-Voigt series viscosity during restretch
%                                    -> may shape peak1_dSL WITHOUT the kSE<->ktr weld
%                                    (it is vel>0-gated, so it doesn't touch redevelopment).
%   c_titin_visc (UseDynamicPassive already on) - velocity-dependent titin force on restretch.
%   UseStretchActivation (k_SA)    - boost ka during restretch (stretch recruitment).
%   UseA2AttachmentShift           - the built-in strained-p2 hop (relative d_actin shift).
%   UseWDetachment / UseVernierVelocity - FV-shape mechanisms.
%
% Writes Analyses/FeatureSensitivity/unusedMechScreen_report.txt
% Run:  cd(root); addpath(genpath('.')); unusedMechScreen
% =========================================================================
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));

SNAP   = 'params/params_2state_a2hop.m';
REPORT = fullfile('Analyses','FeatureSensitivity','unusedMechScreen_report.txt');
MRT    = 90;

params0 = getParams(); run(SNAP); FN = params0.fn;

L = {};
L{end+1} = {'BASE a2hop',      struct()};
L{end+1} = {'viscoSE .02',     struct('UseViscoelasticSE',1,'c_SE_visc',0.02)};
L{end+1} = {'viscoSE .10',     struct('UseViscoelasticSE',1,'c_SE_visc',0.10)};
L{end+1} = {'titinvisc 2',     struct('c_titin_visc',2)};
L{end+1} = {'titinvisc 10',    struct('c_titin_visc',10)};
L{end+1} = {'stretchAct .3',   struct('UseStretchActivation',1,'k_SA',0.3)};
L{end+1} = {'stretchAct 1',    struct('UseStretchActivation',1,'k_SA',1)};
L{end+1} = {'A2shift',         struct('UseA2AttachmentShift',1)};
L{end+1} = {'Wdetach',         struct('UseWDetachment',1)};
L{end+1} = {'vernier',         struct('UseVernierVelocity',1)};

fid = fopen(REPORT,'w');
fprintf(fid, '==== unusedMechScreen on %s ====\n', SNAP);
fprintf(fid, 'Data: peak1_y~96, peak1_dSL~0.025, vall_y~71, peak2~81, ovrsht~1.3, steady~[80..64], FV~[.92 .66 .32 .11], ktr~49\n\n');
fprintf(fid, '%-15s %7s %8s %7s %7s %7s %7s %6s | %-22s %7s\n', ...
    'mechanism','peak1y','pk1dSL','vall_y','peak2','ovrsht','steady','ktr','FV(.5,1,2,4)','COST');

for i = 1:numel(L)
    label = L{i}{1}; extra = L{i}{2};
    try
        [tc,~,fm,~] = costOfSnap(SNAP, FN, MRT, extra);
        g = @(n) lm(fm,n); fv = lfv(fm);
        fprintf(fid, '%-15s %7.1f %8.4f %7.1f %7.1f %7.2f %7.1f %6.1f | %-22s %7.3f\n', ...
            label, g('peak1_y'), g('peak1_dSL'), g('vall_y'), g('peak2'), g('ovrsht_dy'), ...
            g('steady'), g('ktr'), sprintf('%.2f %.2f %.2f %.2f',fv(1),fv(2),fv(3),fv(4)), tc);
    catch e
        fprintf(fid, '%-15s  ERROR: %s\n', label, e.message);
    end
    fprintf('done %d/%d: %s\n', i, numel(L), label);
end
fclose(fid); type(REPORT); disp('DONE unusedMechScreen');

function v = lm(fm,n); if isfield(fm,n)&&~isempty(fm.(n)); v=mean(fm.(n),'omitnan'); else; v=NaN; end; end
function fv = lfv(fm); fv=[NaN NaN NaN NaN]; if isfield(fm,'FV_fnorm')&&numel(fm.FV_fnorm)>=5; x=fm.FV_fnorm(:)'; fv=x(2:5); end; end
