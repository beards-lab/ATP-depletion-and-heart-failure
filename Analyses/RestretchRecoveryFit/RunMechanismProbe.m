% RunMechanismProbe.m — two questions the anatomy raised.
%
% Q1. The post-restretch window recovers +14 kPa NET as the difference of a
%     +22..+26 kPa active rise and a -8..-13 kPa PASSIVE decay. The post-slack
%     window has almost no passive component (-2.6 kPa against +70). Is that
%     large passive transient justified by the passive experiment, or is it a
%     Maxwell-element artefact contaminating the rate?
%     -> switch the Maxwell dashpot off and see what rsK does.
%
% Q2. Is the recovery rate-limited by the ATP-bound reservoir refilling? In
%     the post-restretch window PD (PuATP) moves the WRONG WAY first
%     (-0.84..-1.65 of its excursion) -- heads dump into it and are consumed
%     by re-attachment faster than it refills. If an obligatory slow
%     re-priming step gated re-attachment, the recovery would be clamped.
%     -> slow the hydrolysis/priming step kah and see whether rsK follows
%        while ktr does not.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

SNAP = 'params/rskR2_w025_opt.m';
base = getParams(loadParams(SNAP), [], true, false);
base.MaxRunTime = 600;
base.RunForceVelocity = 1; base.RunSlack = 1; base.RunSlackPassive = 1;
base.RunKtr = 0; base.RunStairs = 0; base.RunForceVelocityTime = 0;

variants = { ...
  'baseline',                struct(); ...
  'Maxwell dashpot OFF',     struct('UseMaxwellDashpot', 0); ...
  'eta_M x0.2',              struct('eta_M', base.eta_M*0.2); ...
  'eta_M x5',                struct('eta_M', base.eta_M*5); ...
  'kah x0.4  (slow priming)',struct('kah',  base.kah*0.4); ...
  'kah x0.6',                struct('kah',  base.kah*0.6); ...
  'kah x2',                  struct('kah',  base.kah*2); ...
  'ka  x0.5',                struct('ka',   base.ka*0.5); ...
  };

fprintf('%-26s %8s %8s %8s %8s %8s %9s %9s\n', 'variant', 'rsK', 'rsK x', ...
        'ktr', 'ktr x', 'steady', 'featL2', 'PS_hold');
res = struct();
for iv = 1:size(variants,1)
    nm = variants{iv,1};
    p  = base;
    f  = fieldnames(variants{iv,2});
    for j = 1:numel(f); p.(f{j}) = variants{iv,2}.(f{j}); end
    p = getParams(p, [], true, false);

    r = evalRs(p, struct('slackOnly', false));
    if ~r.ok
        fprintf('%-26s FAILED: %s\n', nm, r.err); continue;
    end
    [ct2,~,~] = evalFeatureCost(r.features_data, r.features_model, p.fn, 2);
    psh = NaN;
    if isfield(r.features_model,'PS_holdDecayRMSE')
        psh = mean(r.features_model.PS_holdDecayRMSE, 'omitnan');
    end
    fprintf('%-26s %8.1f %8.2f %8.1f %8.2f %8.1f %9.3f %9.2f\n', nm, ...
        mean(r.rsK_m,'omitnan'), r.rsK_ratio, mean(r.ktr_m,'omitnan'), ...
        r.ktr_ratio, mean(r.steady_m,'omitnan'), sum(ct2), psh);
    res.(matlab.lang.makeValidName(nm)) = r;
end

fprintf(['\ndata targets: rsK 43.7 (x1.00), ktr ~49 (x1.00), steady 77.3 kPa\n' ...
         'READ Q1: if "Maxwell OFF" barely moves rsK, the passive transient is a\n' ...
         '  passenger and the rate defect is kinetic. If rsK drops toward data,\n' ...
         '  the model is recovering force by the wrong physics.\n' ...
         'READ Q2: if kah-down slows rsK MORE than it slows ktr, an obligatory\n' ...
         '  slow re-priming step is a live candidate for the missing mechanism.\n']);

save(fullfile(here,'results','mechanism_probe.mat'), 'res', 'SNAP');
