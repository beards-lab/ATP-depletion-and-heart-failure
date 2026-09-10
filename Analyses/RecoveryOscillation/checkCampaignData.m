function ok = checkCampaignData(seedTag)
%CHECKCAMPAIGNDATA  Does this machine have the data the campaign will open?
%
%   ok = checkCampaignData()              % params/rigaware_seed.m
%   ok = checkCampaignData('pooled8_seed')
%
% Run this on the TARGET machine before launching. `data/` is gitignored
% (`.gitignore` has both `*.mat` and `data/`), so a fresh clone has the code
% and the fit target but none of the 2026 recordings, and the failure mode is
% a 30-hour job that dies on its first evaluation.
%
% The required list is derived from the seed's own `*onfile` fields, not
% hardcoded, so it follows the snapshot if the protocol ever changes.
%
% WHAT THE CAMPAIGN ACTUALLY OPENS. RunOptimSpreadW sets
%   RunSlack = true, RunSlackPassive = true, RunForceVelocity = true,
%   RunKtr = false, RunStairs = false, RunForceVelocityTime = false
% (Analyses/FeatureWeighting/RunOptimSpreadW.m, the block above the
% coverage check), so only the slack and passive files are read. The F-V and
% k_tr REFERENCE data are hardcoded arrays in Model/LoadData.m -- a script,
% no file I/O -- and the fit target itself is frozen in the committed
% Model/pooledTarget8mM.m. Neither needs anything from data/.
%
% The k_tr and staircase protocol files are listed as OPTIONAL because they
% are predictions: nothing in the objective touches them, but
% Analyses/OptimResultReport/RunFitReport.m plots them afterwards.
%
% See also pooledTarget8mM, RunOptimRigAware, RunRigAwareReport

if nargin < 1 || isempty(seedTag); seedTag = 'rigaware_seed'; end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
p    = getParams(loadParams(seedTag), [], false, false);

% field name -> what reads it, and whether the campaign needs it
SPEC = { 'velocitytableonfile'         , 'slack protocol + trace + curated features', true
         'passive_velocitytableonfile' , 'PNB+Mava passive run'                     , true
         'ktr_velocitytableonfile'     , 'k_tr protocol (prediction only)'           , false
         'stairs_velocitytableonfile'  , 'staircase protocol (prediction only)'      , false };

fprintf('\ncheckCampaignData -- seed %s\n', seedTag);
fprintf('root %s\n\n', root);

ok = true; missReq = {}; needReq = {}; needOpt = {};
fprintf('%-9s %-48s %10s  %s\n', 'role', 'file', 'bytes', 'what reads it');
fprintf('%s\n', repmat('-', 1, 112));
for i = 1:size(SPEC,1)
    f = SPEC{i,1};
    if ~isfield(p, f) || isempty(p.(f)); continue; end
    rel  = fullfile('data', p.(f));
    full = fullfile(root, rel);
    req  = SPEC{i,3};
    d    = dir(full);
    if isempty(d)
        if req; ok = false; missReq{end+1} = rel; end          %#ok<AGROW>
        fprintf('%-9s %-48s %10s  %s   <-- MISSING\n', ...
            ternary(req,'REQUIRED','optional'), rel, '-', SPEC{i,2});
    else
        if req; needReq{end+1} = rel; else; needOpt{end+1} = rel; end  %#ok<AGROW>
        fprintf('%-9s %-48s %10d  %s\n', ...
            ternary(req,'REQUIRED','optional'), rel, d.bytes, SPEC{i,2});
    end
end

% The target is code, not data -- but check it, because a stale clone is the
% other way this goes wrong, and its header records which preparations and
% which scored segments it was built from.
fprintf('\n');
if isempty(which('pooledTarget8mM'))
    fprintf('TARGET    Model/pooledTarget8mM.m           MISSING -- clone is incomplete\n');
    ok = false;
else
    PT = pooledTarget8mM();
    seg = 'all 5 cycles';
    if isfield(PT,'seg')
        k = find(strcmp(PT.name,'rsK80'), 1);
        if ~isempty(k); seg = sprintf('rsK80 spread on cycles %s', mat2str(PT.seg{k})); end
    end
    fprintf('TARGET    Model/pooledTarget8mM.m (committed)  %d preps %s, %s\n', ...
        PT.nP, strjoin(PT.days,'/'), seg);
end
fprintf('REFERENCE Model/LoadData.m  F-V and k_tr arrays are inline -- no data file\n');

fprintf('\n');
if ok
    fprintf('PASS -- the campaign can run here.\n\n');
    fprintf('Must travel for the CAMPAIGN (%d):\n', numel(needReq));
    fprintf('  %s\n', needReq{:});
    if ~isempty(needOpt)
        fprintf('Add these only if you also want RunFitReport''s prediction panels (%d):\n', ...
                numel(needOpt));
        fprintf('  %s\n', needOpt{:});
    end
else
    fprintf('FAIL -- %d required file(s) missing:\n', numel(missReq));
    fprintf('  %s\n', missReq{:});
    fprintf('\nCopy them from the acquisition machine into %s\n', fullfile(root,'data'));
end
fprintf('\n');
end

function o = ternary(c, a, b); if c; o = a; else; o = b; end; end
