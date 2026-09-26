function T = pooledKtrTarget(cond)
%POOLEDKTRTARGET 3-preparation pooled target for the ktr (Brenner) protocol.
%
%   T = pooledKtrTarget()       8 mM (default)
%   T = pooledKtrTarget(COND)   COND = 8 or 2 (mM MgATP); the 2 mM files are the
%                               *_velocitytable_ktr_2mM.mat of the same three days
%     T.name     feature names (extractKtrAttributes' outputs)
%     T.perPrep  one row per preparation (03/27, 04/03, 04/10)
%     T.val      pooled mean -- what the model is scored against
%     T.s        relative between-preparation SD, RAW (no moderation): the
%                chi-2 scale of each feature, i.e. w = lambda/s^2 and
%                fl = (1/3)*s^2 for a single-number feature
%
%   Computed from the curated data files through the SAME extractor the model
%   goes through, so re-curation can never leave a stale number behind (this
%   replaced a hard-coded k63 = 43.19). Cached per MATLAB process.
%
%   The model is simulated on one day's velocity table; the three days ran the
%   same manoeuvre, so the pooled mean is the fair target for it.
if nargin < 1 || isempty(cond); cond = 8; end
switch round(cond)
    case 8; suf = '';
    case 2; suf = '_2mM';
    otherwise; error('pooledKtrTarget:cond', 'no ktr target for %g mM', cond);
end
persistent PC
key = sprintf('c%d', round(cond));
if isstruct(PC) && isfield(PC, key); T = PC.(key); return; end
root = fileparts(fileparts(mfilename('fullpath')));
DAYS = {'03_27','04_03','04_10'};
F    = {'ktrP_k63','ktrP_k80','ktrP_F20'};
V    = nan(numel(DAYS), numel(F));
for d = 1:numel(DAYS)
    K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr%s.mat', DAYS{d}, suf)));
    f = extractKtrAttributes(K.datatable(:,1), K.datatable(:,3), K.velocitytable);
    for j = 1:numel(F); V(d,j) = f.(F{j}); end
end
P.name    = F;
P.days    = DAYS;
P.perPrep = V;
P.val     = mean(V, 1, 'omitnan');
P.s       = std(V, 0, 1, 'omitnan') ./ abs(P.val);
P.cond    = round(cond);
PC.(key)  = P;
T = P;
end
