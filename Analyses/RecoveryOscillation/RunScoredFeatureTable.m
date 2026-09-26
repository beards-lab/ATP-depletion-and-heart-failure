% RunScoredFeatureTable.m -- every term the 8 mM objective actually scores,
% with what it costs and how far from the data it is.
%
% The deck shows the best and worst features; this is the complete list, so
% "what is in the objective" is answerable without reading an fn string. Three
% kinds of entry are distinguished, because they are not comparable:
%
%   DATA        w = lambda*R/s^2 against the pooled target, with a
%               one-preparation-width deadband. chi-2 is meaningful: 1.0 = one
%               preparation-width from the data.
%   TRACE-LSQ   a raw RMSE with no measured spread, so no floor can be computed
%               and no chi-2 exists. These are the unprincipled ones.
%   GUARD       a boundary entry scored against a fixed range, not against
%               data at all. Cost 0 inside the range.
%
% Output -> results/lean5_features.txt  (and .mat)
%
%   matlab -batch "run('Analyses/RecoveryOscillation/RunScoredFeatureTable.m')"

clearvars -except TAG; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('TAG','var') || isempty(TAG); TAG = 'lean5_a_opt'; end

if isempty(gcp('nocreate')); parpool('local',5); end

p = getParams(); s = loadParams(TAG);
f = fieldnames(s); for i=1:numel(f); p.(f{i}) = s.(f{i}); end
p.RunForceVelocity=true; p.RunSlack=true; p.RunSlackPassive=true;
p.RunKtr=false; p.RunStairs=false; p.RunForceVelocityTime=false;
p.EvalFeatures=true; p.BreakOnODEUnstable=false; p.PlotEachSeparately=0;
p.PlotFeatureFitting=0; p.RunSlackSegments='AllPar'; p.MaxRunTime=300;
p.OptimizeOn='Feats'; p.UsePooledTarget=true; p.mods={}; p.g=ones(1,0);

params0 = p;                                                    %#ok<NASGU>
LoadData; RunBakersExp;
fm = features_model; fd = features_data;
fc = evalFeatureCost(fd, fm, p.fn);
PT = pooledTarget8mM();

fn = p.fn; n = numel(fn);
R = struct('entry',{},'name',{},'kind',{},'nEl',{},'w',{},'fl',{}, ...
           'cost',{},'share',{},'chi2',{});
for i = 1:n
    e  = fn{i};
    tk = strsplit(e,'|');
    nm = regexprep(tk{1}, '\[.*$','');

    isGuard = any(cellfun(@(t) ~isempty(regexp(t,'^-?[\d.]+-{1}-?[\d.]+$','once')), tk));
    fl = NaN; for t=1:numel(tk); if startsWith(tk{t},'fl:'); fl = str2double(tk{t}(4:end)); end; end
    w  = str2double(tk{end}); if isnan(w); w = 1; end

    % What makes a term a chi-2 is that it HAS a floor -- fl = FLOOR_CHI2*nEl*s^2
    % can only be written if some measured spread s exists. Whether that s came
    % from the pooled target or from feature_spread.mat (FV and the passive PS_*
    % terms, which are not pooled yet) changes the provenance, not the kind.
    % Classifying by "is it in pooledTarget8mM" instead would mislabel every
    % floored PS_/FV term as an unfloored RMSE.
    k = find(strcmp(PT.name, nm),1);
    if isGuard
        kind = 'GUARD';                 % boundary range, not scored against data
    elseif isfinite(fl) && ~isempty(k)
        kind = 'CHI2-pool';             % floored, spread from the pooled target
    elseif isfinite(fl)
        kind = 'CHI2-spread';           % floored, spread from feature_spread.mat
    else
        kind = 'UNFLOORED';             % raw RMSE, no measured spread, no floor
    end

    % chi-2 per scored element, from the features directly
    chi2 = NaN; nEl = NaN;
    if ~isempty(k) && isfield(fm, nm)
        m = double([fm.(nm)]); d = PT.val(k,:);
        q = min(numel(m), numel(d)); nEl = q;
        mu = mean(d(1:q),'omitnan');
        if isfinite(mu) && mu ~= 0
            chi2 = mean(((m(1:q)-d(1:q))/mu).^2,'omitnan') / PT.s(k)^2;
        end
    end

    R(end+1) = struct('entry',e,'name',nm,'kind',kind,'nEl',nEl,'w',w, ...
                      'fl',fl,'cost',fc(i),'share',0,'chi2',chi2);  %#ok<SAGROW>
end
tot = sum([R.cost]);
for i=1:numel(R); R(i).share = 100*R(i).cost/tot; end

[~,ord] = sort([R.cost],'descend');
L = {};
L{end+1} = sprintf('SCORED FEATURES -- %s   total %.4f   %d entries', TAG, tot, n);
L{end+1} = '';
L{end+1} = sprintf('%-22s %-10s %4s %10s %10s %9s %7s %8s', ...
    'feature','kind','nEl','weight','floor','cost','share%','chi2/el');
L{end+1} = repmat('-',1,92);
for i = ord
    r = R(i);
    L{end+1} = sprintf('%-22s %-10s %4s %10.4g %10s %9.4f %6.1f%% %8s', ...
        r.name, r.kind, iN(r.nEl), r.w, iF(r.fl), r.cost, r.share, iF2(r.chi2)); %#ok<SAGROW>
end
L{end+1} = repmat('-',1,92);
kinds = {'CHI2-pool','CHI2-spread','UNFLOORED','GUARD'};
for kk = 1:numel(kinds)
    sel = strcmp({R.kind}, kinds{kk});
    L{end+1} = sprintf('  %-12s  %2d entries   %7.4f   %5.1f %% of total', ...
        kinds{kk}, nnz(sel), sum([R(sel).cost]), 100*sum([R(sel).cost])/tot); %#ok<SAGROW>
end
sU = sum([R(strcmp({R.kind},'UNFLOORED')).cost]);
L{end+1} = '';
L{end+1} = sprintf(['WATCH: the UNFLOORED terms are %.1f %% of the objective. They carry no\n' ...
    'measured spread, so no deadband can be computed for them -- which means they\n' ...
    'are the only terms left that a search can keep buying after everything else\n' ...
    'has reached its floor. Their share RISES as the fit improves: they were ~18 %%\n' ...
    'at the lean5 seed and are %.0f %% here.'], 100*sU/tot, 100*sU/tot);

txt = strjoin(L, newline); disp(txt);
fid = fopen(fullfile(resDir,'lean5_features.txt'),'w');
fprintf(fid,'%s\n',txt); fclose(fid);
save(fullfile(resDir,'lean5_features.mat'),'R','tot','TAG');
fprintf('\nwrote %s\n', fullfile(resDir,'lean5_features.txt'));

function s = iN(v);  if isnan(v); s = '-'; else; s = sprintf('%d',v); end; end
function s = iF(v);  if isnan(v); s = '-'; else; s = sprintf('%.4g',v); end; end
function s = iF2(v); if isnan(v); s = '-'; else; s = sprintf('%.2f',v); end; end
