% RunLowAtpDiagnosis.m -- why is the 2 mM fit bad, and why did k_1 stop paying?
%
% Three questions, answered with numbers rather than argument:
%
%   Q1  WHERE is the 2 mM cost? Per-feature, with chi-2 against the pooled
%       2 mM target, so "the target is far" becomes a list of which features
%       are far and by how many preparation-widths.
%
%   Q2  HOW MANY 2 mM terms are FLOORED (paying exactly 0)? If most of the
%       objective is already inside its deadband there is nothing for an extra
%       gate to buy, and a set like P1 collapsing onto C0 is a property of the
%       OBJECTIVE, not evidence that k_1 is mechanistically idle.
%
%   Q3  Does k_1 have leverage at all? Scan it at fixed k2 and watch the
%       PER-FEATURE costs, not just the total. On the old core P1 (k2,k_1)
%       beat C0 (k2 alone) by 4.0; here they came out identical to four
%       decimals, and those two facts need reconciling before a campaign is
%       launched to re-find something that may not be findable.
%
%   matlab -batch "run('Analyses/AtpEffectHypotheses/RunLowAtpDiagnosis.m')"
%
% Output -> results/lowatp_diagnosis.txt

clearvars -except CORE PATP_K2; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
if ~exist('CORE','var');    CORE    = 'lean5_a_opt'; end
if ~exist('PATP_K2','var'); PATP_K2 = 0.44;          end

if isempty(gcp('nocreate')); parpool('local',5); end
LoadData;

p0 = getParams(loadParams(CORE), [], true, false);
p0.EvalFeatures=true; p0.OptimizeOn='Feats'; p0.BreakOnODEUnstable=false;
p0.RunKtr=false; p0.RunStairs=false; p0.RunForceVelocityTime=false;
p0.RunSlackPassive=false; p0.RunSlack=true; p0.RunForceVelocity=true;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=400;
p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.UsePooledTarget=true; p0.mods={}; p0.g=ones(1,0);
p0.UseAtpPowerGates=true; p0.MgATP_ref=8;

prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn = p0.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn = fn(~startsWith(prim(fn),'PS_'));          % passive not pooled at 2 mM
fn = rescaleFnToCondition(fn, 2);              % put the chi-2 on the 2 mM scale

p2 = p0; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8; p2.fn = fn;
p2.velocitytableonfile = 'protocol_03_27_2026_2mM_slack.mat';
r = checkFeatureCoverage(p2);
if ~r.ok
    bad = regexprep(unique([r.missingModel(:);r.missingData(:)])', '\s*\(in\s.*$','');
    fn = fn(~ismember(prim(fn), bad)); p2.fn = fn;
    fprintf(2,'dropped unscorable: %s\n', strjoin(bad,', '));
end

PT2 = pooledTarget2mM();
L = {};

%% ---- Q1/Q2 -------------------------------------------------------------
q = p2; q.pAtp_k2 = PATP_K2;
[fm, fd, fc] = iEval(q, fn);
tot = sum(fc); nm = prim(fn);

L{end+1} = sprintf('LOW-ATP DIAGNOSIS   core %s   pAtp_k2 %.3f', CORE, PATP_K2);
L{end+1} = sprintf('2 mM objective (pooled 2 mM target, 2 mM scale) = %.4f over %d terms', tot, numel(fn));
L{end+1} = '';
L{end+1} = sprintf('%-22s %10s %7s %10s %10s %9s', 'feature','cost','share','model','target','chi2/el');
L{end+1} = repmat('-',1,76);
[~,ord] = sort(fc,'descend');
nPay = 0;
for i = ord
    k = find(strcmp(PT2.name, nm{i}),1);
    mv = NaN; dv = NaN; c2 = NaN;
    if ~isempty(k) && isfield(fm, nm{i})
        m = double([fm.(nm{i})]); d = PT2.val(k,:);
        n = min(numel(m),numel(d)); mv = mean(m(1:n),'omitnan'); dv = mean(d(1:n),'omitnan');
        if isfinite(dv) && dv ~= 0
            c2 = mean(((m(1:n)-d(1:n))/dv).^2,'omitnan') / PT2.s(k)^2;
        end
    end
    if fc(i) > 1e-6; nPay = nPay + 1; end
    L{end+1} = sprintf('%-22s %10.4f %6.1f%% %10.4g %10.4g %9s', ...
        nm{i}, fc(i), 100*fc(i)/tot, mv, dv, iF(c2));                %#ok<SAGROW>
end
L{end+1} = repmat('-',1,76);
L{end+1} = sprintf('PAYING: %d of %d terms. %d are inside their deadband and contribute exactly 0.', ...
                   nPay, numel(fn), numel(fn)-nPay);

%% ---- Q3: does k_1 have leverage? ---------------------------------------
L{end+1} = '';
L{end+1} = 'Q3  k_1 scan at fixed k2 -- total AND the features it moves';
L{end+1} = sprintf('%8s %10s %10s', 'pAtp_k_1', 'total', 'delta');
grid = [-0.4 -0.2 0.2 0.4];   % 0 is the incumbent, already evaluated above
best = struct('v',0,'c',tot,'fc',fc);
for v = grid
    qq = q; qq.pAtp_k_1 = v;
    [~,~,fcv] = iEval(qq, fn);
    L{end+1} = sprintf('%8.2f %10.4f %+10.4f', v, sum(fcv), sum(fcv)-tot);  %#ok<SAGROW>
    if sum(fcv) < best.c; best = struct('v',v,'c',sum(fcv),'fc',fcv); end
end
L{end+1} = '';
if best.v ~= 0
    L{end+1} = sprintf('BEST k_1 = %+.2f  ->  %.4f  (%+.4f). Movers:', best.v, best.c, best.c-tot);
    d = best.fc - fc; [~,o] = sort(abs(d),'descend');
    for i = o(1:min(6,numel(o)))
        if abs(d(i)) < 1e-6; continue; end
        L{end+1} = sprintf('    %-22s %+8.4f  (%.4f -> %.4f)', nm{i}, d(i), fc(i), best.fc(i)); %#ok<SAGROW>
    end
else
    L{end+1} = 'k_1 buys NOTHING anywhere on this grid -- not a search failure, no leverage.';
end

txt = strjoin(L, newline); disp(txt);
fid = fopen(fullfile(resDir,'lowatp_diagnosis.txt'),'w');
fprintf(fid,'%s\n',txt); fclose(fid);
fprintf('\nwrote %s\n', fullfile(resDir,'lowatp_diagnosis.txt'));

%% ---- helpers -----------------------------------------------------------
function [fm, fd, fc] = iEval(p, fn)
params0 = p; params0.fn = fn;                                   %#ok<NASGU>
LoadData; RunBakersExp;
fm = features_model; fd = features_data;
fc = evalFeatureCost(fd, fm, fn);
end

function s = iF(v); if isnan(v); s = '-'; else; s = sprintf('%.2f',v); end; end
