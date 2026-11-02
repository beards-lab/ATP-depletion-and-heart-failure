% RunAtpLeverSensitivity.m -- on top of the k2 gate (R2 x0.54 at 2 mM), which
% SECOND rate could carry the rest of the ATP effect? One lever at a time,
% +-10 % at 2 mM only, scored exactly like RunRefit9K2Sweep.
%
%   matlab -batch "run('Analyses/LowATP_Refit9K2Sweep/RunAtpLeverSensitivity.m')"
%
% Candidate levers, ordered by chemical plausibility (tiers from
% Analyses/AtpEffectHypotheses/RunPlausibleSets.m):
%   A   ksrd2sr  SRX-ADP -> SRX-ATP, nucleotide EXCHANGE: free ATP is a substrate
%   A   ksr0     entry into SRX from the ATP-bound state (supply scales with ATP)
%   AB  kmsrd    SRX-ADP exit to the primed pool (weakly arguable)
%   B   A0       registration-availability floor (thick-filament activation;
%                ATP/SRX-IHM coupling is a hypothesis, not chemistry)
%   B   tau_reg  registration-availability recovery time (same argument)
%   C   k_1      reverse stroke = Pi rebinding (Pi contamination of the 2 mM mix)
%   C   k1       power stroke = Pi release (same)
%   X   kah      hydrolysis of BOUND ATP -- control (see conclusions)
%   X   ka       actin attachment, no nucleotide exchange -- control
%   X   kd       pre-stroke detachment -- control
% Each lever is applied through the merged power gates (exactly 1 at 8 mM):
% factor f at 2 mM  <=>  pAtp_<lever> = log(f)/log(2/8).

if ~exist('CORE','var');   CORE   = 'robust1_seed'; end
if ~exist('G_K2','var');   G_K2   = 0.54; end
if ~exist('LEVERS','var')
    LEVERS = { 'ksrd2sr','A'; 'ksr0','A'; 'kmsrd','AB'; 'A0','B'; 'tau_reg','B'; ...
               'k_1','C'; 'k1','C'; 'kah','X'; 'ka','X'; 'kd','X' };
end
if ~exist('FACT','var');   FACT   = [0.9 1.1]; end
if ~exist('BASE','var');   BASE   = struct(); end   % extra 2 mM gates on top of k2 (e.g. pAtp_ksr0)
if ~exist('TAGL','var');   TAGL   = 'levers'; end    % output file prefix
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resL = fullfile(root,'Analyses','LowATP_Refit9K2Sweep','results');
diary(fullfile(resL,[TAGL '.log'])); diary on;
fprintf('lever sensitivity on %s, k2 gate g = %.2f, started %s\n', CORE, G_K2, datestr(now));

clL = parcluster('local');
jslL = fullfile(tempdir, sprintf('mljob_lev_%d', feature('getpid')));
if ~isfolder(jslL); mkdir(jslL); end
clL.JobStorageLocation = jslL;
if isempty(gcp('nocreate')); parpool(clL, 5); end

%% ---- the two legs: identical set-up to RunRefit9K2Sweep ---------------------
pL = getParams(loadParams(CORE), [], true, false);
pL.RunForceVelocity=true; pL.RunSlack=true; pL.RunStairs=false; pL.RunKtr=true;
pL.RunForceVelocityTime=false; pL.RunSlackPassive=false; pL.EvalFeatures=true;
pL.BreakOnODEUnstable=false; pL.PlotEachSeparately=0; pL.PlotFeatureFitting=0;
pL.RunSlackSegments='AllPar'; pL.MaxRunTime=400; pL.OptimizeOn='Feats';
pL.ghostSave=''; pL.ghostLoad=''; pL.SaveBest=0; pL.mods={}; pL.g=ones(1,0);
pL.recalculateDataFeats=false; pL.UsePooledTarget=true;
pL.UseAtpPowerGates=true; pL.MgATP_ref=8; pL.UseAtpDetachR2=0; pL.UseAtpKmsrd=0;

primL = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn8L = pL.fn; [~,iu] = unique(fn8L,'stable'); fn8L = fn8L(iu);
fn8L = fn8L(~startsWith(primL(fn8L),'PS_'));
fn2L = rescaleFnToCondition(fn8L, 2);
fn2L = fn2L(~(startsWith(fn2L,'FV_normpowerAvg') | startsWith(fn2L,'steady[1]')));

p8L = pL; p8L.MgATP = 8; p8L.fn = fn8L;
p8L.velocitytableonfile     = 'protocol_03_27_2026_8mM_slack.mat';
p8L.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr.mat';
p2L = p8L; p2L.MgATP = 2; p2L.fn = fn2L; p2L.RunForceVelocity = false;
p2L.velocitytableonfile     = 'protocol_03_27_2026_2mM_slack.mat';
p2L.ktr_velocitytableonfile = 'protocol_03_27_2026_velocitytable_ktr_2mM.mat';
p2L.pAtp_k2 = log(G_K2)/log(1/4);
fbL = fieldnames(BASE); for jb = 1:numel(fbL); p2L.(fbL{jb}) = BASE.(fbL{jb}); end

% relevance weights, as in the sweep
PT8L = pooledTarget8mM(); KT8L = pooledKtrTarget(8); KT2L = pooledKtrTarget(2);
nmL = primL(fn2L); aL = ones(1,numel(fn2L));
for i = 1:numel(fn2L)
    if contains(fn2L{i}, ',') || ~isempty(regexp(fn2L{i}, '\|-?[\d.]+-{1}-?[\d.]+(\||$)', 'once')); continue; end
    nm = nmL{i}; if strcmp(nm,'t0_crossing'); nm = 't0'; end
    j = find(strcmp(PT8L.name, nm), 1); z = NaN;
    if ~isempty(j); z = PT8L.z(j);
    else
        j = find(strcmp(KT8L.name, nm), 1);
        if ~isempty(j); r = KT2L.val(j)/KT8L.val(j); sr = r*sqrt(KT2L.s(j)^2+KT8L.s(j)^2); z = abs(r-1)/sqrt(KT8L.s(j)^2+sr^2); end
    end
    if isfinite(z); aL(i) = min(z,3)/3; end
end
FV2L = fv2Target(); zFVL = abs(FV2L.r-1)/sqrt(FV2L.s8r^2 + FV2L.sr^2);
scoreL = @(c, fv) deal(sum(c) + ((fv-FV2L.m2)/FV2L.s2)^2, ...
                       sum(c(:)'.*aL) + min(zFVL,3)/3*((fv-FV2L.m2)/FV2L.s2)^2);

%% ---- reference runs -----------------------------------------------------------
t0L = tic;
[fm8L, fd8L] = iRunL(p8L); fv8L = fv2Model(p8L);
[fmB, fdB]   = iRunL(p2L); fvB = fv2Model(p2L);
cB = evalFeatureCost(fdB, fmB, fn2L); [sPB, sRB] = scoreL(cB, fvB);
fprintf('reference g = %.2f: plain %.4f  rel %.4f  FV2 %.3f  (%.0f s)\n', G_K2, sPB, sRB, fvB, toc(t0L));
save(fullfile(resL,[TAGL '_ref.mat']), 'fm8L','fd8L','fv8L','fmB','fdB','fvB','cB','sPB','sRB','fn2L','aL','FV2L');

%% ---- one lever at a time ------------------------------------------------------
RL = struct('lever',{},'tier',{},'f',{},'fm',{},'cost',{},'fv2',{},'score_plain',{},'score_rel',{},'sec',{});
for iL = 1:size(LEVERS,1)
    for fL = FACT
        tL = tic;
        q = p2L; q.(['pAtp_' LEVERS{iL,1}]) = log(fL)/log(1/4);
        [fmL, fdL] = iRunL(q); fvL = fv2Model(q);
        cL = evalFeatureCost(fdL, fmL, fn2L); [sP, sR] = scoreL(cL, fvL);
        RL(end+1) = struct('lever',LEVERS{iL,1},'tier',LEVERS{iL,2},'f',fL,'fm',fmL,'cost',cL, ...
                           'fv2',fvL,'score_plain',sP,'score_rel',sR,'sec',toc(tL)); %#ok<SAGROW>
        fprintf('%-8s x%.1f: rel %.4f (%+.4f)  plain %.4f (%+.4f)  FV2 %.3f  (%.0f s)\n', LEVERS{iL,1}, fL, ...
                sR, sR-sRB, sP, sP-sPB, fvL, RL(end).sec);
        save(fullfile(resL,[TAGL '.mat']), 'RL','LEVERS','FACT','G_K2','CORE','BASE');
    end
end
fprintf('levers done in %.0f s\n', toc(t0L));
if strcmp(TAGL,'levers'); SummarizeLevers; end
fprintf('DONE\n'); diary off;

function [fm, fd] = iRunL(p)
    params0 = p;                                               %#ok<NASGU>
    LoadData; RunBakersExp;                                    %#ok<*ASGLU>
    fm = features_model; fd = features_data;
end
