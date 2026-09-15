% RunPlausibleSets.m -- which BIOLOGICALLY DEFENSIBLE transitions can carry the
% low-ATP phenotype, at 3-4 degrees of freedom?
%
% ---------------------------------------------------------------------------
% THE PROBLEM WITH H1
%
% H1 (k2, kah, kamh) fits well but two of its three gates are not defensible.
% Mapping the model's rates onto the actual chemistry:
%
%   k2    p2 -> PT   ADP release + ATP BINDING + dissociation.  Free [ATP] is a
%                    substrate. This is the canonical ATP step and the only one
%                    where [ATP] enters the rate law directly.        PLAUSIBLE
%   kah   PT -> PD   hydrolysis of ALREADY-BOUND ATP. Free [ATP] does not appear;
%                    the nucleotide is already in the pocket.       NOT PLAUSIBLE
%   kamh  PD -> PT   reverse of the same.                           NOT PLAUSIBLE
%   ka    PD -> p1   actin binding by M.ADP.Pi. No nucleotide exchange.   NOT
%   kd    p1 -> PD   detachment of the pre-stroke head.                   NOT
%   k1    p1 -> p2   power stroke, gated by Pi RELEASE.   not ATP; yes Pi
%   k_1   p2 -> p1   reverse stroke, driven by Pi REBINDING. not ATP; yes Pi
%
% SRX branch (UseSuperRelaxedADP is on in this core; kmsr = 0 so the SR->PT exit
% is dead and only the SR-ADP branch is live):
%
%   ksr0     PT -> P_SR    entry into SRX FROM THE ATP-BOUND state. Requires an
%                          ATP-loaded head, so scarcity reduces the supply. PLAUSIBLE
%   ksr2srd  P_SR -> P_SRD hydrolysis inside SRX. Bound nucleotide.        NOT
%   ksrd2sr  P_SRD -> P_SR nucleotide EXCHANGE: ADP out, ATP IN. Free [ATP] is a
%                          substrate.                                  PLAUSIBLE
%   kmsrd    P_SRD -> PD   SRX-ADP exit to the primed pool.        weakly arguable
%
% CONTAMINATION. The 2 mM and 8 mM solutions differ in more than [ATP]:
%   ADP  slows ADP release, so the p2 -> PT step slows. In this model that step
%        is LUMPED into k2, so an ADP effect is mathematically DEGENERATE with a
%        low-ATP effect. It cannot be separated here -- stated, not tested.
%   Pi   rebinds and reverses the stroke: k_1 up, k1 down. Pi PREDICTS LOWER
%        force. The data show force x1.34 HIGHER at 2 mM, so Pi contamination
%        works against the observed phenotype; if present it is being
%        over-compensated by something else. Tested explicitly (sets P1-P3).
%
% ---------------------------------------------------------------------------
% WHAT THIS SCRIPT DOES
% Fits each candidate set on the same 2 mM objective as the hypothesis ladder
% (19 features, force-velocity included), core frozen, and ranks them. 8 mM is
% invariant by construction and asserted once.
%
% Run: matlab -batch "run('Analyses/AtpEffectHypotheses/RunPlausibleSets.m')"
%      matlab -batch "BUDGET_MIN=12; run('...')"     % per-set budget

clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root));

if ~exist('CORE','var');       CORE       = 'spreadw_iter_018'; end
if ~exist('BUDGET_MIN','var'); BUDGET_MIN = 12;   end
if ~exist('ACCEPT','var');     ACCEPT     = 0.32; end
if ~exist('NWORKERS','var');   NWORKERS   = 5;    end
if ~exist('TAG','var');        TAG        = 'plaus'; end
if ~exist('RESUME','var');     RESUME     = true; end

resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end
STATE  = fullfile(resDir,[TAG '_state.mat']);

%% ---- candidate sets -------------------------------------------------------
% col1 id, col2 gates, col3 tier, col4 one-line justification
SETS = {
 'C0'  {'k2'}                                'A'  'canonical ATP binding alone (irreducible minimum)'
 'C1'  {'k2','ksrd2sr'}                      'A'  'both true ATP-binding steps: dissociation + SRX nucleotide exchange'
 'C2'  {'k2','ksr0'}                         'A'  'ATP binding + SRX entry from the ATP-bound state'
 'C3'  {'k2','ksr0','ksrd2sr'}               'A'  'all three steps where free ATP is a substrate'
 'C4'  {'k2','ksrd2sr','kmsrd'}              'AB' 'ATP steps + SRX-ADP exit'
 'C5'  {'k2','ksr0','ksrd2sr','kmsrd'}       'AB' 'full defensible SRX branch (4 DOF)'
 'C6'  {'k2','kmsrd'}                        'AB' 'ATP binding + SRX-ADP exit only'
 'P1'  {'k2','k_1'}                          'C'  'ATP + Pi contamination (reverse stroke)'
 'P2'  {'k2','k1','k_1'}                     'C'  'ATP + Pi both directions of the stroke'
 'P3'  {'k2','ksrd2sr','k_1'}                'AC' 'ATP steps + Pi reverse stroke'
 'H1'  {'k2','kah','kamh'}                   'X'  'REFERENCE: the incumbent, hydrolysis gates (not defensible)'
 'X1'  {'k2','ka'}                           'X'  'CONTROL: attachment gate (no nucleotide exchange)'
 'X2'  {'k2','eta_M','kSE_M'}                'X'  'CONTROL: what H4 actually used -- mechanical, not chemistry'
};
ALLG = unique([SETS{:,2}]);

%% ---- objective ------------------------------------------------------------
LoadData;
if exist('parpool','file')==2 && isempty(gcp('nocreate')); parpool('local',NWORKERS); end

p0 = getParams(loadParams(CORE), [], true, false);
p0.UseAtpJointFit=false; p0.UseAtpPowerGates=true; p0.MgATP_ref=8;
p0.UseAtpDetachR2=0; p0.UseAtpKmsrd=0; p0.UseStrainDepAtpGate=0;
p0.UseR2Ceiling=0; p0.UseD0State=0;
p0.EvalFeatures=true; p0.BreakOnODEUnstable=false; p0.OptimizeOn='Feats';
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=400;
p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0; p0.PlotFullscreen=0;
p0.RunKtr=false; p0.RunStairs=false; p0.RunForceVelocityTime=false;
p0.RunSlackPassive=false; p0.RunSlack=true; p0.RunForceVelocity=true;
p0.mods={}; p0.g=ones(1,0);

prim = @(c) cellfun(@(s) regexprep(strtok(s,'|'),'\[.*$',''), c, 'uni', 0);
fn2 = p0.fn; [~,iu]=unique(fn2,'stable'); fn2=fn2(iu);
fn2 = fn2(~startsWith(prim(fn2),'PS_'));

p2 = p0; p2.MgATP=2; p2.MgADP=0; p2.Pi=0.8;
p2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
p2.fn = fn2;
rep = checkFeatureCoverage(p2);
if ~rep.ok
    bad = unique([rep.missingModel(:);rep.missingData(:)])';
    fprintf(2,'dropping unscorable: %s\n', strjoin(bad,', '));
    fn2 = fn2(~ismember(prim(fn2),bad)); p2.fn = fn2;
    assert(checkFeatureCoverage(p2).ok);
end
for i=1:numel(ALLG); p2.(['pAtp_' ALLG{i}]) = 0; end

cost = @(q) iCost(q, fn2, ATP_c, Data_ATP);
c0 = cost(p2);
fprintf('core %s | %d features | H0 (no gates) = %.4f\n', CORE, numel(fn2), c0);

q8 = p2; q8.MgATP=8; q8.Pi=0;
q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
q8g = q8; for i=1:numel(ALLG); q8g.(['pAtp_' ALLG{i}]) = 0.5; end
d8 = abs(cost(q8g) - cost(q8));
fprintf('8 mM invariance |d| = %.3e  (MUST be 0)\n\n', d8);
assert(d8 < 1e-9, 'gates leak into 8 mM');

%% ---- fit each set ---------------------------------------------------------
OUT = struct('id',{},'gates',{},'tier',{},'why',{},'dof',{},'c',{},'exps',{},'conv',{});
if RESUME && isfile(STATE); L=load(STATE); OUT=L.OUT; fprintf('resumed %d set(s)\n',numel(OUT)); end

for k = 1:size(SETS,1)
    id = SETS{k,1};
    if any(strcmp({OUT.id}, id)); continue; end
    g = SETS{k,2};
    fprintf('=== %-3s [%s] %d DOF: %s ===\n', id, SETS{k,3}, numel(g), strjoin(g,', '));

    % A dead parallel pool makes every eval fail and then hang (one set burned
    % 279 min on a single truncated sweep). Check and rebuild before each set.
    pl = gcp('nocreate');
    if exist('parpool','file')==2 && (isempty(pl) || pl.NumWorkers < 2)
        if ~isempty(pl); delete(pl); end
        fprintf('   [pool was dead -- restarting]\n'); parpool('local', NWORKERS);
    end

    % WARM START from the best already-fitted SUBSET of this set. Seeding every
    % set from all-zeros means a superset that runs out of budget can score WORSE
    % than its own subset, which is not a fit result but an artifact -- exactly
    % what happened to C3 (28.13 against C2's 9.04 on a strict superset).
    q = p2; c = c0; seedFrom = 'H0';
    if ~isempty(OUT)
        sub = find(arrayfun(@(o) all(ismember(o.gates, g)), OUT));
        if ~isempty(sub)
            [~,bi] = min([OUT(sub).c]); b = OUT(sub(bi));
            fb = fieldnames(b.exps);
            for i = 1:numel(fb); q.(['pAtp_' fb{i}]) = b.exps.(fb{i}); end
            c = cost(q); seedFrom = b.id;
        end
    end
    fprintf('   seed %.4f (from %s)\n', c, seedFrom);

    tR=tic; improved=true; sweep=0; full=false;
    while improved && toc(tR) < BUDGET_MIN*60
        improved=false; sweep=sweep+1; full=true;
        steps=[0.15 -0.15 0.30 -0.30]; if sweep>1; steps=[0.07 -0.07]; end
        for gi=1:numel(g)
            if toc(tR)>BUDGET_MIN*60; full=false; break; end
            f=['pAtp_' g{gi}]; base=q.(f); bc=c; bv=base;
            for s=steps
                cand=q; cand.(f)=base+s; cc=cost(cand);
                if cc < bc-ACCEPT; bc=cc; bv=cand.(f); end
                if toc(tR)>BUDGET_MIN*60; full=false; break; end
            end
            if bv~=base
                fprintf('   %-10s %+.3f -> %+.3f   %.4f (%+.4f)\n', g{gi}, base, bv, bc, bc-c);
                q.(f)=bv; c=bc; improved=true;
            end
        end
    end
    e=struct(); for i=1:numel(g); e.(g{i}) = q.(['pAtp_' g{i}]); end
    OUT(end+1)=struct('id',id,'gates',{g},'tier',SETS{k,3},'why',SETS{k,4}, ...
                      'dof',numel(g),'c',c,'exps',e,'conv',~improved&&full); %#ok<SAGROW>
    fprintf('   -> %.4f  (%d sweeps, %.1f min, %s)\n\n', c, sweep, toc(tR)/60, ...
            iIf(~improved&&full,'converged','TRUNCATED'));
    save(STATE,'OUT','CORE','fn2','ACCEPT','c0');
end

%% ---- ranking --------------------------------------------------------------
[~,ord] = sort([OUT.c]);
fprintf('\n============ CANDIDATE SETS, RANKED BY 2 mM FIT ============\n');
fprintf('H0 (no ATP effect) = %.4f\n\n', c0);
fprintf('%-4s %-5s %4s %10s %8s  %-34s %s\n','rank','id','DOF','cost','conv','gates','tier');
for r=1:numel(ord)
    o=OUT(ord(r));
    fprintf('%-4d %-5s %4d %10.4f %8s  %-34s %s\n', r, o.id, o.dof, o.c, ...
            iIf(o.conv,'yes','NO'), strjoin(o.gates,','), o.tier);
end
fprintf('\ntier A = free [ATP] is a substrate of every gated step\n');
fprintf('     B = nucleotide-state dependent, arguable\n');
fprintf('     C = Pi-contamination proxy\n');
fprintf('     X = control / not defensible\n');

fprintf('\nfitted exponents:\n');
for r=1:numel(ord)
    o=OUT(ord(r)); f=fieldnames(o.exps);
    s=strjoin(arrayfun(@(i) sprintf('%s %+.3f',f{i},o.exps.(f{i})), 1:numel(f),'uni',0), '  ');
    fprintf('  %-5s %s\n', o.id, s);
end

txt = fullfile(resDir,[TAG '_ranking.txt']);
fid=fopen(txt,'w');
fprintf(fid,'CANDIDATE SETS RANKED BY 2 mM FIT   core %s   H0 = %.4f\n\n', CORE, c0);
fprintf(fid,'%-4s %-5s %4s %10s %8s  %-34s %-5s %s\n','rank','id','DOF','cost','conv','gates','tier','rationale');
for r=1:numel(ord)
    o=OUT(ord(r));
    fprintf(fid,'%-4d %-5s %4d %10.4f %8s  %-34s %-5s %s\n', r,o.id,o.dof,o.c, ...
            iIf(o.conv,'yes','NO'), strjoin(o.gates,','), o.tier, o.why);
end
fclose(fid);
fprintf('\nwrote %s\nPLAUSIBLE SETS DONE\n', txt);

%% ==========================================================================
function c = iCost(q, fn2, ATP_c, Data_ATP)
try
    [~,~,fm,fd] = runSlackExperiment(q);
    qf=q; qf.RunSlack=false;
    [~,~,fmv,fdv] = runFVExperiment(qf, ATP_c, Data_ATP);
    f=fieldnames(fmv); for i=1:numel(f); fm.(f{i})=fmv.(f{i}); end
    f=fieldnames(fdv); for i=1:numel(f); fd.(f{i})=fdv.(f{i}); end
    c = sum(evalFeatureCost(fd, fm, fn2, 2));
catch ME
    fprintf(2,'   [eval failed: %s -> penalty]\n', ME.message);
    c = 1e6;
end
end

function s = iIf(c,a,b); if c; s=a; else; s=b; end; end
