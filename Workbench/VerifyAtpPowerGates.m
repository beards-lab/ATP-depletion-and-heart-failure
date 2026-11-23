% VerifyAtpPowerGates.m — (a) 8 mM must be BIT-IDENTICAL with arbitrary exponents,
% (b) 2 mM must actually move, (c) time the evaluation paths for the optimiser.
root='C:\home\git\ATP-wt-lowatp'; cd(root); addpath(genpath(root)); refreshPool(6); LoadData;

b=getParams(loadParams('rskR2_w025_opt'),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.MaxRunTime=300; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
b.UseAtpDetachR2=0; b.UseAtpKmsrd=0; b.UseStrainDepAtpGate=0;   % power gates replace them
b.UseR2Ceiling=0; b.UseD0State=0; b.MgATP_ref=8;
b.UseAtpPowerGates=1;

P={'k2','kmsrd','kah','kamh','ka','kd','k1','k_1','ksr0','ksr2srd','ksrd2sr','tau_reg','A0','eta_M'};
setp=@(q,v) setPvec(q,P,v);

% ---- (a) 8 mM invariance with arbitrary exponents ----
q0=b; q0.RunSlackSegments='AllPar'; q0.MgATP=8; q0.Pi=0;
q0.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat';
[~,~,f0,~]=runSlackExperiment(setp(q0, zeros(1,numel(P))));
rng(7); vrand = randn(1,numel(P))*0.8;
[~,~,f1,~]=runSlackExperiment(setp(q0, vrand));
d=max(abs(double(f0.steady(:))-double(f1.steady(:))));
d2=max(abs(double(f0.rsK(:))-double(f1.rsK(:))));
fprintf('\n(a) 8 mM invariance with random exponents: max|dsteady| = %.3g, max|drsK| = %.3g  (MUST be 0)\n', d, d2);

% ---- (b) 2 mM responds ----
q2=b; q2.RunSlackSegments='AllPar'; q2.MgATP=2; q2.Pi=0.8;
q2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
[~,~,g0,gd]=runSlackExperiment(setp(q2, zeros(1,numel(P))));
[~,~,g1,~ ]=runSlackExperiment(setp(q2, vrand));
m=@(x)mean(double(x),'omitnan');
fprintf('(b) 2 mM p=0    : steady %.1f rsK %.1f ktr %.1f\n', m(g0.steady),m(g0.rsK),m(g0.ktr));
fprintf('    2 mM p=rand : steady %.1f rsK %.1f ktr %.1f\n', m(g1.steady),m(g1.rsK),m(g1.ktr));
fprintf('    DATA        : steady %.1f rsK %.1f ktr %.1f\n', m(gd.steady),m(gd.rsK),m(gd.ktr));

% ---- (c) timings ----
qA=setp(q2,zeros(1,numel(P))); qA.RunSlackSegments='AllPar';
qS=setp(q2,zeros(1,numel(P))); qS.RunSlackSegments='All';
t=tic; runSlackExperiment(qA); tA=toc(t);
t=tic; runSlackExperiment(qS); tS=toc(t);
fprintf('(c) one 2 mM eval: AllPar %.1f s | All(serial) %.1f s\n', tA, tS);

% ---- (d) can candidates be evaluated inside parfor? ----
ok=true; tp=tic;
try
    cand=repmat({qS},1,6); cc=zeros(1,6);
    parfor i=1:6
        [~,~,fm,fd]=runSlackExperiment(cand{i});
        cc(i)=sum(evalFeatureCost(fd,fm,{'steady|50','ktr|10'},2));
    end
    fprintf('(d) parfor over 6 candidates OK in %.1f s -> %.1f s/eval effective (costs %s)\n', ...
        toc(tp), toc(tp)/6, num2str(cc,'%.3f '));
catch ME
    ok=false; fprintf('(d) parfor FAILED: %s\n', ME.message);
end
fprintf('\nRECOMMEND: %s\n', ternstr(ok,'parallel (1+lambda)-ES over candidates','sequential AllPar evals'));
disp('VERIFY DONE');

function q = setPvec(q,P,v)
    q.UseAtpPowerGates=1;
    for i=1:numel(P); q.(['pAtp_' P{i}]) = v(i); end
end
function s=ternstr(c,a,b); if c; s=a; else; s=b; end; end
