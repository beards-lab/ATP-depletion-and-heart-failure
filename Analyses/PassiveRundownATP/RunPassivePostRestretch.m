% RunPassivePostRestretch.m
% =========================================================================
% CORRECTED framing + the test that actually matches the hypothesis.
%
% DIRECTION (user, and I now agree). Acquisition order 03/27 was
%     8 mM  ->  2 mM  ->  PNB+Mava (passive)
% and rundown removes force-generating material in parallel, which also carries
% titin. So the PASSIVE force falls monotonically through the session:
%     passive(8 mM) > passive(2 mM) > passive(PNB+Mava)
% The model's passive is fitted to PNB+Mava -- i.e. to the LOWEST of the three.
% Therefore, relative to the model, BOTH active runs should want MORE passive,
% and 8 mM should want MORE than 2 mM. My first test compared the 2 mM optimum
% against the 2-vs-8 DIFFERENCE (x0.973) and called it "wrong direction". That
% reference point was wrong: against the PNB+Mava-fitted baseline the expected
% answer is x>1 for both, which is what the scan actually found (x1.10).
%
% MAGNITUDE. 03/27 loses 3.18 kPa / 3.8 % per activation interval. PNB+Mava sits
% after BOTH activations, so cumulatively ~7.6 % below the 8 mM state:
%     model(=PNB+Mava) x1.000 | 2 mM x~1.04 | 8 mM x~1.08
%
% WHERE IT ACTS (user). Passive should shape the SECOND peak and the post-restretch
% DECAY, hence the redevelopment time constant after restretch -- i.e. peak2,
% vall2_dy, ovrsht_dy, coolDownLS and rsK. My previous scan reported only peak1_y
% and steady, so it could not see this. Fixed here.
% =========================================================================
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir=fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end

% Branch-portable baseline: the low-ATP snapshot exists only on lowatp-2state, so
% fall back to the shared high-ATP baseline elsewhere (e.g. on Further-optim).
BASE='lowATP_powergate2_opt';
if ~isfile(fullfile(root,'params',[BASE '.m'])); BASE='rskR2_w025_opt'; end
fprintf('baseline: %s\n', BASE);
b=getParams(loadParams(BASE),[],true,false);
b.RunForceVelocity=false; b.RunForceVelocityTime=false; b.RunKtr=false; b.RunStairs=false;
b.RunSlackPassive=false; b.RunSlack=true; b.EvalFeatures=true; b.BreakOnODEUnstable=false;
b.RunSlackSegments='AllPar'; b.MaxRunTime=400; b.PlotEachSeparately=0; b.PlotFeatureFitting=0;
m=@(x)mean(double(x),'omitnan');
KP=b.k_pas;

f8d=load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat')).features_data;
f2d=load(fullfile(root,'data','protocol_03_27_2026_2mM_slack.mat')).features_data;

fprintf('DATA   8mM: peak2 %6.2f vall2_dy %6.3f ovrsht_dy %6.3f rsK %6.1f\n', ...
    m(f8d.peak2), m(f8d.vall2_dy), m(f8d.ovrsht_dy), m(f8d.rsK));
fprintf('DATA   2mM: peak2 %6.2f vall2_dy %6.3f ovrsht_dy %6.3f rsK %6.1f  | rsK ratio %.3f\n\n', ...
    m(f2d.peak2), m(f2d.vall2_dy), m(f2d.ovrsht_dy), m(f2d.rsK), m(f2d.rsK)/m(f8d.rsK));

% ---- Q1: does passive move the POST-RESTRETCH observables at all? ----
sc=[0.5 1.0 1.5 2.0 3.0];
fprintf('=== Q1  passive sensitivity of the post-restretch window (2 mM) ===\n');
fprintf('%7s | %7s %7s %8s %9s %7s | %7s %7s\n','k_pas x','peak2','vall_y','vall2_dy','ovrsht_dy','rsK','peak1_y','steady');
for s=sc
    q=b; q.MgATP=2; q.Pi=0.8; q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; q.k_pas=KP*s;
    try
        [~,~,f,~]=runSlackExperiment(q);
        fprintf('%7.2f | %7.2f %7.2f %8.3f %9.3f %7.1f | %7.2f %7.2f\n', s, ...
            m(f.peak2), m(f.vall_y), m(f.vall2_dy), m(f.ovrsht_dy), m(f.rsK), m(f.peak1_y), m(f.steady));
    catch ME; fprintf('%7.2f | FAILED %s\n', s, ME.message); end
end

% ---- Q2: the rundown-predicted DIFFERENTIAL -- does it move the rsK ratio? ----
fprintf('\n=== Q2  rundown-predicted differential: 8 mM x1.08, 2 mM x1.04 ===\n');
cases={ {'model as fitted (both x1.00)',1.00,1.00}, ...
        {'rundown differential 1.08/1.04',1.08,1.04}, ...
        {'exaggerated 1.5/1.2',1.50,1.20}, ...
        {'exaggerated 2.0/1.3',2.00,1.30} };
fn=b.fn; [~,iu]=unique(fn,'stable'); fn=fn(iu);
fn=cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput',false);
FN=fn(~(startsWith(fn,'FV')|startsWith(fn,'PS_')));
fprintf('%-32s | %6s %6s %7s | %7s %7s | %8s\n','case','rsK8','rsK2','rsKrat','peak2_2','ovsh_2','cost2');
for i=1:numel(cases)
    nm=cases{i}{1}; s8=cases{i}{2}; s2=cases{i}{3};
    q8=b; q8.MgATP=8; q8.Pi=0; q8.velocitytableonfile='protocol_03_27_2026_8mM_slack.mat'; q8.k_pas=KP*s8;
    q2=b; q2.MgATP=2; q2.Pi=0.8; q2.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat'; q2.k_pas=KP*s2;
    try
        [~,~,g8,~]=runSlackExperiment(q8);
        [~,~,g2,d2]=runSlackExperiment(q2);
        c2=sum(evalFeatureCost(d2,g2,FN,2));
        fprintf('%-32s | %6.1f %6.1f %7.3f | %7.2f %7.3f | %8.4f\n', nm, ...
            m(g8.rsK), m(g2.rsK), m(g2.rsK)/m(g8.rsK), m(g2.peak2), m(g2.ovrsht_dy), c2);
    catch ME; fprintf('%-32s | FAILED %s\n', nm, ME.message); end
end
fprintf('%-32s | %6.1f %6.1f %7.3f | %7.2f %7.3f |\n','DATA TARGET', ...
    m(f8d.rsK), m(f2d.rsK), m(f2d.rsK)/m(f8d.rsK), m(f2d.peak2), m(f2d.ovrsht_dy));
disp('POST-RESTRETCH PASSIVE TEST DONE');
