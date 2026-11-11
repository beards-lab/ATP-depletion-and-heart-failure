% RunToleranceFloor.m -- is the leftover 2 mM difference the reparameterisation, or
% the ODE tolerance?
%
% The converted gates reproduce every rate factor to 3e-16 and the pipeline is
% bit-reproducible, yet the 2 mM cost moved by 0.016. Only two things can do that:
% a real (if tiny) model difference, or the solver's own truncation error being
% re-rolled when the arithmetic changes in the last bit. The model integrates at
% RelTol 1e-2 / AbsTol 1e-4, so the second is the prime suspect.
%
% Test: sweep the tolerance and watch |cost(hill) - cost(power)|. If it is truncation
% error it must SHRINK with the tolerance; if the two forms really differ it will
% converge to a constant. Runtime is reported too -- tightening is not free.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunToleranceFloor.m')"

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if isempty(gcp('nocreate')); parpool('local',6); end

% RelTol; AbsTol tracks it at RelTol/100. 1e-4 is included but is NOT affordable:
% a single slack run exceeds the 1200 s cap (>20x the 1e-2 cost) and is killed by the
% MaxRunTime event, which is recorded as a failed row rather than aborting the sweep.
% That cost is why the production tolerance is loose in the first place.
TOLS = [1e-2 1e-3 1e-4];
PWR  = 'lowATP_powergate2_opt';
HIL  = 'lowATP_hillgate_opt';

qp = iCfg(getParams(loadParams(PWR), [], true, false));
qh = iCfg(getParams(loadParams(HIL), [], true, false));
FN = qp.fn;

D = nan(numel(TOLS),4);   % [cost_power cost_hill |diff| runtime_s]
KEY = {'rsK','ktr','A','steady'};
KD  = nan(numel(TOLS), numel(KEY));
fprintf('\n%-9s %12s %12s %12s %9s\n','RelTol','cost power','cost hill','|difference|','time (s)');
for i = 1:numel(TOLS)
    qp.odeRelTol = TOLS(i); qp.odeAbsTol = TOLS(i)/100;
    qh.odeRelTol = TOLS(i); qh.odeAbsTol = TOLS(i)/100;
    t = tic;
    try
        [~,~,fmP,fd] = runSlackExperiment(qp);
        [~,~,fmH,~ ] = runSlackExperiment(qh);
        cp = sum(evalFeatureCost(fd,fmP,FN,2)); ch = sum(evalFeatureCost(fd,fmH,FN,2));
        D(i,:) = [cp ch abs(ch-cp) toc(t)];
        for k = 1:numel(KEY)
            if isfield(fmP,KEY{k})
                x = double(fmP.(KEY{k})(:)); y = double(fmH.(KEY{k})(:));
                KD(i,k) = max(abs(x-y)./max(abs(x),eps));
            end
        end
        fprintf('%-9.0e %12.6f %12.6f %12.3g %9.0f\n', TOLS(i), cp, ch, abs(ch-cp), D(i,4));
    catch ME
        % a tighter tolerance can simply be unaffordable -- that is data, not a crash
        D(i,4) = toc(t);
        fprintf('%-9.0e %12s %12s %12s %9.0f   <- %s\n', TOLS(i), '-','-','-', D(i,4), ...
                regexprep(ME.message, '\s+', ' '));
    end
    save(fullfile(resdir,'tolerancefloor.mat'),'TOLS','D','KD','KEY');
end

fprintf('\nlargest relative feature difference between the two forms, by tolerance:\n');
fprintf('%-9s','RelTol'); fprintf('%12s', KEY{:}); fprintf('\n');
for i = 1:numel(TOLS)
    fprintf('%-9.0e', TOLS(i)); fprintf('%12.3g', KD(i,:)); fprintf('\n');
end

fprintf(['\nReading: shrinking with the tolerance => the difference is solver truncation error,\n' ...
         'not a difference between the two gate forms, and its size at the production tolerance\n' ...
         '(1e-2) is the noise floor of every cost this project quotes.\n']);
fprintf('\nTOLERANCE FLOOR DONE\n');

function q = iCfg(q)
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false; q.RunStairs=false;
q.RunSlackPassive=false;  q.RunSlack=true; q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=1200;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0; q.UseR2Ceiling=0; q.UseD0State=0;
q.MgATP=2; q.MgADP=0; q.Pi=0.8;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
fn = q.fn; [~,iu] = unique(fn,'stable'); fn = fn(iu);
fn = cellfun(@(s) regexprep(s,'^rsK(\|.*)?$','rsK|0.1'), fn, 'UniformOutput', false);
q.fn = fn(~(startsWith(fn,'FV') | startsWith(fn,'PS_')));
end
