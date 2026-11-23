% RunAtpDoseResponse.m -- where the power law and the saturable gate actually differ.
%
% The two forms are pinned together at the two measured concentrations (8 mM by
% construction, 2 mM by the closed-form conversion in RunHillGateConversion). So no
% fit can tell them apart. They are NOT the same model, though: the power law is
% scale-free and unbounded, the Hill form saturates. This script quantifies the
% divergence over the [ATP] range the project actually cares about -- depletion below
% 2 mM -- both for the bare gate factors and for the model's emergent observables.
%
% Run: matlab -batch "run('Analyses/AtpGateHillForm/RunAtpDoseResponse.m')"
% Knobs: ATPLIST (concentrations, mM), SKIP_SIM (analytic figure only).

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..','..');
cd(root); addpath(genpath(root)); LoadData;
resdir = fullfile(here,'results'); if ~exist(resdir,'dir'); mkdir(resdir); end
if ~exist('SKIP_SIM','var'); SKIP_SIM = false; end
if ~exist('ATPLIST','var'); ATPLIST = [0.25 0.5 1 2 3 4 6 8 12 20]; end

PWR = 'lowATP_powergate2_opt';   % incumbent power-law gates
HIL = 'lowATP_hillgate_opt';     % converted saturable gates
OBS = {'steady','ktr','rsK','peak1_y','A'};

bp = getParams(loadParams(PWR), [], true, false);
bh = getParams(loadParams(HIL), [], true, false);
Aref = bp.MgATP_ref;

%% ---- 1. bare gate factors over [ATP] ------------------------------------------
GATED = {'k2','kmsrd','kah','kamh','ksr0','ksr2srd','ksrd2sr','kSE','mu_neg'};
x = logspace(log10(0.1), log10(30), 400);
fprintf('\n=== gate factor divergence (relative to the 8 mM baseline) ===\n');
fprintf('%-10s %8s %8s | %s\n','rate','p','K (mM)','power / hill factor at 4, 2, 1, 0.5, 0.25 mM');
probe = [4 2 1 0.5 0.25];
Fp = nan(numel(GATED), numel(x)); Fh = Fp;
for i = 1:numel(GATED)
    p = bp.(['pAtp_' GATED{i}]);
    n = bh.(['pAtp_' GATED{i}]); K = bh.(['KAtp_' GATED{i}]);
    Fp(i,:) = atpGateFactor(x, Aref, Inf, p);
    Fh(i,:) = atpGateFactor(x, Aref, K,   n);
    fprintf('%-10s %8.4f %8.4g |', GATED{i}, p, K);
    for a = probe
        fprintf(' %5.3f/%5.3f', atpGateFactor(a,Aref,Inf,p), atpGateFactor(a,Aref,K,n));
    end
    fprintf('\n');
end

try
    f1 = figure('Position',[40 40 1250 520],'Color','w');
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; box on; grid on; set(gca,'XScale','log');
    co = lines(numel(GATED));
    for i = 1:numel(GATED)
        plot(x, Fp(i,:), '--', 'Color', co(i,:), 'LineWidth', 1.1, 'HandleVisibility','off');
        plot(x, Fh(i,:), '-',  'Color', co(i,:), 'LineWidth', 1.8, 'DisplayName', GATED{i});
    end
    xline(8,'k:','8 mM','LabelVerticalAlignment','bottom','HandleVisibility','off');
    xline(2,'k:','2 mM','LabelVerticalAlignment','bottom','HandleVisibility','off');
    xlabel('[MgATP] (mM)'); ylabel('gate factor (1 at 8 mM)');
    title('dashed = power law   solid = saturable Hill'); legend('Location','southeast','NumColumns',2);
    nexttile; hold on; box on; grid on; set(gca,'XScale','log');
    for i = 1:numel(GATED)
        plot(x, Fh(i,:)./Fp(i,:), '-', 'Color', co(i,:), 'LineWidth', 1.6, 'DisplayName', GATED{i});
    end
    yline(1,'k-','HandleVisibility','off');
    xline(8,'k:','HandleVisibility','off'); xline(2,'k:','HandleVisibility','off');
    xlabel('[MgATP] (mM)'); ylabel('hill / power');
    title('identical at 2 and 8 mM by construction, nowhere else');
    title(tl,'ATP gates: the data cannot distinguish the forms, the extrapolation can');
    exportgraphics(f1, fullfile(resdir,'gate_factors.png'), 'Resolution', 140);
catch ME; fprintf('gate figure skipped: %s\n', ME.message); end

save(fullfile(resdir,'gatefactors.mat'),'x','Fp','Fh','GATED','ATPLIST');
if SKIP_SIM, fprintf('SKIP_SIM -> analytic part only.\n'); return; end

%% ---- 2. model dose-response under both forms ----------------------------------
if isempty(gcp('nocreate')); parpool('local',6); end
R = struct('form',{},'atp',{},'ok',{},'feat',{});
for form = 1:2
    if form == 1, q0 = iCfg(bp); nm = 'power'; else, q0 = iCfg(bh); nm = 'hill'; end
    for a = ATPLIST
        q = q0; q.MgATP = a;
        if a == Aref                      % the 8 mM protocol has its own driving file
            q.Pi = 0; q.velocitytableonfile = 'protocol_03_27_2026_8mM_slack.mat';
        end
        t = tic;
        try
            [~,~,fm,~] = runSlackExperiment(q);
            R(end+1) = struct('form',nm,'atp',a,'ok',true,'feat',fm); %#ok<SAGROW>
            fprintf('  [%-5s] %5.2f mM  ok  (%4.0f s)\n', nm, a, toc(t));
        catch ME
            R(end+1) = struct('form',nm,'atp',a,'ok',false,'feat',[]); %#ok<SAGROW>
            fprintf('  [%-5s] %5.2f mM  FAILED: %s\n', nm, a, ME.message);
        end
        save(fullfile(resdir,'doseresponse.mat'),'R','ATPLIST','OBS','PWR','HIL');
    end
end

%% ---- 3. report ----------------------------------------------------------------
V = nan(numel(OBS), numel(ATPLIST), 2);
for k = 1:numel(R)
    if ~R(k).ok, continue; end
    ia = find(ATPLIST == R(k).atp, 1); jf = 1 + strcmp(R(k).form,'hill');
    for io = 1:numel(OBS)
        if isfield(R(k).feat, OBS{io})
            v = double(R(k).feat.(OBS{io})); v = v(isfinite(v));
            if ~isempty(v), V(io,ia,jf) = mean(v); end
        end
    end
end
i8 = find(ATPLIST == Aref, 1);
fprintf('\n=== model observables vs [ATP], normalised to the 8 mM value ===\n');
for io = 1:numel(OBS)
    fprintf('\n%s  (8 mM absolute: %.4g)\n', OBS{io}, V(io,i8,1));
    fprintf('  %-8s','[ATP]'); fprintf('%8.3g', ATPLIST); fprintf('\n');
    fprintf('  %-8s','power'); fprintf('%8.3f', V(io,:,1)/V(io,i8,1)); fprintf('\n');
    fprintf('  %-8s','hill');  fprintf('%8.3f', V(io,:,2)/V(io,i8,2)); fprintf('\n');
    fprintf('  %-8s','h/p');   fprintf('%8.3f', V(io,:,2)./V(io,:,1)); fprintf('\n');
end

try
    f2 = figure('Position',[40 40 1300 760],'Color','w');
    tl2 = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    for io = 1:numel(OBS)
        nexttile; hold on; box on; grid on; set(gca,'XScale','log');
        plot(ATPLIST, V(io,:,1)/V(io,i8,1), 'o--','LineWidth',1.3,'DisplayName','power law');
        plot(ATPLIST, V(io,:,2)/V(io,i8,2), 's-', 'LineWidth',1.8,'DisplayName','saturable Hill');
        xline(8,'k:','HandleVisibility','off'); xline(2,'k:','HandleVisibility','off');
        xlabel('[MgATP] (mM)'); ylabel(sprintf('%s / %s(8 mM)', OBS{io}, OBS{io}));
        title(OBS{io}); if io==1, legend('Location','best'); end
    end
    title(tl2, 'Model ATP dose-response: the two gate forms are welded at 2 and 8 mM only');
    exportgraphics(f2, fullfile(resdir,'dose_response.png'), 'Resolution', 140);
catch ME; fprintf('dose figure skipped: %s\n', ME.message); end

save(fullfile(resdir,'doseresponse.mat'),'R','V','ATPLIST','OBS','PWR','HIL');
fprintf('\nATP DOSE-RESPONSE DONE\n');

%% ---- helper -------------------------------------------------------------------
function q = iCfg(q)
q.RunForceVelocity=false; q.RunForceVelocityTime=false; q.RunKtr=false; q.RunStairs=false;
q.RunSlackPassive=false;  q.RunSlack=true; q.EvalFeatures=true; q.BreakOnODEUnstable=false;
q.RunSlackSegments='AllPar'; q.MaxRunTime=300;
q.PlotEachSeparately=0; q.PlotFeatureFitting=0;
q.UseAtpDetachR2=0; q.UseAtpKmsrd=0; q.UseStrainDepAtpGate=0; q.UseR2Ceiling=0; q.UseD0State=0;
q.MgADP=0; q.Pi=0.8;
q.velocitytableonfile='protocol_03_27_2026_2mM_slack.mat';
end
