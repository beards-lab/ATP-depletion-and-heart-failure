% CheckTimecourses.m -- is the refit actually closer to the data, or has it
% found a way to hit the feature values with a broken trace?
%
% Compares params/ktrfence_seed.m against the campaign's current best
% params/ktrfence_opt.m on every protocol, and runs a pathology sweep over the
% state vector that the feature objective cannot see.
%
% CAVEAT ON "THE POOLED DATA". pooledTargetFeatures swaps the FEATURE VALUES,
% not the trace: the datatable in the protocol file is still 03/27. So the black
% trace below is one preparation and the model is NOT fitted to it. The pooled
% targets -- what the fit actually chases -- are drawn as markers.
%
%   matlab -batch "run('Analyses/RefitStrategy/CheckTimecourses.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RefitStrategy','results');
if ~exist(resdir,'dir'); mkdir(resdir); end

NW = 3;                                   % leave cores for the running campaign
c = parcluster('local');
jsl = fullfile(tempdir, sprintf('mljob_tc_%d', feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end
c.JobStorageLocation = jsl;
if isempty(gcp('nocreate')); parpool(c, NW); end

TAGS = {'ktrfence_seed','ktrfence_opt'};
R = struct();
for k = 1:2
    p = getParams(loadParams(TAGS{k}), [], true, false);
    p.PlotEachSeparately=0; p.ghostSave=''; p.ghostLoad=''; p.SaveBest=0;
    p.EvalFeatures=1; p.BreakOnODEUnstable=0; p.MaxRunTime=600;   % diagnostic, not a fit
    p = getParams(p, p.g, true);
    params0 = p; LoadData; RunBakersExp;          %#ok<*NASGU>
    R(k).tag=TAGS{k}; R(k).p=p; R(k).fm=features_model; R(k).fd=features_data;
    R(k).cost=sum(evalFeatureCost(features_data,features_model,p.fn));
    R(k).os=out_slack; R(k).opas=out_pas; R(k).ok=out_ktr; R(k).outs=outs;
    fprintf('%-16s cost %.4f\n', TAGS{k}, R(k).cost);
end

svt = load(fullfile(root,'data',R(1).p.velocitytableonfile));
SVT = svt.velocitytable; SDT = svt.datatable;
kvt = load(fullfile(root,'data',R(1).p.ktr_velocitytableonfile));
KVT = kvt.velocitytable; KDT = kvt.datatable;
iS = find(SVT(:,2)<-1); iR = find(SVT(:,2)>1);

%% ================= PATHOLOGY SWEEP =================
fprintf('\n================ PATHOLOGY SWEEP (what the features cannot see) ================\n');
fprintf('%-34s %14s %14s\n','check','seed','opt');
chk = {};
for k=1:2
    o = R(k).os; p = R(k).p; Ns=p.NumberOfStates; ss=p.ss;
    tot = o.p1_0 + o.p2_0 + o.PuATP + o.PuR + o.SR + o.SRD;
    iA  = Ns*ss+8;  iAd = Ns*ss + 8 + p.UseRegistrationAvailability + p.UseD0State;
    Areg = o.PU(:,iA); Adef = zeros(size(Areg));
    if p.UseLoadRealign; Adef = o.PU(:,iAd); end
    F = o.Force(:); t = o.t(:);
    dt = diff(t); dF = diff(F);
    % sign changes of dF per second -- a proxy for solver ringing
    sc = sum(diff(sign(dF))~=0) / (t(end)-t(1));
    v(k).probLeak   = max(abs(tot-1));
    v(k).Fmin       = min(F);
    v(k).Fmax       = max(F);
    v(k).negFrac    = mean(F < -0.5);
    v(k).SRmin      = min([o.SR o.SRD]);   v(k).SRmax = max(o.SR+o.SRD);
    v(k).PDmin      = min(o.PU(:,Ns*ss+5)); v(k).PDmax = max(o.PU(:,Ns*ss+5));
    v(k).attMax     = max(o.p1_0+o.p2_0);
    v(k).AregMin    = min(Areg); v(k).AregMax = max(Areg);
    v(k).AdefMax    = max(Adef);
    v(k).edgeMass   = max(max(abs(o.PU(:,[1 ss ss+1 2*ss]))));
    v(k).signChg    = sc;
    v(k).nSamp      = numel(t);
    v(k).SLerr      = max(abs(interp1(SVT(:,1),SVT(:,4)*2,t,'previous','extrap') - o.SL(:)));
end
pr = @(n,f,fmt) fprintf(['%-34s ' fmt ' ' fmt '\n'], n, v(1).(f), v(2).(f));
pr('max |sum(prob) - 1|','probLeak','%14.2e');
pr('min Force (kPa)','Fmin','%14.3f');
pr('max Force (kPa)','Fmax','%14.2f');
pr('fraction of trace below -0.5 kPa','negFrac','%14.4f');
pr('min SR / SRD','SRmin','%14.4f');
pr('max SR+SRD','SRmax','%14.4f');
pr('min PD','PDmin','%14.4f');
pr('max PD','PDmax','%14.4f');
pr('max attached p1+p2','attMax','%14.4f');
pr('min A_reg','AregMin','%14.4f');
pr('max A_reg','AregMax','%14.4f');
pr('max A_def','AdefMax','%14.4f');
pr('max |prob| at a strain-grid edge','edgeMass','%14.3e');
pr('dF sign changes per second','signChg','%14.1f');
pr('max |SL - commanded| (um)','SLerr','%14.4f');
pr('samples in the slack run','nSamp','%14.0f');
fprintf('\nA_0 = %.4f (A_reg floor), Adef_max = %.2f, P_bound_max = %.3f\n', ...
        R(2).p.A0, R(2).p.Adef_max, R(2).p.P_bound_max);

%% ================= FEATURE vs POOLED TARGET =================
fprintf('\n================ MODEL vs the POOLED TARGET it is fitted to ================\n');
FLD = {'A','peak1_y','vall_y','peak2','steady','vall2_dy','ktr','rsK80','peak1_dSL','vall_t'};
for f = FLD
    if ~isfield(R(2).fd,f{1}); continue; end
    d = R(2).fd.(f{1})(:)'; m1 = R(1).fm.(f{1})(:)'; m2 = R(2).fm.(f{1})(:)';
    n = min([numel(d) numel(m1) numel(m2)]);
    fprintf('%-12s pooled %s\n', f{1}, mat2str(round(d(1:n),3)));
    fprintf('%-12s  seed  %s\n', '',    mat2str(round(m1(1:n),3)));
    fprintf('%-12s  OPT   %s   worst %+.1f%%\n', '', mat2str(round(m2(1:n),3)), ...
            100*max(abs(m2(1:n)-d(1:n))./abs(d(1:n))));
end

%% ================= FIGURE =================
f1=figure('Position',[20 20 1700 950],'Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
CS=[0.85 0.33 0.10]; CO=[0 0.45 0.74]; CD=[0 0 0];

% --- A: whole slack protocol ---
nexttile; hold on; box on; grid on
plot(SDT(:,1),SDT(:,3),'-','Color',CD,'LineWidth',1.0);
plot(R(1).os.t,R(1).os.Force,'-','Color',CS,'LineWidth',1.0);
plot(R(2).os.t,R(2).os.Force,'-','Color',CO,'LineWidth',1.4);
xlim([SVT(iS(1),1)-0.3 SVT(end,1)]); xlabel('t (s)'); ylabel('F (kPa)');
title('A   slack protocol  (black = 03/27 trace; the FIT chases POOLED features, not this trace)');
legend({'data 03/27','seed','opt'},'Location','southwest','Box','off');

% --- B: one cycle, with the pooled targets it is actually fitted to ---
nexttile; hold on; box on; grid on
t0=SVT(iS(2),1); t1=SVT(iR(2)+1,1)+0.30;
m=SDT(:,1)>=t0-0.02 & SDT(:,1)<=t1;
plot((SDT(m,1)-t0)*1e3,SDT(m,3),'-','Color',CD,'LineWidth',1.0);
for k=1:2
    mm=R(k).os.t>=t0-0.02 & R(k).os.t<=t1;
    plot((R(k).os.t(mm)-t0)*1e3, R(k).os.Force(mm),'-', ...
        'Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.2+0.3*(k==2));
end
fdp=R(2).fd;
for nm = {'peak1_y','vall_y','peak2','steady'}
    if isfield(fdp,nm{1}); yline(fdp.(nm{1})(2),':','Color',[0.4 0.4 0.4]); end
end
xlabel('t since slack #2 release (ms)'); ylabel('F (kPa)');
title('B   slack cycle 2  (dotted = pooled targets peak1\_y / vall\_y / peak2 / steady)');

% --- C: force-velocity ---
nexttile; hold on; box on; grid on
AllV=[0 0.5 1 2 3 4 5 6]; D1=[56.40,51.8120,37.4459,17.8025,11.4430,6.2643,3.2759,2.2120];
plot(AllV,D1/D1(1),'ks-','MarkerFaceColor','k','LineWidth',1.2);
for k=1:2
    Fa = arrayfun(@(o) o.Force(end), R(k).outs);
    vv = -R(k).p.FV_velocities;
    [vv,ix]=sort(vv); Fa=Fa(ix);
    plot(vv,Fa/Fa(1),'o-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.4);
end
xlabel('shortening velocity (ML/s)'); ylabel('F / F_{iso}');
title('C   force-velocity (scored points only)'); legend({'data','seed','opt'},'Box','off');

% --- D: passive ---
nexttile; hold on; box on; grid on
for k=1:2
    if isempty(R(k).opas); continue; end
    plot(R(k).opas.t, R(k).opas.Force,'-','Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.2);
end
xlabel('t (s)'); ylabel('F (kPa)'); title('D   passive protocol (PNB/Mava, active off)');
legend({'seed','opt'},'Box','off');

% --- E: ktr protocol, held out ---
nexttile; hold on; box on; grid on
tK=KVT(end-1,1);
mk=KDT(:,1)>=tK-0.02;
plot((KDT(mk,1)-tK)*1e3,KDT(mk,3),'-','Color',CD,'LineWidth',1.6);
for k=1:2
    FisoK=interp1(R(k).ok.t,R(k).ok.Force,KVT(2,1),'linear','extrap');
    mm=R(k).ok.t>=tK-0.02;
    plot((R(k).ok.t(mm)-tK)*1e3, R(k).ok.Force(mm)/FisoK,'-', ...
        'Color',(k==1)*CS+(k==2)*CO,'LineWidth',1.2+0.3*(k==2));
end
xlim([-20 59]); xlabel('t since the final release (ms)'); ylabel('F / F_{iso}');
title('E   ktr protocol  — HELD OUT (fenced only)'); legend({'data','seed','opt'},'Location','southeast','Box','off');

% --- F: state sanity over the slack run ---
nexttile; hold on; box on; grid on
o=R(2).os; p=R(2).p; Ns=p.NumberOfStates; ss=p.ss;
tot=o.p1_0+o.p2_0+o.PuATP+o.PuR+o.SR+o.SRD;
plot(o.t,o.p1_0+o.p2_0,'-','LineWidth',1.2);
plot(o.t,o.SR+o.SRD,'-','LineWidth',1.2);
plot(o.t,o.PU(:,Ns*ss+5),'-','LineWidth',1.2);
plot(o.t,o.PU(:,Ns*ss+8),'-','LineWidth',1.2);
plot(o.t,tot-1,'k-','LineWidth',1.0);
xlim([SVT(iS(1),1)-0.3 SVT(end,1)]); xlabel('t (s)'); ylabel('probability');
title('F   opt: state occupancies over the slack protocol');
legend({'attached p1+p2','SR+SRD','PD','A\_reg','\Sigmaprob - 1'},'Location','east','Box','off');

exportgraphics(f1, fullfile(resdir,'timecourse_check.png'),'Resolution',110);
fprintf('\nwrote %s\n', fullfile(resdir,'timecourse_check.png'));
save(fullfile(resdir,'timecourse_check.mat'),'v','TAGS');
fprintf('\nDONE\n');
