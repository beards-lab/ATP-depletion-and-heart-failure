% RunVerify.m -- price the tau_reg fix properly.
%  (a) per-feature cost, incumbent vs fast gate
%  (b) held-out: the ktr protocol of ALL THREE days
%  (c) held-out: the stairs protocol (RunStairs=0 in this paramset)
%  (d) the slack features the fix is supposed to trade against
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

base = getParams(loadParams('estiff125_opt'), [], true, false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=1; base.BreakOnODEUnstable=0;
svt = load(fullfile(root,'data',base.velocitytableonfile)).velocitytable;
iS  = find(svt(:,2)<-1);

VAR = {'V0 incumbent tau_reg=26.6ms', base.tau_reg
       'VF fast gate  tau_reg=0.2ms', 2e-4};
DAYS = {'03_27','04_03','04_10'};

Q = struct();
for v=1:2
    p = base; p.tau_reg = VAR{v,2}; p = getParams(p,p.g,true);

    [Et,Eraw,fc,fn] = evaluateBakersExp([],p,false,false);
    [~, os, fm, fd] = runSlackExperiment(p);
    FisoM = median(os.Force(os.t>=svt(iS(2),1)-0.05 & os.t<=svt(iS(2),1)-0.005));
    ws = recoveryWindows(os.t,os.Force,svt,'slack',FisoM); ws=ws(strcmp({ws.type},'postSlack'));

    % --- ktr protocol, all three days ---
    kk = nan(1,3); ek = nan(1,3);
    for d=1:3
        kvt = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d}))).velocitytable;
        Kd  = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
        pk=p; pk.SL0=2*kvt(1,4); pk=getParams(pk,pk.g,true); pk.Velocity=kvt(:,2);
        [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk);
        FisoK = interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');
        w = recoveryWindows(ok.t,ok.Force,kvt,'ktr',FisoK); kk(d)=w(1).k63;
        Fi = interp1(ok.t,ok.Force/FisoK,Kd.datatable(:,1),'linear','extrap');
        ek(d) = mean((Kd.datatable(:,3)-Fi).^2)*10;
        if d==1; Q(v).ok=ok; Q(v).kvt=kvt; Q(v).FisoK=FisoK; end
    end

    % --- held-out: stairs ---
    ps = p; ps.RunStairs=1;
    try
        [Estairs, ~] = runStairsExperiment(ps);
    catch ME
        Estairs = NaN; fprintf('stairs failed: %s\n', ME.message);
    end

    Q(v).name=VAR{v,1}; Q(v).Et=Et; Q(v).fc=fc; Q(v).fn=fn; Q(v).Eraw=Eraw;
    Q(v).kS=[ws.k63]; Q(v).kk=kk; Q(v).ek=ek; Q(v).Estairs=Estairs;
    Q(v).fm=fm; Q(v).fd=fd; Q(v).os=os; Q(v).FisoM=FisoM;
end

fprintf('\n=========== PER-FEATURE COST ===========\n');
fn=Q(1).fn(:)'; d=Q(2).fc-Q(1).fc; [~,ord]=sort(abs(d),'descend');
fprintf('%-26s %9s %9s %9s\n','feature','incumbent','fast','delta');
for i=ord; fprintf('%-26s %9.4f %9.4f %+9.4f\n', fn{i}, Q(1).fc(i), Q(2).fc(i), d(i)); end
fprintf('%-26s %9.4f %9.4f %+9.4f\n','TOTAL Et',Q(1).Et,Q(2).Et,Q(2).Et-Q(1).Et);

fprintf('\n=========== HELD OUT ===========\n');
fprintf('%-26s %8s %8s %8s | %8s\n','','03/27','04/03','04/10','stairs E');
fprintf('%-26s %8.1f %8.1f %8.1f | %8s\n','DATA ktr k63',49.4,42.1,38.1,'--');
for v=1:2
    fprintf('%-26s %8.1f %8.1f %8.1f | %8.4f\n', ['MODEL ' VAR{v,1}(1:2) ' ktr k63'], Q(v).kk, Q(v).Estairs);
end
for v=1:2
    fprintf('%-26s %8.4f %8.4f %8.4f |\n', ['MODEL ' VAR{v,1}(1:2) ' ktr traceMSE'], Q(v).ek);
end

fprintf('\n=========== SLACK: the features the fix should trade against ===========\n');
g=@(s,f) s.(f)(:)';
flds={'ktr','A','t0','peak1_y','vall_y','peak2','steady','vall2_dy','rsK'};
for k=1:numel(flds)
    if ~isfield(Q(1).fm,flds{k}); continue; end
    fprintf('%-11s data %s\n', flds{k}, mat2str(round(g(Q(1).fd,flds{k}),3)));
    fprintf('%-11s  V0  %s\n', '', mat2str(round(g(Q(1).fm,flds{k}),3)));
    fprintf('%-11s  VF  %s\n', '', mat2str(round(g(Q(2).fm,flds{k}),3)));
end
fprintf('\npost-slack k63 per cycle:\n  V0 %s\n  VF %s\n', mat2str(round(Q(1).kS,1)), mat2str(round(Q(2).kS,1)));
save(fullfile(resdir,'verify.mat'),'Q','fn','-v7.3');
