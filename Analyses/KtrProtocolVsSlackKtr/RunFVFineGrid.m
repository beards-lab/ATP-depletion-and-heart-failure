% RunFVFineGrid.m -- does the fast gate break force-velocity?
% The scored FV uses only FIVE velocities (0, 0.5, 1, 2, 4 ML/s). This resolves
% the curve at 0.25 ML/s up to 4 and 1 ML/s from 4 to the end of the protocol,
% for the incumbent and the fast gate, and puts the two reference datasets on it.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');

VGRID = -[0, 0.25:0.25:4, 5, 6];            % 0.25 steps to 4, then 1 to the end
base=getParams(loadParams('estiff125_opt'),[],true,false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=0; base.BreakOnODEUnstable=0; base.SaveBest=0;

LoadData;                                    % ATP_c, Data_ATP
G={'incumbent tau_reg=26.6ms', base.tau_reg; 'fast gate tau_reg=0.2ms', 2e-4};
R=struct('name',{},'v',{},'F',{},'fn',{},'pw',{});
for g=1:2
    p=base; p.tau_reg=G{g,2}; p.FV_velocities=VGRID; p=getParams(p,p.g,true);
    [~, outs] = runFVExperiment(p, ATP_c, Data_ATP);
    F = arrayfun(@(o) o.Force(end), outs);
    R(g).name=G{g,1}; R(g).v=-VGRID; R(g).F=F;
    R(g).fn=F/F(1); R(g).pw=(F/F(1)).*(-VGRID);
end

% reference data (the two sets runFVExperiment averages)
AllV=[0 0.5 1 2 3 4 5 6];
D1=[56.40,51.8120,37.4459,17.8025,11.4430,6.2643,3.2759,2.2120];   % Baker2022
D2=[67.3942,60.9885,42.7059,18.8539,12.5956,6.4426,3.5136,1.8556]; % IsovelocityForFilip2021
SCORED=[0 0.5 1 2 4];

fprintf('\n============ FORCE-VELOCITY, FINE GRID (8 mM) ============\n');
fprintf('F_iso: incumbent %.2f kPa   fast %.2f kPa   (data %.1f / %.1f)\n', R(1).F(1), R(2).F(1), D1(1), D2(1));
fprintf('\n%5s %2s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s %8s\n', ...
  'v','sc','F inc','F/F0','norm pw','F fast','F/F0','norm pw','d1 F/F0','d2 F/F0','d1 pw','d2 pw');
for i=1:numel(VGRID)
    v=-VGRID(i); j=find(AllV==v,1);
    if isempty(j); d1='       -'; d2='       -'; p1='       -'; p2='       -';
    else
        d1=sprintf('%8.3f',D1(j)/D1(1)); d2=sprintf('%8.3f',D2(j)/D2(1));
        p1=sprintf('%8.3f',D1(j)/D1(1)*v); p2=sprintf('%8.3f',D2(j)/D2(1)*v);
    end
    fprintf('%5.2f %2s | %8.2f %8.3f %8.3f | %8.2f %8.3f %8.3f | %s %s %s %s\n', ...
        v, repmat('*',1,any(abs(SCORED-v)<1e-9)), ...
        R(1).F(i), R(1).fn(i), R(1).pw(i), R(2).F(i), R(2).fn(i), R(2).pw(i), d1,d2,p1,p2);
end

% --- summary numbers ---
fprintf('\n---- peak normalised power and Vmax ----\n');
for g=1:2
    [pm,im]=max(R(g).pw);
    iz=find(R(g).fn<=0.02,1);
    if isempty(iz); vmax=NaN; else; vmax=interp1(R(g).fn(iz-1:iz),R(g).v(iz-1:iz),0.02); end
    fprintf('%-24s peak norm power %.3f at %.2f ML/s   F/F0=0.02 at %s ML/s\n', ...
        R(g).name, pm, R(g).v(im), num2str(vmax,'%.2f'));
end
for q=1:2
    Dq=(q==1)*D1+(q==2)*D2; P=Dq/Dq(1).*AllV; [pm,im]=max(P);
    fprintf('%-24s peak norm power %.3f at %.2f ML/s\n', sprintf('DATA set %d',q), pm, AllV(im));
end

% --- the scored FV cost, unchanged grid ---
fprintf('\n---- the SCORED FV terms (only v = 0, 0.5, 1, 2, 4) ----\n');
for g=1:2
    p=base; p.tau_reg=G{g,2}; p.EvalFeatures=1; p=getParams(p,p.g,true);
    [Efv,~,fm,fd]=runFVExperiment(p, ATP_c, Data_ATP);
    fprintf('%-24s E_fv=%.5f  model F/F0 %s\n', G{g,1}, Efv, mat2str(round(fm.FV_fnorm(:)',3)));
    if g==1; fprintf('%-24s            data  F/F0 %s\n','',mat2str(round(fd.FV_fnorm(:)',3))); end
    fprintf('%-24s            model normpw %s\n','',mat2str(round(fm.FV_normpowerAvg(:)',3)));
    if g==1; fprintf('%-24s            data  normpw %s +/- %s\n','',mat2str(round(fd.FV_normpowerAvg(:)',3)),mat2str(round(sqrt(fd.FV_normpowerVar(:))',3))); end
end
save(fullfile(resdir,'fv_finegrid.mat'),'R','VGRID','AllV','D1','D2','SCORED');
