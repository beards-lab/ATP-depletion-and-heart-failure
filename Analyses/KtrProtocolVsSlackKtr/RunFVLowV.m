% RunFVLowV.m -- the fine grid puts F/F0 = 1.12 at 0.25 ML/s, i.e. the model's
% FV curve is NON-MONOTONIC between the two lowest scored points (0 and 0.5).
% Resolve it, and test whether A_reg's velocity dependence is the cause.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
VLOW = -[0, 0.005, 0.01, 0.02, 0.04, 0.05:0.05:0.60, 0.7, 0.8, 1.0];
base=getParams(loadParams('estiff125_opt'),[],true,false);
base.PlotEachSeparately=0; base.ghostSave=''; base.ghostLoad='';
base.EvalFeatures=0; base.BreakOnODEUnstable=0; base.SaveBest=0;
LoadData;

V={ 'incumbent  tau 26.6ms', @(p) p
    'fast gate  tau 0.2ms',  @(p) setfield(p,'tau_reg',2e-4)
    'A0=1 (gate velocity-blind, diagnostic)', @(p) setfield(setfield(p,'tau_reg',2e-4),'A0',1) };
L=struct('name',{},'F',{},'fn',{},'Ainf',{});
for g=1:size(V,1)
    p=V{g,2}(base); p.FV_velocities=VLOW; p=getParams(p,p.g,true);
    [~,outs]=runFVExperiment(p, ATP_c, Data_ATP);
    F=arrayfun(@(o) o.Force(end), outs);
    L(g).name=V{g,1}; L(g).F=F; L(g).fn=F/F(1);
    L(g).Ainf = p.A0 + (1-p.A0)*(2*-VLOW)./((2*-VLOW)+p.v_ref_reg);  % vel in um/s
end
fprintf('\n=========== FV BELOW THE FIRST SCORED SHORTENING POINT ===========\n');
fprintf('scored velocities are 0 and 0.5 ML/s -- everything between is unconstrained\n\n');
fprintf('%6s %8s | %8s %8s | %8s %8s | %8s %8s\n','v','A_inf', ...
   'F inc','F/F0','F fast','F/F0','F A0=1','F/F0');
for i=1:numel(VLOW)
    fprintf('%6.2f %8.3f | %8.2f %8.3f | %8.2f %8.3f | %8.2f %8.3f\n', -VLOW(i), L(1).Ainf(i), ...
        L(1).F(i),L(1).fn(i), L(2).F(i),L(2).fn(i), L(3).F(i),L(3).fn(i));
end
fprintf('\n---- peak of the hump ----\n');
for g=1:3
    [mx,im]=max(L(g).fn);
    fprintf('%-42s max F/F0 = %.3f at %.2f ML/s  (+%.1f%% over isometric)\n', ...
        L(g).name, mx, -VLOW(im), 100*(mx-1));
end
fprintf(['\nA_inf = A0 + (1-A0)*v/(v+v_ref_reg) with A0=%.3f, v_ref_reg=%.3f um/s.\n' ...
         'At 0.25 ML/s (=0.5 um/s) the gate is already %.2fx its isometric value,\n' ...
         'so attachment is boosted faster than sliding strips the bound population.\n'], ...
         base.A0, base.v_ref_reg, L(1).Ainf(6)/base.A0);
save(fullfile(resdir,'fv_lowv.mat'),'L','VLOW');
