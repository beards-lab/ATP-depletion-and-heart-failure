% StateDiag.m -- what the state vector looks like at the start of each
% redevelopment window, and how the gates move during it.
% Needs results/baseline*.mat from RunKtrVsSlack.m.
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
L = load(fullfile(resdir,'baseline_outs.mat'));
os = L.os; ok = L.ok; svt = L.svt; kvt = L.kvt;

p  = os.params;  Ns = p.NumberOfStates; ss = p.ss;
iA = Ns*ss + 8;                       % A_reg index
fprintf('Ns=%d ss=%d  size(PU)=%s  A_reg idx=%d  A0=%.4f v_ref_reg=%.4f tau_reg=%.4f ms\n', ...
        Ns, ss, mat2str(size(os.PU)), iA, p.A0, p.v_ref_reg, p.tau_reg*1e3);
fprintf('slack out t: %.3f..%.3f (n=%d)   ktr out t: %.3f..%.3f (n=%d)\n', ...
        os.t(1), os.t(end), numel(os.t), ok.t(1), ok.t(end), numel(ok.t));

% ---- window onsets -------------------------------------------------------
iS = find(svt(:,2) < -1);
tS = svt(iS+1,1);          % end of each release ramp = start of the hold
tK = kvt(end-1,1);         % end of the final 1.05->1.00 release = start of redevelopment
tKpre = kvt(2,1);          % onset of the ktr shortening

% true F_iso: last isometric sample before the first perturbation
FisoM = interp1(os.t, os.Force, svt(iS(1),1)-0.002,'linear','extrap');
FisoK = interp1(ok.t, ok.Force, tKpre-0.002,'linear','extrap');
fprintf('\nF_iso  slack = %.2f kPa    ktr = %.2f kPa\n', FisoM, FisoK);

names = {'postSlack1','postSlack2','postSlack3','postSlack4','postSlack5','KTR'};
P = [arrayfun(@(k) probe(os,tS(k),p,FisoM), 1:numel(tS)), probe(ok,tK,p,FisoK)];

f = {'F','att','p1','p2','PD','SR','SRD','NP','A_reg','kaeff','RD1','SL'};
fprintf('\n============ STATE AT THE START OF THE REDEVELOPMENT ============\n');
fprintf('%-8s', 'field'); fprintf('%11s', names{:}); fprintf('\n');
for k=1:numel(f)
    fprintf('%-8s', f{k}); fprintf('%11.4f', arrayfun(@(s) s.(f{k}), P)); fprintf('\n');
end

% ---- gate trajectories during the window --------------------------------
WIN = 0.060;
tq  = (0:2e-4:WIN)';
Ps  = arrayfun(@(t) probe(os,t,p,FisoM), tS(2)+tq);
Pk  = arrayfun(@(t) probe(ok,t,p,FisoK), tK   +tq);
fprintf('\n===== TRAJECTORY  (postSlack2 | KTR) during redevelopment =====\n');
fprintf('%6s %s\n','t[ms]','   F      A_reg    kaeff     PD       SR+SRD   att   |    F      A_reg    kaeff     PD       SR+SRD   att');
for i = 1:round(0.004/2e-4):numel(tq)
    fprintf('%6.1f %8.3f %8.4f %8.2f %8.4f %8.4f %7.4f |%8.3f %8.4f %8.2f %8.4f %8.4f %7.4f\n', tq(i)*1e3, ...
      Ps(i).F,Ps(i).A_reg,Ps(i).kaeff,Ps(i).PD,Ps(i).SR+Ps(i).SRD,Ps(i).att, ...
      Pk(i).F,Pk(i).A_reg,Pk(i).kaeff,Pk(i).PD,Pk(i).SR+Pk(i).SRD,Pk(i).att);
end

% ---- what A_reg did during the manoeuvre itself -------------------------
fprintf('\n===== A_reg THROUGH THE MANOEUVRE =====\n');
fprintf('ktr protocol velocity table (t, v ML/s, ML):\n');
for i=2:size(kvt,1); fprintf('  t=%.4f  v=%8.2f  ML=%.4f  A_reg=%.4f\n', kvt(i,1), kvt(i,2), kvt(i,4), ...
        interp1(ok.t,ok.PU(:,iA),kvt(i,1),'linear','extrap')); end
fprintf('slack cycle 2 (t, v ML/s, ML):\n');
for i=iS(2):iS(2)+2; fprintf('  t=%.4f  v=%8.2f  ML=%.4f  A_reg=%.4f\n', svt(i,1), svt(i,2), svt(i,4), ...
        interp1(os.t,os.PU(:,iA),svt(i,1),'linear','extrap')); end
save(fullfile(resdir,'statediag.mat'),'P','Ps','Pk','tq','names','FisoM','FisoK');

% ---- state probe (local fn) ---------------------------------------------
function s = probe(o,t,p,Fiso)
    Ns=p.NumberOfStates; ss=p.ss; iA=Ns*ss+8;
    g=@(v) interp1(o.t(:),v(:),t,'linear','extrap');
    s.F     = g(o.Force)/Fiso;
    s.FXB   = g(o.FXB)/Fiso;
    s.p1    = g(o.p1_0);  s.p2 = g(o.p2_0);
    s.att   = s.p1+s.p2;
    s.PD    = g(o.PU(:,Ns*ss+5));
    s.SR    = g(o.SR);    s.SRD = g(o.SRD);
    s.NP    = g(o.NP);
    s.PuR   = g(o.PuR);   s.PuATP = g(o.PuATP);
    s.A_reg = g(o.PU(:,iA));
    s.RD1   = g(o.RD1);            % attachment flux  (1/s of total prob)
    s.R1D   = g(o.R1D); s.R12=g(o.R12); s.R21=g(o.R21); s.R2T=g(o.R2T);
    s.RPT2SR= g(o.RPT2SR); s.RSR2PT = g(o.RSR2PT);
    s.OV    = g(o.OV); s.flat = g(o.f_lattice); s.fsat = g(o.f_saturation);
    s.SL    = g(o.SL);
    s.kaeff = p.ka * s.A_reg * s.OV * s.flat * s.fsat;   % per-unit-PD attachment rate
end

