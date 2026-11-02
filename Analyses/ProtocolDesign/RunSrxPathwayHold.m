% RunSrxPathwayHold.m -- does the release-hold-restretch series separate SRX ENTRY (ksr0, DT->SRX.T)
% from SRX EXIT (kmsrd, SRX.D->DD)? Uncompensated kicks (x3 entry vs /3 exit) on robust1_seed, 8 mM.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
p0=getParams(loadParams('robust1_seed'),[],true,false); p0.MaxRunTime=300; p0.ghostLoad='';p0.ghostSave='';p0.SaveBest=0;p0.BreakOnODEUnstable=false;
mf=str2func(p0.modelFcn); T=[5 10 20 40 80 160]*1e-3; dL=0.15;
V={'ref',1,1;'entry ksr0 x3',3,1;'exit kmsrd /3',1,1/3};
fid=fopen(fullfile(root,'Analyses','ProtocolDesign','results','srx_pathway.txt'),'w');
for v=1:3
  k=nan(size(T));
  for i=1:numel(T)
    x=T(i); q=p0; q.ksr0=q.ksr0*V{v,2}; q.kmsrd=q.kmsrd*V{v,3};
    vt=[-2 0;0 -dL/1e-3;1e-3 0;1e-3+x dL/2e-3;3e-3+x 0;0.403+x 0];
    q.SL0=2.2; q.Velocity=vt(:,2); q=getParams(q,q.g,true,false);
    try
      [~,o]=evaluateModel(mf,vt(:,1),q); kk=makeMonotonous(o.t(:)); t=o.t(:); t=t(kk); F=o.Force(:); F=F(kk);
      tr=3e-3+x; m=find(t>=tr&t<=tr+15e-3); [Fmin,j]=min(F(m)); tmin=t(m(j)); Fend=interp1(t,F,vt(end,1)-1e-3);
      lvl=Fmin+0.63*(Fend-Fmin); ii=find(F>=lvl&t>tmin,1); k(i)=1/(t(ii)-tmin);
      if i==1, Fiso=interp1(t,F,-1e-3); end
    catch, end
  end
  fprintf(fid,'%-16s Fiso %6.1f | k63 %s | shape k/k(160) %s\n',V{v,1},Fiso,mat2str(k,3),mat2str(k/k(end),2));
end
fclose(fid); type(fullfile(root,'Analyses','ProtocolDesign','results','srx_pathway.txt'));
