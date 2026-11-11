root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
cl=parcluster('local'); jsl=fullfile(tempdir,sprintf('mj_c4_%d',feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end; cl.JobStorageLocation=jsl;
if isempty(gcp('nocreate')); parpool(cl,5); end
for T={'refit4_seed','refit4_opt'}
  if isempty(which([T{1} '.m'])) && ~isfile(fullfile(root,'params',[T{1} '.m'])); fprintf('%s MISSING\n',T{1}); continue; end
  p=getParams(loadParams(T{1}),[],true,false);
  p.PlotEachSeparately=0;p.ghostSave='';p.ghostLoad='';p.SaveBest=0;
  p.EvalFeatures=1;p.BreakOnODEUnstable=0;p.MaxRunTime=600;p=getParams(p,p.g,true);
  a=evaluateBakersExp([],p,false,false); b=evaluateBakersExp([],p,false,false);
  ps=p; ps.RunStairs=1; ps=getParams(ps,ps.g,true);
  fprintf('%-14s %.6f  (repeat %.6f, delta %.1e)  stairs(heldout) %.4f\n',T{1},a,b,abs(a-b),runStairsExperiment(ps));
end
V=validateKtrProtocol('refit4_opt');
fprintf('\nDONE\n');
