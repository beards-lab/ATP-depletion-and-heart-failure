root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
cl=parcluster('local'); jsl=fullfile(tempdir,sprintf('mj_bk_%d',feature('getpid')));
if ~isfolder(jsl); mkdir(jsl); end; cl.JobStorageLocation=jsl;
if isempty(gcp('nocreate')); parpool(cl,5); end
T={'refit3_seed','refit3_opt'}; C=cell(1,2);
for k=1:2
  p=getParams(loadParams(T{k}),[],true,false);
  p.PlotEachSeparately=0;p.ghostSave='';p.ghostLoad='';p.SaveBest=0;
  p.EvalFeatures=1;p.BreakOnODEUnstable=0;p.MaxRunTime=600;p=getParams(p,p.g,true);
  [c,~,fc,fn]=evaluateBakersExp([],p,false,false); C{k}=fc; FN=fn;
  ps=p; ps.RunStairs=1; ps=getParams(ps,ps.g,true);
  fprintf('%-14s cost %.4f  stairs(heldout) %.4f\n',T{k},c,runStairsExperiment(ps));
end
lab=cellfun(@(s) regexprep(s,'\|.*$',''),FN,'UniformOutput',false);
d=C{2}-C{1}; [~,o]=sort(C{2},'descend');
fprintf('\n%-24s %9s %9s %9s\n','feature','seed','opt','delta');
for i=o(:)'; if C{1}(i)<1e-6 && C{2}(i)<1e-6; continue; end
  fprintf('%-24s %9.4f %9.4f %+9.4f\n',lab{i},C{1}(i),C{2}(i),d(i)); end
fprintf('%-24s %9.4f %9.4f %+9.4f\n','TOTAL',sum(C{1}),sum(C{2}),sum(d));
V=validateKtrProtocol('refit3_opt');
fprintf('\nDONE\n');
