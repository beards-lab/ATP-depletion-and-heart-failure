% RunSpreadReconcile.m -- my re-extraction says peak1_dSL's between-prep spread
% is 4.3 %; pooledTarget8mM says 52 %. Which, and why?
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
DAYS={'03_27','04_03','04_10'}; PT=pooledTarget8mM();
F={'peak1_dSL','vall_t','peak1_y','A','steady','ktr'};
fprintf('\n%-12s %10s %10s   %s\n','feature','PT.s','PT.sraw','per-prep values from each file''s stored features_data');
V=struct();
for q=F
    j=find(strcmp(PT.name,q{1}),1);
    vals=nan(3,5);
    for d=1:3
        fd=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d}))).features_data;
        if isfield(fd,q{1}); v=fd.(q{1})(:)'; vals(d,1:min(5,numel(v)))=v(1:min(5,numel(v))); end
    end
    V.(q{1})=vals;
    s_mine = mean(std(vals,0,1)./abs(mean(vals,1)),'omitnan');
    if isempty(j); sPT=NaN; sRaw=NaN; else
        sPT=PT.s(j);
        if isfield(PT,'sraw'); sRaw=PT.sraw(j); else; sRaw=NaN; end
    end
    fprintf('%-12s %10.3f %10.3f   stored-feature spread %.3f\n', q{1}, sPT, sRaw, s_mine);
    fprintf('%-12s %10s %10s   03/27 %s\n','','','',mat2str(round(vals(1,:),4)));
    fprintf('%-12s %10s %10s   04/03 %s\n','','','',mat2str(round(vals(2,:),4)));
    fprintf('%-12s %10s %10s   04/10 %s\n','','','',mat2str(round(vals(3,:),4)));
end
fprintf('\nPT.seg for peak1_dSL (which cycles s was measured on): ');
j=find(strcmp(PT.name,'peak1_dSL'),1);
if isfield(PT,'seg'); disp(PT.seg{j}); else; fprintf('(no seg field)\n'); end
fprintf('\nDONE\n');
