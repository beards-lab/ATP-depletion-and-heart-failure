% RunInstrument0410.m -- which threshold misfires, exactly.
%
% Reproduces extractSlackAttributes' peak detection verbatim:
%   window   velocity_segment(3) .. velocity_segment(4)   = the restretch ramp
%   detector findpeaks(raw y, MinPeakDistance 1 ms, MinPeakProminence 0.1 kPa)
%   taken    peak1_t(1) -- the FIRST peak in time, not the largest
%   feature  peak1_dSL = SL(t_peak) - SLslack   (SLslack from the hold window)
%
% Lists EVERY peak the detector returns, with its prominence, so it is visible
% which one is being taken and why.
%
%   matlab -batch "run('Analyses/RestretchPeakAudit/RunInstrument0410.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','RestretchPeakAudit','results');
if ~exist(resdir,'dir'); mkdir(resdir); end
DAYS={'03_27','04_03','04_10'};
MIND=0.001; PROM=0.1;                       % extractSlackAttributes.m:30-31

S=struct();
for d=1:3
    L=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
    S(d).dt=L.datatable; S(d).vt=L.velocitytable; S(d).fd=L.features_data;
    v=S(d).vt; iS=find(v(:,2)<-1);
    q=S(d).dt(S(d).dt(:,1)>v(iS(1),1)-0.30 & S(d).dt(:,1)<v(iS(1),1)-0.01,3);
    S(d).sigma=std(q-movmean(q,9));
    S(d).F0=median(q);
end
fprintf('\nnoise sigma vs the prominence threshold:\n');
for d=1:3
    fprintf('  %-6s  F_iso %5.1f kPa   sigma %.3f kPa   threshold %.2f kPa = %.2f sigma\n', ...
            DAYS{d}, S(d).F0, S(d).sigma, PROM, PROM/S(d).sigma);
end
fprintf('\n  -> the threshold is BELOW the noise on every preparation, so findpeaks\n');
fprintf('     returns noise ripples and peak1_t(1) takes the earliest one.\n');

%% ---- every peak the detector sees ----
fprintf('\n=============== ALL PEAKS RETURNED, segment by segment ===============\n');
for d=1:3
    v=S(d).vt; dt=S(d).dt; iS=find(v(:,2)<-1);
    fprintf('\n--- %s (sigma %.2f kPa) ---\n', DAYS{d}, S(d).sigma);
    for s=1:5
        vs = v(iS(s):iS(s)+4, :);
        wh = dt(:,1) >= vs(2,1) & dt(:,1) <= vs(3,1);           % the slack hold
        SLslack = dt(find(wh,1,'last'),2);
        w  = dt(:,1) >= vs(3,1) & dt(:,1) <= vs(4,1);           % the restretch ramp
        t=dt(w,1)-vs(1,1); y=dt(w,3); sl=dt(w,2);
        [py,pt,~,pp]=findpeaks(y,t,'MinPeakDistance',MIND,'MinPeakProminence',PROM);
        if isempty(py); fprintf('  seg%d: NO PEAKS\n',s); continue; end
        i=find(t>=pt(1),1);
        dsl = sl(i)-SLslack;
        [~,ibig]=max(py);
        fprintf('  seg%d: %d peaks | TAKEN #1 at %+.2f ms prom %.2f -> dSL %.4f | LARGEST #%d at %+.2f ms prom %.2f -> dSL %.4f\n', ...
            s, numel(py), 1e3*(pt(1)-(vs(3,1)-vs(1,1))), pp(1), dsl, ibig, ...
            1e3*(pt(ibig)-(vs(3,1)-vs(1,1))), pp(ibig), sl(find(t>=pt(ibig),1))-SLslack);
    end
end

%% ---- sweep the prominence threshold ----
fprintf('\n=============== SWEEP: prominence in units of each prep''s sigma ===============\n');
KS=[0.07 0.5 1 2 3 4 6];
P=nan(3,5,numel(KS));
for ik=1:numel(KS)
    for d=1:3
        v=S(d).vt; dt=S(d).dt; iS=find(v(:,2)<-1);
        for s=1:5
            vs=v(iS(s):iS(s)+4,:);
            wh=dt(:,1)>=vs(2,1)&dt(:,1)<=vs(3,1); SLslack=dt(find(wh,1,'last'),2);
            w=dt(:,1)>=vs(3,1)&dt(:,1)<=vs(4,1);
            t=dt(w,1)-vs(1,1); y=dt(w,3); sl=dt(w,2);
            [py,pt]=findpeaks(y,t,'MinPeakDistance',MIND,'MinPeakProminence',KS(ik)*S(d).sigma);
            if isempty(py); continue; end
            P(d,s,ik)=sl(find(t>=pt(1),1))-SLslack;
        end
    end
end
fprintf('%8s %10s  %s\n','k*sigma','rel SD','peak1_dSL: 03/27 | 04/03 | 04/10   (segment means)');
for ik=1:numel(KS)
    A=P(:,:,ik);
    rs=mean(std(A,0,1)./abs(mean(A,1)),'omitnan');
    fprintf('%8.2f %10.3f  %s | %s | %s\n', KS(ik), rs, ...
        mat2str(round(A(1,:),4)), mat2str(round(A(2,:),4)), mat2str(round(A(3,:),4)));
end
fprintf('\n(k = 0.07 is roughly the shipped 0.1 kPa on 03/27; the shipped value is\n');
fprintf(' a fixed 0.1 kPa, i.e. a DIFFERENT k on every preparation)\n');

%% ---- figure ----
f=figure('Position',[20 20 1600 900],'Color','w');
tiledlayout(3,5,'Padding','compact','TileSpacing','compact');
CC=lines(3);
for d=1:3
    v=S(d).vt; dt=S(d).dt; iS=find(v(:,2)<-1);
    for s=1:5
        nexttile; hold on; box on; grid on
        vs=v(iS(s):iS(s)+4,:);
        wh=dt(:,1)>=vs(2,1)&dt(:,1)<=vs(3,1); SLslack=dt(find(wh,1,'last'),2);
        w=dt(:,1)>=vs(3,1)-0.001 & dt(:,1)<=vs(4,1)+0.002;
        t=(dt(w,1)-vs(3,1))*1e3; y=dt(w,3);
        plot(t,y,'-','Color',CC(d,:),'LineWidth',1.0);
        w2=dt(:,1)>=vs(3,1)&dt(:,1)<=vs(4,1);
        t2=dt(w2,1)-vs(1,1); y2=dt(w2,3); sl2=dt(w2,2);
        [py,pt]=findpeaks(y2,t2,'MinPeakDistance',MIND,'MinPeakProminence',PROM);
        if ~isempty(py)
            plot((pt-(vs(3,1)-vs(1,1)))*1e3, py,'x','Color',[.6 .6 .6],'MarkerSize',7,'LineWidth',1);
            plot((pt(1)-(vs(3,1)-vs(1,1)))*1e3, py(1),'rv','MarkerFaceColor','r','MarkerSize',9);
        end
        [py3,pt3]=findpeaks(y2,t2,'MinPeakDistance',MIND,'MinPeakProminence',3*S(d).sigma);
        if ~isempty(py3)
            plot((pt3(1)-(vs(3,1)-vs(1,1)))*1e3, py3(1),'gs','MarkerFaceColor','g','MarkerSize',9);
        end
        if s==1; ylabel(sprintf('%s\nF (kPa)',strrep(DAYS{d},'_','/'))); end
        if d==3; xlabel('t since restretch (ms)'); end
        title(sprintf('seg %d',s),'FontSize',9);
        if d==1 && s==1
            legend({'force','all peaks','TAKEN (0.1 kPa)','taken at 3\sigma'}, ...
                   'Box','off','FontSize',7,'Location','southeast');
        end
    end
end
sgtitle('extractSlackAttributes peak detection: red = what it takes now, green = at a 3\sigma threshold','FontSize',12);
exportgraphics(f, fullfile(resdir,'instrument_0410.png'),'Resolution',110);
fprintf('\nwrote %s\n', fullfile(resdir,'instrument_0410.png'));
save(fullfile(resdir,'instrument_0410.mat'),'P','KS','S');
fprintf('\nDONE\n');
