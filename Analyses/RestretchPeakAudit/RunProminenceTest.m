% RunProminenceTest.m -- is 04/10's peak1_dSL/vall_t collapse caused by the
% ABSOLUTE prominence threshold?
%
% extractSlackAttributes.m:31 sets PEAK_MIN_PROMINENCE = 0.1 kPa, a fixed
% ABSOLUTE value, and applies it to preparations spanning 1.37x in force
% (Fss22 80.1 / 70.0 / 50.9 kPa). On the weakest prep the real overshoot is
% smaller in kPa and the detector can latch onto an earlier ripple instead.
% Test: rescale the threshold by each preparation's own force and see whether
% 04/10's values join the other two.
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
DAYS={'03_27','04_03','04_10'};
fprintf('\n%-8s %-22s %s\n','day','prominence rule','peak1_dSL per segment');
for scaleIt = [0 1]
    for d=1:3
        L=load(fullfile(root,'data',sprintf('protocol_%s_2026_8mM_slack.mat',DAYS{d})));
        dt=L.datatable; v=L.velocitytable; iS=find(v(:,2)<-1);
        F0=median(dt(dt(:,1)>v(iS(1),1)-0.30 & dt(:,1)<v(iS(1),1)-0.01,3));
        prom = 0.1; rule='absolute 0.1 kPa';
        if scaleIt; prom = 0.1*F0/80.0; rule=sprintf('scaled: %.3f kPa',prom); end
        out=nan(1,5);
        for s=1:5
            iR=find(v(:,2)>1);
            m=dt(:,1)>=v(iR(s),1) & dt(:,1)<=v(iR(s)+1,1);
            t=dt(m,1); y=dt(m,3); sl=dt(m,2);
            [py,pt]=findpeaks(y,t,'MinPeakDistance',0.001,'MinPeakProminence',prom);
            if isempty(py); continue; end
            i=find(t>=pt(1),1); out(s)=sl(i)-sl(end);
        end
        fprintf('%-8s %-22s %s\n', DAYS{d}, rule, mat2str(round(out,4)));
    end
    fprintf('\n');
end
fprintf('NOTE: this uses the ramp window, which is NOT extractSlackAttributes''\n');
fprintf('window, so absolute values differ from the stored features. The point is\n');
fprintf('whether 04/10 behaves like the other two once the threshold is scaled.\n');
fprintf('\nDONE\n');
