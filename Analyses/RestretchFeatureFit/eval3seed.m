% eval3seed.m — single de-risking evaluation of the 3-state seed.
% Checks: does it run (no hang), does FV collapse, what do features look like.
% Run headless via -batch (fresh process; picks up the edited R3 strain-dep).
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));
SNAP   = 'params/params_3state_seedB.m';
REPORT = fullfile('Analyses','RestretchFeatureFit','eval3seed_report.txt');
fid = fopen(REPORT, 'w');
params0 = getParams(); run(SNAP); FN = params0.fn;
fprintf(fid, '3-state seedB: NoS=%d UseKstiff3=%d alpha3=%.0f s3=%.4g k3=%.1f kstiff3=%.0f drp3=%.3f\n', ...
    params0.NumberOfStates, params0.UseKstiff3, params0.alpha3, params0.s3, params0.k3, params0.kstiff3, params0.drp3);
t0 = tic;
try
    [tc, cv, fm, fd] = costOfSnap(SNAP, FN, 60);
    el = toc(t0);
    fprintf(fid, 'RAN OK: cost=%.4f  elapsed=%.0fs\n', tc, el);
    fprintf(fid, 'ktr      = [%s]  (data ~49)\n', num2str(fm.ktr(:)', '%7.1f'));
    fprintf(fid, 'FV_fnorm = [%s]  (data [1 .92 .66 .32 .11]; COLLAPSE if tail ~0)\n', num2str(fm.FV_fnorm(:)', '%6.3f'));
    fprintf(fid, 'peak1_y  = [%s]  (data ~[96 94 90 89 77])\n', num2str(fm.peak1_y(:)', '%7.1f'));
    fprintf(fid, 'vall_y   = [%s]  (data ~[77 75 72 70 61])\n', num2str(fm.vall_y(:)', '%7.1f'));
    fprintf(fid, 'steady   = [%s]  (data ~[80 80 81 81 64])\n', num2str(fm.steady(:)', '%7.1f'));
    fprintf(fid, 'XTOR=%.2f XTOR_vmax=%.2f SRX_ss=%.3f att_ss(p1+p2)=%.3f\n', ...
        mean(fm.XTOR,'omitnan'), mean(fm.XTOR_vmax,'omitnan'), mean(fm.SRX_ss,'omitnan'), mean(fm.attached_ss,'omitnan'));
    [~,ord] = sort(cv,'descend','MissingPlacement','last');
    fprintf(fid, 'top residuals:\n');
    for i=1:6; j=ord(i); fprintf(fid,'  %-42s %.4f\n', FN{j}, cv(j)); end
catch e
    fprintf(fid, 'ERROR after %.0fs: %s\n', toc(t0), e.message);
end
fclose(fid);
type(REPORT);
disp('DONE eval3seed');
