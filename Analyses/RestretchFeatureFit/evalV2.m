% evalV2.m — confirm the expanded-optim result (opt2state_v2_opt.m).
clear; clc;
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root); addpath(genpath(root));
SNAP = 'params/opt2state_v2_opt.m';
params0 = getParams(); run(SNAP); FN = params0.fn;
[tc, cv, fm, fd] = costOfSnap(SNAP, FN, 120);
fprintf('\n==== evalV2: %s ====\n', SNAP);
fprintf('TOTAL cost = %.4f  (was 4.55 at opt2state_opt.m)\n', tc);
fprintf('kmsrd = %.3f\n', params0.kmsrd);
fprintf('ktr   model = [%s]  (data ~49)\n', num2str(fm.ktr(:)', '%7.1f'));
fprintf('steady model = [%s]  (data ~[80 80 81 81 64])\n', num2str(fm.steady(:)', '%7.1f'));
fprintf('peak1_y model = [%s]  (data ~[96 94 90 89 77])\n', num2str(fm.peak1_y(:)', '%7.1f'));
fprintf('FV_fnorm model = [%s]  (data [1 .92 .66 .32 .11])\n', num2str(fm.FV_fnorm(:)', '%6.3f'));
fprintf('XTOR=%.2f XTOR_vmax=%.2f SRX_ss=%.3f att_ss=%.3f\n', ...
    mean(fm.XTOR,'omitnan'), mean(fm.XTOR_vmax,'omitnan'), mean(fm.SRX_ss,'omitnan'), mean(fm.attached_ss,'omitnan'));
% top residual terms
[~,ord] = sort(cv,'descend','MissingPlacement','last');
fprintf('top residual terms:\n');
for i=1:5; j=ord(i); fprintf('  %-40s %.4f\n', FN{j}, cv(j)); end
disp('DONE evalV2');
