% RunDataSpread.m — what the staircase data can actually resolve, and how much
% of it is preparation-to-preparation scatter rather than signal.
%
% Two products:
%
%  (A) MEASUREMENT precision. Per-sample noise sigma is measured inside the
%      isometric holds; every feature is a window mean, so its standard error
%      follows from sigma and the window length. This is the floor a single
%      preparation can be fitted to.
%
%  (B) BETWEEN-PREPARATION spread, over the three 2026 protocol days. The days
%      do NOT share a staircase: 03/27 is 7x1% + 1x3%, while 04/03 and 04/10
%      are 1x1% + 3x3%. So only step-size-matched quantities are pooled:
%      the 1% step response (3 observations), the 3% step response (7), and the
%      2.0 -> 2.2 um envelope (3). Following Analyses/FeatureWeighting, a
%      weight is lambda/s^2 with s the relative between-prep SD.
%
%  (C) A model-free read on the force-strain exponent: the ratio of the 3% to
%      the 1% elastic peak. A Hookean cross-bridge whose resting strain is small
%      against the step gives 3^e; e > 1 needs a ratio above 3, e < 1 below it.
%
% Outputs -> results/data_spread.mat + .png

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resdir = fullfile(here,'results'); if ~isfolder(resdir); mkdir(resdir); end

files = {'protocol_03_27_2026_velocitytable_stairs.mat', ...
         'protocol_04_03_2026_velocitytable_stairs.mat', ...
         'protocol_04_10_2026_velocitytable_stairs.mat'};
days  = {'03/27 M','04/03 F','04/10 M'};

D = struct();
fprintf('=== (A) measurement precision ===\n');
for i = 1:numel(files)
    q  = load(['data/' files{i}]);
    dt = q.datatable; vt = q.velocitytable;
    f  = extractStairsAttributes(dt(:,1), dt(:,3), vt, 'st_', struct(), true);

    % per-sample noise from in-hold detrended residuals
    res = [];
    for r = 1:size(vt,1)-1
        if vt(r,2) ~= 0; continue; end
        m = dt(:,1) >= vt(r,1)+0.005 & dt(:,1) <= vt(r+1,1)-0.002;
        if nnz(m) < 20; continue; end
        p = polyfit(dt(m,1), dt(m,3), 1);
        res = [res; dt(m,3) - polyval(p, dt(m,1))]; %#ok<AGROW>
    end
    sg = std(res);
    fs = 1/median(diff(dt(:,1)));
    % window sample counts: peak 0.30*Tramp, base/plat 18 ms, nadir 10 ms
    Tr   = median(diff(vt(find(vt(:,2)~=0 & (1:size(vt,1))'<size(vt,1),1),:)));
    nPk  = max(round(0.30*0.0195*fs),2); nB = round(0.018*fs); nNd = round(0.010*fs);
    sePk = sg*sqrt(1/nPk + 1/nB);  seNd = sg*sqrt(1/nNd + 1/nB);  sePl = sg*sqrt(2/nB);

    D(i).day=days{i}; D(i).f=f; D(i).sigma=sg; D(i).dSL=f.st_dSL;
    D(i).sePeak=sePk; D(i).seNadir=seNd; D(i).sePlat=sePl;
    fprintf('%s  sigma=%.4f  SE(peak)=%.4f SE(nadir)=%.4f SE(plat)=%.4f  | LSQ floor 10*sigma^2 = %.5f\n', ...
        days{i}, sg, sePk, seNd, sePl, 10*sg^2);
end

%% ---- (B)/(C) step-size matched pooling ---------------------------------
one = []; thr = []; oneD = []; thrD = []; env = []; recO = []; recT = []; nadO = []; nadT = [];
for i = 1:numel(D)
    f = D(i).f; s = f.st_dSL;
    i1 = find(abs(s-0.010) < 0.003);       % 1 % steps
    i3 = find(abs(s-0.030) < 0.006);       % 3 % steps
    one = [one, f.st_peak(i1)];  oneD = [oneD, repmat(i,1,numel(i1))];
    thr = [thr, f.st_peak(i3)];  thrD = [thrD, repmat(i,1,numel(i3))];
    nadO = [nadO, f.st_nadir(i1)]; nadT = [nadT, f.st_nadir(i3)];
    recO = [recO, f.st_rec(i1)];   recT = [recT, f.st_rec(i3)];
    env  = [env, f.st_plat(end)];
end

fprintf('\n=== (B) between-preparation spread (step-size matched) ===\n');
rep = @(nm,v,n) fprintf('%-22s n=%2d  mean=%8.4f  SD=%7.4f  relSD=%5.1f%%  ->  w=lambda/s^2 = %8.1f*lambda\n', ...
        nm, n, mean(v), std(v), 100*std(v)/abs(mean(v)), 1/(std(v)/abs(mean(v)))^2);
rep('peak @1% step',   one, numel(one));
rep('peak @3% step',   thr, numel(thr));
rep('nadir @1% step',  nadO, numel(nadO));
rep('nadir @3% step',  nadT, numel(nadT));
rep('recovery @1% step',recO, numel(recO));
rep('recovery @3% step',recT, numel(recT));
rep('envelope 2.0->2.2', env, numel(env));

fprintf('\n  per-day 1%% peaks: ');  fprintf('%.4f ', one); fprintf('\n  per-day 3%% peaks: '); fprintf('%.4f ', thr); fprintf('\n');

fprintf('\n=== (C) model-free exponent read ===\n');
R = mean(thr)/mean(one);
fprintf('peak(3%%)/peak(1%%) = %.3f   [a linear spring far from rest gives 3.00]\n', R);
fprintf('  -> apparent amplitude exponent  log(R)/log(3) = %.3f\n', log(R)/log(3));
fprintf('  (a LOWER bound on the elastic exponent: the 3%% ramp is 3x faster, so\n');
fprintf('   forced detachment during it pushes this ratio DOWN as well. The purely\n');
fprintf('   elastic ratio is therefore >= 1.14 and the exponent >= 0.12 -- but both\n');
fprintf('   are still far below the 3.00 / 1.00 a Hookean bridge would give.)\n');

% within-day paired ratio for 03/27 (same prep, adjacent steps)
f1 = D(1).f;
fprintf('  03/27 paired, step8(3%%)/step7(1%%) = %.3f\n', f1.st_peak(8)/f1.st_peak(7));
fprintf('  peak SLOPE across the seven 1%% steps = %+.5f /step (%+.1f%% per step)\n', ...
        f1.st_peakSlope, 100*f1.st_peakSlope/f1.st_peakMean);

save(fullfile(resdir,'data_spread.mat'), 'D','one','thr','env','days','R');
fprintf('\nsaved -> %s\n', fullfile(resdir,'data_spread.mat'));
