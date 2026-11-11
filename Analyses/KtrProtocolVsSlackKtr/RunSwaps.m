% RunSwaps.m -- which piece of the state carries the model's ktr/post-slack
% ratio of 1.51 (data: 0.97)?
%
% Both windows are pure isometric holds at SL 2.0 um, so each is fully
% determined by its state vector at the onset. Replay both from their captured
% PU0 (seedWindow.m), then transplant ONE pool at a time from the post-slack
% donor into the ktr recipient and re-measure k63. PT is the slack variable
% (dPUdT_CombinedTransitions.m:222) so every single-pool swap conserves
% probability automatically.
root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','KtrProtocolVsSlackKtr','results');
L  = load(fullfile(resdir,'baseline_outs.mat'));
os = L.os; ok = L.ok; svt = L.svt; kvt = L.kvt;
p  = os.params; Ns = p.NumberOfStates; ss = p.ss;
iSR = Ns*ss+1; iNP = Ns*ss+2; iSL = Ns*ss+3; iLSE = Ns*ss+4;
iPD = Ns*ss+5; iSRD = Ns*ss+6; iXD = Ns*ss+7; iA = Ns*ss+8;
iAd = Ns*ss + 8 + p.UseRegistrationAvailability + p.UseD0State;   % A_def (UseLoadRealign)

iS  = find(svt(:,2) < -1);
tSlk = svt(iS(2)+1,1);      % start of the post-slack #2 hold
tKtr = kvt(end-1,1);        % start of the ktr redevelopment

grab = @(o,t) deal_(o,t);
[PUs, pivS, tS] = grab(os, tSlk);
[PUk, pivK, tK] = grab(ok, tKtr);
fprintf('donor  postSlack2 @ t=%.5f  SL=%.4f  pivot=%.5f\n', tS, PUs(iSL), pivS);
fprintf('recip. ktr        @ t=%.5f  SL=%.4f  pivot=%.5f\n', tK, PUk(iSL), pivK);

DUR = 0.0585;                              % the ktr window's own span
% plateaus: use each window's own measured plateau from the full protocol
Pk = 1.0 * interp1(ok.t,ok.Force,kvt(2,1),'linear','extrap');   % ktr: pre-release F_iso
tsl = os.t(:); Fsl = os.Force(:);
Ps  = median(Fsl(tsl>=tSlk+0.23 & tsl<=tSlk+0.28));

[~,k_S] = seedWindow(p, PUs, pivS, tS, 0.275, 0, Ps);
[~,k_K] = seedWindow(p, PUk, pivK, tK, DUR,   0, Pk);
fprintf('\nreplay control: postSlack2 %.1f /s   ktr %.1f /s   ratio %.2f\n', k_S, k_K, k_K/k_S);

% --------- single-pool transplants donor -> ktr recipient ---------
SW = { 'H0  control',            @(v) v
       'H1  A_reg   <- slack',   @(v) setix(v, iA,  PUs(iA))
       'H2  p1,p2   <- slack',   @(v) setix(v, 1:2*ss, PUs(1:2*ss))
       'H3  SR,SRD  <- slack',   @(v) setix(v, [iSR iSRD], PUs([iSR iSRD]))
       'H4  PD      <- slack',   @(v) setix(v, iPD, PUs(iPD))
       'H5  A_def   <- slack',   @(v) setix(v, iAd, PUs(iAd))
       'H6  A_reg+A_def',        @(v) setix(setix(v,iA,PUs(iA)),iAd,PUs(iAd))
       'H7  A_reg+A_def+SRX',    @(v) setix(setix(setix(v,iA,PUs(iA)),iAd,PUs(iAd)),[iSR iSRD],PUs([iSR iSRD]))
       'H8  everything but SL/LSE', @(v) setix(PUs, [iSL iLSE iXD], v([iSL iLSE iXD])) };
fprintf('\n%-28s %8s %8s   %s\n','swap into the ktr window','k63','vs H0','reading');
kH0 = NaN;
R = struct('name',{},'k',{});
for i=1:size(SW,1)
    v = SW{i,2}(PUk);
    [~,k] = seedWindow(p, v, pivK, tK, DUR, 0, Pk);
    if i==1; kH0 = k; end
    fprintf('%-28s %8.1f %8.1f\n', SW{i,1}, k, k-kH0);
    R(i).name = SW{i,1}; R(i).k = k;
end
fprintf('\ntarget: bring the ktr window down to the donor''s %.1f /s\n', k_S);
save(fullfile(resdir,'swaps.mat'),'R','k_S','k_K','kH0','PUs','PUk');

% ---------------- helpers ----------------
function v = setix(v, ix, val); v(ix) = val; end
function [PU0, piv, tt] = deal_(o, t)
    i = find(o.t(:) >= t, 1);
    PU0 = o.PU(i,:); piv = o.LXBPivot(i); tt = o.t(i);
end
