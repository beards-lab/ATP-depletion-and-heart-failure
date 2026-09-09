% RunOscillationID.m -- is the ripple on every force-recovery a property of
% the preparation, or of the rig?
%
% Every recovery window in the 2026 protocol -- post-slack, post-re-stretch,
% and ktr -- rises like a single exponential with a small ripple on top. The
% ripple is what holds rsR2 at 0.75 while the model, which has no ripple, sits
% at 0.999. Before any weight is changed on that account, the ripple has to be
% shown to be an artefact and not biology.
%
% Five signatures separate a rig resonance from a cross-bridge process:
%
%   S1  ONE frequency.  A resonance has a fixed f0 set by the mass and the
%       compliance of the lever/transducer. A kinetic process moves with the
%       rate it belongs to.
%   S2  ATP-INVARIANT.  2 mM slows every measured rate by 2-3x (ktr x0.54,
%       rsK80 x0.38). A resonance does not know about nucleotide.
%   S3  PRESENT WITHOUT CROSS-BRIDGES.  The PNB+Mava run is the same rig and
%       the same manoeuvre with the cross-bridges switched off.
%   S4  AMPLITUDE TRACKS THE EXCITATION, not the force level.
%   S5  IN THE LENGTH CHANNEL.  If the motor position rings at f0 the
%       argument is over.
%
% Output -> results/oscillation_id.mat / .txt / .png

clearvars; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

CFG = { ...
 '03/27','8mM','protocol_03_27_2026_8mM_slack.mat'
 '03/27','2mM','protocol_03_27_2026_2mM_slack.mat'
 '03/27','PNB','protocol_03_27_2026_ActivePNBMava_slack.mat'
 '04/03','8mM','protocol_04_03_2026_8mM_slack.mat'
 '04/03','2mM','protocol_04_03_2026_2mM_slack.mat'
 '04/03','PNB','protocol_04_03_2026_ActivePNBMava_slack.mat'
 '04/10','8mM','protocol_04_10_2026_8mM_slack.mat'
 '04/10','2mM','protocol_04_10_2026_2mM_slack.mat'
 '04/10','PNB','protocol_04_10_2026_ActivePNBMava_slack.mat'};

W = {};
for i = 1:size(CFG,1)
    S  = load(fullfile('data', CFG{i,3}));
    dtb = S.datatable;  vt = S.velocitytable;
    t = dtb(:,1); L = dtb(:,2); F = dtb(:,3);
    keep = [true; diff(t) > 0];  t=t(keep); L=L(keep); F=F(keep);
    fs = 1/median(diff(t));

    iS = find(vt(:,2) < -1);  iR = find(vt(:,2) > 1);
    nC = min(numel(iS), numel(iR));
    for c = 1:nC
        % (A) post-slack rise: release-ramp END -> re-stretch START
        W{end+1} = oscWindow(CFG{i,1}, CFG{i,2}, 'postSlack', c, t, L, F, ...
                             vt(iS(c)+1,1), vt(iR(c),1), vt(iS(c),2), fs); %#ok<SAGROW>
        % (B) post-re-stretch rise: re-stretch ramp END -> next release
        ta = vt(iR(c)+1,1);
        if c < nC; tb = vt(iS(c+1),1); else; tb = min(max(t), ta+0.30); end
        W{end+1} = oscWindow(CFG{i,1}, CFG{i,2}, 'postRestretch', c, t, L, F, ...
                             ta, tb, vt(iR(c),2), fs); %#ok<SAGROW>
    end
end
W = [W{~cellfun(@isempty,W)}];
fprintf('RunOscillationID: %d recovery windows\n', numel(W));

save(fullfile(resDir,'oscillation_id.mat'), 'W', '-v7.3');
reportOscillation(W, resDir);
plotOscillation(W, resDir);
