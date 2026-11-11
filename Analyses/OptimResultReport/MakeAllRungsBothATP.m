% MakeAllRungsBothATP.m -- the whole both-[ATP] figure set, once per rung of the
% ATP hypothesis ladder, each into its own subfolder.
%
% The 8 mM side never moves: every ATP gate is exactly 1 at MgATP_ref = 8, so
% all four rungs share the ladder's frozen core (spreadw_iter_018 = spreadw_opt).
% Only the 2 mM model changes, which is precisely what makes the four sets worth
% putting side by side.
%
%   H1   3 gates   chemistry only            k2 kah kamh
%   H2   7 gates   + SRX nucleotide states   + kmsrd ksr0 ksr2srd ksrd2sr
%   H3  10 gates   + mechanical              + eta_M kSE mu_neg
%   H4  21 gates   any rate                  + kstiff2 kSE_M dr dr2 P_bound_max
%                                              ka kd k1 k_1 tau_reg A0
%
% That is the FREE set. The FITTED set is much smaller, and H1 and H2 landed on
% exactly the same point -- so their figures are identical, which is itself the
% result: the SRX nucleotide gates bought nothing.
%
% Prerequisite: Analyses/AtpEffectHypotheses/RunBothAtpSlideRuns.m in the
% ATP-wt-lowatp worktree, which writes results/bothatp/<rung>/lowatp_runs.mat.
% The data envelope is model-independent and shared, built once by
% RunEnvelopeBothATP.m into results/bothatp/envelope_bothatp.mat.
%
% Run: matlab -batch "run('Analyses/OptimResultReport/MakeAllRungsBothATP.m')"

clearvars -except RUNGS TAG_NEW; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));

if ~exist('RUNGS','var');   RUNGS   = {'H1','H2','H3','H4'}; end
if ~exist('TAG_NEW','var'); TAG_NEW = 'spreadw_opt'; end

for iR = 1:numel(RUNGS)
    RUNG = RUNGS{iR};
    f = fullfile(here,'results','bothatp',RUNG,'lowatp_runs.mat');
    if ~isfile(f)
        fprintf(2,'\n== %s SKIPPED -- no %s\n', RUNG, f); continue;
    end
    fprintf('\n================ %s ================\n', RUNG);
    run(fullfile(here,'MakeSlideFiguresBothATP.m'));
    run(fullfile(here,'MakeFeatureFiguresBothATP.m'));
    % Force-velocity comes from its own sweep (RunFVFineAllRungs.m); skip
    % rather than fail if only the slack/staircase export is present.
    if isfile(fullfile(here,'results','bothatp',RUNG,'fv_fine.mat'))
        run(fullfile(here,'MakeFVFiguresBothATP.m'));
    else
        fprintf(2,'  no fv_fine.mat for %s -- F-V panels skipped\n', RUNG);
    end
end
fprintf('\nALL RUNGS DONE\n');
