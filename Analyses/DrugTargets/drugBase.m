function [base, th, tb, Fd, D] = drugBase(core)
%DRUGBASE  Base human-EHT twitch model for the drug-target tests.
%   [BASE, TH, TB, FD, D] = DRUGBASE(CORE) returns the params struct BASE and the Ca/thin-filament
%   options TH of a HumanTwitch tier optimum (CORE = 't1_isoT' realistic core [default], 's3_free'
%   possible core, 't2_allc' corridor corner), the data grid TB (s from the stimulus, 0..0.95), the
%   Plate 6 control twitch FD (uN) on it and the full data struct D. The Ca transient is the fixed
%   hiPSC-EHT set (de Lange 2021 anchor) of Analyses/HumanTwitch/literature_ca_transient.md.
%   Both ATP gates (UseAtpDetachR2 on R2, UseAtpKmsrd on the SRX-ADP return) are switched ON with the
%   snapshot's K values so that MgATP can be lowered; at MgATP_ref = 8 mM both gates are exactly 1 and
%   the fitted twitch is unchanged.
if nargin < 1, core = 't1_isoT'; end
ht = fullfile(fileparts(mfilename('fullpath')), '..', 'HumanTwitch');
addpath(ht);
cfg = tierConfigs(core);
R = load(fullfile(ht, 'results', [core '_fit.mat']));  R = R.R;
[base, th] = buildModel(cfg, R.x);
base.UseAtpDetachR2 = 1;  base.UseAtpKmsrd = 1;          % inert at 8 mM (normalised to MgATP_ref)
th.sigma = R.sigma;
L = load(fullfile(ht, '..', 'CalciumTwitch', 'results', 'eht_controls.mat'));  D = L.D;
sel = D.tb >= 0;  tb = D.tb(sel);  Fd = D.meanF(sel);
end
