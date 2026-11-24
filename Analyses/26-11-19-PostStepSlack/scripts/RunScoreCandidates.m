%% RunScoreCandidates  What the series-element candidates of #77 do to the 8 mM objective.
%
% Scores robust1_seed (must give 0.424047) and the candidates with a toe in the
% series element (exponent n, same force and same stiffness at the isometric
% force of the 1.10 Lo hold) on the pooled 8 mM battery, as
% Analyses/LowATP_Refit9K2Sweep/RunBaseline8mM.m does. Parameter change only.
% Runs in the main checkout, where the data are; needs a pool of 5.
%   CAND = [3 4] (default), ROOT = checkout with data/

if ~exist('ROOT', 'var'), ROOT = 'C:\home\git\ATP-depletion-and-heart-failure'; end
if ~exist('CAND', 'var'), CAND = [3 4]; end
here = fileparts(mfilename('fullpath'));  ana = fileparts(here);
cd(ROOT);  addpath(genpath(ROOT));
wt = strsplit(genpath(fullfile(ROOT, '.claude')), pathsep);  wt = wt(~cellfun(@isempty, wt));
if ~isempty(wt), rmpath(wt{:}); end
fprintf('model: %s\n', which('dPUdT_CombinedTransitions'));
refreshPool(5);
Fss = 71.03;  L0 = 0.03278;                 % kPa and um at the 1.10 Lo hold, 8 mM (RunStepScan)
p0 = getParams(loadParams('robust1_seed'), [], true, false);
n0 = p0.ekSE;
S = struct('name', {'incumbent'}, 'ekSE', {n0}, 'kSE', {p0.kSE});
for n = CAND, S(end+1) = struct('name', sprintf('n%g', n), 'ekSE', n, 'kSE', Fss^(1/n)/(L0*n/n0)); end %#ok<SAGROW>
outF = fullfile(ana, 'results', 'score.mat');
for i = 1:numel(S)
    t0 = tic;
    p = p0;  p.ekSE = S(i).ekSE;  p.kSE = S(i).kSE;
    p.RunForceVelocity = true;  p.RunSlack = true;  p.RunStairs = false;  p.RunKtr = true;
    p.RunForceVelocityTime = false;  p.RunSlackPassive = true;  p.EvalFeatures = true;
    p.BreakOnODEUnstable = false;  p.PlotEachSeparately = 0;  p.PlotFeatureFitting = 0;
    p.RunSlackSegments = 'AllPar';  p.MaxRunTime = 300;  p.OptimizeOn = 'Feats';
    p.ghostSave = '';  p.ghostLoad = '';  p.SaveBest = 0;  p.mods = {};  p.g = ones(1, 0);
    try
        [cost, ~, fc, fn] = evaluateBakersExp(1, p, true, false);
        S(i).cost = cost;  S(i).featCost = fc(:)';  S(i).featNames = fn(:)';
    catch ME
        S(i).cost = NaN;  S(i).err = ME.message;
    end
    S(i).sec = toc(t0);
    fprintf('%-10s ekSE %.4g kSE %.5g : Et = %.6f (%.0f s)\n', S(i).name, S(i).ekSE, S(i).kSE, S(i).cost, S(i).sec);
    save(outF, 'S');
end
% the ten features that moved most, per candidate
for i = 2:numel(S)
    if isnan(S(i).cost), continue; end
    d = S(i).featCost - S(1).featCost;  [~, k] = sort(abs(d), 'descend');
    fprintf('%s: largest feature-cost changes:', S(i).name);
    for j = k(1:min(10, numel(k))), fprintf(' %s %+.3f;', S(1).featNames{j}, d(j)); end
    fprintf('\n');
end
