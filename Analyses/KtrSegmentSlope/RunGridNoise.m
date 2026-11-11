% RunGridNoise.m -- is the objective's numerical roughness the strain-grid
% shift events, and does a window wide enough to never shift remove it?
%
% WHY. refit8 has been stalled for 9 rounds: re-evaluations of its incumbent
% at 1e-6 perturbations span 0.156-0.182, more than every gain it banked.
% The memory (objective-noise-floor) blames discrete decisions, event-driven
% strain-grid shifts first. In evaluateModel the strain window follows the
% cross-bridge pivot; when the pivot comes within A1AttachmentWidth of an edge,
% an event shifts the grid by nds bins and ZEROES the bins that fall off the
% end (PU0(1:nds) = 0). How much probability mass that discards depends on
% the exact event time, which moves with any parameter -- a jump, not a slope.
%
% TEST, on refit8_opt, alongside the running campaign (3 workers):
%   A  the current window   Slim 1.70-2.21 um (255 bins per state)
%   B  twice as wide         Slim 1.45-2.46 um
% For each: the number of shift events in one slack run, the cost and rate
% features at 5 points 1e-6 apart (SD = the numerical noise), and wall time.
%
%   matlab -batch "run('Analyses/KtrSegmentSlope/RunGridNoise.m')"
root='C:\home\git\ATP-depletion-and-heart-failure'; cd(root); addpath(genpath(root));
resdir=fullfile(root,'Analyses','KtrSegmentSlope','results');
cl=parcluster('local'); jsl=fullfile(tempdir,sprintf('mljob_gn_%d',feature('getpid'))); if ~isfolder(jsl); mkdir(jsl); end
cl.JobStorageLocation=jsl; if isempty(gcp('nocreate')); parpool(cl,3); end

p0 = getParams(loadParams('refit8_opt'), [], true, false);
p0.RunForceVelocity=true; p0.RunSlack=true; p0.RunStairs=false; p0.RunKtr=true;
p0.RunForceVelocityTime=false; p0.RunSlackPassive=true; p0.EvalFeatures=true;
p0.BreakOnODEUnstable=false; p0.PlotEachSeparately=0; p0.PlotFeatureFitting=0;
p0.RunSlackSegments='AllPar'; p0.MaxRunTime=300; p0.OptimizeOn='Feats';
p0.ghostSave=''; p0.ghostLoad=''; p0.SaveBest=0; p0.mods={}; p0.g=ones(1,0);

CFG = { 'A current  cap 60',        1.70, 2.21,  60
        'B cap 128 (Slim range)',   1.70, 2.21, 128
        'C wide 1.45-2.46, 253',    1.45, 2.46, 253 };
NREP = 5; FEAT = {'ktr80','ktrfit','rsK80','rsT20'};
OUT = struct();
for c = 1:size(CFG,1)
    pc = p0; pc.Slim_l = CFG{c,2}; pc.Slim_r = CFG{c,3}; pc.MaxStrainArraySize = CFG{c,4};
    pc = getParams(pc, pc.g, true, false);                 % rebuilds s and PU0
    % how many grid shifts in one (serial) slack run
    pe = pc; pe.RunSlackSegments = 'All'; pe.ShowArrayShiftWarnings = 1;
    txt = evalc('[~,~,fmE] = runSlackExperiment(pe);');
    nEv = numel(strfind(txt, 'Hovna took'));
    % numerical noise
    rng(3); C = nan(NREP,1); X = struct(); tw = nan(NREP,1);
    for r = 1:NREP
        pv = pc;
        if r > 1
            for i = 1:numel(pv.leanPool)
                nm = pv.leanPool{i};
                if isempty(strfind(nm,'__')) && isfield(pv,nm) && isscalar(pv.(nm)) && pv.(nm) ~= 0 %#ok<STREMP>
                    pv.(nm) = pv.(nm) * (1 + 1e-6*randn);
                end
            end
            pv = getParams(pv, pv.g, true, false);
        end
        t0 = tic;
        params0 = pv; LoadData; RunBakersExp;                  %#ok<NASGU>
        [fc,~,~] = evalFeatureCost(features_data, features_model, pv.fn);
        tw(r) = toc(t0); C(r) = sum(fc);
        for q = FEAT; X.(q{1})(r,:) = features_model.(q{1})(:)'; end
    end
    OUT(c).name = CFG{c,1}; OUT(c).nEv = nEv; OUT(c).C = C; OUT(c).X = X; OUT(c).tw = tw;
    OUT(c).bins = numel(pc.s);
    fprintf('\n=== %s  (%d bins per state) ===\n', CFG{c,1}, numel(pc.s));
    fprintf('  grid-shift events in one slack run: %d\n', nEv);
    fprintf('  cost: %s  -> mean %.4f  SD %.4f (%.1f %%)\n', mat2str(round(C',4)), mean(C), std(C), 100*std(C)/mean(C));
    fprintf('  wall per evaluation: %.0f s\n', mean(tw));
    for q = FEAT
        M = X.(q{1}); se = std(M,0,1) ./ abs(mean(M,1));
        fprintf('  %-7s mean %s   jitter %s %%\n', q{1}, mat2str(round(mean(M,1),2)), mat2str(round(100*se,2)));
    end
end
save(fullfile(resdir,'grid_noise.mat'),'OUT','CFG');
fprintf('\nDONE\n');
