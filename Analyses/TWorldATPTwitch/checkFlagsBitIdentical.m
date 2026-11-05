function checkFlagsBitIdentical(oldDir)
% CHECKFLAGSBITIDENTICAL  Bit-identity of the new default-off flags: RHS of dPUdT_CombinedTransitions with the branch code
% vs the HumanTwitch code, (a) cross-bridge only at MgATP 2 mM with the R2 ATP gate on (UseAdpCompR2
% absent/off), (b) with the embedded T-World Ca layer (CaInput 'tworld', coupling 'land').
% OLDDIR holds the HumanTwitch versions of Model/{dPUdT_CombinedTransitions,calciumActivation,tworldCaLayer,getParams}.m,
% e.g. for f in ...; do git show HumanTwitch:Model/$f.m > old/$f.m; done. Result 2026-11-05: 0 and 0 -> PASS.
root = fullfile(fileparts(mfilename('fullpath')), '..', '..');  addpath(genpath(root));
cases = {'xb', 'tworld'};  vn = {'new', 'old'};  pf = {'PASS', 'FAIL'};
for c = 1:2
    d = cell(1, 2);
    for v = 1:2
        if v == 2, addpath(oldDir, '-begin'); end
        clear functions; %#ok<CLFUNC>
        P = loadParams('human_twitch_t1_isoT');  P.UseAtpDetachR2 = 1;  P.UseAtpKmsrd = 1;  P.MgATP = 2;
        if c == 2
            P.CaModel = 'twoState';  P.CaInput = 'tworld';  P.CaPeriod = 1;  P.CaFirstBeat = 0;
            P.CaTworld = struct('coupling', 'land');
        end
        P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
        rng(1);  y = P.PU0(:) .* (1 + 0.05*rand(numel(P.PU0), 1));
        d{v} = dPUdT_CombinedTransitions(0.0123, y, P);
        fprintf('%s %s: %s | %s\n', cases{c}, vn{v}, which('dPUdT_CombinedTransitions'), which('tworldCaLayer'));
        if v == 2, rmpath(oldDir); end
    end
    fprintf('CHECK %-6s max|dRHS| = %.3g (n = %d) -> %s\n', cases{c}, max(abs(d{1} - d{2})), numel(d{1}), pf{1 + ~isequal(d{1}, d{2})});
end
end
