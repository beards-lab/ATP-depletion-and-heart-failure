function C = standaloneConditions(setName)
%STANDALONECONDITIONS  Standalone T-World runs (1 Hz, endo, 37 C, normoxic) of each set.
%   'twatp'  TW_ATP model_TWorld_shock_ATP: grid ATP 8/6/4/2 x pH 7.2/7.0/6.8/6.6 (GKATP 0),
%            ATP 8/6/4/2 at pH 7.2 with GKATP 0.05 (Michailova), its own reference
%            (7.9875 / 0.0125 mM) and TW_ATP's plain model_TWorld.
%   'lowatp' TWorld@lowAtp: ATP 8/6/4/2 with I_KATP on (default) and off.
[ATP, ADP] = atpLadder();  pHs = [7.2 7.0 6.8 6.6];  B = 150;
C = struct('label', {}, 'flavor', {}, 'model', {}, 'ATP', {}, 'ADP', {}, 'pH', {}, 'gK', {}, 'beats', {});
add = @(C, lab, fl, mdl, a, d, p, g) [C, struct('label', lab, 'flavor', fl, 'model', mdl, 'ATP', a, 'ADP', d, 'pH', p, 'gK', g, 'beats', B)];
switch setName
    case 'twatp'
        for p = pHs, for i = 1:4
            C = add(C, sprintf('twatp_A%g_pH%.1f', ATP(i), p), 'twatp', 'shock_ATP', ATP(i), ADP(i), p, 0); end, end
        for i = 1:4
            C = add(C, sprintf('twatp_A%g_pH7.2_gK', ATP(i)), 'twatp', 'shock_ATP', ATP(i), ADP(i), 7.2, 0.05); end
        C = add(C, 'twatp_ref', 'twatp', 'shock_ATP', 7.9875, 0.0125, 7.2, 0);
        C = add(C, 'twatp_plain', 'twatp', 'plain', NaN, NaN, NaN, NaN);
    case 'lowatp'
        for g = [1 0], for i = 1:4
            C = add(C, sprintf('lowatp_A%g_K%d', ATP(i), g), 'lowatp', 'lowAtp', ATP(i), ADP(i), 7.2, g); end, end
    case 'adp'                    % ADP study (adpCases): TW_ATP at CK-equilibrium ADP, GKATP 0
        [~, S] = adpCases();
        for i = 1:size(S, 1)
            C = add(C, extractAfter(saLabel(S(i, 1), S(i, 2), S(i, 3)), 'sa_'), 'twatp', 'shock_ATP', S(i, 1), S(i, 2), S(i, 3), 0); end
    otherwise, error('unknown set %s', setName);
end
end
