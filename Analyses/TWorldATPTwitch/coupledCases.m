function K = coupledCases()
%COUPLEDCASES  Coupled Ca -> XB twitch cases: the Ca side (TW_ATP at ATP aCa, pH p) and the
%   XB side (MgATP aXB in both ATP gates) are set independently, so every ATP effect splits into
%   a Ca-side part, an XB-side part and their interaction:
%     full    (aCa = a, aXB = a)   all 4 ATP x 4 pH
%     caOnly  (aCa = a, aXB = 8)   a < 8
%     xbOnly  (aCa = 8, aXB = a)   a < 8
%   plus the alternate XB core s3_free, full, at 8 and 2 mM, pH 7.2 (robustness).
ATP = atpLadder();  pHs = [7.2 7.0 6.8 6.6];
K = struct('label', {}, 'core', {}, 'aCa', {}, 'aXB', {}, 'pH', {}, 'kind', {});
add = @(K, core, aCa, aXB, p, kind) [K, struct('label', sprintf('cp_%s_Ca%g_XB%g_pH%.1f', core, aCa, aXB, p), ...
    'core', core, 'aCa', aCa, 'aXB', aXB, 'pH', p, 'kind', kind)];
for p = pHs
    for a = ATP
        K = add(K, 't1_isoT', a, a, p, 'full');
        if a < 8
            K = add(K, 't1_isoT', a, 8, p, 'caOnly');
            K = add(K, 't1_isoT', 8, a, p, 'xbOnly');
        end
    end
end
K = add(K, 's3_free', 8, 8, 7.2, 'robust');
K = add(K, 's3_free', 2, 2, 7.2, 'robust');
end
