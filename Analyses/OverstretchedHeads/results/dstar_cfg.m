rip = @(K) struct('UseD0State', 1, 'k2rip', 200, 'alphaRip', 100, 'dr3', -0.015, 'k_D0T', 500, 'K_D0T', K);
VAR = { 'as fitted @ 2 mM',                        struct(),   2
        'D* rip >15 nm, K_D0T 0.15 mM',            rip(0.15),  8
        'D* rip >15 nm, K_D0T 0.15 mM @ 2 mM',     rip(0.15),  2
        'D* rip >15 nm, K_D0T 2 mM',               rip(2),     8
        'D* rip >15 nm, K_D0T 2 mM @ 2 mM',        rip(2),     2 };
PARAMSET = 'refit9_opt';
