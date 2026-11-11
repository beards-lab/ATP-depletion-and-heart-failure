% load gs
g = load("gopt.csv");

% set the function
fcn = @dPUdTCa;
g_names = {"ka", "kd", "k1", "k_1", "k2", "ksr", "sigma_0", "kmsr", "\alpha_3", "k3", "K_{T1}", "s3", "k_{stiff1}", "k_{stiff2}", "K_{T3}", "\alpha_1", "\alpha_2" ,"A_{max0}", "\mu_{v0}", "\mu_h0", "k_pas"};

% run the command
tic
[Etot, E1] = evaluateProblem(fcn, g, true, [1 1 0 0 0 1])
toc
E1