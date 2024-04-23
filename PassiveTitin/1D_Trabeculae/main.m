%%  Define some parameters ...
    clear
%%  Import matlaws ...
    import matlaws.*
    import frac.*

    % Material Parameters for myocardium
    pars    = [2.46, 10.02, 0.5];
    % Caputo Fractional Derivative Parameters
    delta  = 0.023;
    alpha  = 0.187;
    % Prony Approximation terms
    np     = 9;
    Tf     = 100.0;
    % bundle fractional parameters into struct
    frac_pars = frac_parameters(alpha, delta, Tf, np);

rd = .1;
    % %%   Run simulation
figure(3232);clf;
nexttile;
        nt   = 10^(3);
        dt   = Tf/nt;

        time = linspace(0, Tf, nt + 1);
        displacement = 0.2*sin(2* pi *time);
        displacement = 0.2/rd*(time-1+rd).*(time < 1 & time > 1-rd) + 0.2*(time>=1);

plot(time, displacement);hold on;
nexttile;
% Generate Arguments
t_cpu = [];
orders = [2,2.5, 3, 3.5, 4]

for order = orders % Number of time steps
        % paramaters of the simulation
        nt   = 10^(order);
        dt   = Tf/nt;
        time = linspace(0, Tf, nt + 1);
        displacement = 0.2*sin(2* pi *time);
        displacement = 0.2/rd*(time-1+rd).*(time < 1 & time > 1-rd) + 0.2*(time>=1);

        % run the simulation
        tic;
        force = diffeq_sim(@trabeculae3D, pars, displacement, dt, 3, frac_pars);
        t2 = toc;
        t_cpu = [t_cpu t2];
        disp(['  - finished running in ' num2str(t_cpu(end))])
        loglog(time, force);hold on;
        % %   exporting
        % if ~exist(fullfile('output'), 'dir')
        %    mkdir(fullfile('output'))
        % end
        % dlmwrite(fullfile('output', sprintf('force_%d.txt', order)), force, 'delimiter', '\t', 'precision', '%18.12f')
end
nexttile;
plot(orders, t_cpu, 's-')



