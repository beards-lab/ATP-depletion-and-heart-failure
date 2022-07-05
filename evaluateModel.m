function F_active = evaluateModel(vel,T,MgATP,Pi,MgADP,g0)

    if size(T, 2) == 1
        % T is an endpoint, thus output F_active is just a
        Tspan = [0 T];
        vector_output = false;
    else
        % T is a vector, thus output is  F_active per time points
        Tspan = T;
        vector_output = true;
    end;

    F_active = zeros(size(vel));
    N = 20; % space (strain) discretization--number of grid points in half domain
    Slim = 0.075; 
    dS = Slim/N;
    s = (-N:1:0)*dS; % strain 

    % Initial variables for Force-velocity experiment
    p1 = zeros(N+1,1);
    p2 = zeros(N+1,1);
    p3 = zeros(N+1,1);
    U_NR = 1;
    % State variable vector concatenates p1, p2, p2, and U_NR
    PU0 = [p1; p2; p3; U_NR];

    % moments and force
    dr = 0.01; % Power-stroke Size; Units: um
    kstiff1 = g0(13)*1500; 
    kstiff2 = g0(14)*10000;     


    if vel == 0
      % Zero velocity:
      [t,PU] = ode15s(@dPUdT,Tspan,PU0,[],N,dS,MgATP,Pi,MgADP,g0);
    else
        dt = dS/abs(vel);
        tend = 0.20/abs(vel); % ending time of simulation
        Nstep = round(tend/dt);
        % simulate kinetics for 1/2 timestep
        [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP,Pi,MgADP,g0);
        PU = PU(end,:); 
        for i = 1:(Nstep-1)
          % advection (sliding step)
          PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
          PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
          PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
          % simulate kinetics for full step
          [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0);
          PU = PU(end,:); 
        end
        % final advection (sliding step)
        PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
        PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
        PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
        % final 1/2 timestep for kinetics
        [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP,Pi,MgADP,g0);
    end

        if ~ vector_output 
            PU = PU(end,:);
            p1 = PU(1:1*N+1);
            p2 = PU(1*N+2:2*N+2);
            p3 = PU(2*N+3:3*N+3);
        else
            % vector output, at each timestep
            p1 = PU(:,1:1*N+1)';
            p2 = PU(:,1*N+2:2*N+2)';
            p3 = PU(:,2*N+3:3*N+3)';
            s = s';
        end

        p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
        p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
        p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);       
        F_active = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 );      
