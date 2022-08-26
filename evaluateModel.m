function Force = evaluateModel(fcn, vel,T,MgATP,Pi,MgADP,g0, opts)

    if ~exist('opts')
        opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0);
    end
    

    if size(T, 2) == 1
        % T is an endpoint, thus output F_active is just a
        Tspan = [0 T];
        vector_output = false;
    else
        % T is a vector, thus output is  F_active per time points
        Tspan = T;
        vector_output = true;
    end;

    Force = zeros(size(vel));
    N = opts.N; % space (strain) discretization--number of grid points in half domain
    Slim = opts.Slim; 
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
    dr = g0(12)*0.01; % Power-stroke Size; Units: um
    kstiff1 = g0(13)*2500; 
    kstiff2 = g0(14)*200;
    mu = g0(19)*0.5; % viscosity
    


    if vel == 0
      % Zero velocity:
      [t,PU] = ode15s(fcn,Tspan,PU0,[],N,dS,MgATP,Pi,MgADP,g0);
    else
        dt = dS/abs(vel);
        tend = Tspan(end)/abs(vel); % ending time of simulation
        Nstep = round(tend/dt);% = Tspan(end)/dS

        
        % simulate kinetics for 1/2 timestep
        [~,PU] = ode15s(fcn,[0 dt/2],PU0,[],N,dS,MgATP,Pi,MgADP,g0);
        PU = PU(end,:); 
          
        for i = 1:(Nstep-1)
          
%           if opts.ValuesInTime
%             PU = PU(end,:);
%             p1 = PU(1:1*N+1);
%             p2 = PU(1*N+2:2*N+2);
%             p3 = PU(2*N+3:3*N+3);
%             p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
%             p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
%             %p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);       
%             %Force = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 );      
% 
%             p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);       
%             F(i) = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1), 0) + mu*vel;
%             if i > 1
%                 t(i) = t(i-1) + dt;
%                 SL(i) = SL(i-1) +vel*dt;
%             end
%           end
            
          % advection (sliding step)
          PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
          PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
          PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
          % simulate kinetics for full step
          [~,PU] = ode15s(fcn,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0);
          PU = PU(end,:); 
        end
        % final advection (sliding step)
        PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
        PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
        PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
        % final 1/2 timestep for kinetics
        [~,PU] = ode15s(fcn,[0 dt/2],PU,[],N,dS,MgATP,Pi,MgADP,g0);
                
    end

        if ~ vector_output
            PU = PU(end,:);
            p1 = PU(1:1*N+1);
            p2 = PU(1*N+2:2*N+2);
            p3 = PU(2*N+3:3*N+3);
        elseif vel == 0
            % vector output, at each timestep for zero velocity only
            p1 = PU(:,1:1*N+1)';
            p2 = PU(:,1*N+2:2*N+2)';
            p3 = PU(:,2*N+3:3*N+3)';
            s = s';
        else
           error("Non-zero velocity vector output not implemented") 
           
        end

        p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
        p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
        %p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);       
        %Force = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1 );      
        
        p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);       
        Force = kstiff2*p3_0 + kstiff1*( p2_1 + p3_1 ) + mu*vel;
        
        if ~opts.PlotProbsOnFig
            return
        end
        
        figure(opts.PlotProbsOnFig);clf;hold on;
        
%         if opts.ValuesInTime
%             subplot(122);
%             plot(t, F);
%             xlabel('t');
%             ylabel('Force');
%             
%             subplot(121);
%         end
        
        plot(s,p1,s,p2,s,p3,'x-', 'linewidth',1.5);
        ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
        xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14);
        set(gca,'xlim',[-Slim 0]);
        legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');
        
        
        
        
