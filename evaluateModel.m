function [Force, out] = evaluateModel(fcn, vel,T,MgATP,Pi,MgADP,g0, opts)

    if ~exist('opts')
        opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
    end
    

%     Force = zeros(size(vel));
    N = opts.N; % space (strain) discretization--number of grid points in half domain
    Slim = opts.Slim; 
    dS = Slim/N;
    s = (-N:1:0)*dS; % strain 

    % Initial variables for Force-velocity experiment
    p1 = zeros(N+1,1);
    p2 = zeros(N+1,1);
    p3 = zeros(N+1,1);
    out = struct();
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
        [t,PU] = ode15s(fcn,T,PU0,[],N,dS,MgATP,Pi,MgADP,g0);
        if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
            % vector output, at each timestep for zero velocity only
            p1 = PU(:,1:1*N+1)';
            p2 = PU(:,1*N+2:2*N+2)';
            p3 = PU(:,2*N+3:3*N+3)';
            s = s';
        else
            PU = PU(end,:);
            p1 = PU(1:1*N+1);
            p2 = PU(1*N+2:2*N+2);
            p3 = PU(2*N+3:3*N+3);
        end
        
        p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
        p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    %     p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
        % TODO explain the difference
        p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);   
        
        Force = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1), 0).^g0(20) + mu*vel;
        
        if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
            out.F = Force; out.t = t;
            
            out.SL = zeros(length(out.t), 1);% TODO later if needed
            out.p1_0 = p1_0;out.p2_0 = p2_0;out.p3_0 = p3_0;
            out.p1_1 = p1_0;out.p2_1 = p2_0;out.p3_1 = p3_1;
        end

    else
        dt = dS/abs(vel);
        tend = T(end)/abs(vel); % ending time of simulation
        Nstep = round(tend/dt);% = Tspan(end)/dS
        
        if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
            zer = zeros(Nstep,1);
            out.F = zer;
            out.t = zer;
            out.SL = zer;
            % dsitribution in time
            out.p1_0 = zer;
            out.p2_0 = zer;
            out.p3_0 = zer;
            out.p1_1 = zer;
            out.p2_1 = zer;
            out.p3_1 = zer;            
%             % peaks in time: to plot against the means
%             out.p1p_t = zer;
%             out.p2p_t = zer;
%             out.p3p_t = zer;
            % minimal value from the right - meaning we overflow our span
            out.ps0_t = zer;
        end

        
        % simulate kinetics for 1/2 timestep
        [~,PU] = ode15s(fcn,[0 dt/2],PU0,[],N,dS,MgATP,Pi,MgADP,g0);
        PU = PU(end,:); 
          
        for i = 1:(Nstep-1)         
            
          % advection (sliding step)
          PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
          PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
          PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
          
          % simulate kinetics for full step
          [~,PU] = ode15s(fcn,[0 dt],PU,[],N,dS,MgATP,Pi,MgADP,g0);
          PU = PU(end,:);           
          
          if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
            p1 = PU(1:1*N+1); p2 = PU(1*N+2:2*N+2); p3 = PU(2*N+3:3*N+3);
            out.p1_0(i) = dS*sum(p1); out.p1_1(i) = dS*sum(s.*p1);
            out.p2_0(i) = dS*sum(p2); out.p2_1(i) = dS*sum(s.*p2);
            out.p3_0(i) = dS*sum(p3); out.p3_1(i) = dS*sum((s+dr).*p3); 
                   
            out.F(i) = kstiff2*out.p3_0(i) ...
                - max(-kstiff1*(out.p2_1(i) + out.p3_1(i)), 0)^g0(20) + mu*vel;
            if i > 1
                out.t(i) = out.t(i-1) + dt;
                out.SL(i) = out.SL(i-1) +vel*dt;
            end
            % check the overflow
            out.ps0_t(i) = max([p1(1), p2(1), p3(1)]);
          end
          
        end
        % final advection (sliding step)
        PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
        PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
        PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
        % final 1/2 timestep for kinetics
        [~,PU] = ode15s(fcn,[0 dt/2],PU,[],N,dS,MgATP,Pi,MgADP,g0);
        
        PU = PU(end,:);
        p1 = PU(1:1*N+1);
        p2 = PU(1*N+2:2*N+2);
        p3 = PU(2*N+3:3*N+3);

        p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
        p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    %     p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
        % TODO explain the difference
        p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);       

        Force = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1), 0).^g0(20) + mu*vel;

        if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
            out.F(end) = Force; out.t(end) = tend;out.SL(end) = out.SL(i-1) + vel*dt;
            out.p1_0(end) = p1_0;out.p2_0(end) = p2_0;out.p3_0(end) = p3_0;
            out.p1_1(end) = p1_1;out.p2_1(end) = p2_1;out.p3_1(end) = p3_1;
            
            if max(out.ps0_t) > 1e-3
                warning("Boundary broken at vel " + num2str(vel) + ...
                    " Extend the Slim from " + num2str(Slim) );
            end
        end
    
    end % end the velocity dependent condition


    
    if ~opts.PlotProbsOnFig
        return
    end

%%
    figure(opts.PlotProbsOnFig);hold on;


%         subplot(122);
%         plot(t, F);
%         xlabel('t');
%         ylabel('Force');
% 
%         subplot(121);

    plot(s,p1,s,p2,s,p3,'x-', 'linewidth',1.5);
    ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
    xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
    set(gca,'fontsize',14);
    set(gca,'xlim',[-Slim 0]);
    legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');
        
        
        
        
        
        
