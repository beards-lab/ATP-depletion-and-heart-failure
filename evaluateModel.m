function [Force, out] = evaluateModel(fcn, T, params, g0, opts)
% params: model parameter structure (required)s
% default: params = struct('Pi;, 0,'MgADP', 0, 'velocity', -1);
% opts: optional simulation options, otherwise reverting to default
% default: opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
% T must be a vector [start end] TODO remove correction for velocity at this point
% if Velocity in params needs to be vector too

    if ~exist('opts')
        opts = struct()
    end
    defs = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, ...
        'ValuesInTime', 0, 'MatchTimeSegments', 0, ...
        'SL0', 2.2, ... % initial SL length
        'ML', 2.2 ... % muscle length (um)(for calculating velocity)
        );
    
    opts = fillInDefaults(opts, defs);
    
    vel = params.Velocity;
    
%     Force = zeros(size(vel));
    N = opts.N; % space (strain) discretization--number of grid points in half domain
%     Slim = opts.Slim; 
    opts.dS = opts.Slim/opts.N;
    opts.s = (-opts.N:1:0)*opts.dS; % strain 
    opts.T = T;

    % Initial variables for Force-velocity experiment
    p1 = zeros(opts.N+1,1);
    p2 = zeros(opts.N+1,1);
    p3 = zeros(opts.N+1,1);
%     out = struct();
    U_NR = 1;
    NP = 0;
    SL0 = opts.SL0;
    % State variable vector concatenates p1, p2, p2, and U_NR
    PU = [p1; p2; p3; U_NR;NP;SL0];
%     PU = PU0;

    % moments and force
    dr = g0(12)*0.01; % Power-stroke Size; Units: um
    kstiff1 = g0(13)*2500; 
    kstiff2 = g0(14)*200;
    mu = g0(19)*0.5; % viscosity
    
%     if all(vel == 0)
%         % Zero velocity:
%         [t,PU] = ode15s(fcn,T,PU,[],N,dS,params,g0);
%         if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
%             % vector output, at each timestep for zero velocity only
%             p1 = PU(:,1:1*N+1)';
%             p2 = PU(:,1*N+2:2*N+2)';
%             p3 = PU(:,2*N+3:3*N+3)';
%             s = s';
%         else
%             PU = PU(end,:);
%             p1 = PU(1:1*N+1);
%             p2 = PU(1*N+2:2*N+2);
%             p3 = PU(2*N+3:3*N+3);
%         end
%         
%         p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
%         p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
%     %     p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
%         % the p3 is the post ratcheted state, since the heads are moved by
%         % the dr lengths
%         p3_0 = dS*sum(p3); p3_1 = dS*sum((s+dr).*p3);   
%         
%         Force = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1), 0).^g0(20) + mu*vel;
%         
%         if isfield(opts, 'ValuesInTime') && opts.ValuesInTime
%             out.F = Force; out.t = t;
%             
%             out.SL = zeros(length(out.t), 1);% TODO later if needed
%             out.p1_0 = p1_0;out.p2_0 = p2_0;out.p3_0 = p3_0;
%             out.p1_1 = p1_1;out.p2_1 = p2_1;out.p3_1 = p3_1;
%             out.NR = PU(:, 3*N+4);
%             out.NP = PU(:, 3*N+5);
%             out.SL = PU(:, 3*N+6);
%             
%             % calculate XB turnover rate
%             K_T1 = g0(11)*1; % (mM) ATP binding for detachment
%             K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
%             g2 = (params.MgATP/K_T1)/(1 + params.MgADP/K_D + params.MgATP/K_T1);
%             k3  = g2*g0(10)*25;%;
%             alpha3 = g0(9)*5000;
%             alpha3o = dS*sum((exp(alpha3*(s + dr) .^ 2) .*  p3));
%             out.XB_TOR = k3*alpha3o; % XB turnover rate
%         end
% 
%     else
        
        
        if opts.ValuesInTime
            out = struct('F', [], ...
                't', [] , ...
                'SL', [], ...
                'p1_0', [], ...
                'p2_0', [], ...
                'p3_0', [], ...
                'p1_1', [], ...
                'p2_1', [], ...
                'p3_1', [], ...
                'v', [],... % velocity in ML/s
                'NR', [], ...
                'NP', [], ...
                'ps0_t', [], ...
                'dr', dr);
%             % just an estim TODO
%             Nstep = 100;
%             zer = zeros(Nstep,1);
%             out.F = zer;
%             out.t = zer;
%             out.SL = zer;
%             % dsitribution in time
%             out.p1_0 = zer;
%             out.p2_0 = zer;
%             out.p3_0 = zer;
%             out.p1_1 = zer;
%             out.p2_1 = zer;
%             out.p3_1 = zer;
%             out.NR = zer;
%             out.NP = zer;
% %             % peaks in time: to plot against the means
% %             out.p1p_t = zer;
% %             out.p2p_t = zer;
% %             out.p3p_t = zer;
%             % minimal value from the right - meaning we overflow our span
%             out.ps0_t = zer;
        end
        
        % vs for VelocitySegment
        for vs = 1:length(vel)
            ts = T(vs);
%             et = 0; %elapsed time
            tend = T(vs+1); % ending time of simulation in the current segment
            
            % TODO account for SL = 1 correction
            if abs(vel(vs)) == 0
                v = 0;
                dt = tend - ts;
                Nstep = 0;
            else
                v = vel(vs);
                dt = opts.dS/abs(v);
                Nstep = round((tend-ts)/dt);% = Tspan(end)/dS
            end
        
            % simulate kinetics for 1/2 timestep
            [t,PU] = ode15s(fcn,[ts ts+dt/2],PU,[],opts.N,opts.dS,params,g0, v*opts.ML);
            out = storeOutputs(out, PU, opts, t, v);
%             et = dt/2;
            PU = PU(end,:); 

            for j = 1:(Nstep-1)         

              % advection (sliding step)
              PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
              PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
              PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;

              % simulate kinetics for full step
              [t,PU] = ode15s(fcn,[t(end) t(end) + dt],PU,[],N,opts.dS,params,g0, v*opts.ML);
              PU = PU(end,:);
%               et = et + dt;

              if opts.ValuesInTime
                out = storeOutputs(out, PU, opts, t(end), v);
              end

            end
%             if opts.MatchTimeSegments
%                 PU_1 = PU; % save for later
% %                 t_1 = t;
%             end
            
            if Nstep > 0
                % final advection (sliding step) only for non-zero velocities
                PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
                PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
                PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
            end
            % final 1/2 timestep for kinetics
            [t,PU] = ode15s(fcn,[t(end) t(end)+ dt/2],PU,[],N,opts.dS,params,g0,v*opts.ML);
            
            
            if opts.MatchTimeSegments
                % extrapolate to match the exact end time
                tis = (t(end)-t(1))/length(t);
                % TODO fit previous data points instead of making up new positions
                ti = [t(1):tis:(tend)];
                PU = interp1(t, PU, ti, 'linear', 'extrap');
                t = ti;
%                 t = interp1([t_1(end) t(end)], [PU_1;PU], tend-ts, 'linear', 'extrap');
                tf = tend-ts - t(1); % final timestep length
                if tf > dt || tf < 0 
                    warning('Extrapolattion out of reasonable bounds');
                end
%             else
%                 [t,PU] = ode15s(fcn,[0 dt/2],PU,[],N,opts.dS,params,g0,v*opts.ML);
%                 tf = dt/2;
            end
            
            if opts.ValuesInTime
                out = storeOutputs(out, PU, opts, t, v);
                
                % reconstruct Force
                out.F = kstiff2*out.p3_0 ...
                    - max(-kstiff1*(out.p2_1 + out.p3_1), 0).^g0(20) ...
                    - abs(mu*out.v);
                
                if max(out.ps0_t) > 1e-3
                    warning("Boundary broken at vel " + num2str(v) + ...
                        " Extend the Slim from " + num2str(Slim) );
                end
            end
            

            PU = PU(end, :);
            p1 = PU(1:1*N+1);
            p2 = PU(1*N+2:2*N+2);
            p3 = PU(2*N+3:3*N+3);

%             p1_0 = opts.dS*sum(p1); p1_1 = opts.dS*sum(s.*p1);
%             p2_0 = opts.dS*sum(p2); 
            p2_1 = opts.dS*sum(opts.s.*p2);
        %     p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
            % TODO explain the difference
            p3_0 = opts.dS*sum(p3); p3_1 = opts.dS*sum((opts.s + dr).*p3);       

            Force = kstiff2*p3_0 - max(-kstiff1*(p2_1 + p3_1), 0).^g0(20) - abs(mu*v);

        end % end the velocity segment
%     end % end the velocity dependent condition


    
    if ~opts.PlotProbsOnFig
        return
    end

%%
    figure(opts.PlotProbsOnFig);hold on;

    plot(s,p1,s,p2,s,p3,'x-', 'linewidth',1.5);
    ylabel('Probability density ($\mu$m$^{-1}$)','interpreter','latex','fontsize',16);
    xlabel('strain, $s$ ($\mu$m)','interpreter','latex','fontsize',16);
    set(gca,'fontsize',14);
    set(gca,'xlim',[-Slim 0]);
    legend('$p_1(s)$','$p_2(s)$','$p_3(s)$','interpreter','latex','fontsize',16,'location','northwest');
        
        
end

function out = storeOutputs(out, PU, opts, T, v)
    if ~opts.ValuesInTime
        return;
    end
        
    % extend the curent size
    for j = 1:length(T)
%         dt = T(j);
        i = length(out.t) + 1;
        p1 = PU(j, 1:1*opts.N+1); p2 = PU(j, 1*opts.N+2:2*opts.N+2); p3 = PU(j, 2*opts.N+3:3*opts.N+3);
        out.p1_0(i) = opts.dS*sum(p1); out.p1_1(i) = opts.dS*sum(opts.s.*p1);
        out.p2_0(i) = opts.dS*sum(p2); out.p2_1(i) = opts.dS*sum(opts.s.*p2);
        out.p3_0(i) = opts.dS*sum(p3); out.p3_1(i) = opts.dS*sum((opts.s+out.dr).*p3); 

        % calculated post-process
        %     out.F(i) = kstiff2*out.p3_0(i) ...
        %         - max(-kstiff1*(out.p2_1(i) + out.p3_1(i)), 0)^g0(20) + mu*v;
        out.v(i) = v;
        out.t(i) = T(j);
%         if i > 1
%             out.t(i) = out.t(i-1) + dt;
%         else
%             out.t(1) = opts.T(1) + dt;
%         end
        out.NR(i) = PU(j, 3*opts.N+4);
        out.NP(i) = PU(j, 3*opts.N+5);
        out.SL(i) = PU(j, 3*opts.N+6);

        % check the overflow
        out.ps0_t(i) = max([p1(1), p2(1), p3(1)]);
    end
end

function opts = fillInDefaults(opts, defs)
    % thanks to Adam Danz from MatlabCentral
    % List fields in both structs
    optsfn = fieldnames(opts); 
    defsfn = fieldnames(defs); 
    % List fields missing in s
    missingIdx = find(~ismember(defsfn,optsfn));
    % Assign missing fields to s
    for i = 1:length(missingIdx)
        opts.(defsfn{missingIdx(i)}) = defs.(defsfn{missingIdx(i)}); 
    end
end