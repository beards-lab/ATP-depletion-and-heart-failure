function [Force, out] = evaluateModel(fcn, T, params)
% params: model parameter structure (required)s
% default: params = struct('Pi;, 0,'MgADP', 0, 'velocity', -1);
% opts: optional simulation options, otherwise reverting to default
% default: opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
% T must be a vector [start end] TODO remove correction for velocity at this point
% if Velocity in params needs to be vector too
    
    params = getParams(params, params.g,false); % update the init vectors
    PU = params.PU0;
    out = [];

    % vs for VelocitySegment
    for vs = 1:length(T) - 1
        ts = T(vs);
%             et = 0; %elapsed time
        tend = T(vs+1); % ending time of simulation in the current segment

        params.v = params.Velocity(vs);
        params.Vums = params.v*params.ML; % velocity in um/s

        opts = odeset();%odeset('AbsTol',1e-4, 'RelTol', 1e-2);
        
        % test odess
%         tic
%         [t,PU] = ode45(fcn,[ts tend],PU(end,:), opts, params);
%         save('ode45', 'PU');
%         disp(['Ode45: ' num2str(length(t))])
%         toc
%         tic
%         [t,PU] = ode23(fcn,[ts tend],PU(end,:), opts, params);
%         save('ode23', 'PU');
%         disp(['Ode45: ' num2str(length(t))])
%         toc
%         tic        
%         [t,PU] = ode23s(fcn,[ts tend],PU(end,:), opts, params);
%         save('ode23s', 'PU');
%         disp(['Ode45: ' num2str(length(t))])
%         toc
%         tic        
%         [t,PU] = ode23t(fcn,[ts tend],PU(end,:), opts, params);
%         save('ode23t', 'PU');
%         disp(['Ode45: ' num2str(length(t))])
%         toc
%         tic        
        [t,PU] = ode15s(fcn,[ts tend],PU(end,:), opts, params);
%         save('ode15s', 'PU');
%         disp(['Ode45: ' num2str(length(t))])
%         toc
%         tic        
        
        if params.RescaleOutputDt
             %%
%              params.RescaleOutputDt = 1e-5;
             i = find(diff(t) > params.RescaleOutputDt);
             if isempty(i)
                 i = length(t);
             end
%              clf;hold on;plot(t);plot(i, t(i), '*');
%              plot(diff(i), 'd-');
%              clf;
%              plot(out.t(i), out.SL(i), '|')
%             t_n = t(1):params.RescaleOutputDt:t(end);
%             if t_n(end) < t(end)
%                 t_n = [t_n t(end)];
%             end
%             PU = interp1(t, PU, t_n);
%             t = t_n;
        t = t(i);PU = PU(i, :);
        end
        %%
        out = storeOutputs(out, PU, params, t);

        if params.ValuesInTime                
            % reconstruct Force
%                 out.F =  out.LSE*params.kSE;
            is = find(out.t >= ts, 1);
            if max(out.ps0_t(is:end)) > 1e-3
                warning("Boundary broken at vel " + num2str(params.v) + ...
                    "( " + num2str(max(out.ps0_t)) + ")" + ...
                    " Extend the Slim from " + num2str(params.Slim) );
            end
        end      

    end % end the velocity segment
    
    %% Check for the length crossing IN THE LAST SEGMENT ONLY
    if params.OutputAtSL < Inf
        SL = PU(:, 3*params.ss+3);
        ma = max(SL);
        mi = min(SL);
        if ma > params.OutputAtSL && mi < params.OutputAtSL
            % there is a crossing - check the directio frist
            if params.OutputAtSL > SL0
                % growing
                i = find(SL > opts.OutputAtSL, 1);
            else
                % shrinking
                i = find(SL < opts.OutputAtSL, 1);
            end
        else
            error(['Sarcomere did not cross required length of ' ...
                num2str(opts.OutputAtSL) ...
                ', ranged from ' num2str(mi) ' to ' num2str(ma)])
        end
        % find the exact time
        v = (SL(i) - SL(i-1))/(t(i)-t(i-1));                s = opts.OutputAtSL - SL(i-1); 
        tc = s/v + t(i-1); % [t(i-1) tc t(i)]
%         tc = (SL(i) - (opts.OutputAtSL))./((SL(i) - SL(i-1))/(t(i)-t(i-1))) + t(i-1);
        PUi = interp1(t, PU, tc);        
%         SLc = PUi(:, 3*params.ss+3) % control value 
%         Force = PUi(3*params.ss+4)*params.kSE;
        % importance of interp
%         [PU(i - 1, 3*params.ss+4)*params.kSE Force PU(i, 3*params.ss+4)*params.kSE]
    else
%         Force = PU(end, 3*params.ss+4)*params.kSE;
        PUi = PU(end, :)';
    end
    
    % Use the dpudt func to get the actual force (depends on config)
    [~, outputs] = fcn(0, PUi, params);
    Force = outputs(1);
    
    
    if ~params.PlotProbsOnFig
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

function out = storeOutputs(out, PU, params, T)
    if isempty(out)
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
            'ps0_t', []);
    end  

    if ~params.ValuesInTime
        out.PU = PU(end, :);
        T = T(end);
%         return;
    end  

    % extend the curent size
%     The first point of the simulation overlaps with last point of the
%     previous one. Lets cut the frist point then
%%
if length(T) > 1
    fp = 2;% skip the first point to seamless stitch the velocity segments together
else
    fp = 1; % do not skip, we have just one datapoint!
end

    for j = fp:length(T)
%         dt = T(j);
        i = length(out.t) + 1;
        out.PU(i, :) = PU(j, :);
        p1 = PU(j, 1:params.ss); p2 = PU(j, 1*params.ss+1:2*params.ss); p3 = PU(j, 2*params.ss+1:3*params.ss);

        % first moments invalid due to shifting in strain s        
%         out.p1_0(i) = params.dS*sum(p1); out.p1_1(i) = params.dS*sum(params.s.*p1);
%         out.p2_0(i) = params.dS*sum(p2); out.p2_1(i) = params.dS*sum(params.s.*p2);
%         out.p3_0(i) = params.dS*sum(p3); out.p3_1(i) = params.dS*sum((params.s+params.dr).*p3); 

        % calculated post-process
        %     out.F(i) = kstiff2*out.p3_0(i) ...
        %         - max(-kstiff1*(out.p2_1(i) + out.p3_1(i)), 0)^g0(20) + mu*v;
        out.v(i) = params.v;
        out.t(i) = T(j);

        out.NR(i) = PU(j, 3*params.ss+1);
        out.NP(i) = PU(j, 3*params.ss+2);
        out.SL(i) = PU(j, 3*params.ss+3);
        out.LSE(i) = PU(j, 3*params.ss+4);
            
        
        
        % get the XB force from the dpudt directly        
        [~, outputs] = dPUdTCa(0, PU(j, :)', params); 
        out.Force(i) = outputs(1);
        out.FXB(i) = outputs(2);
        out.FXBPassive(i) = outputs(3);
        out.OV(i) = outputs(4);
        out.XB_TOR(i, :) = outputs(5:end);
        out.XB_TORs(i) = params.dS*sum(outputs(5:end));
%         params.kstiff2*out.p3_0(i) - max(-params.kstiff1*(out.p2_1(i) + out.p3_1(i)), 0);
        
        out.LXB = out.SL - out.LSE;
        if i > 1
            out.Vxb(i) = (out.LXB(i) - out.LXB(i-1))/(out.t(i) - out.t(i-1));            
        else
            out.Vxb(i) = 0;
        end

        % check the overflow
        % TODO repair the overflow for both directions
        out.ps0_t(i) = 0;
%         if params.s_i0 == 1 
%             % positive velocities, right side only
%             out.ps0_t(i) = max([p1(end), p2(end), p3(end)]);
%         elseif params.s_i0 == params.ss
%             % negative velocities, left side only
%             out.ps0_t(i) = max([p1(1), p2(1), p3(1)]);
%         else
%             % whole space, mixed velocities, better check both sides
%             out.ps0_t(i) = max([[p1(1), p2(1), p3(1)], p1(end), p2(end), p3(end)]);
%         end
    end
end