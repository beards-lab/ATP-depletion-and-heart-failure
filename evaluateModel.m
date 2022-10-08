function [Force, out] = evaluateModel(fcn, T, params, g0, opts)
% params: model parameter structure (required)s
% default: params = struct('Pi;, 0,'MgADP', 0, 'velocity', -1);
% opts: optional simulation options, otherwise reverting to default
% default: opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
% T must be a vector [start end] TODO remove correction for velocity at this point
% if Velocity in params needs to be vector too

    if ~exist('opts')
        opts = struct();
    end
    defs = struct('N', 20, 'Slim', 0.04, 'PlotProbsOnFig', 0, ...
        'ValuesInTime', 0, 'MatchTimeSegments', 0, ...
        'SL0', 2.2, ... % initial SL length
        'ML', 2.2, ... % muscle length (um)(for calculating velocity)
        'PlotProbsOnStep', false, ...
        'ReduceSpace', false, ...
        'SLmax', Inf, ...
        'OutputAtSL', Inf, ...
        'LSE0', 0 ...
        );
    
    opts = fillInDefaults(opts, defs);
    
    vel = params.Velocity;
    params.PlotProbsOnStep = opts.PlotProbsOnStep;
    
%     Force = zeros(size(vel));
    params.N = opts.N; % space (strain) discretization--number of grid points in half domain
    N = opts.N;
%     Slim = opts.Slim; 
    params.dS = opts.Slim/opts.N;
%     params.Slim = opts.Slim;
    if opts.ReduceSpace && all(params.Velocity == 0)
        params.s = [-opts.N 0 opts.N]*params.dS; % strain space
        params.s_i0 = 1; % index of the origin zero strain    
    elseif opts.ReduceSpace && all(params.Velocity >= 0)
        params.s = (0:opts.N)*params.dS; % strain space
        params.s_i0 = 1; % index of the origin zero strain    
    elseif opts.ReduceSpace && all(params.Velocity <= 0)
        params.s = (-opts.N:0)*params.dS; % strain space
        params.s_i0 = length(params.s); % index of the origin zero strain    
    else
        params.s = (-opts.N:opts.N).*params.dS; % strain space
        params.s_i0 = opts.N + 1; % index of the origin zero strain    
    end

    params.ss = length(params.s); % strain step (number of Ns in one set)

    % Initial variables for Force-velocity experiment
%     p1 = 1*ones(1, params.ss);
%     p2 = 2*ones(1, params.ss);
%     p3 = 3*ones(1, params.ss);   
% 
%     p1([3, params.ss-2]) = 4;
%     p2([3, params.ss-2]) = 4;
%     p3([3, params.ss-2]) = 4;
if isfield(params,'PU0')
    PU = params.PU0;
else
    p1 = zeros(1, params.ss);
    p2 = zeros(1, params.ss);
    p3 = zeros(1, params.ss);
    U_NR = 0;
    NP = 0;
    SL0 = opts.SL0;
    LSE = opts.LSE0;
    % State variable vector concatenates p1, p2, p2, and U_NR
    PU = [p1, p2, p3, U_NR,NP,SL0,LSE];
end

    % moments and force
    dr = g0(12)*1; % Power-stroke Size; Units: um
    params.kstiff1 = g0(13)*2500; 
    params.kstiff2 = g0(14)*20000;
    params.mu = g0(19)*0.1; % viscosity
    params.kSE = g0(20)*5000;
        
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
        else
            out = struct();
        end
        
        % vs for VelocitySegment
        for vs = 1:length(vel)
            ts = T(vs);
%             et = 0; %elapsed time
            tend = T(vs+1); % ending time of simulation in the current segment
            
            % TODO account for SL = 1 correction
            params.v = vel(vs);
            params.Vums = params.v*opts.ML; % velocity in um/s
            params.g0 = g0;
        
            [t,PU] = ode15s(fcn,[ts tend],PU(end,:),[], params,g0);
            out = storeOutputs(out, PU, params, t, opts.ValuesInTime);
            
            if opts.ValuesInTime                
                % reconstruct Force
%                 out.F =  out.LSE*params.kSE;                
                if max(out.ps0_t) > 1e-3
                    warning("Boundary broken at vel " + num2str(params.v) + ...
                        " Extend the Slim from " + num2str(opts.Slim) );
                end
            end      
        
        end % end the velocity segment
    
    %% Check for the length crossing IN THE LAST SEGMENT ONLY
    if opts.OutputAtSL < Inf
        SL = PU(:, 3*params.ss+3);
        ma = max(SL);
        mi = min(SL);
        if ma > opts.OutputAtSL && mi < opts.OutputAtSL
            % there is a crossing - check the directio frist
            if opts.OutputAtSL > SL0
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
        Force = PUi(3*params.ss+4)*params.kSE;
        % importance of interp
%         [PU(i - 1, 3*params.ss+4)*params.kSE Force PU(i, 3*params.ss+4)*params.kSE]
    else
        Force = PU(end, 3*params.ss+4)*params.kSE;
    end
    
    
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

function out = storeOutputs(out, PU, params, T, store)
    if ~store
        out.PU = PU(end, :);
        return;
    end
        
    % extend the curent size
%     The first point of the simulation overlaps with last point of the
%     previous one. Lets cut the frist point then
%%
if length(T) > 1
    fp = 2;% skip the first point
else
    fp = 1; % do not skip, we have just one datapoint!
end

    for j = fp:length(T)
%         dt = T(j);
        i = length(out.t) + 1;
        out.PU(i, :) = PU(j, :);
        p1 = PU(j, 1:params.ss); p2 = PU(j, 1*params.ss+1:2*params.ss); p3 = PU(j, 2*params.ss+1:3*params.ss);
        out.p1_0(i) = params.dS*sum(p1); out.p1_1(i) = params.dS*sum(params.s.*p1);
        out.p2_0(i) = params.dS*sum(p2); out.p2_1(i) = params.dS*sum(params.s.*p2);
        out.p3_0(i) = params.dS*sum(p3); out.p3_1(i) = params.dS*sum((params.s+out.dr).*p3); 

        % calculated post-process
        %     out.F(i) = kstiff2*out.p3_0(i) ...
        %         - max(-kstiff1*(out.p2_1(i) + out.p3_1(i)), 0)^g0(20) + mu*v;
        out.v(i) = params.v;
        out.t(i) = T(j);

        out.NR(i) = PU(j, 3*params.ss+1);
        out.NP(i) = PU(j, 3*params.ss+2);
        out.SL(i) = PU(j, 3*params.ss+3);
        out.LSE(i) = PU(j, 3*params.ss+4);
            
        out.Force(i) = out.LSE(i)*params.kSE;
        
        % get the XB force from the dpudt directly        
        [f out.FXB(i)] = dPUdTCa(0, PU(j, :)', params, params.g0); 
%         params.kstiff2*out.p3_0(i) - max(-params.kstiff1*(out.p2_1(i) + out.p3_1(i)), 0);
        
        out.LXB = out.SL - out.LSE;
        if i > 1
            out.Vxb(i) = (out.LXB(i) - out.LXB(i-1))/(out.t(i) - out.t(i-1));            
        else
            out.Vxb(i) = 0;
        end

        % check the overflow
        % TODO repair the overflow for both directions
        if params.s_i0 == 1 
            % positive velocities, right side only
            out.ps0_t(i) = max([p1(end), p2(end), p3(end)]);
        elseif params.s_i0 == params.ss
            % negative velocities, left side only
            out.ps0_t(i) = max([p1(1), p2(1), p3(1)]);
        else
            % whole space, mixed velocities, better check both sides
            out.ps0_t(i) = max([[p1(1), p2(1), p3(1)], p1(end), p2(end), p3(end)]);
        end
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