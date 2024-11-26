function [Force, out] = evaluateModel(fcn, T, params)
% params: model parameter structure (required)s
% default: params = struct('Pi;, 0,'MgADP', 0, 'velocity', -1);
% opts: optional simulation options, otherwise reverting to default
% default: opts = struct('N', 50, 'Slim', 0.05, 'PlotProbsOnFig', 0, 'ValuesInTime', 0);
% T must be a vector [start end] TODO remove correction for velocity at this point
% if Velocity in params needs to be vector too

params = getParams(params, params.g,false, true); % update the init vectors
PU0 = params.PU0;
out = [];
ss = params.ss;

params = updateRates(params);

% inactive switches
params = rmfield(params, 'UseSlack');
params = rmfield(params, 'UseP31Shift');
params = rmfield(params, 'F_act_UseP31');
params = rmfield(params, 'UseMutualPairingAttachment');
params = rmfield(params, 'UseSpaceDiscretization');
params = rmfield(params, 'UseSpaceInterpolation');
params = rmfield(params, 'UseKstiff3');
params = rmfield(params, 'EvalAtp');
params = rmfield(params, 'vmax');
params = rmfield(params, 'UseSerialStiffness');

% params = rmfield(params, 'UseTitinModel');


if isfield(params, 'ksrm')
    error('Change ksrm to kmsr!!')
end

ticId = tic;

    function  [value, isterminal, direction] = movingWindow(t, y, ~)
        if params.NumberOfStates == 2
            SL = y(2*ss + 3);
            LSE = y(2*ss + 4); % length of the serial stiffness
        elseif params.NumberOfStates == 3
            SL = y(3*ss + 3);
            LSE = y(3*ss + 4); % length of the serial stiffness
        end            
        s = params.s(1) + (-(SL - LSE) + params.LXBpivot)/2;
        s_p0 = 1 + round(-s(1)/params.dS, 6);
        value(1) = floor(s_p0) - 1;
        direction(1) = -1;

        value(2) = ceil(s_p0) - params.ss;
        direction(2) = 1;
        isterminal = [true,true, true];
        % if any(abs(value)< 1e-3)
        %     a = 3;
        % end

        elapsed = toc(ticId); % counting current time
        value(3) = elapsed - params.MaxRunTime; %
        isterminal(3) = 1; % stop
        direction(3) = 1; % find all direction 0
        
        if params.UseForceOnsetShift
            %% calculate using XB force
            % % register force onset
            % s = params.s' + (-(SL - LSE) + params.LXBpivot)/2;
            % s = flipud(-s);
            % ss = params.ss; % space size (length of the s for each of p1-p3)
            % p1_ = y(1:ss);
            % p2_ = y(ss+1:2*ss);
            % 
            % % sum of all probabilities
            % p1_1 = params.dS*sum(s.*p1_);
            % p2_1 = params.dS*sum((s+params.dr).*p2_);
            % 
            % F_active = params.kstiff2*(p2_1) + params.kstiff1*(p1_1);
            % fp = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
            % F_passive = fp(0.4, -0.4, 7.9, (SL-LSE)/2);
            % % Ftotal > eps
            % value(4) = F_passive + F_active - 1e-3;

            %% calculate using LSE force
            LSE = y(2*ss + 4); % length of the serial stiffness
            value(4) = params.kSE*LSE - 1e-3;
            if t > 0
                isterminal(4) = true;
            else
                % do not break before zero
                isterminal(4) = false;
            end
            direction(4) = 1;
        end
    end


if params.UseTitinInterpolation
    % addpath(genpath('PassiveTitin'));
    params.TitinTable = load("PassiveTitin\titin-slack.mat").tit;
end

% vs for VelocitySegment
for vs = 1:length(T) - 1
    ts = T(vs);
    %             et = 0; %elapsed time
    tend = T(vs+1); % ending time of simulation in the current segment

    params.v = params.Velocity(vs);
    params.Vums = params.v*params.ML; % velocity in um/s


    opts = odeset('Events', @movingWindow, 'AbsTol',1e-4, 'RelTol', 1e-2);

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
    t = ts;
    % no event for initialization
    te = [];
    imax = 4;
    while t < tend

        if ~isempty(te)
            % we have an event from previous run
            ts = t(end);
            nds = params.WindowsOverflowStepCount;
            if ie == 1
                % move right
                PU0(1:ss-nds) = PU0(1+nds:ss);
                PU0(ss+1:2*ss-nds) = PU0(ss+1+nds:2*ss);
                % zero the new space
                PU0(ss-nds:ss) = 0; PU0(2*ss-nds:2*ss) = 0;
                if params.NumberOfStates > 2
                    PU0(2*ss +1:3*ss-nds) = PU0(2*ss +1+nds:3*ss);
                    PU0(3*ss-nds:3*ss) = 0;
                end
                params.LXBpivot = params.LXBpivot - nds*params.dS*2;
                fprintf('Hovna took %d steps right\n', length(t))
                imax = imax - 1;
            elseif ie == 2
                % move left
                PU0(1+nds:ss) = PU0(1:ss-nds);
                PU0(ss+1+nds:2*ss) = PU0(ss+1:2*ss-nds);
                
                PU0(1:nds) = 0; PU0(ss+1:ss+nds) = 0; 
                if params.NumberOfStates > 2
                    PU0(2*ss +1+nds:3*ss) = PU0(2*ss +1:3*ss-nds);
                    PU0(2*ss +1:2*ss+nds) = 0;
                end

                % dS is in half-sarcomere space, converting to sarcomere space by 2
                params.LXBpivot = params.LXBpivot + nds*params.dS*2;
                fprintf('Hovna took %d steps left\n', length(t))
                imax = imax - 1;
            elseif ie == 3
                try
                    error('ODEslower took longer than MaxRunTime %d s', params.MaxRunTime);
                catch e
                    handleAndRethrowCostException(e, T(end) - t(end))
                end                
            elseif ie == 4
                % time at which the data start to generate force
                t_f0 = [1.1616    1.4637    1.8162    2.2691    2.7724];
                % t0 = interpts = t(end)
                % i_tf0 = interp1(t_f0, ts, 'nearest');
                ts = interp1(t_f0, t_f0, t(end), 'nearest', 'extrap');
                sprintf('Resetting %f to %f', t(end), ts)
                % ts = t_f0(i_tf0);
                % ts = i_tf0;
                % ts = t_f0
            
                %% unfinished business
                % set states
                % We start from state dsitribution before the slack, this
                % should be in PU0
                % ss = params.ss; 
                % dS = params.dS;
                % 
                % 
                % % p1
                % p10 = sum(PU0(1:ss))*dS;
                % PU0(1:ss) = zeros(1, ss);
                % % p2
                % p20 = sum(PU0(ss+1:2*ss))*dS;
                % PU0(ss+1:2*ss) = zeros(1, ss);
                % 
                % if params.UseSuperRelaxed
                %     SR0 = PU0(2*ss+1);
                %     % all p2 states to SRT states
                %     PU0(2*ss+1) = SR0 + p20;
                % end
                % 
                % % all p1 states back to UD
                % PD0 = PU0(2*ss + 5);
                % PU0(2*ss + 5) = PD0 + p20;
                % 
                % if params.UseSuperRelaxedADP
                %     SRD0 = PU0(2*ss+6);
                % end
                % 
                % % test: braking force is in the p1 only
                % y = @(k_pas, x0, gamma, x) k_pas.*(x-x0).^gamma - 4*0 - x0*0 + 0.5e9.*(x-0.95).^13;    
                % F_passive = y(0.4, -0.4, 7.9, (SL-LSE)/2)
                % 
                % p1_1 = F_passive / params.kstiff1; 
                % 
                % % ok... how is this?
                % s = params.s' + (-(SL - LSE) + params.LXBpivot)/2;
                % 
                % s_i0 = floor(s_p0);
                % s_i1 = ceil(s_p0);
                % s_i0k = s_i1 - s_p0;
                % 
                % p1(s_i0) = p1_1/dS;
                % p1(s_i1) = p1_1/dS;
            end            
        end
        
        lastwarn('', ''); 
        if isfield(params, 'ResetSRat') && ~isempty(params.ResetSRat) && ~isempty(find(params.ResetSRat(:, 1) == ts, 1))
            SR = params.ResetSRat(find(params.ResetSRat(:, 1) == ts, 1), 2);
            PU0(2*ss+1) = SR;
        end
        [t,PU, te, ye, ie] = ode15s(fcn,[ts tend],PU0, opts, params);

        % te contains the times when events occurred
        % ye contains the solutions at the times when events occurred
        % ie contains the indices of the triggered events
        out = storeOutputs(fcn,out, PU, params, t);

        if ~isempty(lastwarn) || imax < 0 || (~params.UseSpaceExtension && ~isempty(te))
            warning('ODEslower is not stable')
            if params.BreakOnODEUnstable
                try
                    error('ODEslower is not stable');
                catch e
                    handleAndRethrowCostException(e, T(end) - t(end))
                end
            end
            % break;
        end        

        PU0 = PU(end,:);
        %%
        % if params.UseTitinModel
        %     if ~exist('x0', 'var')
        %         % for the first run, when it does not exist. Then it shuold be reused
        %         x0 = params.SL0/2 - 0.95;
        %         titin.Time = [];titin.Length = [];titin.Force = [];
        %     end
        % 
        %     % identified separately
        %     mod = [468, 3.83e+04, 2.3, 9, 2.33, 8.36e+06, 4.98, 84.9, 1.73e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 1, 0.175, NaN, NaN, 5.04e+04, 0, ];
        %     [Time, L_t, F_t, ~, x0] = evaluateTitinModel(mod, x0, [ts tend], {params.Velocity(vs)}, 4.4, []);
        %     titin.Time = [titin.Time;Time];
        %     titin.Length = [titin.Length;L_t'];
        %     titin.Force = [titin.Force;F_t'];
        % end
        %%
    end
end % end the velocity segment

%% Check for the length crossing IN THE LAST SEGMENT ONLY
if params.OutputAtSL < Inf
    SL = PU(:, params.NumberOfStates*params.ss+3);
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

if params.UseTitinModel && length(titin.Time) > 2
    FL_i = interp1(titin.Time, [titin.Force,titin.Length], out.t, 'linear', 'extrap');
    out.TitinPassive = FL_i(:, 1)';
    out.TitinLength = FL_i(:, 2)';
    out.Force = out.Force - out.FXBPassive + FL_i(:, 1)';
end


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

function out = storeOutputs(fcn, out, PU, params, T)
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
nS = params.NumberOfStates; % number of strain-dependent states
    for j = fp:length(T)
%         dt = T(j);
        i = length(out.t) + 1;
        out.PU(i, :) = PU(j, :);
        % p1 = PU(j, 1:params.ss); p2 = PU(j, 1*params.ss+1:2*params.ss); p3 = PU(j, 2*params.ss+1:3*params.ss);

        % first moments invalid due to shifting in strain s        
%         out.p1_0(i) = params.dS*sum(p1); out.p1_1(i) = params.dS*sum(params.s.*p1);
%         out.p2_0(i) = params.dS*sum(p2); out.p2_1(i) = params.dS*sum(params.s.*p2);
%         out.p3_0(i) = params.dS*sum(p3); out.p3_1(i) = params.dS*sum((params.s+params.dr).*p3); 

        % calculated post-process
        %     out.F(i) = kstiff2*out.p3_0(i) ...
        %         - max(-kstiff1*(out.p2_1(i) + out.p3_1(i)), 0)^g0(20) + mu*v;
        out.v(i) = params.v;
        out.t(i) = T(j);

        out.SR(i) = PU(j, nS*params.ss+1);
        out.NR(i) = 1- PU(j, nS*params.ss+1);
        out.NP(i) = PU(j, nS*params.ss+2);
        out.SL(i) = PU(j, nS*params.ss+3);
        out.LSE(i) = PU(j, nS*params.ss+4);
        out.PuR(i) = PU(j, nS*params.ss+5);
        out.SRD(i) = PU(j, nS*params.ss+6);
        
        % get the XB force from the dpudt directly        
        [~, outputs, rates] = fcn(T(j), PU(j, :)', params); 
        out.Force(i) = outputs(1);
        out.FXB(i) = outputs(2);
        out.FXBPassive(i) = outputs(3);
        out.OV(i) = outputs(4);
        % out.XB_TOR(i, :) = outputs(5:params.ss+4);
        % out.XB_TORs(i) = params.dS*sum(out.XB_TOR(i, :));
        out.LXBPivot(i) = params.LXBpivot;

        % rates = [RTD, RD1, sum([R1D; R12;R21;XB_Ripped])*dS];
        out.RTD(i) = rates(1);
        out.RD1(i) = rates(2);
        out.R1D(i) = rates(3);
        out.R12(i) = rates(4);
        out.R21(i) = rates(5);
        out.R2T(i) = rates(6);
        out.XB_Ripped(i) = rates(7);
        out.RSR2PT(i) = rates(8);
        out.RPT2SR(i) = rates(9);
        out.RSRD2PD(i) = rates(10);
        out.RPD2SRD(i) = rates(11);
        out.RSRD2SR(i) = rates(12);

        % first moments invalid due to shifting in strain s        
        % p1_0, p2_0, p3_0, p2_1, p3_1_stroke
        if nS == 2
            out.p1_0(i) = outputs(5);
            out.p2_0(i) = outputs(6);
            out.p1_1(i) = outputs(7);
            out.p2_1(i) = outputs(8);
            out.PuATP(i) = outputs(9);
        elseif nS == 3
            out.p1_0(i) = outputs(5);
            out.p2_0(i) = outputs(6);
            out.p3_0(i) = outputs(7);
            out.p1_1(i) = outputs(8);
            out.p2_1(i) = outputs(9);
            out.p3_1(i) = outputs(10);
            out.PuATP(i) = outputs(11);
        end
        
        

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
