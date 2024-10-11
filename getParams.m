function params = getParams(params, g, updateInit, updateModifiers)
% Updated from
% https://github.com/beards-lab/Cardiac-Crossbridge-Explicit-Space-Discretization
% Beard, Daniel A., et al. "Reduced cardiac muscle power with low ATP simulating heart failure." Biophysical Journal 121.17 (2022): 3213-3223.
% https://www.sciencedirect.com/science/article/pii/S0006349522006026
if nargin == 0 || isempty(params)
    params = struct();
end

if nargin < 2 
    g = ones(1, 30); % better longer than sorry    
end

if nargin < 3
    updateInit = true;
end

if nargin < 4
    updateModifiers = false;
end

%% Build default params0
    ML = 2.0; % reference muscle length

    SL0 = 2.0*1.0;
    
    params0 = struct(...
        'N', 30, ... % number of bins. Could be overwritten when UseCalculatedN = 1
        'Slim', 0.2, ... % left and right distance. Overridden by Slim_l/r when UseCalculatedN = 1
        'LXBpivot', 2.1, ... % reference starting point for the probability distribution dicretization (um)
        'dS', 50e-3, ... % default distance between strain bins (um)
        'Slim_l', 1.6, ... % minimal XB length (left bound)
        'Slim_r', 2.2, ... % maximal XB length (right bound)
        'SL0', SL0, ... % initial SL length
        'LSE0', 0, ... % initial length of the spring    
        'SLmax', Inf, ...
        'OutputAtSL', Inf, ...    
        'ML', ML, ... % muscle length (um)(for calculating velocity)
        'Pi', 0,  ...
        'MgATP', 8, ...
        'MgADP', 0, ...
        'Ca', 1000,...
        'Velocity', 0,...
        'UseCalculatedN', true, ... % Use dS and Slim_m / Slim_p
        'UseCa', false,...
        'UseOverlap', false, ...
        'UsePassive', false, ... % parallel passive stiffness
        'PlotProbsOnFig', 0, ... % 0 - false, any number: figure to plot on
        'ValuesInTime', true, ... % export values in time. Outputs just last value otherwise.
        'MatchTimeSegments', true, ... % interpolate for exactly given last time point
        'ReduceSpace', false, ... % use only half- to no- of the discretized space
        'UseSerialStiffness', true, ... % serial stiffness used with dashpot viscosity
        'UseSlack', true, ... % Enable XB slacking
        'UseKtrProtocol', true, ... % reproduce the protocol for acquiring Ktr
        'PlotEachSeparately', true , ... % show each plot on separate figure
        'PlotFullscreen', false, ... % Each plot is in fullscreen instead on common figure
        'UseSLInput', false, ... % Use SL as a driving instead of velocities, provide input in datatable
        'RescaleOutputDt', 0,... % downsamples unnecessary complex output vector. False or value (e.g. 1e-5)
        'UseP31Shift', false, ... % Shifts the s by dr in p3_1
        'F_act_UseP31', false, ... Use kstiff2*p3_1 instead of p3_0*dr
        'UseAtpOnUNR', false, ... Enables ATP effect via g4 from SR to NR
        'UseTORNegShift', false, ... XB TOR uses s - s3 instead of s + s3
        'UseMutualPairingAttachment', false, ... % Pu to P1 state transient relative to Pu^2
        'UseSpaceDiscretization', false, ...
        'UseSpaceInterpolation', true, ...
        'UseKstiff3', false, ... % uses additional parameter kstiff3 for overstroke stifness (=kstiff2 otherwise)
        'EvalAtp', [1],... % which ATPs should be evaluated from the range [8 4 2] mM - so [1 3] evals 8mM and 2mM and [1] evals 8mM only 
        'SaveBest', true, ... % save g on each iter, if better than previous
        'ghostSave', '', 'ghostLoad', '', ...
        'UseTitinModel', false, ...
        'RunForceVelocity', true, ...
        'RunKtr', true, ...
        'RunStairs', true, ...
        'RunSlack', true, ...
        'WindowsOverflowStepCount', 1, ... % number of dS extensions of the array in case we hit the boundary with the moving window
        'UseSuperRelaxed', false, ...
        'UseSpaceExtension', false,...
        'MaxRunTime', Inf, ... % maximal model runtime in seconds before is killed. Used to unstuck the optimization
        'NumberOfStates', 2, ... % number of strain-dependent states
        'UseTitinInterpolation', true, ... % interpolating titin pre-simulated force at each model eval
        'MaxSlackNegativeForce', 0, ... 
        'justPlotStateTransitionsFlag', false, ...
        'ShowStatePlots', false, ...
        'ShowResidualPlots', false, ...
        'UseDirectSRXTransition', false, ...
        'UsePassiveForSR', false, ...   
        'UseSuperRelaxedADP', false, ...
        'RunSlackSegments', 'All' ...
        );
 
    params0.mods = {}; % names of the modifiers in the cell array. First is modified by g(1), second g(2) etc    

%% SARCOMERE PARAMETERS
    params0.g = g;
    % transition from NP to P, only when UseCa = true
    % g0 params from cross bridge model identrification
    g0 = [ 1.5*0.3977    2.0478    1.4903    0.3765    0.5219    0.2726    1.25  1.0471    0.2382    0.9342];
    params0.K_coop = 5.7;
    params0.k_on   = g0(1)*100;
    params0.k_off  = g0(2)*1.5*100;
    params0.datatable = [];
    
    
    
    params0.vmax = 10;
    
    % rate constants
    params0.kah = 80; % rate of ATP hydrolysis state change
    params0.kadh = 20; % % rate of ATP de-hydrolysis state change
    params0.ka  = 373.23;
    params0.kd  = 102.94; 
    params0.k1  = 40.116;%
    params0.k_1 = 17.103;%
    params0.k2  = 419.39;
    params0.k_2 = 2.7901; 
    params0.k3  = 44.255;%;

    % transitions between super relaxed state and non relaxed state
    params0.ksr0   = 9.0853; % 
    params0.sigma0 = 33.125;
    params0.kmsr   = 250; % 

    % dissociation constants
    params0.K_Pi = 4.007;
    params0.K_T1 = 0.5; % (mM) ATP binding for detachment
%     K_T2 = 0.05; % (mM) ATP binding to P state


    % moments and force
    params0.dr = 0.01; % Power-stroke Size; Units: um
    params0.kstiff1 = 1393.2; 
    params0.kstiff2 = 13275;
    params0.kstiff3 = params0.kstiff2;

    params0.K_T3 = 4; % (mM)
    params0.K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
    
    % strain-associated parameters
    params0.alpha0 = 0;
    params0.alpha1 = 15.14;
    params0.alpha_1 = 0;
    params0.alpha2 = 10.06;
    params0.alpha3 = 276.67;
    params0.s3     = 0.0099383;
    params0.alphaRip = 0;
    params0.k2rip = 0;

    % offsets
    params0.dr0 = 0;
    params0.dr_1 = 0;
    params0.dr2 = 0;
    params0.dr3 = 0;

    % % linear ripping stuff
    % params0.TK0 = 0;
    % params0.TK = 0;
    % 
    % 

    params0.Amax = 1;
    
    params0.mu = 1e-3; % viscosity
    params0.kSE = 500;

    % passive force coeff
    params0.k_pas = 200; % From Kim Salla et al.

    % other
    params0.k2_L = 0;
            
    %% Fill in the missing input params
    
    params = fillInDefaults(params, params0);
    
    %% MODIFIERS

    if updateModifiers
        %     mods = {'kstiff1', 'kstiff2'};
        for i = 1:length(params.mods)
            params.(params.mods{i}) = params.(params.mods{i})*g(i);
        end

        % ensure we delete the modifiers right after
        params.mods = {};
        params.g = [];
    end
        

    %% SIMULATION PARAMETERS
    if params.UseCalculatedN
        % params.N = ceil((params.Slim_r - params.Slim_l)/params.dS/2);
        % params.LXBpivot = params.SL0;
        % params.LXBpivot = params.Slim_l + (params.Slim_r - params.Slim_l)/2;
        params.LXBpivot = params.Slim_l;
        % params.ss = params.N;
        params.s = ((params.Slim_l:2*params.dS:params.Slim_r) - params.LXBpivot)/2;
        params.ss = length(params.s);
%         params.s_i0 = 0; % not used in this context, searched for in each iteration
    else
    
        % refresh these
        params.dS = params.Slim/(params.N+1);
    %     if params.ReduceSpace && all(params.Velocity == 0)
    %         params.s = [-params.N 0 params.N]*params.dS; % strain space
    %         params.s_i0 = 1; % index of the origin zero strain    
        if params.ReduceSpace && all(params.Velocity > 0)
            params.s = (0:params.N)*params.dS; % strain space
            params.s_i0 = 1; % index of the origin zero strain    
        elseif params.ReduceSpace && all(params.Velocity <= 0)
            params.s = (-params.N:0)*params.dS; % strain space
            params.s_i0 = length(params.s); % index of the origin zero strain    
        else
            params.s = (-params.N:params.N).*params.dS; % strain space
            params.s_i0 = params.N + 1; % index of the origin zero strain    
        end
        params.ss = length(params.s); % strain step (number of Ns in one set)
    end
    
    
    % Build the initialization
    if ~isfield(params, 'PU0') || updateInit
        p0 = zeros(1, params.ss);
        U_SR = 0;
        U_SRD = 0;
        NP = 0;
        PuATP = 0;
        SL0 = params.SL0;
        LSE = params.LSE0;
        % State variable vector concatenates p1, p2, p2, and U_NR
        if params.NumberOfStates == 2
            params.PU0 = [p0, p0, U_SR,NP,SL0,LSE, PuATP, U_SRD];
        elseif params.NumberOfStates == 3
            params.PU0 = [p0, p0, p0,U_SR,NP,SL0,LSE, PuATP, U_SRD];
        end
    end
    
end

%%
function params = fillInDefaults(params, defaults)
    % thanks to Adam Danz from MatlabCentral
    % List fields in both structs
    paramsfn = fieldnames(params); 
    defsfn = fieldnames(defaults); 
    % List fields missing in s
    missingIdx = find(~ismember(defsfn,paramsfn));
    % Assign missing fields to s
    for i = 1:length(missingIdx)
        params.(defsfn{missingIdx(i)}) = defaults.(defsfn{missingIdx(i)}); 
    end
end   