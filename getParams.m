function params = getParams(params, g)

if nargin == 0 || isempty(params)
    params = struct();
end

if nargin < 2 && ~isfield(params, 'g')
    g = ones(1, 30); % better longer than sorry
elseif isfield(params, 'g')
    g = params.g;
end

%% Build default params0
    ML = 2.2; % reference muscle length

    SL0 = 2.2*1.1;
    
    params0 = struct(...
        'N', 30, ...
        'Slim', 0.06, ...    
        'SL0', SL0, ... % initial SL length
        'LSE0', 0, ... % initial length of the spring    
        'SLmax', Inf, ...
        'OutputAtSL', Inf, ...    
        'ML', ML, ... % muscle length (um)(for calculating velocity)
        'Pi', 0,  ...
        'MgATP', 2, ...
        'MgADP', 0, ...
        'Ca', 1000,...
        'Velocity', 0,...
        'UseCa', false,...
        'UseOverlap', true, ...
        'UsePassive', false, ... % parallel passive stiffness
        'PlotProbsOnFig', 0, ... % 0 - false, any number: figure to plot on
        'ValuesInTime', true, ...
        'MatchTimeSegments', true, ...
        'ReduceSpace', false, ...
        'UseSerialStiffness', true, ... % serial stiffness used with dashpot viscosity
        'UseKtrProtocol', true, ... % reproduce the protocol for acquiring Ktr
        'Terminator', false);

    % transition from NP to P, only when UseCa = true
    % g0 params from cross bridge model identrification
    g0 = [ 1.5*0.3977    2.0478    1.4903    0.3765    0.5219    0.2726    1.25  1.0471    0.2382    0.9342];
    params0.K_coop = 5.7;
    params0.k_on   = g0(1)*100;
    params0.k_off  = g0(2)*1.5*100;    
    

    % rate constants
    params0.ka  = g(1)*50 ;
    params0.kd  = g(2)*5; 
    params0.k1  = g(3)*1000;%
    params0.k_1 = g(4)*10;%
    params0.k2  = g(5)*1000;
    params0.k_2 = 10; % not identified
    params0.k3  = g(10)*100;%;

    % transitions between super relaxed state and non relaxed state
    params0.ksr0   = g(6)*10; % 
    params0.sigma0 = g(7)*40;
    params0.kmsr   = g(8)*20; % 
    % kmsr   = g(8)*250*(1-g3); % 
    % dissociation constants
    params0.K_Pi = 15;
    params0.K_T1 = g(11)*1; % (mM) ATP binding for detachment
    % K_T2 = 0.05; % (mM) ATP binding to P state


    % moments and force
    params0.dr = +g(12)*0.01; % Power-stroke Size; Units: um
    params0.kstiff1 = g(13)*2500; 
    params0.kstiff2 = g(14)*20000;

    params0.K_T3 = g(15)*4; % (mM)
    params0.K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
    
    % strain-associated parameters
    params0.alpha1 = g(16)*50;
    params0.alpha2 = g(17)*50;
    params0.alpha3 = g(9)*10000;
    params0.s3     = 0.0025;
    
    params0.Amax = g(18)*1;
    
    params0.mu = g(19)*1; % viscosity
    params0.kSE = g(20)*5000;

    % passive force coeff
    params0.k_pas = 230*g(21); % From Kim Salla et al.
        
    %% Fill in the missing input params
    
    params = fillInDefaults(params, params0);
    
    % refresh these
    params.dS = params.Slim/params.N;
    if params.ReduceSpace && all(params.Velocity == 0)
        params.s = [-params.N 0 params.N]*params.dS; % strain space
        params.s_i0 = 1; % index of the origin zero strain    
    elseif params.ReduceSpace && all(params.Velocity >= 0)
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
    
    
    % Reset the initialization
    p1 = zeros(1, params.ss);
    p2 = zeros(1, params.ss);
    p3 = zeros(1, params.ss);
    U_NR = 0;
    NP = 0;
    SL0 = params.SL0;
    LSE = params.LSE0;
    % State variable vector concatenates p1, p2, p2, and U_NR
    params.PU0 = [p1, p2, p3, U_NR,NP,SL0,LSE];
    
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