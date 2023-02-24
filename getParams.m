function params = getParams(params, g, updateInit)
% Updated from
% https://github.com/beards-lab/Cardiac-Crossbridge-Explicit-Space-Discretization
% Beard, Daniel A., et al. "Reduced cardiac muscle power with low ATP simulating heart failure." Biophysical Journal 121.17 (2022): 3213-3223.
% https://wwwhttps://www.sciencedirect.com/science/article/pii/S0006349522006026.sciencedirect.com/science/article/pii/S0006349522006026
if nargin == 0 || isempty(params)
    params = struct();
end

if nargin < 2 && ~isfield(params, 'g')
    g = ones(1, 30); % better longer than sorry
elseif isfield(params, 'g')
    g = params.g;
end

if nargin < 3
    updateInit = true;
end

%% Build default params0
    ML = 2.0; % reference muscle length

    SL0 = 2.0*1.0;
    
    params0 = struct(...
        'N', 30, ...
        'Slim', 0.06, ...    
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
        'UseCa', false,...
        'UseOverlap', false, ...
        'UsePassive', false, ... % parallel passive stiffness
        'PlotProbsOnFig', 0, ... % 0 - false, any number: figure to plot on
        'ValuesInTime', true, ... % export values in time. Outputs just last value otherwise.
        'MatchTimeSegments', true, ... % interpolate for exactly given last time point
        'ReduceSpace', false, ... % use only half- to no- of the discretized space
        'UseSerialStiffness', false, ... % serial stiffness used with dashpot viscosity
        'UseSlack', false, ... % Enable XB slacking
        'UseKtrProtocol', true, ... % reproduce the protocol for acquiring Ktr
        'PlotEachSeparately', false , ... % show each plot on separate figure
        'UseSLInput', false, ... % Use SL as a driving instead of velocities, provide input in datatable
        'RescaleOutputDt', 1e-5,... % downsamples unnecessary complex output vector. False or value (e.g. 1e-5)
        'UseP31Shift', false, ... % Shifts the s by dr in p3_1
        'F_act_UseP31', false, ... Use kstiff2*p3_1 instead of p3_0*dr
        'UseAtpOnUNR', false, ... Enables ATP effect via g4 from SR to NR
        'UseTORNegShift', false, ... XB TOR uses s - s3 instead of s + s3
        'UseMutualPairingAttachment', false, ... % Pu to P1 state transient relative to Pu^2
        'Terminator', false);

    % transition from NP to P, only when UseCa = true
    % g0 params from cross bridge model identrification
    g0 = [ 1.5*0.3977    2.0478    1.4903    0.3765    0.5219    0.2726    1.25  1.0471    0.2382    0.9342];
    params0.K_coop = 5.7;
    params0.k_on   = g0(1)*100;
    params0.k_off  = g0(2)*1.5*100;
    params0.datatable = [];
    
    
    
    params0.vmax = g(22)*10;
    
    % rate constants
    params0.ka  = g(1)*373.23;
    params0.kd  = g(2)*102.94; 
    params0.k1  = g(3)*40.116;%
    params0.k_1 = g(4)*17.103;%
    params0.k2  = g(5)*419.39;
    params0.k_2 = 2.7901; 
    params0.k3  = g(10)*44.255;%;

    % transitions between super relaxed state and non relaxed state
    params0.ksr0   = g(6)*9.0853; % 
    params0.sigma0 = g(7)*33.125;
    params0.kmsr   = g(8)*250; % 

    % dissociation constants
    params0.K_Pi = 4.007;
    params0.K_T1 = g(11)*0.5; % (mM) ATP binding for detachment
%     K_T2 = 0.05; % (mM) ATP binding to P state


    % moments and force
    params0.dr = +g(12)*0.01; % Power-stroke Size; Units: um
    params0.kstiff1 = g(13)*1393.2; 
    params0.kstiff2 = g(14)*13275;

    params0.K_T3 = g(15)*4; % (mM)
    params0.K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
    
    % strain-associated parameters
    params0.alpha1 = g(16)*15.14;
    params0.alpha2 = g(17)*10.06;
    params0.alpha3 = g(9)*276.67;
    params0.s3     = 0.0099383;
    
    params0.Amax = g(18)*1;
    
    params0.mu = g(19)*1e-3; % viscosity
    params0.kSE = g(20)*500;

    % passive force coeff
    params0.k_pas = 200*g(21); % From Kim Salla et al.
        
    %% Fill in the missing input params
    
    params = fillInDefaults(params, params0);
    
    % refresh these
    params.dS = params.Slim/(params.N+1);
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
    
    
    % Build the initialization
    if ~isfield(params, 'PU0') || updateInit
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