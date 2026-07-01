function params = getParams(params, g, updateInit, updateModifiers)
% GETPARAMS  Build or update the model parameter struct.
%
%   PARAMS = GETPARAMS(PARAMS, G, UPDATEINIT, UPDATEMODIFIERS) is the single
%   source of truth for all model parameters. It performs four steps:
%     1. Fill any missing fields from hard-coded defaults.
%     2. Apply multiplicative modifiers G to the fields listed in PARAMS.MODS.
%     3. Resolve linked fields whose value is a string starting with '=' by
%        evaluating it as a MATLAB expression (see resolveParams).
%        Example: params.kamh = '=0.1*kah'  => kamh = 0.1 * kah
%     4. Reconstruct array fields from fieldname__index notation, compute
%        the strain grid PARAMS.S, and initialize the state vector PARAMS.PU0.
%
%   Inputs:
%     PARAMS         - Existing params struct to update, or [] / omitted to
%                      start from defaults.
%     G              - Modifier vector (optional). When non-empty, G(i)
%                      multiplies the field named PARAMS.MODS{i}. May be
%                      empty [] to apply no modifiers.
%     UPDATEINIT     - If true (default), recompute PU0 and strain grid.
%                      Pass false to skip reinitialisation (e.g. when only
%                      rate constants need refreshing).
%     UPDATEMODIFIERS- If true, re-apply G to the current field values.
%                      Default false (applied once at construction time).
%
%   Outputs:
%     PARAMS  - Complete parameter struct ready to pass to evaluateModel.
%
%   Important: always call GETPARAMS again after changing SL0 — it affects
%   the strain grid and the initial state vector PU0.
%
%   Reference: Beard et al. 2022, Biophysical Journal 121.17
%   https://www.sciencedirect.com/science/article/pii/S0006349522006026
%   Adapted from: https://github.com/beards-lab/Cardiac-Crossbridge-Explicit-Space-Discretization
%
%   See also: evaluateModel, resolveParams, dPUdT_CombinedTransitions
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

%% Input validation
% g may be empty (no modifiers applied); only validate size when non-empty
if ~isempty(g) && isfield(params, 'mods') && ~isempty(params.mods) && numel(g) ~= numel(params.mods)
    warning('getParams: g has %d elements but params.mods has %d entries — must match when non-empty', ...
        numel(g), numel(params.mods));
end

%% Build default params0
    ML = 2.0; % reference muscle length

    SL0 = 2.0*1.0;
    
    params0 = struct(...
        'N', 30, ... % number of bins. Could be overwritten when UseCalculatedN = 1
        'Slim', 0.2, ... % left and right distance. Overridden by Slim_l/r when UseCalculatedN = 1
        'LXBpivot', SL0, ... % reference starting point for the probability distribution dicretization (um)
        'dS', 50e-3, ... % default distance between strain bins (um)
        'Slim_l', 1.6, ... % minimal XB length (left bound)
        'Slim_r', 2.2, ... % maximal XB length (right bound)
        'SL0', SL0, ... % initial SL length
        'rsl0', 1, ... % default extra stretch given by serial spring element
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
        'UseOverlapFactor', false, ...
        'UsePassive', false, ... % parallel passive stiffness
        'PlotProbsOnFig', 0, ... % 0 - false, any number: figure to plot on
        'ValuesInTime', true, ... % export values in time. Outputs just last value otherwise.
        'MatchTimeSegments', true, ... % interpolate for exactly given last time point
        'ReduceSpace', false, ... % use only half- to no- of the discretized space
        'UseSerialStiffness', true, ... % serial stiffness used with dashpot viscosity
        'UseViscoelasticSE', false, ... % Kelvin-Voigt dashpot in parallel with kSE spring (one-sided: only active during restretch vel>0)
        'UseDynamicPassive', false, ... % velocity-proportional titin viscoelastic force during restretch
        'UseCatchBond', false, ...      % reduce R1D+R2 detachment during rapid stretch (catch bond / force enhancement)
        'UseCatchBondR1DOnly', false, ... % catch bond on R1D only (not R2) — preserves valley timing; test of p1-burst hypothesis
        'CatchBondStrainMax', inf, ...    % upper strain limit for R2 catch bond — bridges above this detach normally (slip region)
        'UseStretchActivation', false, ... % boost ka during rapid restretch (stretch-activated attachment)
        'UseSlack', true, ... % Enable XB slacking
        'UseKtrProtocol', true, ... % reproduce the protocol for acquiring Ktr
        'PlotEachSeparately', true , ... % show each plot on separate figure
        'PlotFullscreen', false, ... % Each plot is in fullscreen instead on common figure
        'UseSLInput', false, ... % Use SL as a driving instead of velocities, provide input in datatable
        'RescaleOutputDt', 0,... % downsamples unnecessary complex output vector. False or value (e.g. 1e-5)
        'UseP31Shift', false, ... % Shifts the s by dr in p3_1
        'F_act_UseP31', false, ... Use kstiff2*p3_1 instead of p3_0*dr
        'UseAtpOnUNR', false, ... Enables ATP effect via g4 from SR to NR
        'UseAtpK2', false, ...    % Direct [ATP] dependence on k2: k2_eff = k2*[ATP]/([ATP]+K_T1). Low ATP -> slower detachment -> higher F0, lower Vmax.
        'UseAtpKah', false, ...   % Direct [ATP] dependence on kah: kah_eff = kah*[ATP]/([ATP]+K_T_kah). Low ATP -> fewer PD heads ready to attach -> lower Ktr.
        'K_T_kah', 4.0, ...       % Km(ATP) for hydrolysis step [mM]. Larger than K_T1: kah more sensitive -> decouples Ktr from Vmax.
        'UseAdpTrap', false, ...   % ADP-trapping (active path): R2 (ADP release P2->P3) *= g2 (=1 at MgADP=0). Elevated [ADP] traps force-bearing P2 -> force up, ktr down.
        'UsePiReversal', false, ... % Pi power-stroke inhibition (active path): R12 *= 1/(1+Pi/K_Pi). No effect at Pi=0. Elevated [Pi] lowers force (Cooke&Pate).
        'UsePiReverseStroke', false, ... % if true, ALSO enhance reverse stroke R21 *= (1+Pi/K_Pi) — lowers force more but SPEEDS ktr (spurious). Default off = ktr-neutral Pi.
        'UsePiForce', false, ...   % Pi weakens strong (P2) binding force: kstiff2_eff = kstiff2/(1+Pi/K_Pi). ktr-neutral force reduction; spares P3 rigor. No effect at Pi=0.
        'UseAtpDetach', false, ...  % ATP-limited rigor detachment: R3 (P3->PD) *= MgATP/(MgATP+K_T_detach). Low [ATP] -> rigor (P3) accumulates -> force/stiffness up, ktr down.
        'K_T_detach', 4.0, ...     % Km(ATP) [mM] for rigor detachment; sets the 2mM/8mM ratio of effective k3 (with K_T_detach=4: ratio ~0.5).
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
        'RunForceVelocityTime', false, ...
        'RunKtr', true, ...
        'RunStairs', true, ...
        'RunSlack', true, ...
        'WindowsOverflowStepCount', 1, ... % number of dS extensions of the array in case we hit the boundary with the moving window
        'UseSuperRelaxed', false, ...
        'UseSpaceExtension', false,...
        'MaxRunTime', Inf, ... % maximal model runtime in seconds before is killed. Used to unstuck the optimization
        'NumberOfStates', 2, ... % number of strain-dependent states
        'UseTitinInterpolation', false, ... % interpolating titin pre-simulated force at each model eval
        'MaxSlackNegativeForce', 0, ... 
        'justPlotStateTransitionsFlag', false, ...
        'ShowStatePlots', false, ...
        'ShowResidualPlots', false, ...
        'UseDirectSRXTransition', false, ...
        'UsePassiveForSR', false, ...   
        'UseTitinIdentifiedPassive', false, ...
        'UseNegativeKstiff', false, ...
        'UseWDetachment', false, ...
        'useHalfActiveForSR', false, ...
        'UseA2Reattaching', false, ...
        'UseA2Popping', false, ...
        'UseA2AttachmentShift', false, ... % p2 heads at high strain hop to adjacent actin site (A2 hopping)
        'UseA2MechanicalRecocking', false, ... % p2 heads at high strain ratchet back to p1
        'UseForceOnsetShift', false, ...
        'L_thick', 1.67, ... % Length of thick filament, um
        'L_hbare', 0.10,... % Length of bare region of thick filament, um
        'L_thin', 1.20, ... % Length of thin filament, um
        'UseGlobalOccupancySaturation', false, ... % Global (mean-field) occupancy attenuation of attachment: RD1 *= (1 - P_bound/P_bound_max). Only faithful occupancy form here (strain axis is not a site axis). See Analysis/BindingSiteOccupancy/OccupancySaturation_Report.md
        'P_bound_max', 0.30, ...             % Max total bound fraction P_bound=p1_0+p2_0+p3_0 [-]; rat-cardiac max-Ca isometric ~0.12-0.20, relaxed default (uncertain) up to ~0.3-0.4
        'OccupancyForm', 'linear', ...       % 'linear' => 1-P_bound/P_bound_max ; 'langmuir' => 1/(1+P_bound/P_bound_max)
        'UseTargetZoneSaturation', false, ... % [DEPRECATED as occupancy: per-strain-bin, treats strain as a site axis] phenomenological strain-modifier only
        'max_attached_per_bin', 0.01, ...    % [DEPRECATED, dS-dependent] superseded by rho_attach_max
        'rho_attach_max', 16.0, ...          % Max attached-head density [1/um strain] for dS-invariant per-bin form (Stewart/Baker 2021 K_M ~16-26 heads/um)
        'UseVernierVelocity', false, ...     % Velocity-dependent vernier test (alternative to UseTargetZoneSaturation)
        'alpha_vernier', 0.3, ...            % Max fractional increase in ka at high velocity [-]
        'v_ref_vernier', 1.0, ...            % Half-saturation velocity [um/s]
        'UseVelGaussAttachment', false, ...  % 1-Gaussian dip in ka: congestion at isometric, recovery at shortening
        'v_att_sigma', 2.0, ...             % Gaussian width [um/s]; width of isometric congestion dip
        'v_att_amplitude', 0.4, ...         % Depth of dip at v=0 [0-1]; 0.4 = 40% reduction at isometric
        'v_att_center', 0.0, ...            % Center of Gaussian [um/s]; 0 = dip at isometric, >0 = shifted peak
        'UseLatticeSpacing', false, ... % Lattice spacing correction on attachment rate
        'd10_ref', 37e-3, ...          % Reference (1,0) lattice spacing [um] (37 nm, Irving 2017)
        'SL_ref_lattice', 2.0, ...     % SL at which d10_ref was measured [um]
        'R_thick', 8e-3, ...           % Thick filament radius [um] (8 nm)
        'R_thin', 4e-3, ...            % Thin filament radius [um] (4 nm)
        'd_optimal', 25e-3, ...        % Optimal surface-to-surface distance [um] (= d10_ref - R_thick - R_thin)
        'sigma_lattice', 3e-3, ...     % Gaussian width for lattice factor [um] (3 nm, ~thermal fluctuation amplitude)
        'Lsc0', 1.51, .... % minimal sarcomere length for passive
        'gamma', 1, ... % passive exponent
        'FudgeVmax', false, ...
        'e2L', 2, ...
        'e2R', 2,...
        'estiff', 1, ...
        'drawForceOnset', false, ...
        'UseLimitedSRDTransition', false, ...
        'UseSuperRelaxedADP', false, ...        
        'RunSlackSegments', 'All', ...
        'UseUniformTransitionFunc', false, ...
        'EvalFitSlackOnset', false, ...
        'RunForceLengthEstim', false, ...
        'BreakOnODEUnstable' , true, ...
        'DryRun', false, ... % only constant output from dxdt to save computational time and check the infrastructure
        'SimTitle', 'Baseline sim', ...
        'A1AttachmentWidth', 0,... % probability distribution of attachment around 0
        'velocitytableonfile', 'bakers_slack8mM_all.mat', ...
        'modelFcn', '@dPUdT_CombinedTransitions', ...
        'UseNegativeForceRip', false, ...
        'UseStrictDetachmentAt', 0, ...
        'StrainExp', 2,...
        'MaxStrainArraySize', 0, ...
        'UseSafeSRReset', false, ...
        'ShowArrayShiftWarnings', true, ...
        'LegacyStrainFlipping', false, ...
        'UsePieceWiseStrainDep', false, ...
        'PieceWiseStrainDepParams', [1 1], ...
        'PieceWiseStrainDepX', [-0.1 0.1],...
        'UseA1AttachmentKernel', false,...
        'OptimizeFVInit', true, ...
        'EvalFeatures', false, ...
        'SkipParametersInBounds', true, ... % when true, only out-of-bounds params are auto-added to fn (default: true)
        'EvalPeaks', false, ...
        'recalculateDataFeats', false, ...
        'MaxSpaceExtensionCount', Inf, ...
        'UseRegistrationAvailability', false, ... % target-zone registration availability (FV "shoulder"); appends +1 scalar ODE state A_reg that gates attachment (ka_eff)
        'A0', 0.7, ...             % isometric registration availability floor, A_reg in [A0,1]; A0=1 => mechanism neutral. Used when UseRegistrationAvailability
        'v_ref_reg', 0.8, ...      % sliding speed (um/s, same units as Vums; ~0.4 ML/s) at which availability is half-recovered. Used when UseRegistrationAvailability
        'tau_reg', 0.02, ...       % registration relaxation time (s, ~1/ktr) so ktr stays single-order. Used when UseRegistrationAvailability
        'UseMaxwellDashpot', false, ...
        'mu_neg', '=mu', ...
        'mu2', 0, ...
        'ksrd2sr', 0, ...          % SRD→SR rate (s⁻¹); used when UseSuperRelaxed && UseSuperRelaxedADP
        'slope', 11000, ...        % A2 hopping rate per unit excess strain; used when UseA2AttachmentShift
        's_threshold_R', 0.0046, ... % positive strain threshold for A2 hopping (um)
        's_threshold_L', 5.5, ...  % negative strain threshold for A2 hopping (large = disabled by default)
        'd_actin', 0.0055, ...     % actin site repeat distance (um)
        'kSE_M', 500, ...          % Maxwell dashpot spring stiffness; used when UseMaxwellDashpot
        'eta_M', 1.0, ...          % Maxwell dashpot viscosity; used when UseMaxwellDashpot
        'FudgeA', 0, ...           % quadratic SL coeff for slack velocity; used when FudgeVmax
        'FudgeB', 0, ...           % linear SL coeff for slack velocity
        'FudgeC', 0, ...           % constant SL coeff for slack velocity
        'kstiff2_n', '=kstiff2', ... % P2 stiffness for negative strain; used when UseNegativeKstiff
        'kstiff1_n', '=kstiff1', ... % P1 stiffness for negative strain
        'k2d', 0, ...              % A2 mechanical recocking rate (s⁻¹); used when UseA2MechanicalRecocking
        'drmr', 0, ...             % strain threshold for mechanical recocking (um)
        'pop_s', -Inf, ...         % A2 popping strain threshold (um; -Inf = never triggers); used when UseA2Popping
        'Srd_max', 0, ...          % max SRD→PD transition rate; used when UseLimitedSRDTransition
        'Srd_n', 1, ...
        'UseJPattern', false, ...  % precompute Jacobian sparsity pattern for ode15s (default off)
        'UseFastPPval', true, ... % replace ppval with precomputed table lookup (default)
        'UseCompiledMex', false, ...  % use compiled C MEX for ODE RHS (requires compileMex.m; config-specific)
        'PieceWiseStrainDepR21X_logOffset', 1, ...
        'PieceWiseStrainDepR1DX_logOffset', 1, ...
        'PieceWiseStrainDep2X_logOffset', 1, ...
        'PieceWiseStrainDepX_logOffset', 1, ...
        'LoadPassiveMat', '', ...  % path to preprocessed driving signal .mat; if non-empty, loads sig into params.driveSig
        'xrate', 1, ...
        'PlotFeatureFitting', false, ...
        'FV_dataset', 'Baker2022' ... % or 'IsovelocityForFilip2021'
        );    

% , false, ...
% , false, ...
% , false, ...
 
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
    params0.kamh = '=0.1*kah'; % rate of reverse ATP hydrolysis (M·ADP·Pi → M·ATP)
    params0.kadh = 20; % % rate of ATP de-hydrolysis state change
    params0.ka  = 373.23;
    params0.kd  = 102.94; 
    params0.k1  = 40.116;%
    params0.k_1 = 17.103;%
    params0.k2  = 419.39;
    params0.k_2 = 2.7901; 
    params0.k3  = 44.255;%;
    params0.k3m = 0;
    params0.drp3 = 0.01;

    % transitions between super relaxed state and non relaxed state
    params0.ksr0   = 9.0853; % 
    params0.sigma1 = 33.125;
	params0.sigma2 = 1e6;
    params0.kmsr   = 250; % 
	
	% SRD state
	params0.kmsrd = 0;
    params0.ksrd = 0;
    params0.ksr2srd = 0;
    params0.sigma_srd1 = 33.125;
	params0.sigma_srd2 = 1e6;

    % Initial SRX fractions — probability "parked" in off-state heads.
    % Default 0 preserves prior behaviour. Set non-zero so that turning SRX
    % on/off only gates transitions, not probability mass (no cycling transient).
    params0.SRXT_0 = 0;   % initial P_SR fraction (super-relaxed, ATP-bound)
    params0.SRXD_0 = 0;   % initial P_SRD fraction (super-relaxed, ADP-bound)

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
    params0.alpha0_L = 0;
    params0.alpha0_R = 0;    
    params0.alpha1 = 15.14;
    params0.alpha_1 = 0;
    params0.alpha2 = 10.06;
    params0.alpha2_L = 0;
    params0.alpha2_R = 0;
    params0.alpha3 = 276.67;
    params0.s3     = 0.0099383;
    params0.alphaRip = 0;
    params0.k2rip = 0;

    % offsets
    params0.dr0 = 0;
    params0.dr1 = 0;
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
    params0.ekSE = 1;
    params0.c_SE_visc = 0;    % Kelvin-Voigt SE dashpot coefficient (0 = pure spring)
    params0.c_titin_visc = 0; % titin viscoelastic force coefficient: F += c_titin_visc * max(0,vel)
    params0.k_catch_bond = 0; % catch bond: detachment suppression factor during stretch
    params0.k_SA = 0;         % stretch activation: ka multiplier per unit velocity

    % passive force coeff
    params0.k_pas = 200; % From Kim Salla et al.

    % other
    params0.k2_L = 0;
    params0.k2_R  = 419.39;
    params0.dr2_L = 0;
    params0.dr2_R = 0.01;
            
    %% Step 1: Fill missing fields from defaults

    params = fillInDefaults(params, params0);

    %% Step 3: Apply multiplicative modifiers G for optimization
    if updateModifiers
        %     mods = {'kstiff1', 'kstiff2'};
        for i = 1:length(params.mods)
            try
                % initialize array mods
                if ~isfield(params, params.mods{i}) && contains(params.mods{i}, '__')
                    arr = split(params.mods{i}, '__');
                    if isfield(params, arr(1)) ...
                        && length(params.(arr{1})) > 1 ... % the array exists
                        && ~isempty(str2num(arr{2})) % the index is numeric
                        % ok, this is valid
                        params.(params.mods{i}) = params.(arr{1})(str2num(arr{2}));
                    end
                end

                params.(params.mods{i}) = params.(params.mods{i})*g(i);
            catch e
                    switch e.identifier
                        case 'MATLAB:nonExistentField'
                            error('Mod "%s" unknown. Initialize it first.', params.mods{i});
                        otherwise
                            error('Failed to update "%s": check g(%d). %s', params.mods{i}, i, e.message);
                    end
            end
        end

        % ensure we delete the modifiers right after
        params.mods = {};
        params.g = [];
        
        %% Step 4: Resolve linked parameters (fields starting with '=')
        params = resolveParams(params);
        
        %% Step 3a: Scale rates (if params.xrate is set, multiplies all turnover rates)
        params = updateRates(params);
    else
        %% Step 4: Resolve linked parameters (fields starting with '=')
        params = resolveParams(params);

    end


    

    %% Step 5: Reconstruct arrays: fields like params.arr__2 = 3 -> params.arr(2) = 3
    paramsfn = fieldnames(params);
    for i = 1:length(paramsfn)
        if ~contains(paramsfn{i}, '__')
            continue;
        end
        arr = split(paramsfn{i}, '__');
        if isfield(params, arr(1)) ...
            && length(params.(arr{1})) > 1 ... % the array exists
            && ~isempty(str2num(arr{2})) % the index is numeric
                % apply value to an array
                params.(arr{1})(str2num(arr{2})) = params.(paramsfn{i});
                % remove the field from the struct - it is already applied
                params = rmfield(params, paramsfn{i});
        end
    end

    %% Step 6: Build piecewise strain-dependent rate interpolants
    pwsdX = params.PieceWiseStrainDepX;
    pwsdP = params.PieceWiseStrainDepParams;

    % params.PieceWiseStrainDepFun = @(x)pchip(pwsdX, pwsdP, x);   
    params.PieceWiseStrainDep = pchip(pwsdX, pwsdP);

    % R2: own piecewise if provided, otherwise fall back to PieceWiseStrainDep shape
    if isfield(params, 'PieceWiseStrainDep2X') && isfield(params, 'PieceWiseStrainDep2Params') ...
            && ~isempty(params.PieceWiseStrainDep2X) && ~isempty(params.PieceWiseStrainDep2Params)
        params.PieceWiseStrainDep2 = pchip(params.PieceWiseStrainDep2X + log(params.PieceWiseStrainDep2X_logOffset)/10, params.PieceWiseStrainDep2Params);
    else
        params.PieceWiseStrainDep2 = params.PieceWiseStrainDep;
    end

    % R1D: Gaussian decay, own piecewise, or flat 1 (no strain dependence)
    if isfield(params, 'UseStrainDep4R1D') && params.UseStrainDep4R1D
        params.PieceWiseStrainDepR1D = pchip(pwsdX, exp(-(pwsdX.^2) / (2*params.kd_sigma^2)));
    elseif isfield(params, 'PieceWiseStrainDepR1DX') && isfield(params, 'PieceWiseStrainDepR1DParams') ...
            && ~isempty(params.PieceWiseStrainDepR1DX) && ~isempty(params.PieceWiseStrainDepR1DParams)
        params.PieceWiseStrainDepR1D = pchip(params.PieceWiseStrainDepR1DX + log(params.PieceWiseStrainDepR1DX_logOffset)/10, params.PieceWiseStrainDepR1DParams);
    else
        params.PieceWiseStrainDepR1D = pchip(pwsdX, ones(size(pwsdX)));
    end

    % R21: own piecewise if provided, otherwise zero (no reverse P2→P1 transition)
    if isfield(params, 'PieceWiseStrainDepR21X') && isfield(params, 'PieceWiseStrainDepR21Params') ...
            && ~isempty(params.PieceWiseStrainDepR21X) && ~isempty(params.PieceWiseStrainDepR21Params)
        params.PieceWiseStrainDepR21 = pchip(params.PieceWiseStrainDepR21X + log(params.PieceWiseStrainDepR21X_logOffset)/10, params.PieceWiseStrainDepR21Params);
    else
        params.PieceWiseStrainDepR21 = pchip(pwsdX, zeros(size(pwsdX)));
    end

    %% Step 6b (optional): Pre-unpack pchip structs for inline polynomial evaluation
    % Enabled by params.UseFastPPval = true. Requires UsePieceWiseStrainDep = true.
    % Stores raw breaks + coefficients from unmkpp so the ODE RHS can evaluate
    % the piecewise polynomial directly (no ppval/unmkpp call overhead per step).
    if params.UseFastPPval && isfield(params, 'UsePieceWiseStrainDep') && params.UsePieceWiseStrainDep
        % Store breaks as COLUMN vectors — this is critical for shape-safe indexing:
        % b(si) where b is column and si is (ss×1) always returns (ss×1).
        [brk, cof] = unmkpp(params.PieceWiseStrainDep);
        params.cached.pwsd_brk   = brk(:);   params.cached.pwsd_cof   = cof;

        if ~isempty(params.PieceWiseStrainDep2)
            [brk, cof] = unmkpp(params.PieceWiseStrainDep2);
            params.cached.pwsd2_brk  = brk(:);   params.cached.pwsd2_cof  = cof;
        else
            params.cached.pwsd2_brk = [];  params.cached.pwsd2_cof = [];
        end

        if ~isempty(params.PieceWiseStrainDepR1D)
            [brk, cof] = unmkpp(params.PieceWiseStrainDepR1D);
            params.cached.pwsdR1D_brk = brk(:);  params.cached.pwsdR1D_cof = cof;
        else
            params.cached.pwsdR1D_brk = [];  params.cached.pwsdR1D_cof = [];
        end

        if ~isempty(params.PieceWiseStrainDepR21)
            [brk, cof] = unmkpp(params.PieceWiseStrainDepR21);
            params.cached.pwsdR21_brk = brk(:);  params.cached.pwsdR21_cof = cof;
        else
            params.cached.pwsdR21_brk = [];  params.cached.pwsdR21_cof = [];
        end
    end

    %% Step 6c: Pack flat MEX parameter array (UseCompiledMex)
    % Must run AFTER step 6b (pchip unpacking) and AFTER step 7 (strain grid).
    % Deferred — called at end of this function after params.s is set.

    %% Step 7: Compute strain grid (params.s) from Slim_l, Slim_r, dS, SL0
    % Warn if SL0 is outside the physiological sarcomere length range
    if isfield(params, 'SL0') && (params.SL0 < 1.4 || params.SL0 > 3.0)
        warning('getParams: SL0=%.2f um is outside physiological range [1.4, 3.0]', params.SL0);
    end
    if params.UseCalculatedN
        % params.N = ceil((params.Slim_r - params.Slim_l)/params.dS/2);
        % params.LXBpivot = params.SL0;
        % params.LXBpivot = params.Slim_l + (params.Slim_r - params.Slim_l)/2;
        
        % this causes triubles for reinitializations
        % params.LXBpivot = params.SL0; 
        LXBpivot = params.SL0;

        % we have 2 half-sarcomeres, so the step is double
        SL_strain = (params.Slim_l:2*params.dS:params.Slim_r) - LXBpivot;
        if params.MaxStrainArraySize > 0 && params.MaxStrainArraySize < length(SL_strain)
            % SL_strain = -LXBpivot/2:params.dS:params.MaxStrainArraySize/2;
            SL_strain = ((-(params.MaxStrainArraySize-1)/2 : (params.MaxStrainArraySize-1)/2)) * params.dS*2;


            % % half-sarcomere strain limit - ready for contraction
            % if abs(params.Slim_r - LXBpivot) > abs(params.Slim_l - LXBpivot)
            %     % its closer to slim_r - lets cut off the Slim_l
            %     SL_strain = SL_strain (1:params.MaxStrainArraySize);                    
            % else
            %     % pivot closer to slim_l - lets cut off the slim_r
            %     SL_strain = SL_strain(end-params.MaxStrainArraySize+1:end);            
            % end
        end
        
        if params.LegacyStrainFlipping
            % we have two half-sarcomeres, so divide by 2
            params.s = SL_strain/2;
        else
            % the xontraction direction is positive
            params.s = flipud(-SL_strain')/2;
        end
        params.ss = length(params.s);


        % LSE = 0;
        % s = params.s' + (-(params.SL0 - LSE) + params.LXBpivot)/2;
        % s = flipud(-s);
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
    
    
    %% Step 8: Initialize state vector PU0
    expected_PU0 = params.NumberOfStates * params.ss + 7 + params.UseRegistrationAvailability;
    PU0_size_mismatch = isfield(params, 'PU0') && ~isempty(params.PU0) && numel(params.PU0) ~= expected_PU0;
    if ~isfield(params, 'PU0') || isempty(params.PU0) || updateInit || PU0_size_mismatch
        % only during init - otherwise keep it

        p0 = zeros(1, params.ss);
        U_SR  = params.SRXT_0;
        U_SRD = params.SRXD_0;
        NP = 0;
        PuATP = 0;
        SL0 = params.SL0*params.rsl0;
        params.LXBpivot = SL0;
        x_dash = 0; % X_visc: initial viscous stress deviation is zero
        LSE = params.LSE0;
        % State variable vector concatenates p1, p2, p2, and U_NR
        % Registration-availability: append a single scalar state A_reg at the tail
        % (index Ns*ss+8), initialized to the isometric floor A0. [] when the feature
        % is off => nothing appended => byte-identical to the legacy state vector.
        if params.UseRegistrationAvailability
            A_reg0 = params.A0;
        else
            A_reg0 = [];
        end
        if params.NumberOfStates == 2
            params.PU0 = [p0, p0, U_SR,NP,SL0,LSE, PuATP, U_SRD, x_dash, A_reg0];
        elseif params.NumberOfStates == 3
            params.PU0 = [p0, p0, p0,U_SR,NP,SL0,LSE, PuATP, U_SRD, x_dash, A_reg0];
		end
		
		if params.UseSuperRelaxedADP
			params.PU0 = [params.PU0];
        end
    end
    
%% Step 6c (deferred): Pack flat MEX parameter array
    % Requires: UseFastPPval=true (pchip cache), strain grid (params.s), and
    % all numerical params to be resolved. Called here so params.s is available.
    if params.UseCompiledMex && params.UseFastPPval && isfield(params,'s')
        params.mex_vals = packMexVals(params);
    end

%% Step 6d: Load preprocessed driving signal for passive force lookup
    % If params.LoadPassiveMat is a non-empty path, load the sig struct and
    % store it as params.driveSig. The file is produced by preprocessDrivingSignal.m.
    % Used by dPUdT_CombinedTransitions when params.LoadPassiveMat is set.
    if isfield(params, 'LoadPassiveMat') && ~isempty(params.LoadPassiveMat)
        tmp = load(['data\' params.LoadPassiveMat], 'sig');
        params.driveSig = tmp.sig;
    end

%% Step 9: Precompute attachment kernel (Gaussian/triangular, for UseA1AttachmentKernel)
    % Precompute kernel parameters
    if params.UseA1AttachmentKernel && params.A1AttachmentWidth > 0
        halfSpan = ceil(params.A1AttachmentWidth/params.dS);
        jj = -halfSpan:halfSpan;       % integer offsets
        sigma = params.A1AttachmentWidth/3;                  % Gaussian width
        
        params.A1AttachmentKernel = makeKernel(params.dS, params.A1AttachmentWidth, 'tri');
        % Precompute *centered* Gaussian kernel (unnormalized)
        K0 = exp(-( (jj*params.dS).^2 ) / (2*sigma^2));
        % K0 = max(0, 1 - abs(jj)*params.dS/params.A1AttachmentWidth);
        
        % Normalize so sum(K)*dS = 1
        K0 = K0' / (sum(K0) * params.dS);
        
        % Save constants
        params.cached.halfSpan = halfSpan;
        params.cached.jj       = jj;
        params.cached.K0       = K0;
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

function K = makeKernel(dS, Ax, kernelType)
    arguments
        dS (1,1) double {mustBePositive}
        Ax (1,1) double {mustBePositive}
        kernelType (1,:) char {mustBeMember(kernelType,{'tri','gauss'})} = 'gauss'
    end

    % number of grid points on each side
    halfSpan = ceil(Ax / dS);    % in bins
    
    % make odd length
    L = 2*halfSpan + 1;

    % spatial coordinates (centered)
    x = (-halfSpan:halfSpan) * dS;

    switch kernelType
        case 'tri'
            % triangular kernel: w = max(0, 1 - |x|/Ax)
            K = max(0, 1 - abs(x)/Ax);

        case 'gauss'
            % gaussian: sigma chosen so Ax spans roughly 3σ
            sigma = Ax/3;
            K = exp(-(x.^2)/(2*sigma^2));
    end

    % normalize so sum(K)*dS = 1
    K = K / (sum(K));
end
