function [E_slack, out_slack, features_model, features_data] = runSlackExperiment(params0)
%RUNSLACKEXPERIMENT Run the slack (repeated shortening/re-stretch) protocol.
%
%   [E_slack, out_slack, features_model] = runSlackExperiment(params0)
%
%   Loads Baker slack data, builds the velocity table according to
%   params0.RunSlackSegments, runs the simulation (optionally in parallel
%   chunks for the 'AllPar' case), and returns the cost vector, merged
%   output struct, and extracted feature struct.
%
%   Inputs:
%     params0  - Base parameter struct. Must include params0.modelFcn (string).
%                Key fields used:
%                  RunSlackSegments    - String selecting the protocol segment.
%                  EvalFitSlackOnset   - If true, also fit the force onset.
%                  EvalFeatures        - If true, extract slack features.
%                  PlotEachSeparately  - If true, produce plots.
%                  ShowStatePlots      - If true, add a state-probability tile.
%                  ShowResidualPlots   - If true, plot residuals in figure 1001.
%                  ghostLoad / ghostSave - Ghost management strings.
%
%   Outputs:
%     E_slack        - Row vector of cost values:
%                        E_slack(1) = main force-trace MSE (scaled by 20)
%                        E_slack(2) = slack-onset dt costparam (if EvalFitSlackOnset)
%                        E_slack(3) = slack-onset params.ktr cost (if EvalFitSlackOnset)
%                      Elements 2–3 are omitted when EvalFitSlackOnset is false.
%     out_slack      - Merged output struct from all simulation chunks.
%     features_model - Struct of simulated slack features (populated only when
%                      params0.EvalFeatures is true; otherwise empty struct).
%     features_data  - Struct of reference data slack features (populated only
%                      when params0.EvalFeatures is true; otherwise empty struct).

    modelFcn = str2func(params0.modelFcn);

    % recalculateDataFeats = params0.recalculateDataFeats;
    features_model = struct();
    features_data  = struct();

    params = params0;

    if isfield(params, 'velocitytableonfile')
        datastruct = load(['data/' params.velocitytableonfile]);
    else
        datastruct = load('data/bakers_slack8mM_all.mat');
    end
    datatable = datastruct.datatable;
    if isfield(datastruct, 'features_data')
        features_data = datastruct.features_data;
    end

    

    validZone = datatable(:, 1) > 1;

    datastruct.velocitytable(1, 1) = -20;
    PU0 = [];
    par_velocitytable = [];
    parple = [];

    switch params.RunSlackSegments
        % first slack
        case 'First'
            velocitytable = datastruct.velocitytable(1:7, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1);
        case 'FirstTwo'
            % two slacks
            velocitytable = datastruct.velocitytable(1:11, :);
            validZone = datatable(:, 1) > 1;
        case 'Fourth-rampuponly'
            params.SL0 = 1.9194;
            % only the last slack
            velocitytable = [datastruct.velocitytable(1, :); datastruct.velocitytable(16:19, :)];
            % ramp-up and peak
            validZone = datatable(:, 1) > datastruct.velocitytable(17, 1) & datatable(:, 1) < datastruct.velocitytable(19, 1) - 0.1;
        case 'Fourth'
            params.SL0 = 2.2;
            % only the pre-last slack
            velocitytable = [datastruct.velocitytable(1, :); datastruct.velocitytable(15:19, :)];
            % ramp-up and peak
            validZone = datatable(:, 1) > datastruct.velocitytable(15, 1) - 0.1 & datatable(:, 1) < datastruct.velocitytable(19, 1) - 0.1;
        case 'FirstAndLast'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) | datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(21, 1);
        case 'FirstAndLastShort'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            velocitytable([5:6 9:10], 1) = velocitytable([5:6, 9:10], 1) + 1;
            validZone = datatable(:, 1) > datastruct.velocitytable(3, 1) - 0.1 & datatable(:, 1) < datastruct.velocitytable(3, 1)+0.04 ...
                | datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(19, 1) + 0.04;
        case 'FirstAndLastExtended'
            velocitytable = datastruct.velocitytable([1:6 19:23], :);
            velocitytable([5:6 9:10], 1) = velocitytable([5:6, 9:10], 1) + 1;
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) ...
                | datatable(:, 1) > datastruct.velocitytable(18, 1) & datatable(:, 1) < datastruct.velocitytable(21, 1);
        case 'AllButLast'
            % all but the last
            velocitytable = datastruct.velocitytable(1:19, :);
            validZone = datatable(:, 1) > 1;
        case 'Last'
            % only the last slack
            velocitytable = datastruct.velocitytable(18:end, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(23, 1);
            velocitytable(1, 1) = -2;
        case 'All'
            % all
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1);
        case 'AllPar'
            % all, but run all in parallel
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1);

            % 1. precalculate the init to be reused
            if isfield(params, 'PU0')
                params = rmfield(params, 'PU0');
            end
            params.Velocity = 0;
            [~, out_init] = evaluateModel(modelFcn, [-10 velocitytable(3, 1)], params);
            PU0 = out_init.PU(end, :);
            LXBpivot0 = out_init.params.LXBpivot;   % capture pivot after init run

            wins = 2:4:(size(velocitytable, 1) - 1);

            % Preallocate cell array
            par_velocitytable = cell(length(wins), 1);

            % Fill the cells
            for k = 1:length(wins)
                par_velocitytable{k} = velocitytable(wins(k):wins(k)+4, :);
            end

            if isempty(gcp('nocreate')) || params.justPlotStateTransitionsFlag
                parple = [];
            else
                parple = gcp('nocreate');
                N_par = length(par_velocitytable);
                futures = cell(1, N_par);
            end

        case 'AllNoBump'
            % all, but evaluate only force onset, not the restretch
            velocitytable = datastruct.velocitytable(1:end, :);
            validZone = datatable(:, 1) > datastruct.velocitytable(2, 1) & datatable(:, 1) < datastruct.velocitytable(5, 1) |...
                datatable(:, 1) > datastruct.velocitytable(7, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(9, 1) |...
                datatable(:, 1) > datastruct.velocitytable(11, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(13, 1) |...
                datatable(:, 1) > datastruct.velocitytable(15, 1) - 0.05 & datatable(:, 1) < datastruct.velocitytable(17, 1) |...
                datatable(:, 1) > datastruct.velocitytable(19, 1)-.1 & datatable(:, 1) < datastruct.velocitytable(21, 1);
        case 'ramp-up'
            velocitytable = [-5, 0; 0, 1; .15, 0; 1, 0];
            params.SL0 = 1.9;
            validZone = datatable(:, 1) > 1;
        case 'ramp-down'
            params.SL0 = 2.2; %#ok<NASGU> (assigned but only used via params downstream)
            velocitytable = [-2, 0; 0, -.01; 15, 0; 20, 0];
        case 'stairs-down'
            rampvel  = -0.05;
            rampup   = 0.5;
            ramphold = 1;
            params.SL0 = 2.2;
            velocitytable = [-2, 0; 0, 0];
            for iv = 2:5
                velocitytable = [velocitytable; velocitytable(end, 1) + ramphold, rampvel/rampup; velocitytable(end, 1) + rampup+ramphold, 0]; %#ok<AGROW>
            end
            velocitytable = [velocitytable; velocitytable(end, 1) + ramphold*2, 0];
        case 'stairs-up'
            rampvel  = 0.05;
            rampup   = 0.5;
            ramphold = 10;
            params.SL0 = 1.9;
            velocitytable = [-2, 0; 0, 0];
            for iv = 2:5
                velocitytable = [velocitytable; velocitytable(end, 1) + ramphold, rampvel/rampup; velocitytable(end, 1) + rampup+ramphold, 0]; %#ok<AGROW>
            end
            velocitytable = [velocitytable; velocitytable(end, 1) + ramphold*2, 0];
    end

    %% Run simulation (possibly in parallel chunks)
    if isempty(par_velocitytable)
        par_velocitytable = {velocitytable};
    end
    N = length(par_velocitytable);
    outs = cell(1, N);
    params = getParams(params, params.g, true);
    if exist('LXBpivot0', 'var')            % only set for AllPar branch
        params.LXBpivot = LXBpivot0;        % restore pivot from end of init sim
    end

    for i_par_chunk = 1:N
        velocitytable_chunk = par_velocitytable{i_par_chunk};
        params.Velocity = velocitytable_chunk(:, 2);
        params.PU0 = PU0;

        if isempty(parple)
            [~, out_chunk] = evaluateModel(modelFcn, velocitytable_chunk(:, 1), params);
            outs{i_par_chunk} = out_chunk;
        else
            futures{i_par_chunk} = parfeval(parple, @evaluateModel, 2, modelFcn, velocitytable_chunk(:, 1), params);
        end
    end

    if ~isempty(parple)
        for j = 1:N
            if ~isempty(futures{j})
                [~, out_chunk] = fetchOutputs(futures{j});
                outs{j} = out_chunk;
            end
        end
    end
    out_slack = mergeOutStructs([outs{:}]);
    out_slack.datatable  = datatable(validZone, :);
    % out_slack.validZone  = validZone;

    %% Cost function
    nonrepeating = makeMonotonous(out_slack.t);
    Fi = interp1(out_slack.t(nonrepeating), out_slack.Force(nonrepeating), datatable(validZone, 1));

    e = (datatable(validZone, 3) - max(0, Fi)).^2;
    E_slack = mean(e(~isnan(e))) * 20;
    % plot(datatable(validZone, 1), 10*cumsum(e)/length(e), datatable(validZone, 1), datatable(validZone, 3))

    %% peaks for the passive - TODO CLEAR
    
    if params.EvalPeaks
        % peaks_data = findpeaks(datatable(:, 3), datatable(:, 1),  MinPeakProminence=10, MinPeakHeight=10, Annotate="peaks");
        peaks_data =    [21.4488;   21.1131;    23.4577;   24.4631];
        peaks_model = findpeaks(out_slack.Force(nonrepeating), out_slack.t(nonrepeating), MinPeakProminence=10, MinPeakHeight=10, Annotate="peaks")';
        N = min(length(peaks_data), length(peaks_model));
        E_p = 0.5*sum((peaks_data(1:N) - peaks_model(1:N)).^2);
        E_slack = E_slack + E_p;
    end

    %% Plotting
    if params.PlotEachSeparately
        nexttile;
        yyaxis right;
        pl1 = plot(datatable(:, 1), datatable(:, 2), '-', out_slack.t, out_slack.SL, '--', out_slack.t, out_slack.LXB, ':', 'Linewidth', 2, 'MarkerSize', 3);

        yyaxis left; hold on;

        % manage GHOST
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_slack.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_slack']);
            ghost = ghost.ghost;
            gp = plot(ghost(:, 1), ghost(:, 2), '-', 'Linewidth', 3, 'Color', [0.5843    0.8157    0.9882]);
        else
            clear gp;
        end

        plot(datatable(:, 1), datatable(:, 3), 'k-', 'linewidth', 0.5);
        pl5 = plot(datatable(validZone, 1), datatable(validZone, 3), 'k-', 'linewidth', 1.5);

        pl3 = plot(out_slack.t, out_slack.Force, 'b-', 'linewidth', 1.5);
        pl4 = plot(out_slack.t, out_slack.FXBPassive, 'b--', 'linewidth', 1.5);

        xlabel('$t$ (sec.)','interpreter','latex','fontsize',16);
        ylabel('Force (rel.)','interpreter','latex','fontsize',16);
        set(gca,'fontsize',14); box on;
        title('Slack');
        xlim([par_velocitytable{1}(2, 1)-0.02, velocitytable(end, 1)])
        yl = ylim;
        ylim(yl)
        xl = xlim();

        if exist('gp', 'var') && isvalid(gp)
            legend([gp; pl1(1:3); pl5; pl3; pl4], ['Ghost ' params.ghostLoad], 'MLx2.0', 'SL*', 'LXB*', 'Force (data)', 'Force (sim)', 'Passive (sim)', 'Location', 'best');
        else
            legend([pl1(1:3); pl5; pl3; pl4], 'MLx2.0', 'SL*', 'LXB*', 'Force (data)', 'Force (sim)', 'Passive (sim)', 'Location', 'best');
        end

        if params.ShowStatePlots
            nexttile;
            if params.NumberOfStates == 2
                plot(out_slack.t, out_slack.p1_0, '-', out_slack.t, out_slack.p2_0, '-', ...
                    out_slack.t, out_slack.PuATP, '-', out_slack.t, out_slack.PuR, '-', out_slack.t, out_slack.SR, out_slack.t, out_slack.SRD, 'LineWidth', 1.5, 'LineStyle', '-')
                legend('A1','A2','UT','UD','SR','SRD')
            elseif params.NumberOfStates == 3
                plot(out_slack.t, out_slack.p1_0, '-', out_slack.t, out_slack.p2_0, '-', out_slack.t, out_slack.p3_0, '-', ...
                    out_slack.t, out_slack.PuATP, '-', out_slack.t, out_slack.PuR, '-', out_slack.t, out_slack.SR, out_slack.t, out_slack.SRD, 'LineWidth', 1.5, 'LineStyle', '-')
                legend('A1','A2','A3','UT','UD','SR','SRD')
            end
            xlim(xl);
        end
    end

    if params.ShowResidualPlots
        figure(1001);
        if ~isempty(params.ghostLoad) && exist(['Ghost_' params.ghostLoad '_slack.mat'],'file')
            ghost = load(['Ghost_' params.ghostLoad '_slack']);
            ghost = ghost.ghost;
            Gi = interp1(ghost(:, 1), ghost(:, 2), datatable(validZone, 1));
            plot(datatable(validZone, 1), (Fi - Gi), 'linewidth', 2); hold on;
            plot(datatable(validZone, 1), (Gi - datatable(validZone, 3)), 'k');
        else
            plot(datatable(validZone, 1), (Fi - datatable(validZone, 1)));
        end
    end

    if ~isempty(params.ghostSave)
        ghost = [out_slack.t; out_slack.Force]';
        % do not save ghostSave command!
        gs = params0.ghostSave; params0.ghostSave = '';
        save(['Ghost_' params.ghostSave '_slack'], 'ghost', 'params0');
        params0.ghostSave = gs;
    end

    %% SAVE FIG
    if params.PlotEachSeparately
        fig = gcf;
        if ~isempty(params.ghostSave)
            saveas(fig, ['XBBakersDataFit_' params.ghostSave '.png']);
        end
    end

    %% Optional: fit slack onset timing
    if params.EvalFitSlackOnset
        try
            [e_dt, e_ktr] = fitSlackForceOnset(datatable, velocitytable, out_slack.t, out_slack.SL, out_slack.Force, params.PlotEachSeparately & params.drawForceOnset);
            E_slack(2) = 0*5e6*e_dt;
            E_slack(3) = 0*e_ktr;
        catch
            E_slack(2) = 1e3;
            warning('Eval slack onset failed');
        end
    end

    %% Optional: extract features if missing

    if params.EvalFeatures
        if isempty(features_data) && recalculateDataFeats
            features_data = extractSlackAttributes(datatable(:, 1), datatable(:, 3), datatable(:, 2), velocitytable, features_data, [], false);
            % Print extracted values to console so they can be pasted into the else branch below
            fieldNames = fieldnames(features_data);
            for k = 1:numel(fieldNames)
                fname = fieldNames{k};
                val   = features_data.(fname);
                if ~isnumeric(val), continue; end
                s = mat2str(round(val, 4, 'significant'));
                fprintf('features_data.%s = %s;\n', fname, s);
            end
        elseif isempty(features_data)
            % Hardcoded Baker lab 8 mM slack reference values (extracted 2024-11)
            features_data.ktr                 = [39.76 34.18 31.63 28.36 27.6];
            features_data.A                   = [68.4 65.47 58.92 52.2 43.51];
            features_data.t0                  = [0.002165 0.004643 0.007818 0.01117 0.01589];
            features_data.Am                  = [58.18 55.29 52.1 47.33 41.98];
            features_data.SLslack             = [2.04 2 1.96 1.92 1.88];
            features_data.SLdiff              = [0.162 0.202 0.242 0.282 0.322];
            features_data.restretchSlopeStart = [1131 1011 1006 907 953.4];
            features_data.restretchSlopeEnd   = [111.3 122.2 108.7 119.7 135.3];
            features_data.v_restretch         = [4 5 6 7 3];
            features_data.peak1_y             = [77.59 76.02 74.2 71.11 56.89];
            features_data.peak1_t             = [0.0035 0.003 0.0025 0.002 0.005];
            features_data.peak1_SL            = [2.068 2.034 1.99 1.954 1.912];
            features_data.peak1_dSL           = [0.028 0.034 0.03 0.034 0.032];
            features_data.vall_y              = [68.97 66.95 65.33 62.49 54.37];
            features_data.vall_t              = [0.0085 0.006 0.006 0.005 0.0095];
            features_data.peak2               = [77.22 78.47 80.01 82.42 61.85];
            features_data.steady              = [76.91 76.48 75.98 75.53 56.75];
            features_data.vall2_dy            = [-13.56 -13.84 -13.6 -13.08 -8.826];
            features_data.vall2_t             = [0.0085 0.0055 0.008 0.0065 0.0095];
            features_data.ovrsht_dy           = [0.431 1.392 1.244 1.552 0.3205];
            features_data.ovrsht_t            = [0.2116 0.225 0.247 0.239 0.2676];
            features_data.XTOR               = [10 10 10 10 10];
        end        
    end
    features_model = extractSlackAttributes(out_slack.t, out_slack.Force, out_slack.SL, velocitytable, features_model, out_slack, params0.PlotFeatureFitting);
end
