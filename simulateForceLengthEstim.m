function Es = simulateForceLengthEstim(g, params, plotThat)

if nargin < 2
    plotThat = false;
end

SL = [2.2 2.0400    2.0000    1.9600    1.9200  1.8800];
df = [76.5 68.3878   65.5438   59.2362   52.5507   43.9655];

% SL = 1.8:0.05:2.1; 
% df = ones(size(SL));

for i = 1:length(SL)
    params.SL0 = SL(i);
    params.Velocity = 0;
    if isfield(params, 'PU0')
        params = rmfield(params, 'PU0');
    end

    % reset the PU0
    params = getParams(params, params.g, true);
    
    % [F out] = evaluateModel(modelFcn, velocitytable(:, 1), params);
    [F out] = evaluateModel(@dPUdTCaSimpleAlternative2State, [0 100], params);
    F_total(i) = F;
    F_passive(i) = out.FXBPassive(end);
p1_0(i) = out.p1_0(end);
p2_0(i) = out.p2_0(end);
PuT(i) = out.PuATP(end);
PuD(i) = out.PuR(end);
SR(i) = out.SR(end);
end
 

if plotThat
    %%
    % figure(2); clf;
    nexttile;hold on;
    plot(SL, df, 'o-', LineWidth=2)
    plot(SL,F_total,'x-', SL,F_passive,'x--', LineWidth=2);
    title('Steady-state tension');xlabel('SL');ylabel('Tension (kPa)');
    legend('Data (estimated)', 'Simulation', 'Passive tension component (model)')
    nexttile();
    plot(SL, p1_0, '|-', SL, p2_0, '|-', SL, PuT, '|-',SL, PuD, '|-', SL, SR, '|-', LineWidth=1.5)
    legend('A1','A2','UT','UD', 'SR')

end

E = (F_total - df).^2;
Es = sum(E);

end