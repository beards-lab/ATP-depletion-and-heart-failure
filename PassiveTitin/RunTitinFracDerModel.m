% function force = evaluateFracModel()
% addpath(genpath("1D_Trabeculae"))

% pars    = [10.46, 15.02, 0.9];

options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

pars = fminsearch(@evalModel, pars, options)
%%
pars = [16.1471    8.0026    0.7662];
% pars = [2.46, 10.02, 0.5];
figure(2);
evalModel(pars)

function cost = evalModel(pars)
import matlaws.*
import frac.*

pCa = 4.4;
rampSet = [1 2 3 4]; % only 100ms

% Material Parameters for myocardium
% pCa 11
% pars    = [2.46, 10.02, 0.5];
% pCa 4.4
% pars    = [16.1471    8.0026    0.7662];

Fss = 4; d = 0; b = 0; c = 2;
% Caputo Fractional Derivative Parameters
delta  = 0.023;
alpha  = 0.187;
% Prony Approximation terms
np     = 9;
Tf     = 160.0;
% bundle fractional parameters into struct
frac_pars = frac_parameters(alpha, delta, Tf, np);

order = 4; % Number of time steps
% paramaters of the simulation
nt   = 10^(order);
dt   = Tf/nt;
time = linspace(0, Tf, nt + 1);
F = [];

%% Load data
rds = fliplr([0.1, 1, 10, 100]);
for i_rd = 1:length(rds)
    if isinf(pCa) || pCa >= 10
        % newest format of experiments
        datatables{i_rd} = readtable(['..\..\Data\AvgRelaxed_' num2str(rds(i_rd)) 's.csv']);
    else
        if exist(['..\..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv'], "file")
            datatables{i_rd} = readtable(['..\..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv']);
        else
            datatables{i_rd} = [];
        end
    end
end
%% init the sim
L0 = 0;

Lmax = 1.175 - 0.95; % Lmax = 0.225; % identified to ramp of this height
Vlist = Lmax./rds; %  half-sarcomere velocity (um/s)
Force = cell(1, 5);
Time = cell(1, 5);
Length = cell(1, 5);
LengthPar = cell(1, 5);

for j = rampSet
    if isempty(datatables{j})
        fprintf('Skipping pCa %0.2f %0.0fs dataset\n', pCa, rds(j))
        continue;
    end

    displacement = Lmax/rds(j)*(time).*(time < rds(j) & time > 0) + Lmax*(time>=rds(j));

    % run the simulation
    Force_v{j} = diffeq_sim(@trabeculae3D, pars, displacement, dt, 3, frac_pars);
    Length{j} = displacement';
    Time{j} = time;
    a = (Fss - d)/((Lmax -b)^c);
    Force_par{j} = a*max(Length{j} - b, 0).^c + d;
    % calc force
    Force{j} = Force_v{j} + Force_par{j};

end

% Get error for the whole ramp-up and decay
t_endFreeware = zeros(1, 5); % time when we start counting the costs
% alternatively, get the error from decay only
% t_endFreeware  = Lmax./Vlist + 2;

%% Evaluating all ramps at once
En = cell(1, length(rampSet));
Es = cell(0);
% PeakModel = zeros(size(rampSet));
for j = rampSet
    datatable_cur = datatables{j};
    if isempty(datatable_cur)
        continue;
    end
    inds = find(datatable_cur.Time >= t_endFreeware(j));
    datatable_cur = datatable_cur(inds, :);
    t_int{j} = datatable_cur.Time - 2;
    % t_int{j} = datatable_cur.Time;
    Ftot_int{j} = interp1(Time{j}, Force{j}, t_int{j}); % total force interpolated
    Es{j} = ((Ftot_int{j} - datatable_cur.F)/max(Ftot_int{j})).^2; % error set
    Es{j} = Es{j}(~isnan(Es{j})); % zero outside bounds
    En{j} = 1e3*sum(Es{j})/length(Es{j}); % normalized error

    PeakData(j, 1) = rds(j);
    PeakData(j, 2) = max(datatable_cur.F);
    % PeakModel(j) = max(Ftot_int{j}(:));
end

PeakModel = nan(1, length(PeakData));
for j = 1:length(PeakModel)
    m = max(Force{j});
    if ~isempty(m)
        PeakModel(j) = m;
    end
end

Ep = nansum(((PeakData(:, 2) - PeakModel(:))/max(Ftot_int{j})).^2);
% discarding peak fit for high Ca's
if pCa < 10
    Ep = 0;
end

cost = Ep*100 + sum([En{1:end}], 'all');

% if ~plotAll
%     return;
% end
%%
% return;
figure(1);clf;nexttile;
for j = 1:length(rds)
    % plot(datatables{j}.Time, datatables{j}.L);hold on;
    % plot(Time{j}, Length{j})
    semilogx(datatables{j}.Time - 2, datatables{j}.F, 'k');hold on;
    semilogx(t_int{j}, Ftot_int{j}, '-', LineWidth=2)
end
title('Data and simulation for 0.1s - 100s ramps in semilog')
legend('Data', 'Fractional derivative simulation')
xlabel("Time (s)")
ylabel('Stress (kPa)')
nexttile;
for j = 1:length(rds)
    % plot(datatables{j}.Time, datatables{j}.L);hold on;
    % plot(Time{j}, Length{j})
    loglog(datatables{j}.Time - 2, datatables{j}.F - Force_par{j}(end), 'k');hold on;
    loglog(Time{j}, Force_v{j}, '-', LineWidth=2)
end
ylim([1 inf])
xlabel("Time (s)")
ylabel('Stress - Stress_{parallel} (kPa)', Interpreter='tex');
title('Loglog view of the dynamic part of the stress only')
end