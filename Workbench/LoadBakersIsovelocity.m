%% ISovelocity table

opts = detectImportOptions('../data/2021 06 15 isovelocity fit Filip.xlsx', 'Sheet', 1); % Specify the sheet to read from
opts.DataRange = 'A5'; % Set the data range starting from the 5th row
dataTable = readtable('../data/2021 06 15 isovelocity fit Filip.xlsx', opts); % Load the data into a table

% Limit to 3 columns and rename variables
dataTable = dataTable(:, 1:3);
dataTable.Properties.VariableNames = {'t', 'L', 'F'};
%%
figure(1); clf; hold on;
plot(dataTable.t, dataTable.L*100, 'DisplayName', 'L');
plot(dataTable.t, dataTable.F, 'DisplayName', 'F');
xlabel('Time (s)');
ylabel('Values');
title('Plot of L and F over Time');
legend show;
grid on;

%% Fitting ktr to force recover from drop after first peak after slack
% ktr8 = load('data/bakers_ktr_8.mat');
% win = ktr8.datatable(:, 1) > 0;
% Subsample the data to reduce density
% subsampleFactor = 10; % Adjust this factor as needed for desired density
% ft = ktr8.datatable(win, 1);
% ff = ktr8.datatable(win, 3);
% ft = ft(1:subsampleFactor:end);
% ff = ff(1:subsampleFactor:end);


threshold = 0.01; % Define a threshold for fast rise
riseVelocities = find(diff(dataTable.L)./diff(dataTable.t) > threshold) + 1; % Find indices of fast rises
riseIndices = riseVelocities([false; diff(riseVelocities) > 1000]);

winstart = 480 - 480;
winsize = 4000;
% Align F to start at these bump-ups
alignedF = zeros(size(dataTable.F)); % Initialize aligned F
% Plot the aligned data
figure(2); clf; 
for i = 1:length(riseIndices)
    win = riseIndices(i) + winstart:min(length(dataTable.L), riseIndices(i)+winsize);
    % nexttile(1);hold on;plot(dataTable.t(win) - dataTable.t(win(1)), dataTable.L(win), 'DisplayName', 'L');   
    nexttile();hold on;
    ft = (dataTable.t(win) - dataTable.t(win(1)) +1)/1000; % conversion to s from ms
    ff = dataTable.F(win);
    plot(ft, ff, 'DisplayName', 'F');   
    
    y_exp = @(A, k, t0, B, x) max(0, A*(1-exp(-(x-t0)*k)) + B);
    
    y_2nd = @(A, zeta, wn, t0, B, x) ...
    (A*(1 - 1./sqrt(1-zeta.^2) .* ...
    exp(-zeta*wn.*(x-t0)) .* ...
    sin(wn*sqrt(1-zeta.^2).*(x-t0) + acos(zeta))) + B);
    
    % (x < t0).*B;

    [ae be] = fit(ft, ff, y_exp, 'StartPoint', [50, 30, 0.01, 50], 'Lower', [10, 1, 0.0, 0], 'Upper', [200 100 0.0, 100]);

    init_vals = [25.6, 0.7, 0.046*1000, 0, 41.5];
    
    [ae2 be2] = fit(ft, ff, y_2nd, 'StartPoint', init_vals, 'Lower', init_vals/10, 'Upper', init_vals*10);

    
    K(i) = ae.k;
    A(i) = ae.A;
    t0(i) = ae.t0;
    B(i) = ae.B;
    
    
    % Store outputs of ae2 in independent arrays
    A2(i) = ae2.A;
    zeta(i) = ae2.zeta;
    wn(i) = ae2.wn;
    t02(i) = ae2.t0;
    B2(i) = ae2.B;

    good(i) = (be.rmse - be2.rmse)/be.rmse*100;

    plot(ft, ae(ft),'--', LineWidth=2);
    plot(ft, ae2(ft), ':', LineWidth=3);
    % plot(ft, y_exp(50, 0.05, 0.01, 50, ft), ':', LineWidth=2)
end

% Quadratic coefficients: k1^2 - (k1+k2)*k1 + k1*k2 = 0


a = 1;
b = -2*zeta.*wn;
c = wn.^2;

discriminant = b.^2 - 4*a.*c;

if discriminant < 0
    warning('Underdamped system: choose k3 as free parameter; rates will remain real and positive.');
end

% Use positive real solution for k1
k1 = (-b + sqrt(discriminant))/2;

% Step 2: Compute k2 and k3
k2 = wn.^2 ./ k1;
k3 = 2*zeta.*wn - k1 - k2;
i = 1:length(K);
nexttile; plot(i, K, '+-', i, A, 'x-', i, B, '|-'); title('In order'); legend('K', 'A', 'offset');
nexttile; plot(i, A2, '+-', i, zeta, 'x-', i, wn, '|-',i, t0, '|-', i, B2, '|-'); title('In order for A2, zeta, wn'); legend('A2', 'zeta', 'wn', 't0', 'b2');
nexttile; plot(good, 'x-')

% i = [15 15 15 15 15 10];
% nexttile; plot(i, K, '+-', i, A, 'x-', i, B, '|-');title('By speed'); legend('K', 'A', 'offset');
% i = [1.1 1.1 1.1 1.1 1.1 1.0];
% nexttile; plot(i, K, '+-', i, A, 'x-', i, B, '|-');title('By height'); legend('K', 'A', 'offset');


%% Plot the aligned data
figure(2); clf; hold on;
plot(dataTable.t, dataTable.L * 100, 'DisplayName', 'L');
plot(dataTable.t, alignedF, 'DisplayName', 'Aligned F');
xlabel('Time (s)');
ylabel('Values');
title('Aligned F with Fast Rise in L');
legend show;
grid on;

% opts = detectImportOptions(filename, 'Sheet', 1, ); % Specify the range starting from row 5 and first three columns
% opts.VariableNamesLine = 3; % Set the variable names from the 3rd row
% dataTable = readtable(filename, opts); % Load the data into a table

%% load isovelocity

isov = load('../data/isovelocity.mat');
datable = isov.datatable;

clf;
t = datable(:, 1);
f = datable(:, 3);
l = datable(:, 2);
v = [0; diff(l)./diff(t)];
N = 1:length(l);
plot(t, v, t, l);hold on;

%% extract manually
isov = load('../data/isovelocity.mat');
datable = isov.datatable;

clf;
t = datable(:, 1);
f = datable(:, 3);
l = datable(:, 2);

% Sample data (replace with your actual signal)
win = t < 4000 & t > 0;
x = t(win)/1000;
y = l(win);
f = f(win);
v = [0;diff(y)./diff(x)]/100;


plot(x, y, x, v); hold on;

% segment 1
velocities1 = [0, -0.15, -6e-3, -0.5, 0, 0.015];
positions1 = [500.25, 500.65, 509.2, 509.58, 519.5, 539.55];

% segment 2
velocities2 = [0, -0.05, -1e-3, -0.5, 0, 0.015];
positions2 = 500 + [500.25, 500.65, 509.2+82, 509.58+82, 519.5+82, 539.55+82];

% segment 3
velocities3 = [0, -0.1, -3e-3, -0.5, 0, 0.015];
positions3 = 1000 + [500.25, 500.65, 509.2+15, 509.58+15, 519.5+15, 539.55+15];

% segment 4
velocities4 = [0, -0.145, -5e-3, -0.5, 0, 0.015];
positions4 = 1500 + [500.25, 500.65, 509.2+1.8, 509.58+1.8, 519.5+2, 539.55+2];

% segment 5
velocities5 = [0, -0.025, -0.5e-3, -0.5, 0, 0.015];
positions5 = 2000 + [500.25, 500.65, 509.2+192, 509.58+192, 519.5+192, 539.55+192];

% segment 6
velocities6 = [0, -0.13, -4e-3, -0.5, 0, 0.015];
positions6 = 2600 + [500.25, 500.65, 509.2 + 7, 509.58+7, 519.0+7, 539.2+7];

% segment 7
velocities7 = [0, -0.1, -2e-3, -0.5, 0, 0.01];
positions7 = 3100 + [500.25, 500.65, 509.2+27, 509.58+27, 519.5+27, 539.55+27];

velocities = [0, velocities1,velocities2,velocities3,velocities4,velocities5,velocities6,velocities7, 0]*1000;
positions = [0, positions1,positions2,positions3,positions4,positions5,positions6,positions7, 4000]/1000;

segs = 1.1;

for i = 2:length(velocities)
    segs(i) = segs(i-1) + velocities(i)*(positions(i) - positions(i-1));   % use i, not 1
end

% segs now contains cumulative values after each segment
% disp(segs);
plot(positions, segs, '-o');
% xlim([positions(2)-1, positions(end)]);
ylim([0.75, 1.2])
xlabel('segment index'); ylabel('cumulative value');
grid on;


datatable = struct();
datatable.t = x;
datatable.L = 2*y;
datatable.F = f;
datatable = struct2array(datatable);
datatable = datatable(1:10:end, :);
velocitytable = [positions; [velocities(2:end), 0]; [velocities(2:end), 0]/2; segs*2]';

% save('../data/bakers_isovelocity.mat', 'datatable', 'velocitytable');
%% align to start
% [positions1(1),positions2(1),positions3(1),positions4(1),positions5(1),positions6(1),positions7(1)]
ta = x - [positions1(1),positions2(1),positions3(1),positions4(1),positions5(1),positions6(1),positions7(1)]/1000;
clf;
plot(ta, repmat(f, [1, 7]))
%%
d = load('../data/bakers_slack8mM_all_fix.mat');


%% 

params0 = getParams();
ModelParamsInitNiceSlack;

params0.velocitytableonfile = '../data/bakers_isovelocity.mat';
RunBakersExp;
