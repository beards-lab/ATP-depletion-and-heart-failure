S1 = dir('../data/PassiveCaSrc2/');

S1 = S1([S1.isdir]);
fmaxdata = {};

for i_s1 = 1:length(S1)

    % fmaxdata = readtable('../data/PassiveCaSrc2/20241220/02_03_Fmax_stiff_ktr.txt');
    fname = ['../data/PassiveCaSrc2/' S1(i_s1).name '/02_03_Fmax_stiff_ktr.txt'];
    if ~exist(fname, "file")
        %isempty(str2num(S1(i_s1).name))
        continue;
    end
    fmaxdata{end+1} = readtable(fname);
end


%%
figure(404);clf;

datatable = table2array([1 2.0 1].*fmaxdata{2}(:, ["x_s_", "x_Lo_", "x_kPa_"]));

for i_fm = 1:length(fmaxdata)

validRng = abs(fmaxdata{i_fm}.x_Lo_ - 1) < 0.002;
t = fmaxdata{i_fm}.x_s_;
F = fmaxdata{i_fm}.x_kPa_;


% 1. Fit a linear model to the data: F = p1*t + p2
% 'polyfit' returns [slope, intercept]
p = polyfit(t(validRng), F(validRng), 1); 

% 2. Calculate the drift values based on the model
drift = polyval(p, t);

% 3. Subtract the drift from the original data
F_corrected = F;
% F_corrected = F - drift;

% Optional: Add the initial offset back if you want to preserve the starting magnitude
% F_corrected = F_corrected + drift(1);

% Plotting the results
% figure(403);clf;
% plot(t, F, 'r', 'DisplayName', 'Original (Drifting)'); hold on;
% plot(t, drift, 'k--', 'DisplayName', 'Linear Trend');
% plot(t, F_corrected, 'b', 'DisplayName', 'Corrected');
% legend; xlabel('Time (t)'); ylabel('Force (F)');

    normf = fmaxdata{i_fm}.x_kPa_(end);
    normf = F_corrected(end);
    normf = 1;
    nexttile(1);hold on;
    % plot(fmaxdata{i_fm}.x_s_, fmaxdata{i_fm}.x_kPa_/normf);
    plot(t, F_corrected/normf);
    nexttile(2);hold on;
    % plot(fmaxdata{i_fm}.x_s_, fmaxdata{i_fm}.x_Lo_*2);
    % plot(fmaxdata{i_fm}.x_Lo_,fmaxdata{i_fm}.x_um_);
    plot(fmaxdata{i_fm}.x_s_, fmaxdata{i_fm}.x_um_);


%% 1. Definice modelu: f(t) = 1 - A * exp(-k * t)
% 'a' odpovídá parametru A, 'k' odpovídá parametru k
t0 = 1.6;
modelfun = @(A, k, t0, x) A*(1 - exp(-(x-t0)*k));
% modelfun = fittype(mf, 'independent', 't', 'coefficients', {'A', 'k', 't0'});

% 2. Odhad počátečních hodnot (StartPoints)
% Je důležité trefit se blízko, aby algoritmus konvergoval.
% Předpokládáme A = 1 a k = 0.1 jako startovní výstřel.
start_A = 50;
start_k = 50;

rng = t > t0;
% 3. Provedení fitu
[fitresult, gof] = fit(t(rng) - t0, F(rng), modelfun, 'StartPoint', [start_A, start_k, t0]);

% 4. Získání výsledných parametrů
A_fit(i_fm) = fitresult.A;
k_fit(i_fm) = fitresult.k;
err_fit(i_fm) = gof.rsquare;

% Zobrazení výsledků v konzoli
fprintf('Nalezené parametry:\n A = %.4f\n k = %.4f\n', A_fit(i_fm), k_fit(i_fm));
fprintf('Kvalita fitu (R^2): %.4f\n', err_fit(i_fm));

%% 5. Vykreslení
% clf;
% plot(fitresult);hold on;
% plot(t(rng)-t0, F(rng));
% legend('Naměřená data', 'Fit: 1 - A*exp(-kt)');
% xlabel('t'); ylabel('F');
% grid on;    

end

nexttile(1);legend;
%%

plot(fmaxdata{i_fm}.x_s_, [0;diff(fmaxdata{i_fm}.x_Lo_)]);