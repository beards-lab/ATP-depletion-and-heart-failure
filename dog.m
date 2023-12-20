% get degrees of freedom

% Generate a time series (replace this with your own time series data)
data = randn(100, 1);
data = tab_rmpAvg.F;
% clf;
figure(22);
% Set the maximum lag for ACF
ai = 2:200;
for i = ai
    maxlag = i;
    
    % Calculate autocorrelation function (ACF)
    lags = 0:maxlag;
    acf = autocorr(data, maxlag);
    
    % Effective Sample Size Calculation
    neff = length(data) / (1 + 2 * sum(acf));
    
    % Number of Parameters in Your Model (replace this with the actual number)
    num_parameters = 19;
    
    % Degrees of Freedom Calculation
    df = neff - num_parameters;
    
    % Display Results
    fprintf('%d: Effective Sample Size: %.2f\n', i, neff);
    fprintf('%d: Degrees of Freedom: %.2f\n', i, df);
    aneff(i) = neff;
    adf(i) = df;
end

plot(1:ai(end), aneff);hold on;
figure(4)

% Function to calculate autocorrelation
function acf = autocorr(x, maxlag)
    n = length(x);
    acf = zeros(1, maxlag + 1);
    
    % Calculate sample mean and sample variance
    xbar = mean(x);
    Sxx = sum((x - xbar).^2);
    
    % Calculate autocorrelation for each lag
    for k = 0:maxlag
        gamma_k = sum((x(1:end-k) - xbar) .* (x(1+k:end) - xbar));
        acf(k+1) = gamma_k / Sxx;
    end
end