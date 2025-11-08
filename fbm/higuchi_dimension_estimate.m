% MATLAB script to estimate the fractal dimension using Higuchi's method
% on the fractional derivative D^alpha W(x) from values.csv.
% Corrected formula: D = 1 - slope; refined normalization per Higuchi (1988).

clear; clc; close all;

% Load data from values.csv (columns: x, Dw); use Dw for derivative
data = readmatrix('values.csv');
x = data(:, 1);
y = data(:, 2);  % y = Dw (or W for alpha=0)

% Ensure sorted (should be)
[x, sortIdx] = sort(x);
y = y(sortIdx);

N = length(y);
kmax = min(200, floor(N/4));  % Max lag: ~N/4 for stability, cap at 200
Lk = zeros(1, kmax);
for k = 1:kmax
    Lmk = zeros(1, k);
    for m = 1:k
        idx_start = m;
        idx = idx_start:k:N;
        if length(idx) > 1
            diffs = abs(y(idx(2:end)) - y(idx(1:end-1)));
            M = length(idx);  % Subseries length
            Lm = (sum(diffs) * (N - 1)) / (k * (M - 1));  % Standard Higuchi normalization
            Lmk(m) = Lm;
        end
    end
    Lk(k) = mean(Lmk(Lmk > 0));  % Average positive lengths
end

% Filter valid points
valid = (Lk > 0) & isfinite(Lk);
k_valid = (1:kmax)';
logk = log(k_valid(valid));
logLk = log(Lk(valid));
if length(logk) < 10
    error('Too few valid points; increase N or reduce kmax.');
end

% Linear fit: slope s = -(D - 1)
p = polyfit(logk, logLk, 1);
s = p(1);
D_higuchi = 1 - s;  % Corrected: D = 1 - s
residuals = logLk - polyval(p, logk);
SS_res = sum(residuals.^2);
SS_tot = sum((logLk - mean(logLk)).^2);
R2 = 1 - (SS_res / SS_tot);

% Display results
fprintf('Higuchi fractal dimension: %.4f\n', D_higuchi);
fprintf('R-squared of fit: %.4f\n', R2);
fprintf('Used %d valid lags out of %d\n', sum(valid), kmax);

% Plot log-log for visualization
f = figure;
semilogx(k_valid, Lk, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;  % semilogx for log k
semilogx(k_valid(valid), Lk(valid), 'go', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', 'g');
fitk = logspace(log10(min(k_valid(valid))), log10(max(k_valid(valid))), 100);
plot(fitk, polyval(p, log(fitk)), 'g-', 'LineWidth', 3);
xlabel('k (lag)');
ylabel('L(k) (curve length)');
title('Higuchi Dimension Estimate');
legend('Full Data', 'Fitted Lags', sprintf('Fit (D = %.3f)', D_higuchi), 'Location', 'best');
grid on;

% Set figure size to 2400x1800 pixels
f.Position(3:4) = [2400 1800];

% Save plot
exportgraphics(f, 'higuchi_fit.png', 'Units', 'pixels', 'Width', 2400, 'Height', 1800);

disp('Script completed. Dimension printed above; plot saved as higuchi_fit.png.');