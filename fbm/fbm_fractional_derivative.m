% MATLAB script to generate fractional Brownian motion (fBm) and compute its
% Riemann-Liouville fractional derivative numerically (shifted Grünwald-Letnikov approx for higher accuracy using FFT convolution).
% Uses zilany2014_ffgn for fBm generation (included as subfunction).
% Saves plot.png and values.csv (t, D^alpha B_H(t)).

clear; clc; close all;

% Parameters
H = 0.5;          % Hurst parameter (0 < H < 1; D0 = 2 - H = 1.5 for H=0.5)
alpha = 0.2;      % Fractional order 0 < alpha < H
num_points = 1e6 + 1; % Number of time points (adjusted for B(0)=0)
t = linspace(0, 1, num_points);  % Domain [0,1]
dt = t(2) - t(1);

% Generate fBm path using zilany2014_ffgn
Hinput = H;  % For fGn, Hinput = H
noiseType = 1;  % Variable fGn
mu = 0;
sigma = 1;
fGn = zilany2014_ffgn(num_points - 1, dt, Hinput, noiseType, mu, sigma);
B_H = [0, cumsum(fGn) * dt^H];

% Validate fBm: Check increment variance
lags = [10, 100, 1000];
for lag = lags
    diffs = B_H(lag+1:end) - B_H(1:end-lag);
    var_diff = var(diffs);
    expected_var = (lag * dt)^(2*H);
    fprintf('Lag %d: Variance = %.4e, Expected ~%.4e\n', lag, var_diff, expected_var);
end
% Validate power spectrum (expect P(f) ~ f^-(2H+1) at low frequencies)
psd = abs(fft(B_H)).^2 / num_points;
freq = (0:num_points/2-1) / (num_points * dt);
% Fit over low frequencies to avoid high-frequency discretization effects
fit_max_freq = 10000;  % Adjust based on resolution; covers several decades for num_points=1e6
idx = find(freq(2:end) <= fit_max_freq);
log_freq_fit = log(freq(1 + idx));
log_psd_fit = log(psd(1 + idx));
p = polyfit(log_freq_fit, log_psd_fit, 1);
fprintf('Power spectrum slope (low freq): %.2f (Expected ~%.2f for H=%.1f)\n', p(1), -(2*H+1), H);

% Alternative (requires Wavelet Toolbox):
% B_H = wfbm(H, num_points);  % Wavelet-based synthesis
% B_H = B_H / std(B_H);  % Normalize

% Compute Riemann-Liouville fractional derivative (shifted Grünwald-Letnikov approx for higher accuracy using FFT convolution)
Db = zeros(size(t));
if alpha == 0
    Db = B_H;  % Zeroth derivative is the original function
else
    % Increase truncation for better approximation (recursive avoids overflow)
    max_j = min(1e5, num_points-1);  % Adjust as needed for runtime/accuracy trade-off; full is num_points-1
    coeffs = zeros(1, max_j + 1);
    coeffs(1) = 1;
    for j = 1:max_j
        coeffs(j+1) = coeffs(j) * (j - 1 - alpha) / j;
    end
    % FFT-based linear convolution for speed
    pad = num_points + max_j + 1;  % Sufficient padding for linear conv
    padded_B = [0, B_H(1:end), zeros(1, pad - num_points - 1)];  % Pad with 0 at start for shift
    padded_coeffs = [coeffs, zeros(1, pad - length(coeffs))];
    conv_result = ifft(fft(padded_B) .* fft(padded_coeffs));
    Db = real(conv_result(1:num_points)) / dt^alpha;  % Take real to discard numerical imag parts
end
if any(isnan(Db) | isinf(Db))
    error('NaN or Inf detected in derivative Db');
end

% Save values of fractional derivative to CSV (full range)
data = [t(:), Db(:)];
writematrix(data, 'values.csv');

% Plot
f = figure('Visible', 'off');  % Create invisible figure to avoid display issues
plot(t, B_H, 'r-', 'LineWidth', 1.5); hold on;
plot(t, Db, 'b-', 'LineWidth', 1.5);
xlabel('t');
grid on;

% Input for zoom range
zoom_str = input('Enter zoom range as [x_min x_max] or leave empty for full range: ', 's');
if isempty(zoom_str)
    xlim([0 1]);
else
    % Clean and parse without str2num to avoid performance warning
    clean_str = strrep(strrep(zoom_str, '[', ''), ']', '');
    zoom_range = sscanf(clean_str, '%f %f');
    if length(zoom_range) ~= 2 || zoom_range(1) >= zoom_range(2) || zoom_range(1) < 0 || zoom_range(2) > 1
        error('Invalid zoom range. Must be [x_min x_max] with 0 <= x_min < x_max <= 1.');
    end
    xlim(zoom_range);
end

% Set figure size to 2400x1800 pixels
f.Units = 'pixels';
f.Position(3:4) = [2400 1800];

% Save plot with exact pixel size (will capture the zoomed view if set)
exportgraphics(f, 'plot.png', 'Resolution', 300);

disp('Script completed. Plot saved as plot.png (2400x1800 pixels) and values saved as values.csv.');
disp('Note: values.csv contains the full range [0,1]; plot.png reflects the chosen zoom.');

% Subfunction: zilany2014_ffgn (fast fGn generator, modified for exact Davies-Harte)
% Derived from the AMT toolbox implementation (GPLv3), see readme.
function y = zilany2014_ffgn(N, tdres, Hinput, noiseType, mu, sigma)
    if nargin < 6
        sigma = 1;
    end
    if nargin < 5
        mu = 0;
    end
    if nargin < 4
        noiseType = 1;
    end
    if nargin < 3
        error('Requires at least three input arguments.');
    end
    if N <= 0
        error('Length of the return vector must be positive.');
    end
    if tdres > 1
        error('Original sampling rate should be checked.');
    end 
    if (Hinput < 0) || (Hinput > 2)
        error('The Hurst parameter must be in the interval (0,2].')
    end
    eigVal = zeros(1,2*N);  % Eigenvalues of the circulant matrix
    if (noiseType == 0)
        S = zeros(1,2*N);  % Power Spectral Density
        S(1) = 1;
        for k = 2:2*N
            S(k) = (k-1)^(2*Hinput - 2) * (1 - (2*Hinput - 1)/ (k-1) + (2*Hinput - 2)/(2*(k-1)^2));
        end
        eigVal(1) = S(1);
        for k = 2:N+1
            eigVal(k) = S(k) + S(2*N - k + 2);
        end
        for k = N+2:2*N
            eigVal(k) = eigVal(2*N - k + 2);
        end
    else
        % Exact Davies-Harte method for variable fGn
        rho = zeros(1, N+1);
        for j = 0:N
            rho(j+1) = 0.5 * ((j+1).^(2*Hinput) + abs(j-1).^(2*Hinput) - 2 * j.^(2*Hinput));
        end
        acf_ext = [rho(1:N), rho(N+1), rho(N:-1:2)];
        eigVal = real(fft(acf_ext));
        if any(eigVal < 0)
            warning('Negative eigenvalues encountered, setting to zero for approximation.');
            eigVal(eigVal < 0) = 0;
        end
    end
    Z = randn(1,2*N) + 1i * randn(1,2*N);
    Y = fft(Z);
    Y = Y .* sqrt(eigVal / (2*N));
    y = ifft(Y);
    y = real(y(1:N));
    y = sigma * y / std(y);
    y = mu + y;
    if (Hinput > 1)
        y = cumsum(y) * tdres^(Hinput - 1);
    end
end