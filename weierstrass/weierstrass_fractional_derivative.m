% MATLAB script to compute the Riemann-Liouville fractional derivative
% of the Weierstrass function using the analytic term-by-term expression
% with phase shift. Approximates with finite sum.

clear; clc; close all;

% Parameters for Weierstrass function
a = 0.5;          % 0 < a < 1
b = 3;            % odd integer, ab > 1 + 3*pi/2
N = 20;           % number of terms for approximation
alpha = 0.2;      % fractional order 0 < alpha < 1

% Domain
x = linspace(0, 2, 500000);

% Compute Weierstrass function W(x)
W = zeros(size(x));
for n = 0:N-1
    W = W + a^n * cos(b^n * pi * x);
end

% Compute fractional derivative D^alpha W(x)
Dw = zeros(size(x));
for n = 0:N-1
    freq = b^n * pi;
    Dw = Dw + a^n * freq^alpha * cos(freq * x + alpha * pi / 2);
end

% Plot
f = figure;
plot(x, W, 'r-', 'LineWidth', 1.5); hold on;
plot(x, Dw, 'b-', 'LineWidth', 1.5);
xlabel('x');
grid on;

% Set figure size to 2400x1800 pixels
f.Position(3:4) = [2400 1800];

% Save plot with exact pixel size using exportgraphics (requires MATLAB R2020a or later)
exportgraphics(f, 'plot.png', 'Units', 'pixels', 'Width', 2400, 'Height', 1800);

% Save values of fractional derivative to CSV
data = [x(:), Dw(:)];
writematrix(data, 'values.csv');

disp('Script completed. Plot saved as plot.png (2400x1800 pixels) and values saved as values.csv.');