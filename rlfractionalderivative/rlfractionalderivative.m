% MATLAB script to read function values from CSV, compute left or right Riemann-Liouville fractional derivative, and plot it
% Improvements for numerical accuracy:
% - Recursive computation of binomial coefficients to avoid gamma function overflow for large N
% - Option for shifted Grünwald-Letnikov approximation for higher-order accuracy (O(h^2) instead of O(h))
% - Boundary handling assumes f=0 outside the interval
% - Higher-order central difference for first derivative to improve precision

% Read data from CSV
data = readmatrix('values.csv');
x = data(:,1);
y = data(:,2);

% Assume uniform spacing
h = x(2) - x(1);
N = length(x);

% Compute first derivative using 4th-order central difference for higher accuracy
dy = zeros(size(y));
if N >= 5
    % Forward differences at start (O(h^2))
    dy(1) = (-3*y(1) + 4*y(2) - y(3)) / (2*h);
    dy(2) = (-y(1) + y(3)) / (2*h);
    
    % 4th-order central for interior
    for i = 3:N-2
        dy(i) = (-y(i+2) + 8*y(i+1) - 8*y(i-1) + y(i-2)) / (12*h);
    end
    
    % Backward differences at end (O(h^2))
    dy(N-1) = (y(N-2) - y(N)) / (2*h);
    dy(N) = (y(N-2) - 4*y(N-1) + 3*y(N)) / (2*h);
else
    % Fallback to gradient if too few points
    dy = gradient(y, x);
end

% Input fractional order
alpha = input('Enter fractional order alpha (0<alpha<1): ');

% Input choice
choice = input('Enter "left" or "right" for the derivative: ', 's');

% Option for shifted (higher-order) approximation
use_shifted = input('Use shifted Grünwald-Letnikov for higher order accuracy (yes/no): ', 's');
use_shifted = strcmp(use_shifted, 'yes');

% Compute coefficients recursively (length N+1 to handle shifted case)
binom = zeros(1, N+1);
binom(1) = 1;
for j = 1:N
    binom(j+1) = binom(j) * (alpha - j + 1) / j;
end
coeff = (-1).^(0:N) .* binom;

% Compute the selected derivative
der = zeros(N,1);
if strcmp(choice, 'left')
    for i = 1:N
        s = 0;
        if use_shifted
            max_k = i;  % Sum up to k = i for shifted
            for k = 0:max_k
                pos = i - k + 1;
                if pos < 1 || pos > N
                    fval = 0;
                else
                    fval = y(pos);
                end
                s = s + coeff(k+1) * fval;
            end
        else
            max_k = i-1;  % Standard sum up to k = i-1
            for k = 0:max_k
                pos = i - k;
                fval = y(pos);
                s = s + coeff(k+1) * fval;
            end
        end
        der(i) = s / h^alpha;
    end
elseif strcmp(choice, 'right')
    for i = 1:N
        s = 0;
        if use_shifted
            max_k = N - i + 1;  % Sum up to k = N-i+1 for shifted
            for k = 0:max_k
                pos = i + k - 1;
                if pos < 1 || pos > N
                    fval = 0;
                else
                    fval = y(pos);
                end
                s = s + coeff(k+1) * fval;
            end
        else
            max_k = N - i;  % Standard sum up to k = N-i
            for k = 0:max_k
                pos = i + k;
                fval = y(pos);
                s = s + coeff(k+1) * fval;
            end
        end
        der(i) = s / h^alpha;
    end
else
    error('Invalid choice. Must be "left" or "right".');
end

% Save to CSV: x, der
writematrix([x, der], 'values.csv');

% Input for zoom range
zoom_str = input('Enter zoom range as [min_x max_x] or leave empty for full range: ', 's');
if isempty(zoom_str)
    zoom_min = min(x);
    zoom_max = max(x);
else
    zoom_range = str2num(zoom_str);
    if length(zoom_range) ~= 2 || zoom_range(1) >= zoom_range(2)
        error('Invalid zoom range. Must be [min_x max_x] with min_x < max_x.');
    end
    zoom_min = zoom_range(1);
    zoom_max = zoom_range(2);
end

% Desired resolution
resolution = 300;  % DPI

% Create figure with specified size in inches (2400x1800 pixels at 300 DPI -> 8x6 inches)
fig = figure('Units', 'inches', 'Position', [0 0 8 6]);

% Plot the original (green solid), first derivative (red solid), and fractional derivative (blue solid) without labels
hold on;
plot(x, y, 'g-');
plot(x, dy, 'r-');
plot(x, der, 'b-');
hold off;
grid on;
xlim([zoom_min, zoom_max]);

% Set paper properties
set(fig, 'PaperUnits', 'inches', 'PaperSize', [8 6]);
set(fig, 'PaperPosition', [0 0 8 6]);
set(fig, 'PaperPositionMode', 'manual');

% Save the image with specified resolution
print(fig, 'plot.png', '-dpng', ['-r' num2str(resolution)]);