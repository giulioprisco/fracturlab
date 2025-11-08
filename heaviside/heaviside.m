% MATLAB script to illustrate how the fractional derivative "rotates" the tangent to the Heaviside function,
% depending on past history (x0 and A).
% Plots Heaviside functions for 4 different A between x_min and x_max.
% At 2 different values of x > x0, displays the fractional derivative as rotated tangent lines on each function.
% Uses analytic formula for Riemann-Liouville fractional derivative: A * (x - x0)^(-alpha) / gamma(1 - alpha)
% Saves a high-resolution plot of size 2400 x 1800 pixels.

% Parameters
alpha = 0.9;
x0 = 0;
gamma_val = gamma(1 - alpha);
A_vals = [1, 2, 3, 4];  % 4 different values of A (varying history amplitude)
x_tang = [1, 1.5];  % 2 different values of x > x0 (varying distance from history)
x_min = 1;  % Plot between x_min and x_max
x_max = 2;
x_plot = linspace(x_min, x_max, 1000);
colors = {'green', 'red', 'magenta', 'blue'};  % Colors for tangents: green (A=1), red (A=2), magenta (A=3), blue (A=4)
arrow_len = 0.3;  % Fixed line length for visualization

% Create figure and set size for export (8x6 inches for 2400x1800 at 300 dpi)
fig = figure;
set(fig, 'Units', 'inches');
set(fig, 'Position', [0 0 8 6]);
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPosition', [0 0 8 6]);

hold on;

for i = 1:length(A_vals)
    A = A_vals(i);
    col = colors{i};
    
    % Compute Heaviside function (constant A for x >= x0 in the plot range)
    f = A * ones(size(x_plot));
    
    % Plot the function in black (representing the standard horizontal tangent)
    plot(x_plot, f, 'Color', 'k', 'LineWidth', 2);
    
    % Display fractional derivative as rotated tangent lines at specified x
    for j = 1:length(x_tang)
        x = x_tang(j);
        y = A;  % f(x) = A for x >= x0
        
        % Analytic fractional derivative (generalized slope, showing rotation from standard 0 slope)
        if x > x0
            deriv = A * (x - x0)^(-alpha) / gamma_val;
        else
            deriv = 0;
        end
        
        % Direction vector (1, deriv) normalized and scaled
        vec = [1, deriv];
        nrm = norm(vec);
        if nrm == 0
            u = arrow_len;
            v = 0;
        else
            unit = vec / nrm;
            u = unit(1) * arrow_len;
            v = unit(2) * arrow_len;
        end
        
        % Plot the rotated tangent line segment (no arrow)
        plot([x, x + u], [y, y + v], 'Color', col, 'LineWidth', 2);
    end
end

% Labels (no ylabel)
xlabel('x');
grid on;
hold off;

% Save the plot as high-resolution PNG (2400 x 1800 pixels at 300 dpi)
exportgraphics(fig, 'plot.png', 'Resolution', 300);