% Matlab script to compute and plot a function% Input: number of points
n = input('Enter the number of equally spaced points n: ');% Define the function here (easy to change)
f = @(x) x.^2;% Generate x values
x = linspace(0, 5, n);% Compute y
y = f(x);% Save the values to a CSV file
data = [x', y'];  % Transpose to column vectors
writematrix(data,'values.csv');% Desired resolution
resolution = 300;  % DPI% Create figure with specified size in inches (2400x1800 pixels at 300 DPI -> 8x6 inches)
fig = figure('Units', 'inches', 'Position', [0 0 8 6]);% Plot the function without labels
plot(x, y);
grid on;% Set paper properties
set(fig, 'PaperUnits', 'inches', 'PaperSize', [8 6]);
set(fig, 'PaperPosition', [0 0 8 6]);
set(fig, 'PaperPositionMode', 'manual');% Save the image with specified resolution
print(fig, 'plot.png', '-dpng', ['-r' num2str(resolution)]);