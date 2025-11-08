% Complete MATLAB program for computing RL fractional integral, constructing fractal region,
% approximating Hausdorff content via self-similar measure integral, verifying scaling,
% normalizing to match fractional integral, and plotting/saving results.

% Note: Removed 'clear all;' to avoid performance warning; script assumes clean start or manual clear if needed.
close all;
clc;

% Parameters
alpha = 0.85;          % Fractional order (0 < alpha < 1)
lambda = 2;           % Scaling factor for verification
x_max = 1;            % Domain [0, x_max]
N_samples = 1e6;      % Number of samples for fractal approximation
iterations = 30;      % Max depth for sampling (to avoid infinite loop)
slope = 1;            % Slope b for linear f(x) = b*x
f = @(x) slope * x;   % Linear function

% Step 1: Compute RL fractional integral analytically for linear f(x) = b*x
gamma_alpha2 = gamma(alpha + 2);
I_alpha_original = slope * x_max^(alpha + 1) / gamma_alpha2;
I_alpha_scaled = slope * (lambda * x_max)^(alpha + 1) / gamma_alpha2;
scaling_ratio_I = I_alpha_scaled / I_alpha_original;  % Should be lambda^{1+alpha}

% Step 2 & 3: Construct generalized Cantor set with dim alpha, approximate "fractal content"
% as integral_C f(t) d mu(t) ≈ average f(t_samples) where t ~ mu (self-similar measure)
r = 2^(-1 / alpha);  % Scaling ratio for Cantor branches

% Function to sample from Cantor measure on [0, scale]
function t = sample_cantor(r, scale, iterations)
    t = 0;
    current_scale = scale;
    for d = 1:iterations
        if rand < 0.5
            % Left branch: add 0
        else
            % Right branch: add (1 - r) * current_scale
            t = t + (1 - r) * current_scale;
        end
        current_scale = r * current_scale;
    end
end

% Sample for original (scale=1)
t_samples_orig = zeros(N_samples, 1);
for i = 1:N_samples
    t_samples_orig(i) = sample_cantor(r, x_max, iterations);
end
mean_f_orig = mean(f(t_samples_orig));
fractal_content_orig = x_max^alpha * mean_f_orig;

% Sample for scaled (scale=lambda, f_scaled(y) = lambda * f(y / lambda))
t_samples_scaled = zeros(N_samples, 1);
g = @(y) lambda * f(y / lambda);
for i = 1:N_samples
    t_samples_scaled(i) = sample_cantor(r, lambda * x_max, iterations);
end
mean_f_scaled = mean(g(t_samples_scaled));
fractal_content_scaled = (lambda * x_max)^alpha * mean_f_scaled;
scaling_ratio_fractal = fractal_content_scaled / fractal_content_orig;  % Should ≈ lambda^{1+alpha}

% Step 3.5: Normalization to match RL integral (principled via scaling alignment, approximated as constant factor for this f)
factor = I_alpha_original / fractal_content_orig;
normalized_content_orig = factor * fractal_content_orig;
normalized_content_scaled = factor * fractal_content_scaled;

% Step 4: Verify scaling (both should be lambda^{1+alpha})
expected_scaling = lambda^(1 + alpha);
disp('--- Console Diagnostics ---');
disp(['Alpha: ', num2str(alpha)]);
disp(['Lambda: ', num2str(lambda)]);
disp(['RL Integral at x=1: ', num2str(I_alpha_original)]);
disp(['RL Integral at x=lambda: ', num2str(I_alpha_scaled)]);
disp(['RL Scaling Ratio: ', num2str(scaling_ratio_I), ' (expected: ', num2str(expected_scaling), ')']);
disp(['Fractal Content Original: ', num2str(fractal_content_orig)]);
disp(['Fractal Content Scaled: ', num2str(fractal_content_scaled)]);
disp(['Fractal Scaling Ratio: ', num2str(scaling_ratio_fractal), ' (expected: ', num2str(expected_scaling), ')']);
disp(['Normalization Factor: ', num2str(factor)]);
disp(['Normalized Fractal Original: ', num2str(normalized_content_orig), ' (matches RL: ', num2str(I_alpha_original), ')']);
disp(['Normalized Fractal Scaled: ', num2str(normalized_content_scaled), ' (matches RL: ', num2str(I_alpha_scaled), ')']);

% Step 5: Merged single plot
% Generate grid for f(x) and I^alpha
x_grid = linspace(0, x_max, 1000);
f_grid = f(x_grid);
I_alpha_grid = slope * x_grid.^(alpha + 1) / gamma_alpha2;

% Create single figure with merged elements
figure('Position', [100, 100, 2400, 1800], 'PaperUnits', 'inches', 'PaperSize', [8 6], 'PaperPosition', [0 0 8 6]);
hold on;

% Plot Cantor bar set (fuzzy shaded lines for fractal suggestion)
subsample = 5000;  % Larger subsample for more bars/fractal look
idx = randperm(N_samples, subsample);
t_sub = t_samples_orig(idx);
for k = 1:subsample
    tt = t_sub(k);
    plot([tt tt], [0 f(tt)], 'b-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5 0.05], 'HandleVisibility', 'off');  % Gray color, low-alpha for fuzzy shading when dense
end

% Dummy line for Cantor in legend
plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Cantor');

% Overlay function f(x)
f_legend = ['f(x) = ' num2str(slope) 'x'];
plot(x_grid, f_grid, 'b-', 'LineWidth', 2, 'DisplayName', f_legend);

% Overlay fractional integral on secondary y-axis (since scales differ)
yyaxis right;
plot(x_grid, I_alpha_grid, 'r-', 'LineWidth', 2, 'DisplayName', ['I^{\alpha} f(x), \alpha=' num2str(alpha)]);
ax = gca;
ax.YAxis(2).Color = 'r';  % Explicitly set right y-axis color to red
ax.YAxis(1).Color = 'b';  % Set left y-axis color to blue for consistency
ax.Toolbar.Visible = 'off';  % Disable axes toolbar to avoid export warning

% Finalize plot
yyaxis left;
title('');
legend('Location', 'northwest');
xlim([0, x_max]);
ylim([0, max(f_grid)]);
grid on;
hold off;

% Step 6: Save high-resolution plot
print('-dpng', '-r300', 'plot.png');
disp('High-resolution plot saved as plot.png');