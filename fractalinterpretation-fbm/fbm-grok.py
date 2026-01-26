import numpy as np
from scipy.linalg import cholesky
from math import sqrt, pi, gamma
import matplotlib.pyplot as plt

def approximate_local_time_integral(t, path, f, eps):
    dt = np.diff(t)  # Support non-uniform grid
    dt = np.append(dt, dt[-1])  # Approximate last dt
    indicator = (1 / (eps * sqrt(2 * pi))) * np.exp(-path**2 / (2 * eps**2))
    integral = np.sum(f(t) * indicator * dt)
    return integral

def average_over_runs(num_runs, t, L, eps, f, base_seed=42):
    total = 0.0
    N = len(t)
    np.random.seed(base_seed)
    for run in range(num_runs):
        z = np.random.normal(size=N)
        path = L @ z
        integ = approximate_local_time_integral(t, path, f, eps)
        total += integ
        if (run + 1) % (num_runs // 10) == 0:  # Simple console progress every 10%
            print(f"Progress: {(run + 1) / num_runs * 100:.0f}%")
    return total / num_runs

# Parameters
T = 1.0
H = 0.15
alpha = 1 - H
N = 16384  # Increased for better accuracy
eps = 0.0025  # Decreased for sharper approximation
num_runs = 4096  # Increased for lower variance
a = 1.0
b = 1.0
f = lambda t: a + b * t
p = 0.4  # For non-uniform grid denser near 0 (optional, comment out for uniform)

# Generate t (non-uniform if p defined, else uniform)
u = np.linspace(0, 1, N)
t = T * u**(1/p)  # Non-uniform; comment this and uncomment next for uniform
# t = np.linspace(0, T, N)  # Uniform grid

# Covariance on grid
cov = 0.5 * (np.power(np.abs(t[:, None]), 2*H) + np.power(np.abs(t[None, :]), 2*H) - np.power(np.abs(t[:, None] - t[None, :]), 2*H))
L = cholesky(cov + 1e-15 * np.eye(N), lower=True)

# Numerical average
avg_numerical = average_over_runs(num_runs, t, L, eps, f)

# Analytical expectation
analytical = (1 / sqrt(2 * pi)) * (a * (T ** (1 - H) / (1 - H)) + b * (T ** (2 - H) / (2 - H)))

# Fixed I_alpha_g
I_alpha_g = (a + b * T) * T**alpha / (gamma(alpha) * alpha) - b * T**(alpha + 1) / (gamma(alpha) * (alpha + 1))

# Full proportional expression
proportional_fractional = gamma(alpha) / sqrt(2 * pi) * I_alpha_g

relative_error = abs(avg_numerical - analytical) / analytical * 100
print(f"Numerical average: {avg_numerical}")
print(f"Analytical value: {analytical}")
print(f"Proportional fractional integral: {proportional_fractional}")  # Should match analytical
print(f"Relative error (numerical vs analytical): {relative_error:.2f}%")

# Plotting code
num_paths = 3  # Reduced to 3 realizations as requested
colors = ['blue', 'green', 'red']  # Three distinct colors

bar_width = np.median(np.diff(t)) * 1.0  # Adjusted for non-uniform

# Generate paths and zeros
paths = []
zero_times_list = []  # List of zero times per path
np.random.seed(42)  # Reset seed for reproducible plots
for i in range(num_paths):
    z = np.random.normal(size=N)
    path = L @ z
    paths.append(path)
    
    # Approximate zero crossings via interpolation
    sign_changes = np.diff(np.sign(path)) != 0
    zero_indices = np.where(sign_changes)[0]
    zero_times = []
    for idx in zero_indices:
        t1, t2 = t[idx], t[idx+1]
        p1, p2 = path[idx], path[idx+1]
        zero_t = t1 - p1 * (t2 - t1) / (p2 - p1)
        zero_times.append(zero_t)
    zero_times_list.append(np.array(zero_times))

# Plot in repo style with increased size and resolution
fig, ax = plt.subplots(figsize=(16, 10))  # Larger figure size (16x10 instead of 10x6)

# For Gaussian blur simulation
num_sub_bars = 31  # Increased for smoother blur
sigma = bar_width * 2  # Doubled sigma for stronger blur with higher N
offsets = np.linspace(-4 * sigma, 4 * sigma, num_sub_bars)  # Wider range for more blur
alpha_scale = np.exp(-offsets**2 / (2 * sigma**2))
alpha_scale /= alpha_scale.max()  # Normalize to 1 at center
sub_bar_width = bar_width / (num_sub_bars / 8)  # Adjusted for more overlap and fuzziness at high N

for i in range(num_paths):
    zero_times = zero_times_list[i]
    for zt in zero_times:
        for j, off in enumerate(offsets):
            ax.bar(zt + off, f(zt), width=sub_bar_width, color=colors[i], alpha=0.08 * alpha_scale[j], align='center')  # Lower base alpha for softer, fuzzier effect

# Overlay proportional fractional integral (for general linear f, running analytical expectation)
tt = np.linspace(0, T, 100)
I_alpha = (1 / sqrt(2 * pi)) * (a * (tt ** (1 - H) / (1 - H)) + b * (tt ** (2 - H) / (2 - H)))  # Running integral for overlay
ax.plot(tt, I_alpha, color='black', linewidth=3, label='Proportional I^α g(t)')  # Thicker line for visibility

# Overlay the function f(t) = a + b t
ax.plot(t, f(t), color='gray', linestyle='--', linewidth=2, label='f(t) = {:.1f} + {:.1f} t'.format(a, b))

# Add dummy plots for H, alpha, and relative error in legend
ax.plot([], [], ' ', label='H = {:.2f}'.format(H))
ax.plot([], [], ' ', label='α = {:.2f}'.format(alpha))
ax.plot([], [], ' ', label='Relative Error = {:.2f}%'.format(relative_error))

ax.set_xlabel('Time t', fontsize=14)
ax.set_ylabel('Height', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, linestyle='--', alpha=0.5)  # Added light grid for better resolution feel

# For higher resolution when saving (optional, for local high-res output)
# plt.savefig('fbm_plot_high_res.png', dpi=300)  # Uncomment to save with 300 DPI

plt.show()  # Uncomment to display the plot when running locally
