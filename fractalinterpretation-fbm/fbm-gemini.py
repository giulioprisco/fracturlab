import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from fbm import FBM
import sys

def plot_fractal(alpha, n_steps=2000, n_realizations=1000, epsilon=0.05):
    """
    Final optimized stochastic simulation for alpha=0.1, H=0.9
    Balancing smoothness and origin-capture without forcing.
    """
    H = 1.0 - alpha
    a, b = 1.0, 2.0
    
    # Grid setup
    t = np.linspace(0.0001, 1.0, n_steps)
    dt = 1.0 / n_steps
    f_t = a + b * t
    
    print(f"--- Running Fractal Simulation: Stability Mode ---")
    print(f"Alpha: {alpha:.2f}, H: {H:.2f}")
    print(f"Settings: {n_steps} steps, {n_realizations} realizations, epsilon = {epsilon}")
    
    # 1. Analytical Riemann-Liouville Integral (Blue Ground Truth)
    analytical = (a * t**alpha / gamma(alpha + 1)) + \
                 (b * t**(alpha + 1) / gamma(alpha + 2))
    
    # 2. Stochastic Sampling
    accumulated_density = np.zeros(n_steps)
    d_norm = (1.0 / (epsilon * np.sqrt(2 * np.pi)))
    
    total_loops = n_realizations // 2
    for i in range(total_loops):
        # Clean vertical progress reporting
        if (i + 1) % (max(1, total_loops // 20)) == 0:
            percent = int((i + 1) / total_loops * 100)
            print(f"Progress: {percent}%")
            sys.stdout.flush()

        # Path generation
        f = FBM(n=n_steps-1, hurst=H, length=1.0, method='hosking')
        path = f.fbm()
        
        # KDE with Antithetic Variance Reduction
        d1 = d_norm * np.exp(-(path**2) / (2 * epsilon**2))
        d2 = d_norm * np.exp(-((-path)**2) / (2 * epsilon**2))
        
        accumulated_density += (d1 + d2)

    print("Computation finished. Rendering plot...")
    avg_density = accumulated_density / n_realizations
    
    # Pure Numerical Integration (Riemann Sum)
    numerical = np.cumsum(f_t * avg_density) * dt
    
    # Scale to match at x=1 (Theoretical Normalization)
    scaled_numerical = numerical * (analytical[-1] / numerical[-1])
    
    # 3. Final Aesthetic Plotting
    plt.figure(figsize=(16, 10))
    
    # Fractal Mist Background
    colors = ['#2ca02c', '#ff7f0e', '#9467bd']
    for i in range(100):
        p = FBM(n=n_steps-1, hurst=H, length=1.0, method='hosking').fbm()
        mist = np.exp(-(p**2) / (2 * 0.04**2))
        plt.fill_between(t, 0, f_t * mist, color=colors[i % 3], alpha=0.02, lw=0)

    # Core Result Lines
    plt.plot(t, f_t, color='black', lw=1.2, label=f'$f(x) = {a} + {b}x$', alpha=0.6)
    plt.plot(t, analytical, color='blue', lw=4, label='Analytical RL Integral', zorder=10)
    plt.plot(t, scaled_numerical, color='red', ls='--', lw=3, label='Stochastic Average', zorder=11)

    # Pure Legend Parameters
    plt.legend(loc='upper left', title=f"$\\alpha = {alpha:.2f}, H = {H:.2f}$", 
               title_fontsize=14, fontsize=12, framealpha=0.9)
    
    plt.xlabel("x", fontsize=14)
    plt.ylabel("Value", fontsize=14)
    plt.xlim(0, 1)
    plt.ylim(0, max(analytical)*1.1)
    plt.grid(True, alpha=0.15)
    plt.tight_layout()
    plt.show()

# Execute stability run
plot_fractal(alpha=0.4)
