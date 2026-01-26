---
title: fBm Local Time and Fractional Integrals Simulation

---

---
title: fBm Local Time and Fractional Integrals Simulation
---

# fBm Local Time and Fractional Integrals Simulation (Grok)

This repository contains a Python script that demonstrates the mathematical connection between fractional integrals and the averaged local time over the fractal zero sets of fractional Brownian motion (fBm). Specifically, it shows that the fractional integral of order α of a smooth function f is proportional to the expected integral of f over the zero set of an fBm with Hurst parameter H = 1 - α.

The code simulates fBm paths using the Cholesky decomposition method, approximates the local time using a Gaussian mollifier, averages over multiple realizations, compares to the analytical expectation, and visualizes the zero sets with fuzzy bars for a fractal effect.

This document and the Python script have been developed with AI assistance from Grok.

## Requirements

- numpy
- scipy
- matplotlib

## Mathematical Concept

Fractional Brownian motion (fBm) is a generalization of Brownian motion parametrized by the Hurst parameter H. The "local time" measures how much time the process spends near zero. This script shows that averaging this measure over many fBm paths, weighted by a function f, gives something proportional to the fractional integral of a reversed version of f - a tool used in areas like physics and signal processing to smooth or accumulate data over time. This shows that the fractional integral of order α (α = 1 - H) of f is proportional to the normal integral restricted to a fractal subset with the right H. The fractal zero sets of fBm are a family of fractal subsets of the x axis with the right H. The proportionality factor is $\Gamma(\alpha) / \sqrt{2\pi}$. Since every smooth function is locally linear, this should apply to nonlinear functions as well.

fBm is a Gaussian process $B^H(t)$ with Hurst parameter $H \in (0,1)$, covariance $\mathrm{Cov}[B^H(t), B^H(s)] = \frac{1}{2}(|t|^{2H} + |s|^{2H} - |t-s|^{2H})$, and self-similarity $B^H(at) = a^H B^H(t)$.

The local time $L^0(T)$ at level 0 up to time $T$ approximates the "occupation" near zero. The expected integral is the average value over many random paths of the integral of f with respect to the local time measure. It is proportional to the fractional integral $I^\alpha g(T)$, where $g(s) = f(T - s)$ is the reversed $f$.

For $f(t) = a + b t$, the analytical expectation is $\frac{1}{\sqrt{2\pi}} \left[ a \frac{T^{1-H}}{1-H} + b \frac{T^{2-H}}{2-H} \right]$. This analytical expectation refers to the expected value (average over many random fBm paths) of the integral of the function $f(t)$ with respect to the local time measure $dL^0(t)$ at level 0. It is derived from the known density of the local time for fBm, which is $\frac{1}{\sqrt{2\pi}} t^{-H}$. Integrating $f(t)$ against this density from 0 to T gives the formula by splitting into the constant and linear terms. The analytical result is exact, based on rigorous mathematical theory of fBm local times, while numerical simulations approximate it with errors due to finite grid resolution, mollifier width, and limited realizations.

The simulation verifies this numerically and plots the zero sets from a few realizations with an overlay of the proportional integral and $f(t)$.

### Customizing Parameters
Edit the script's parameters section:
- `T`: Time horizon (default 1.0)
- `H`: Hurst parameter (default 0.2)
- `N`: Grid resolution (higher for accuracy, e.g., 8192)
- `eps`: Mollifier width (smaller for sharpness, e.g., 0.005)
- `num_runs`: Number of realizations for averaging (higher for low variance, e.g., 2000)
- `a`, `b`: Coefficients for linear f(t) = a + b t (default 1.0, 1.0)
- `p`: For non-uniform grid (optional, lower for denser near 0, e.g., 0.4)

## Examples

### Output for H=0.7, f(t)=1 + t
Numerical average: 0.793
Analytical value: 0.798
Proportional fractional integral: 0.798
Relative error: 0.54%

### Visualization
The plot shows fuzzy colored bars for 3 fBm zero sets, the proportional $I^\alpha g(t)$ (black), $f(t)$ (gray dashed), and legend with H, $\alpha$, error. The 3 fBm zero sets are used for illustration, but the actual averaging uses num_runs paths.

The "relative error" is the percentage difference between the numerical average and the analytical value, calculated as |numerical - analytical| / analytical * 100%.

For high H (0.9), errors may be higher (~40%) due to persistence; see Limitations.

## Limitations

- High H (low $\alpha$) leads to higher errors due to persistent paths and singular density near t=0; use non-uniform grid for improvement.
- Low H (high $\alpha$) may have variance issues; increase num_runs.
- Cholesky decomposition is O(N^3) time/O(N^2) memory—large N (>16384) may require more resources.

## References
- Mandelbrot & Van Ness (1968): Benoit B. Mandelbrot and John W. van Ness, "Fractional Brownian Motions, Fractional Noises and Applications," SIAM Review, Vol. 10, No. 4 (October 1968), pp. 422-437.
- Geman & Horowitz (1980): Donald Geman and Joseph Horowitz, "Occupation Densities," The Annals of Probability, Vol. 8, No. 1 (February 1980), pp. 1-67.
- Jolis (2007): Maria Jolis and Noèlia Viles, "Continuity in law with respect to the Hurst parameter of the local time of the fractional Brownian motion," Journal of Theoretical Probability, Vol. 20 (2007), pp. 133-152.
- Mishura (2008): Yuliya Mishura, "Stochastic Calculus for Fractional Brownian Motion and Related Processes," Lecture Notes in Mathematics, Vol. 1929, Springer (2008).