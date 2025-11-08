# README: Fractional Integral and Fractal Analogy MATLAB Script

## Summary

The value of a smooth function f(x), where both x and f(x) are lengths, has dimension 1, intuitively and formally because locally (in the vicinity of a point) it scales linearly when the scale is changed. The integral of the function (the area under the function) has dimension 2, intuitively and formally because the local contribution to the integral scales quadratically when the scale is changed. Changing the scale means zooming (expanding or shrinking) by a factor lambda. This multiplies lengths by lambda (that is, by lambda^1), areas by lambda^2, and volumes by lambda^3, formalizing the intuitive meaning of integer dimensions (1, 2, 3...).

The fractional integral of order alpha of the function scales as lambda^(1 + alpha), reflecting a blend of local linear (dimension 1) and accumulative quadratic (dimension 2) behaviors. That is, changing the scale by a factor lambda multiplies the value of the fractional integral by lambda^(1 + alpha). This recovers dimension 1 when alpha tends to zero and dimension 2 when alpha tends to 1. But what about intermediate values of alpha?

There is a similar scaling behavior in fractals: the size of a fractal, or more precisely its Hausdorff content (a simpler version of the Hausdorff measure), scales as lambda^D, where D is the fractal (Hausdorff) dimension. This suggests interpreting the fractional integral as analogous to the size (Hausdorff content) of a fractal with dimension D = 1 + alpha.

This Matlab script starts with a simple linear function f(x) = bx and derives its Riemann-Liouville fractional integral. Then the script constructs a fractal with dimension D = 1 + alpha: the region under the function restricted to a Cantor set-like fractal with dimension alpha on the x-axis, Then the script confirms the matching scaling behaviors under zooming in/out, and shows that with appropriate normalization the size (D-dimensional Hausdorff content) of the fractal is equal to the fractional integral.

This suggests interpreting the fractional integral of order alpha of a function f(x) as the integral of f over a fractal subset of the x axis with fractal dimension alpha. This interpretation is suggested in the literature.

## Overview

This MATLAB script demonstrates the analogy between the Riemann-Liouville (RL) fractional integral of order alpha and the Hausdorff content of a fractal set with dimension D = 1 + alpha. It uses a linear function f(x) = b x (where both x and f(x) have dimension of length) on [0, 1], computes the RL fractional integral analytically, constructs a generalized Cantor set with dimension alpha on the x-axis, approximates the fractal's content under the curve, verifies scaling under zoom factor lambda, and applies normalization to match the fractional integral value.

The script generates console diagnostics for verification and a high-resolution plot (`plot.png`) overlaying the function, its fractional integral, and a visualization of the fractal region (using fuzzy shaded bars for the Cantor set).

## Prerequisites

- MATLAB (tested on R2020a and later).
- No additional toolboxes required (uses built-in functions like `gamma`, `linspace`, etc.).

## How to Run the Script

1. **Open MATLAB**: Launch MATLAB and navigate to the directory containing the script.

2. **Customize Parameters**: Edit the parameters at the top:
   - `alpha`: Fractional order (0 < alpha < 1; e.g., 0.5 for dimension 1.5).
   - `lambda`: Scaling factor for verification (e.g., 2).
   - `x_max`: Upper bound of domain [0, x_max] (default 1).
   - `N_samples`: Number of Monte Carlo samples for fractal approximation (e.g., 1e6 for accuracy).
   - `iterations`: Depth of recursive sampling for the Cantor set (e.g., 30 for fine resolution).
   - `slope`: Slope b for f(x) = b x (e.g., 1 or 5; assumes no intercept for pure scaling).

3. **Run the Script**: In the MATLAB Command Window, type the script name and press Enter, or use the Run button in the Editor.

4. **Output**:
   - **Console Diagnostics**: Displays values like RL integrals, fractal contents, scaling ratios, and normalization factor.
   - **Plot**: Generates `plot.png` (2400x1800 pixels) with:
     - Blue: f(x) and left y-axis (blue ticks).
     - Red: I^alpha f(x) and right y-axis (red ticks).
     - Gray fuzzy bars: Visualization of the fractal region under the curve.
   - Example diagnostics for alpha=0.5, lambda=2, slope=1:
     - RL Integral at x=1: ~0.7523
     - Scaling Ratio: ~2.8284 (matches lambda^(1+alpha))
     - Normalization Factor: ~2 / Gamma(alpha + 2) ≈ 1.5045

The script runs in seconds on standard hardware; increase `N_samples` or `iterations` for better accuracy at the cost of time.

## Fractal Construction

The fractal is constructed as a generalized Cantor set on the x-axis [0, x_max] with Hausdorff dimension alpha, then extended to the region under f(x) (effective dimension 1 + alpha).

### Construction Details:
- **Scaling Ratio r**: Set to r = 2^(-1/alpha) to achieve dimension alpha = log(2) / log(1/r), ensuring two branches per iteration with contraction r.
- **Recursive Sampling**: The `sample_cantor` function generates points via an iterative self-similar process:
  - Start at t=0, current_scale = x_max.
  - For each iteration (up to `iterations`):
    - With 50% probability, stay in the left branch (add 0).
    - Otherwise, jump to the right branch (add (1 - r) * current_scale).
    - Shrink scale by r.
- This produces samples distributed according to the invariant measure mu (normalized Hausdorff measure) on the Cantor-like set, clustering self-similarly with gaps.

For alpha near 0, the set is sparse (point-like); near 1, it's dense (line-like). Finite iterations approximate the limit set; higher values improve resolution but increase computation.

### Evaluation of Fractal Content
- **Approximation**: The "content" (Hausdorff-like measure at dimension 1 + alpha) is approximated as the integral of f(t) over the fractal base with respect to mu, scaled by x_max^alpha:
  - Compute mean E[f(t)] (expected value) over samples (Monte Carlo average).
  - Multiply by x_max^alpha to account for the measure's total mass (since mu is normalized to 1).
- **Normalization**: Factor = exact RL integral / approximated content, typically ~2 / Gamma(alpha + 2) for this construction, making normalized content match I^alpha f(x).
  - **Why 2 / Gamma(alpha + 2)?**: For f(x) = b x, the exact RL is b / Gamma(alpha + 2). The fractal content approximates b * (1/2) due to symmetric mu with E[t] = 1/2. The factor is thus [b / Gamma(alpha + 2)] / (b / 2) = 2 / Gamma(alpha + 2), independent of b (exact in limit; approximate with finite samples/iterations).
- **Scaling Verification**: Repeat for scaled domain [0, lambda x_max] and scaled function g(y) = lambda f(y / lambda), confirming ratio ~ lambda^(1+alpha).

## Conceptual Points

- **Dimensional Scaling**: Functions with length dimension scale as lambda^1; their integrals as lambda^2. Fractional integrals blend this via lambda^(1+alpha), mirroring fractal scaling lambda^D with D = 1 + alpha.
- **Local Linearity**: Smooth functions are locally linear (f(x) ≈ b (x - x_0)), so the analogy extends via Taylor approximation, focusing on the slope-driven term for primary scaling.
- **Fractal Interpretation**: The RL fractional integral can be viewed as integrating f over a fractal subset of the x-axis with dimension alpha, capturing the "blended" dimension. This geometric link aids understanding intermediate alpha.
- **Literature Connection**: This interpretation draws from works like Nigmatullin's, who established relationships between fractional integrals and Cantor sets, showing fractional calculus can model fractal geometries physically. Further explorations include precise mappings between fractals and fractional operators.

For extensions (e.g., non-linear f, asymmetric fractals), consult the literature and modify the script.