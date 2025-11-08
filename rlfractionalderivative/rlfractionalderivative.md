# MATLAB Script for Computing Riemann-Liouville Fractional Derivatives

This MATLAB script is designed to numerically compute and visualize the left or right Riemann-Liouville (RL) fractional derivative of a function provided as discrete, uniformly spaced points in a CSV file. It uses the Grünwald-Letnikov approximation method, with enhancements for improved accuracy, and generates a plot comparing the original function, its first derivative, and the fractional derivative. The results are saved back to a CSV file and as a high-resolution PNG image.

## Key Features
- **Input Data Handling**: Reads x and y values from `values.csv` (assuming two columns: x in the first, y in the second).
- **User Inputs**:
  - Fractional order α (0 < α < 1).
  - Choice of derivative: "left" or "right".
  - Option to use the shifted Grünwald-Letnikov approximation for higher-order accuracy (O(h²) vs. standard O(h)).
  - Optional zoom range for the plot ([minimum and maximum x] or full range if left empty).
- **Numerical Improvements**:
  - Recursive computation of generalized binomial coefficients to prevent gamma function overflow for large datasets.
  - Shifted approximation option reduces truncation errors.
  - Higher-order (4th-order central difference) approximation for the first integer derivative to enhance precision and reduce oscillations.
  - Boundary handling assumes the function is zero outside the provided interval to mitigate edge effects.
- **Computations**:
  - Approximates the first integer derivative using a 4th-order central difference scheme (with forward/backward at boundaries) for plotting.
  - Computes the RL fractional derivative via finite sums in loops (backward for left, forward for right).
- **Output**:
  - Overwrites `values.csv` with x and the computed fractional derivative values.
  - Generates a plot without labels: original function (green solid line), first derivative (red solid line), and fractional derivative (blue solid line), with a grid.
  - Saves the plot as `plot.png` at 300 DPI, sized for 2400x1800 pixels (8x6 inches).
- **Assumptions and Limitations**:
  - Uniform grid spacing (h calculated from x(2) - x(1)).
  - Suitable for smooth functions; boundary singularities may appear for non-zero f at limits.
  - Computationally intensive (O(N²)) for large N; vectorization could optimize further.
  - No support for high-precision variables (e.g., vpa) in the current version.

## Summary of Observations on Numerical Accuracy
Tests reveal that accuracy doesn't always improve with the number of points N. If N is too large, the results can exhibit degradation due to  round-off errors dominating over truncation errors as step size h decreases. This instability stems from subtractive cancellation in alternating sums (condition number O(1/h^α)), amplified in double-precision arithmetic. The method's convergence is also limited by boundary effects.

## Usage Workflow
1. Prepare `values.csv` with your function data (e.g., from a prior script like `inputfunction.m`).
2. Run the script in MATLAB.
3. Respond to prompts for α, derivative side, shifted option, and zoom range.
4. View the generated `plot.png` and updated `values.csv`.

This script is part of a modular workflow for fractional calculus explorations. For best accuracy use large intervals and many points to approximate infinite limits, but monitor for round-off issues as described.