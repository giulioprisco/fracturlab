# Explanation of the MATLAB Script

This MATLAB script generates a visualization to illustrate the effect of the Riemann-Liouville fractional derivative on the Heaviside step function. The Heaviside function is defined as 0 for (x < x_0) and A for (x >= x0), where (x0 = 0) in this case. The script plots the constant portions of multiple Heaviside functions (for different amplitudes \(A\)) and overlays "rotated" tangent lines at specific points, where the rotation (slope) is determined by the fractional derivative value, highlighting how fractional calculus incorporates memory of the function's past.

This script has been developed and tested with AI assistance from Grok.

## Key Parameters and Setup
- **Fractional Order (alpha)**: Set to 0.9, which is close to 1 (the standard integer derivative).
- **Gamma Value**: Computed as gamma(1 - alpha) for use in the fractional derivative formula.
- **Amplitudes (A)**: [1, 2, 3, 4], representing different "past histories" via the step height.
- **Tangent Points (x_tang)**: [1, 1.5], points where the fractional derivative is evaluated and visualized.
- **Plot Range**: From (x = 1 to (x = 2), focusing on the region after the step at x_0 = 0.
- **Colors**: Green, red, magenta, blue for the tangent lines corresponding to increasing (A).
- **Arrow Length**: 0.3 units, the fixed length of the tangent line segments.

## What the Script Does Step by Step
1. **Figure Setup**: Creates a figure with specific dimensions (8x6 inches) for high-resolution export at 300 DPI.

2. **Loop Over Amplitudes \(A\)**:
   - For each (A), computes the constant Heaviside value across the plot range (a horizontal line at y = A).
   - Plots this horizontal line in black with a line width of 2, representing the standard derivative (slope 0 for (x > x_0)).

3. **Tangent Line Calculation and Plotting**:
   - For each tangent point \(x\) in [1, 1.5]:
     - Calculates the fractional derivative using the analytic formula: (deriv = A cdot (x - x_0)^{-alpha} / gamma(1 - alpha)).
     - Normalizes and scales a direction vector [1, deriv] to a fixed length (0.3 units).
     - Plots a short line segment starting at (x, A) in the corresponding color, with slope equal to the fractional derivative value. This visualizes the "rotation" from the horizontal standard tangent.

4. **Final Plot Elements**:
   - Adds an x-label ('x').
   - Enables grid lines.
   - Exports the plot as 'plot.png' at high resolution.

## Interpretation
The script demonstrates the non-local nature of fractional derivatives. The standard (integer-order) derivative of the Heaviside function for (x > x_0) is zero, resulting in horizontal tangents. However, the fractional derivative introduces a non-zero slope that depends on the distance from the discontinuity at (x_0) and the step amplitude (A), effectively "rotating" the tangent away from horizontal. This rotation diminishes as (alpha) approaches 1 or as (x) increases far from (x_0), but persists due to the memory effect in fractional calculus.

For the simple case of the Heaviside function, this shows how the fractional derivative of order alpha of a function, for alpha close to 1, rotates from the standard derivative of the function in a way that depends on the past of the function, represented in this case by x0 and A.