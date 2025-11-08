# Readme

This MATLAB script is designed to numerically compute and visualize the left or right Riemann-Liouville (RL) fractional derivative of a function provided as discrete, uniformly spaced points in a CSV file. It uses the Grünwald-Letnikov approximation method, with enhancements for improved accuracy, and generates a plot comparing the original function, its first derivative, and the fractional derivative. The results are saved back to a CSV file and as a high-resolution PNG image.

See rlfractionalderivative.md for details.

The best introduction to fractional derivatives and fractional calculus that I know of is in ["*Fractional Derivative Modeling in Mechanics and Engineering*](https://link.springer.com/book/10.1007/978-981-16-8802-7)" (2022). 

I have written this script with AI assistance from Grok 4.

Plots for simple test functions and other functions are in the directory Plots-August-2025.

## Simple test functions

I have tested the script with functions for which the RL fractional derivative is known analytically. There's good numeric agreement, so the script works. However, I've realized that this computation is often sensitive to numeric errors (see rlfractionalderivative.md), so more sophisticated methods like those described in "[*Fractional Calculus: High-Precision Algorithms and Numerical Implementations*](https://link.springer.com/book/10.1007/978-981-99-2070-9)" (2024) are likely needed.

### lx(02)(0--5).png

Left RL fractional derivative of f(x) = x. Order 0.2. Integration limits 0 -- 5.

### lx(08)(0--5).png

Left RL fractional derivative of f(x) = x. Order 0.8. Integration limits 0 -- 5.

### lx2(02)(0--5).png

Left RL fractional derivative of f(x) = x^2. Order 0.2. Integration limits 0 -- 5.

### lx2(08)(0--5).png

Left RL fractional derivative of f(x) = x^2. Order 0.8. Integration limits 0 -- 5.

### lsinx(02)(-4pi--4pi).png

Left RL fractional derivative of f(x) = sin(x). Order 0.2. Integration limits -4pi -- 4pi.

### lsinx(08)(-4pi--4pi).png

Left RL fractional derivative of f(x) = sin(x). Order 0.8. Integration limits -4pi -- 4pi.

### lconst(02)(0--5).png

Left RL fractional derivative of f(x) = 1. Order 0.2. Integration limits 0 -- 5.

### lconst(08)(0--5).png

Left RL fractional derivative of f(x) = 1. Order 0.8. Integration limits 0 -- 5.

## Other tests

There's no simple analytic expression of the RL fractional derivative of the functions below, but these plots seem plausible.

### ltanhx(02)(-10--10).png

Left RL fractional derivative of f(x) = tanh(x). Order 0.2. Integration limits -10 -- 10.

### ltanhx(08)(-10--10).png

Left RL fractional derivative of f(x) = tanh(x). Order 0.8. Integration limits -10 -- 10.

### latanx(02)(-10--10).png

Left RL fractional derivative of f(x) = 2atan(4x)/pi. Order 0.2. Integration limits -10 -- 10.

### latanx(08)(-10--10).png

Left RL fractional derivative of f(x) = 2atan(4x)/pi. Order 0.8. Integration limits -10 -- 10.

## Left and right derivatives

All plots above are for the left RL fractional derivative. The plots below are for the right RL fractional derivative and illustrate the relationship between the two. Compare these right derivatives to the left ones.

### rtanhx(08)(-10--10).png

Right RL fractional derivative of f(x) = tanh(x). Order 0.8. Integration limits -10 -- 10.

### ratanx(08)(-10--10).png

Right RL fractional derivative of f(x) = 2atan(4x)/pi. Order 0.8. Integration limits -10 -- 10.